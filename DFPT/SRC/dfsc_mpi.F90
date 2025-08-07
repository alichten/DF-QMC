program dfpt_mpi

  ! ==========================================================
  ! DF-QMC for d-wave Supterconductivity in Square Lattice
  ! Dual Fermion Perturbation Theory DQMC-code fo d-wave SC 
  ! developed by
  ! 
  ! Alexander Lichtenstein 
  ! alexander.lichten@uni-hamburg.de
  ! 
  ! with help for QUEST-version by Edin Kapetanović (PhD)
  ! ==========================================================

  use dqmc_cfg
  use dqmc_geom_wrap
  use dqmc_hubbard
  use dqmc_tdm1
  ! 
  use dfpt_tools
  use mpi

  implicit none
  
  ! Original declarations
  real                :: t1, t2
  type(config)        :: cfg
  type(Hubbard)       :: Hub
  type(GeomWrap)      :: Gwrap
  type(tdm1)          :: tm
  type(Gtau)          :: tau
  character(len=slen) :: gfile
  logical             :: tformat
  integer             :: ii, jj, kk, slice, nhist, comp_tdm
  integer             :: nBin, nIter
  character(len=60)   :: ofile  
  integer             :: OPT
  !integer             :: HSF_output_file_unit
  integer             :: symmetries_output_file_unit
  integer             :: FLD_UNIT, TDM_UNIT
  real(wp)            :: randn(1)
  ! 
  ! =======================================
  
  ! 
  ! My declarations
  ! 
  
  ! MPI Variables
  integer :: mpi_size, mpi_rank, ierr, rc, master, iran
  
  
  
  integer             :: tmp, num_d, nMeas, N, L
  
  ! GF-Matrix, R- & K-vector
  complex(wp), allocatable :: GF_up_full(:,:), GF_dn_full(:,:)
  real(wp)   , allocatable :: k_vecs(:,:), r_vecs(:,:)
  
  ! F1 will contain Fourier coefficients, F2 will be its 
  ! hermitian conjugate ( for (i,j) -> (k,k') FT )
  complex(wp), allocatable :: F1(:,:), F2(:,:)
  ! 
  ! Reference GF's from "build_ref"-run, real and k-space
  real(wp)   , allocatable :: ref_GF_d_tau(:,:)
  real(wp)   , allocatable :: ref_GF_ij_tau(:,:,:)
  complex(wp), allocatable :: ref_GF_kkp_tau(:,:,:)
  real(wp)   , allocatable :: ref_GF_full(:,:)
  
  ! Bare dual GF and (t',mu) for perturbation
  complex(wp), allocatable :: dual_G0(:,:), dual_G0_full(:)
  complex(wp), allocatable :: dual_G0sc(:,:), dual_G0sc_full(:)
  complex(wp), allocatable :: Sdual(:), SdualMPI(:), lattG(:)
  complex(wp), allocatable :: Sdualsc(:), SdualscMPI(:), lattGsc(:)
  real(wp) :: tp, mu, hsc
  
  ! Number of (positive) Matsubara frequencies
  integer  :: Nom
  
  ! Array for "minus" k
  integer, allocatable :: mk(:)
  
  ! ===============================
  ! 
  ! Get in contact with MPI
  ! 
  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS .and. mpi_rank .eq. master) then
      write(*,*) 'MPI_INIT returned an error. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if
  ! 
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_size,ierr)
  master = 0
   
  ! ===============================
  
  
  ! =========== Original Initialization stuff from QUEST-"ggeom", no changes ==================
 
  ! Timer
  call cpu_time(t1)
  !Read input
  call DQMC_Read_Config(cfg)
  !Get output file name header
  call CFG_Get(cfg, "ofile", ofile)
  !Get general geometry input
  call CFG_Get(cfg, "gfile", gfile)
  !Save whether to use refinement for G used in measurements.
  call CFG_Get(cfg, "nhist", nhist)
  !if (nhist > 0) then
  !   call DQMC_open_file(adjustl(trim(ofile))//'.HSF.stream','unknown', HSF_output_file_unit)
  !endif
  call DQMC_open_file(adjustl(trim(ofile))//'.geometry','unknown', symmetries_output_file_unit)
  !Determines type of geometry file
  call DQMC_Geom_Read_Def(Hub%S, gfile, tformat)
  if (.not.tformat) then
     !If free format fill gwrap
     call DQMC_Geom_Fill(Gwrap, gfile, cfg, symmetries_output_file_unit)
     !Transfer info in Hub%S
     call DQMC_Geom_Init(Gwrap,Hub%S,cfg)
  endif
  call DQMC_Geom_Print(Hub%S, symmetries_output_file_unit)

  ! Initialize the rest data
  call DQMC_Hub_Config(Hub, cfg)

  ! Perform input parameter checks
  if (Hub%nTry >= Gwrap%Lattice%nSites .and. mpi_rank .eq. master) then
    write(*,*)
    write(*,"('  number of lattice sites =',i5)") Gwrap%Lattice%nSites
    write(*,"('  ntry =',i5)") Hub%nTry
    write(*,*) " Input 'ntry' exceeds the number of lattice sites."
    write(*,*) " Please reset 'ntry' such that it is less than"
    write(*,*) " the number of lattice sites."
    write(*,*) " Program stopped."
    stop
  end if

  ! Initialize time dependent properties if comp_tdm > 0
  call CFG_Get(cfg, "tdm", comp_tdm)
  if (comp_tdm > 0) then
     !call DQMC_open_file(adjustl(trim(ofile))//'.tdm.out','unknown', TDM_UNIT)
     call DQMC_Gtau_Init(Hub, tau)
     call DQMC_TDM1_Init(Hub%L, Hub%dtau, tm, Hub%P0%nbin, Hub%S, Gwrap)
  endif
  ! ======================== Original  Initialization stuff end =============================
  
  
  
  
  ! =========================================================================================
  ! 
  ! DFPT initializations
  ! 
      
  N = tau%n
  L = tm%L

  ! Arrays which store possible R-& K-vectors
  allocate(k_vecs(N,2))
  allocate(r_vecs(N,2))
  
  ! Matrices for the spatial FT
  allocate(F1(N,N))
  allocate(F2(N,N))
    
  ! Get R- & K-vectors
  call build_r_vecs(tau%n,r_vecs)
  call build_k_vecs(tau%n,k_vecs)
  
  if (mpi_rank .eq. master) then
    write(*,*) "Lattice sites and coordinates:"
    do tmp = 1, N
        write(*,'(I6,3X,F6.3,3X,F6.3)') tmp, r_vecs(tmp,1), r_vecs(tmp,2)
    end do
  
    write(*,*) "K-Vectors:"
    do tmp = 1, N
      write(*,'(I6,3X,F6.3,3X,F6.3)') tmp, k_vecs(tmp,1), k_vecs(tmp,2)
    end do
  end if
  
  ! Build Matrix for FT from (R,R') to (k,k')
  call build_FT_matrix(F1,F2,r_vecs,k_vecs,N)
   
  ! Initialize Reference G and Read from Files
  num_d = tm%properties(1)%nclass
  allocate(ref_GF_d_tau(num_d, 0:L-1))
  call DFPT_Read_Ref_G(tm, ref_GF_d_tau)
  
  ! Extend Reference G from (D,tau)- to full spatial (i,j,tau)-dependence
  allocate(ref_GF_ij_tau(N,N,0:L-1))
  call DFPT_Extend_Ref_G_ij(ref_GF_d_tau,ref_GF_ij_tau,tm)
  
  ! Extend to full (tau,tau')-dependence (used for the substraction)
  allocate(ref_GF_full(N*L,N*L))
  call DFPT_Extend_Ref_G_ttp(ref_GF_ij_tau,tm,ref_GF_full)
  
  ! Initialize the array for measured GFs (and g_tilde) in the same shape
  allocate(GF_up_full(N*L,N*L))
  allocate(GF_dn_full(N*L,N*L))
  GF_up_full = ZERO
  GF_dn_full = ZERO
  
  ! Get indexing for "minus" k
  allocate(mk(N*L))
  call DFPT_minusk(mk, N, L)
  
  ! FT from (i,j,tau) -> (k,k',tau) for Reference GF
  allocate(ref_GF_kkp_tau(N,N,0:L-1))
  call GF_R_to_K(ref_GF_ij_tau,ref_GF_kkp_tau,F1,F2,L)
  
  ! 
  ! Build "bare" dual G0 in (k,iw)-space for fixed perturbation (t',mu)
  ! 
  ! Note: Nom should stay L/2 here!
  ! 
  Nom = L / 2
  ! 
  call DFPT_Read_Pert_Params(tp,mu,hsc)
  ! 
  allocate(dual_G0(N,0:Nom))
  allocate(dual_G0_full(N*L))
  allocate(Sdual(N*L))
  allocate(SdualMPI(N*L))
  ! For Nambu-Gor'kov part
  allocate(dual_G0sc(N,0:Nom))
  allocate(dual_G0sc_full(N*L))
  allocate(Sdualsc(N*L))
  allocate(SdualscMPI(N*L))
  ! 
  dual_G0(:,:)    = cmplx(1.0,0.0,wp)
  dual_G0_full(:) = cmplx(1.0,0.0,wp)
  Sdual        = ZERO
  SdualMPI     = ZERO
  !
  dual_G0sc(:,:)    = cmplx(1.0,0.0,wp)
  dual_G0sc_full(:) = cmplx(1.0,0.0,wp)
  Sdualsc        = ZERO
  SdualscMPI     = ZERO
  !  
  call DFPT_build_dual_G0(tp,mu,hsc,N,L,k_vecs,ref_GF_kkp_tau,dual_G0,dual_G0_full,dual_G0sc,dual_G0sc_full,Nom,tm%dtau)
  
  ! =========================================================================
  ! 
  ! New Seed for each MPI process!
  ! (In the most trivial way for now, change later)
  ! 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Hub%seed = Hub%seed + mpi_rank
  Hub%idum = Hub%idum + mpi_rank
  write(*,*) "Seed ran1", Hub%seed, "PC", mpi_rank
  !write(*,*) "Seed ran2", Hub%idum, "PC", mpi_rank
  ! 
  ! =========================================================================
  
  
 
 
  ! If no sweeps, just stop the program
  if (Hub%nWarm + Hub%nPass == 0) then
    write(*,*) "Number of Sweeps is 0!"
    stop
  endif
  
  ! Warmup sweep
  do ii = 1, Hub%nWarm
     if (mpi_rank .eq. master) then
        if (mod(ii, 10)==0.and.mpi_rank .eq. master) write(*,'(A,I6,1X,I3)')' Warmup Sweep, nwrap  : ', ii, Hub%G_up%nwrap
     end if
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)
     call DQMC_Hub_Sweep2(Hub, Hub%nTry)
  end do
  
  
  ! We divide all the measurement into nBin,
  ! each having nPass/nBin pass.
  nBin   = Hub%P0%nBin
  nIter  = Hub%nPass / Hub%tausk / nBin
  nMeas  = Hub%nPass / Hub%tausk
  if (nIter > 0) then
     do ii = 1, nBin
        do jj = 1, nIter
           do kk = 1, Hub%tausk
              call DQMC_Hub_Sweep(Hub, NO_MEAS0)
              call DQMC_Hub_Sweep2(Hub, Hub%nTry)
           enddo

           ! Fetch a random slice for (equal time)-measurement 
           call ran0(1, randn, Hub%seed)
           slice = ceiling(randn(1)*Hub%L)
           ! Output of Master
           if (mpi_rank .eq. master) then
               write(*,'(A,3I6)') ' Measurement Sweep, QMC, bin, slice : ', jj, ii,  slice
           end if

           if (comp_tdm > 0) then
              ! 
              ! Compute A-Matrix, which contains full GF for specific (l,l').
              ! 
              ! Note: - Contains all (l,l') if nOrth = 1!
              ! 
              call DQMC_Gtau_LoadA(tau, TAU_UP, slice, Hub%G_up%sgn)
              call DQMC_Gtau_LoadA(tau, TAU_DN, slice, Hub%G_dn%sgn)
              ! Measure equal-time properties
              call DQMC_Hub_FullMeas(Hub, tau%nnb, tau%A_up, tau%A_dn, tau%sgnup, tau%sgndn)
              ! 
              ! =================
              ! Newest version
              ! =================
              ! 
              ! From loaded A-Matrix: Measure full G(i,j,t,t') from current HS-Config; substract
              ! G_ref to obtain G_tilde, perform 6D-FFT and evaluate the dual Self-Energy contribution
              ! 
              call DFPT_Meas_Dual_Sigma(tm,tau,mk,ref_GF_full,dual_G0_full,dual_G0sc_full,GF_up_full,GF_dn_full,Sdual,Sdualsc)
              ! 
           else if (comp_tdm == 0 .and. mpi_rank .eq. master) then
              write(*,*) "TDM parameter is zero! Stopping program..."
              stop
           endif

           !Write fields 
           !if (nhist > 0) call DQMC_Hub_Output_HSF(Hub, .false., slice, HSF_output_file_unit)
        end do

        ! Stuff here belongs to original "ggeom"
        ! Accumulate results for each bin
        call DQMC_Phy0_Avg(Hub%P0)
        call DQMC_tdm1_Avg(tm)

  
        if (Hub%meas2) then
           if(Hub%P2%diagonalize)then
             call DQMC_Phy2_Avg(Hub%P2, Hub%S)
           else
             call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
           endif
        end if

     end do
  endif
  
  ! ===================
  ! 
  ! MPI: Collect and reduce everything to master process
  ! 
  if (mpi_rank .eq. master) write(*,*) "Sweeps over! From Rank: ", mpi_rank
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  call MPI_REDUCE(Sdual, SdualMPI, N*L ,MPI_DOUBLE_COMPLEX, &
               &  MPI_SUM, master, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(Sdualsc, SdualscMPI, N*L ,MPI_DOUBLE_COMPLEX, &
               &  MPI_SUM, master, MPI_COMM_WORLD, ierr)
               
  if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'MPI_REDUCE returned an error.'
  end if
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  ! ===================
  
  if (mpi_rank .eq. master) then
    ! Normalize by Number of Measurements and MPI-Ranks
    Sdual(:) = SdualMPI(:) / nMeas / mpi_size
    Sdualsc(:) = SdualscMPI(:) / nMeas / mpi_size
    ! Get Lattice Green's function
    allocate(lattG(N*L))
    allocate(lattGsc(N*L))    
    lattG = ZERO
    call DFPT_Get_Lattice_G(tp,mu,hsc,N,L,k_vecs,ref_GF_kkp_tau,Sdual,Sdualsc,Nom,tm%dtau,lattG,lattGsc)
    ! Write to File
    call DFPT_Write_Lattice_G(N,L,tm%dtau,lattG,lattGsc,Sdual,Sdualsc)
  end if
  
  

  ! Clean up the used storage
  call DQMC_TDM1_Free(tm)
  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)
  
  ! ===================
  ! 
  ! End MPI
  ! 
  call MPI_FINALIZE(ierr)
  
  ! ===================
  
  call cpu_time(t2)
  if (mpi_rank .eq. master) then
  write(STDOUT,*) "Running time:",  t2-t1, "(second)"
  end if !master


end program dfpt_mpi

