module dfpt_tools

  use dqmc_cfg
  use dqmc_geom_wrap
  use dqmc_hubbard
  use dqmc_tdm1
  use dqmc_util
  use dqmc_struct
  use dqmc_gtau

  implicit none

contains
 
  subroutine get_coords(ii,N,vec_i)
    ! 
    ! Purpose
    ! =======
    !    Returns the (x,y) Coordinates for site "ii". N is the total
    !    number of sites. For now, only works for square lattice.
    ! 
    integer, intent(in) :: ii, N
    real(wp), intent(out) :: vec_i(2)
    ! 
    integer :: mod_num, tmp
    
    ! Evaluate
    mod_num = int(sqrt(float(N)))
    tmp = floor(float( (ii-1) / mod_num))
    vec_i(1) = float((ii-1) - tmp * mod_num)
    vec_i(2) = float(tmp)
  
  end subroutine get_coords
  
  subroutine get_latt_vec(ii,jj,N,vec_ij)
    ! 
    ! Purpose
    ! =======
    !    Returns the vector (xx,yy) from site ii to site jj on a periodic
    !    N-site square lattice
    ! 
    integer, intent(in)   :: ii, jj, N
    real(wp), intent(out) :: vec_ij(2)
    ! 
    real(wp) :: vec_i(2), vec_j(2)
    integer  :: length, half_len
    
    ! Get Lattice vectors and get difference
    call get_coords(ii,N,vec_i)
    call get_coords(jj,N,vec_j)
    vec_ij(:) = vec_j(:) - vec_i(:)
    
    ! Periodicity
    length   = int(sqrt(real(N)))
    half_len = length / 2
    ! X component
    if (vec_ij(1) > half_len) then
      vec_ij(1) = vec_ij(1) - length
    endif
    ! Y component
    if (vec_ij(2) > half_len) then
      vec_ij(2) = vec_ij(2) - length
    endif
  
  end subroutine get_latt_vec
  
  
  subroutine linspace(x, a, b, n, endpoint)
    ! 
    ! Basically the same as numpy.linspace() in Python,
    ! might be useful
    ! 
    real(wp), intent(out) :: x(n)
    real(wp), intent(in) :: a, b
    integer, intent(in) :: n
    logical, intent(in) :: endpoint

    real(wp) :: dx
    integer :: i

    if ( .not. endpoint ) then
      dx = (b-a) / real(n, wp)
    else
      dx = (b-a) / real(n-1, wp)
    end if

    x = [(i*dx+a, i = 0, n-1)]

  end subroutine linspace 
  
  
  subroutine build_k_vecs(N,k_vecs)
    ! 
    ! Builds the possible K vectors (square lattice, full BZ)
    ! 
    integer, intent(in) :: N
    real(wp), allocatable, intent(inout) :: k_vecs(:,:)
    ! 
    real(wp), allocatable :: k_tmp(:)
    real(wp) :: ONEPI = 3.14159265359_wp
    integer :: ii, tmp, length
    
    ! Allocate for N Sites and x/y component
    length = int(sqrt(real(N)))
    allocate(k_tmp(length))
    
    ! Build K-Vectors
    !call linspace(k_tmp,-ONEPI,ONEPI,length,.false.)
    call linspace(k_tmp, 0.0_wp, 2.0_wp * ONEPI,length,.false.)
    tmp = 0
    !write(*,*) "K-Vectors:"
    do ii = 1, N
      k_vecs(ii,1) = k_tmp(ii-tmp*length)
      k_vecs(ii,2) = k_tmp(tmp+1)
      !write(*,'(F6.3,3X,F6.3)') k_vecs(ii,1), k_vecs(ii,2)
      if (mod(ii,length)==0) tmp = tmp + 1
    end do
    
  end subroutine build_k_vecs
  
  subroutine build_r_vecs(N,r_vecs)
    ! 
    ! Builds the possible Site vectors and collects them in one array
    ! 
    integer, intent(in) :: N
    real(wp), allocatable, intent(inout) :: r_vecs(:,:)
    ! 
    integer :: ii
    ! 
    !write(*,*) "R-Vectors:"
    do ii = 1, N
      call get_coords(ii,N,r_vecs(ii,:))
      !write(*,'(I6,3X,F6.3,3X,F6.3)') ii, r_vecs(ii,1), r_vecs(ii,2)
    end do
    
  end subroutine build_r_vecs
  
  
  subroutine build_FT_matrix(F1,F2,r_vecs,k_vecs,N)
    ! 
    ! Builds the (k,i)-dependent matrix with Fourier-coefficients exp(-i*R_i*k) / sqrt(N) in F1.
    ! F2 is its conjugate transpose.
    ! 
    integer, intent(in)                     :: N
    real(wp), intent(in)                    :: r_vecs(:,:), k_vecs(:,:)
    complex(wp), allocatable, intent(inout) :: F1(:,:), F2(:,:)
    ! 
    integer :: ii, kk
    complex(wp) :: tmp(N,N)
    ! 
    do kk = 1, N
      do ii = 1, N
        F1(kk,ii) = exp( - dcmplx(0.d0,1.d0) * dot_product(k_vecs(kk,:),r_vecs(ii,:)) )
      end do
    end do
    F1 = F1 / sqrt(float(N))
    ! 
    tmp(:,:) = conjg(F1)
    F2(:,:)  = transpose(tmp)
  
  end subroutine build_FT_matrix
  
  
  
  subroutine GF_R_to_K(gf_ij_tau,gf_kkp_tau,F1,F2,L)
  ! 
  ! Transforms G_(i,j)(tau) to G_(k,k')(tau)
  ! 
  real(wp)   , intent(in)   , allocatable :: gf_ij_tau(:,:,:)
  complex(wp), intent(inout), allocatable :: gf_kkp_tau(:,:,:)
  ! 
  complex(wp), intent(in) :: F1(:,:), F2(:,:)
  integer, intent(in) :: L
  ! 
  integer :: it
  ! 
  do it = 0, L-1
    gf_kkp_tau(:,:,it) = matmul( matmul( F1, gf_ij_tau(:,:,it) ), F2 )
  end do
  
  end subroutine GF_R_to_K
  
  
  
  subroutine DFPT_Meas_Ref_G(T1, tau)
    ! 
    ! Purpose
    ! =======
    ! 
    !    Modified Measurement Routine, adapted from the original in DQMC_TDM1,
    !    but only measures G(tau) so it can be saved as the "reference" system
    !    
    !    - Assumes that "Gtau_LoadA" has already been done.
    ! 
    
    
    ! Original arguments
    type(TDM1), intent(inout)   :: t1
    type(Gtau), intent(inout)   :: tau
    ! 
    ! My arguments
    real(wp), allocatable :: gf_tmp(:,:,:)
    
    ! ... Local var ...
    integer  :: i, k, m, L, cnt, dt, i0, it, j0, jt, dtau, iprop
    real(wp) :: sgn, factor
    real(wp), pointer :: up0t(:,:)
    real(wp), pointer :: upt0(:,:)
    real(wp), pointer :: dn0t(:,:)
    real(wp), pointer :: dnt0(:,:)
    real(wp), pointer :: up00(:,:)
    real(wp), pointer :: uptt(:,:)
    real(wp), pointer :: dn00(:,:)
    real(wp), pointer :: dntt(:,:)
    real(wp), pointer :: values(:,:,:)
    
    
    ! Allocate temporary GF tensor for measurements

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(tau%north,2)
    m   = (tau%north-k) / 2

    ! Aliases
    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t
    up00 => tau%up00
    uptt => tau%uptt
    dn00 => tau%dn00
    dntt => tau%dntt
    


    blocks: do i0 = 1, tau%nb
       do dtau = 0, tau%nb-1
          it = mod(i0+dtau-1,tau%nb) + 1
          ! Stored value
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up
          
          ! Get GF
          call GF_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
          
          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             ! Get GF
             call GF_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
             
             
          enddo
          
          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif
          

          ! Increment index tau%it
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif
             
             jt = tau%it_up; j0 = tau%i0_up
             ! Get GF
             call GF_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
             
          enddo
          

       enddo

       cnt = cnt + 1
    enddo blocks
 
    if (i0 .ne. tau%nb+1) then
       write(*,*) "Up and down time slices are mismatched. Stop"
       stop
    endif

    ! Normalize by "bond number" F(i) and number "cnt" of GF-Computations
    ! and put result in correct bin
    sgn = tau%sgnup * tau%sgndn
    ! Only spin-average GF
    iprop = 1
       values => T1%properties(iprop)%values
       do it = 0, L-1
          do i = 1, T1%properties(iprop)%nClass
             factor = sgn/(T1%properties(iprop)%F(i)*cnt)
             values(i,it,T1%idx)   = values(i,it,T1%idx)   + factor*values(i,it,T1%tmp)
          end do
       end do
       values(:,:,T1%tmp)   = ZERO
    
    ! Put sign into bin and increment measurement count within bin
    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

  end subroutine DFPT_Meas_Ref_G
  
  
  
  
  subroutine DFPT_Test_Time(T1, tau)
    ! Original arguments
    type(TDM1), intent(inout)   :: t1
    type(Gtau), intent(inout)   :: tau
    
    ! ... Local var ...
    integer  :: i, k, m, L, dt, d0, i0, it, j0, jt, dtau, iprop
    real(wp) :: sgn, factor
    real(wp), pointer :: up0t(:,:)
    real(wp), pointer :: upt0(:,:)
    real(wp), pointer :: dn0t(:,:)
    real(wp), pointer :: dnt0(:,:)
    real(wp), pointer :: up00(:,:)
    real(wp), pointer :: uptt(:,:)
    real(wp), pointer :: dn00(:,:)
    real(wp), pointer :: dntt(:,:)
    real(wp), pointer :: values(:,:,:)
    
    ! 
    integer :: nOrth, cnt
    nOrth = tau%north

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(north,2)
    m   = ( north - k ) / 2

    ! Aliases
    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t
    up00 => tau%up00
    uptt => tau%uptt
    dn00 => tau%dn00
    dntt => tau%dntt
    
    do i0 = 1, tau%nb
       do dtau = 0, tau%nb-1
          it = mod(i0+dtau-1,tau%nb) + 1
          ! Current stored value, starting point
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up
          ! Compute GF here
          !write(*,*) jt, j0
          
          ! ===============
          ! First: i0 fixed, go through all it
          ! ===============
          ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             ! Compute GF here
             !write(*,*) jt, j0
             
          enddo
          
          ! Go back to starting point
          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif         

          ! Decrement index tau%it (via TMINUS)
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif
             
             jt = tau%it_up; j0 = tau%i0_up
             ! Compute GF here
             !write(*,*) jt, j0
          enddo
          
          ! ===============
          ! Same as above, but shift i0!
          ! Note: Come up with a more compact way of writing this
          ! ===============
          
          ! =================================================================================
          ! Increment index tau%i0 (ZPLUS) - If nOrth is even do only north/2-1 increments
          do d0 = 1, m-1+k
          
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0! (d0 times!)
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_UP)
            enddo
            
            jt = tau%it_up; j0 = tau%i0_up
            ! Compute GF here
            !write(*,*) jt, j0
            
            ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
            do dt = 1, m-1+k
              call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              ! Compute GF here
              !write(*,*) jt, j0
            enddo
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0 back to correct place!
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_UP)
            enddo
             
            ! Decrement index tau%it (via TMINUS)
            do dt = 1, m
              call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              ! Compute GF here
              !write(*,*) jt, j0
            enddo
            
          enddo
          
          ! ================================================================================
          ! Decrement index tau%i0 (ZMINUS))
          do d0 = 1, m
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0! (d0 times!)
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_UP)
            enddo
            
            jt = tau%it_up; j0 = tau%i0_up
            ! Compute GF here
            !write(*,*) jt, j0
            
            ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
            do dt = 1, m-1+k
              call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              ! Compute GF here
              !write(*,*) jt, j0
            enddo
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0 back to correct place!
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_UP)
            enddo
             
            ! Decrement index tau%it (via TMINUS)
            do dt = 1, m
              call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              ! Compute GF here
              !write(*,*) jt, j0
            enddo
          
          end do
          ! ================================================================================
          

       enddo
    enddo
    !write(*,*) "Program ends!"
    !stop
  
  
  end subroutine DFPT_Test_Time
  
  
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DF-QMC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  subroutine DFPT_Meas_Dual_Sigma(T1,tau,mk,ref_GF_full,dual_G0_full,dual_G0sc_full,GF_up_full,GF_dn_full,Sdual,Sdualsc)
    ! 
    ! Purpose
    ! =======
    !    
    !    Measures current GF, obtains g_tilde (k,k',iw,iw') through 6D-FFT
    !    and evaluates the dual Self-Energy
    !    
    ! 
    ! 
    ! Original arguments
    type(TDM1), intent(inout)   :: T1
    type(Gtau), intent(inout)   :: tau
    ! 
    ! My arguments
    real(wp)   , allocatable, intent(in)    :: ref_GF_full(:,:)
    complex(wp), allocatable, intent(inout) :: GF_up_full(:,:), GF_dn_full(:,:)
    complex(wp), allocatable, intent(inout) :: Sdual(:)
    complex(wp), allocatable, intent(inout) :: Sdualsc(:)
    complex(wp), allocatable, intent(in)    :: dual_G0_full(:)
    complex(wp), allocatable, intent(in)    :: dual_G0sc_full(:)
    integer    , allocatable, intent(in) :: mk(:)
    real(wp)    :: ONEPI, dtau, beta, tau0, tau1
    complex(wp) :: Xi, zf, Hc, Sk, Sd, sdsc, Eq(T1%properties(1)%n*tau%L)
    integer :: N, L, V, Nx, iks, jks, i1, i2, k1, k2, km1, km2
    integer :: NN(6), Ndim
    
    ! ... Local var ...
    integer  :: i0, it, idt, j0, jt, dt, d0, k, m, north, cnt
    integer  :: it1, it0
    real(wp) :: sgn
    real(wp), pointer :: up0t(:,:)
    real(wp), pointer :: upt0(:,:)
    real(wp), pointer :: dn0t(:,:)
    real(wp), pointer :: dnt0(:,:)
    real(wp), pointer :: up00(:,:)
    real(wp), pointer :: uptt(:,:)
    real(wp), pointer :: dn00(:,:)
    real(wp), pointer :: dntt(:,:)
    real(wp), pointer :: values(:,:,:)
    
    ! Signs
    real(wp), pointer :: sgnup, sgndn
    sgnup => tau%sgnup
    sgndn => tau%sgndn

    ! 
    ! Some notes on notation:
    ! L: Number of time slices
    ! N: Number of atoms (Nat in Sasha's code)
    ! V: "Volume" of the space-time lattice
    ! Nx: Number of sites in one direction (sqrt(N) for square 
    !     lattice, called "n" in Sasha's code)
    ! 
    L   = tau%L
    N   = T1%properties(1)%n
    V   = N*L
    Nx  = NINT( sqrt(real(N)) )
    dtau = T1%dtau
    Beta = L * dtau
    ! 
    north = tau%north
    k   = mod(north,2)
    m   = ( north - k ) / 2
    
    
    ! Pi and imag. unit
    ONEPI = 3.14159265359_wp
    Xi    = dcmplx(0.0_wp,1.0_wp)
    
    Ndim = 6
    NN(1) = Nx
    NN(2) = Nx
    NN(3) = L
    NN(4) = Nx
    NN(5) = Nx
    NN(6) = L

    ! Aliases
    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t
    up00 => tau%up00
    uptt => tau%uptt
    dn00 => tau%dn00
    dntt => tau%dntt
    
    
    if ( tau%north == 1 ) then
      ! 
      ! Loop over both time-indices
      ! 
      do i0 = 1, L
          do idt = 0, L-1  
              it = mod(i0+idt-1,L) + 1
              tau0 = (i0 - 1) * dtau
              tau1 = (it - 1) * dtau
              zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
            
              ! Get the Green's functions for (it,i0) from full A-matrix
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              ! (if necessary, for Spin-Dn as well)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)              
              ! Note: Important! The measured GF's (e.g. upt0) are defined in QMC 
              !       notation (with a "+"-sign), while the reference GF already has
              !       the correct sign. Hence, the addition in the brackets and the
              !       "-"-sign before 

              ! Loops over lattice sites
              !do i1 = 1, N
              !    iks = i1 + (it - 1) * N
              !    do i2 = 1, N
              !        jks = i2 + (i0 - 1) * N
                      ! Measure current HS-QMC-g and substract reference g 
              !        GF_up_full(iks,jks) = - zf * ( upt0(i1,i2) + ref_GF_full(iks,jks) ) 
              !        GF_dn_full(iks,jks) = - zf * ( dnt0(i1,i2) + ref_GF_full(iks,jks) ) 
              !    end do
              !end do
              
              ! Cast (i.j)-block directly into big array, avoid double loop over sites
              it1 = (it - 1) * N + 1
              it0 = (i0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              
          end do
      end do
      
    else
      ! ==================================
      ! 
      ! Block for nOrth != 1 - evaluate the whole g(i,j,l,l') through use of T/Z PLUS/MINUS
      ! arguments in DQMC_change_gtau_time function
      ! 
      !===================================
      
      do i0 = 1, tau%nb
        do idt = 0, tau%nb-1
          it = mod(i0+idt-1,tau%nb) + 1
          ! Current stored value, starting point
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up
          tau0 = (j0 - 1) * dtau
          tau1 = (jt - 1) * dtau
          zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
          
          ! Get GF
          it1 = (jt - 1) * N + 1
          it0 = (j0 - 1) * N + 1
          GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
          GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
          
          
          ! ===============
          ! First: i0 fixed, go through all it
          ! ===============
          ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             tau0 = (j0 - 1) * dtau
             tau1 = (jt - 1) * dtau
             zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
             
             ! Get GF
             it1 = (jt - 1) * N + 1
             it0 = (j0 - 1) * N + 1
             GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
             GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
             
          enddo
          
          ! Go back to starting point
          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif         

          ! Decrement index tau%it (via TMINUS)
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif
             
             jt = tau%it_up; j0 = tau%i0_up
             tau0 = (j0 - 1) * dtau
             tau1 = (jt - 1) * dtau
             zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
             
             ! Get GF
             it1 = (jt - 1) * N + 1
             it0 = (j0 - 1) * N + 1
             GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
             GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
          enddo
          
          ! ===============
          ! Same as above, but shift i0!
          ! Note: Come up with a more compact way of writing this
          ! ===============
          
          ! =================================================================================
          ! Increment index tau%i0 (ZPLUS) - If nOrth is even do only up to north/2-1 increments
          do d0 = 1, m-1+k
          
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0! (d0 times! Because DumpA always resets to starting point)
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_UP)
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_DN)
            enddo
            
            jt = tau%it_up; j0 = tau%i0_up
            tau0 = (j0 - 1) * dtau
            tau1 = (jt - 1) * dtau
            zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
            
            ! Get GF
            it1 = (jt - 1) * N + 1
            it0 = (j0 - 1) * N + 1
            GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            
            ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
            do dt = 1, m-1+k
              call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              tau0 = (j0 - 1) * dtau
              tau1 = (jt - 1) * dtau
              zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
              
              ! Get GF
              it1 = (jt - 1) * N + 1
              it0 = (j0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            enddo
            
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0 back to correct place!
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_UP)
              call DQMC_change_gtau_time(tau, ZPLUS, TAU_DN)
            enddo
             
            ! Decrement index tau%it (via TMINUS)
            do dt = 1, m
              call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              tau0 = (j0 - 1) * dtau
              tau1 = (jt - 1) * dtau
              zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
              
              ! Get GF
              it1 = (jt - 1) * N + 1
              it0 = (j0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            enddo
            
          enddo
          
          ! ================================================================================
          ! Decrement index tau%i0 (ZMINUS))
          do d0 = 1, m
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0! (d0 times!)
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_UP)
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_DN)
            enddo
            
            jt = tau%it_up; j0 = tau%i0_up
             tau0 = (j0 - 1) * dtau
             tau1 = (jt - 1) * dtau
             zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
             
             ! Get GF
             it1 = (jt - 1) * N + 1
             it0 = (j0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            
            ! Increment index tau%it (TPLUS) - If nOrth is even do only north/2-1 increments
            do dt = 1, m-1+k
              call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              tau0 = (j0 - 1) * dtau
              tau1 = (jt - 1) * dtau
              zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
              
              ! Get GF
              it1 = (jt - 1) * N + 1
              it0 = (j0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            enddo
            ! Go back to starting point
            if (m .gt. 0) then
              call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
              if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
            endif
            
            ! Shift i0 back to correct place!
            do cnt = 1, d0
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_UP)
              call DQMC_change_gtau_time(tau, ZMINUS, TAU_DN)
            enddo
             
            ! Decrement index tau%it (via TMINUS)
            do dt = 1, m
              call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
              if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
              elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
              endif
              jt = tau%it_up; j0 = tau%i0_up
              tau0 = (j0 - 1) * dtau
              tau1 = (jt - 1) * dtau
              zf = exp( Xi * (tau0+tau1) * ONEPI / Beta ) * dtau * dtau
              
              ! Get GF
              it1 = (jt - 1) * N + 1
              it0 = (j0 - 1) * N + 1
              GF_up_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( upt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
              GF_dn_full( it1:(it1+N-1), it0:(it0+N-1) ) = - zf * ( dnt0(:,:) + ref_GF_full( it1:(it1+N-1), it0:(it0+N-1) ) )
            enddo
          
          enddo
        enddo
      enddo
          ! ================================================================================
      
    end if

    !!!!!!!!!!!!!!!!!!!! DF-QMC part !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now: Call the subroutine to get 6D-FFT of G-Tilde!
    call FFT6D(T1, tau, V, GF_up_full, GF_dn_full)

    ! Get dual self-energy, for now: up to first order
    Hc = dcmplx(0.d0,0.d0)
    
    do k1 = 1, V
        km1 = mk(k1)
        Sk = dcmplx(0.d0,0.d0)
        sdsc=dcmplx(0.d0,0.d0)
        Hc = Hc + ( GF_up_full(k1,km1) + GF_dn_full(k1,km1) ) * dual_G0_full(k1)
        ! 
        do k2 = 1, V
           km2 = mk(k2)
           Sk = Sk + ( GF_up_full(k1,km2) * GF_up_full(k2,km1) &
                   & + GF_dn_full(k1,km2) * GF_dn_full(k2,km1) ) * dual_G0_full(k2)
           sdsc=sdsc - GF_up_full(k1,km2) * GF_dn_full(km1,k2) * dual_G0sc_full(k2)
        end do
        Eq(k1) = Sk
        Sdualsc(k1) = Sdualsc(k1) - real(sdsc) * sgnup*sgndn / (beta * N)**2  
    end do
    
    ! Final evaluation of first-order correction
    do k1 = 1, V
        km1 = mk(k1)
          Sd = ( GF_up_full(k1,km1) + GF_dn_full(k1,km1) ) * Hc - Eq(k1)
        ! Normalization: Factor 2 due to Spin-Average
          Sdual(k1) = Sdual(k1) - Sd * sgnup*sgndn / 2.0_wp / (beta * N)**2
    end do
  end subroutine DFPT_Meas_Dual_Sigma
  
  
  subroutine FFT6D(T1, tau, V, GF_up_full, GF_dn_full)
    ! 
    ! Purpose
    ! =======
    !    
    ! 6D complex FFT from FFTW-3 version   
    !    
    use, intrinsic :: iso_c_binding
    implicit none

    include 'fftw3.f03'

    type(C_PTR) :: plan
    ! 
    type(TDM1), intent(in)   :: T1
    type(Gtau), intent(in)   :: tau
    ! 
    integer :: N, L, V, Nx
    integer :: NN(6), Ndim
    complex(wp) :: GF_up_full(V,V), GF_dn_full(V,V)

 
    L   = tau%L
    N   = T1%properties(1)%n
    Nx  = NINT( sqrt(real(N)) )
  
  
    Ndim = 6
    NN(1) = Nx
    NN(2) = Nx
    NN(3) = L
    NN(4) = Nx
    NN(5) = Nx
    NN(6) = L

    ! FFTW interface old-style Fortran: 
    call dfftw_plan_dft(plan,ndim,nn,GF_up_full , GF_up_full,'FFTW_FORWARD','FFTW_ESTIMATE')
    call dfftw_execute_dft(plan, GF_up_full , GF_up_full)
    call dfftw_execute_dft(plan, GF_dn_full , GF_dn_full)
    call dfftw_destroy_plan(plan)

  
  end subroutine FFT6D
  
  
  
  
  
  subroutine DFPT_MinusK(mk, N, L)
    ! 
    ! Calculates "minus" k: G^*(k)=G(N-k)
    ! - Adapted from Sasha's Code
    !  
    integer, intent(in)    :: N, L
    integer, intent(inout) :: mk(:)
    integer :: Nx, V
    integer :: iL, imL, iX, iY, imX, imY, iK, imK
  
    V   = N*L
    Nx  = NINT( sqrt(real(N)) )
  
    do iL = 0, L-1
        imL = L - 1 - iL
        ! 
        do iY = 0, Nx-1
            if (iY .eq. 0) then
                imY = 0
            else
                imY = Nx - iY
            end if
            ! 
            do iX = 0, Nx-1
                if (iX .eq. 0) then
                    imX = 0
                else
                    imX = Nx - iX
                end if
                ! Indices
                iK  = iX + 1 + iY*Nx + iL*N
                imK = imX + 1 + imY*Nx + imL*N
                mk(iK) = imK
                !write(*,*) iK, mk(iK)
            end do
        end do
    end do
  
  end subroutine DFPT_MinusK
  
 
  subroutine GF_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, it, i0)
  ! 
  ! Purpose
  ! =======
  !      Modified TDM1_Compute-routine which only computes the GF, in order to 
  !      save time.
  ! 
  ! - Only used in build_ref to build the reference system
  ! 
  ! 
    ! Original Arguments
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    real(wp), intent(in)         :: up00(:,:), uptt(:,:)
    real(wp), intent(in)         :: dn00(:,:), dntt(:,:)
    integer, intent(in)          :: it, i0
    

    ! ... Local scalar ...

    integer  :: i, j, k, dt, dt1, dt2
    real(wp), pointer :: value1(:)
    real(wp), pointer :: value2(:)
    real(wp) :: factor


    ! ... Executable ...
    if (.not.T1%compute) return

    dt = it - i0
    if (dt .gt. 0) then
       ! it > i0
       dt1  =  dt
       dt2  =  T1%L - dt
       factor  = 0.25d0
    elseif (dt .lt. 0) then
       ! it < i0
       dt1  =  dt + T1%L
       dt2  = -dt
       factor = -0.25d0
    else
       dt1 = 0
       dt2 = 0
       factor = 0.5d0
    endif
    
    if (dt .ne. 0) then

       ! Spin-Averaged GF
       value1  => T1%properties(1)%values(:, dt1, T1%tmp)
       value2  => T1%properties(1)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(1)%n
          do j = 1, T1%properties(1)%n
             ! k is the distance index of site i and site j
             k = T1%properties(1)%D(i,j)
             value1(k)  = value1(k) + factor*(upt0(i,j) + dnt0(i,j))
             value2(k)  = value2(k) - factor*(up0t(i,j) + dn0t(i,j))
          end do
       end do

       ! Spin-Up GF
       value1  => T1%properties(2)%values(:, dt1, T1%tmp)
       value2  => T1%properties(2)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(2)%n
          do j = 1, T1%properties(2)%n
             ! k is the distance index of site i and site j
             k = T1%properties(2)%D(i,j)
             value1(k)  = value1(k) + 2*factor*upt0(i,j)
             value2(k)  = value2(k) - 2*factor*up0t(i,j)
          end do
       end do

       ! Spin-Down GF
       value1  => T1%properties(3)%values(:, dt1, T1%tmp)
       value2  => T1%properties(3)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(3)%n
          do j = 1, T1%properties(3)%n
             ! k is the distance index of site i and site j
             k = T1%properties(3)%D(i,j)
             value1(k)  = value1(k) + 2*factor*dnt0(i,j)
             value2(k)  = value2(k) - 2*factor*dn0t(i,j)
          end do
       end do

    else

       ! Spin-Averaged GF
       value1  => T1%properties(1)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(1)%n
          do j = 1, T1%properties(1)%n
             ! k is the distance index of site i and site j
             k = T1%properties(1)%D(i,j)
             value1(k)  = value1(k) + factor*(upt0(i,j) + dnt0(i,j))
          end do
       end do

       ! Spin-Up GF
       value1  => T1%properties(2)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(2)%n
          do j = 1, T1%properties(2)%n
             ! k is the distance index of site i and site j
             k = T1%properties(2)%D(i,j)
             value1(k)  = value1(k) + 2*factor*upt0(i,j)
          end do
       end do

       ! Spin-Down GF
       value1  => T1%properties(3)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(3)%n
          do j = 1, T1%properties(3)%n
             ! k is the distance index of site i and site j
             k = T1%properties(3)%D(i,j)
             value1(k)  = value1(k) + 2*factor*dnt0(i,j)
          end do
       end do

    endif

  end subroutine GF_Compute
  
  
  
  
  subroutine DFPT_Print_Ref_G(T1,Hub)
    ! =======
    ! 
    ! Prints the averaged Green's functions (only spin-avg)
    ! to files for later use as the "reference system"
    ! 
    ! =========
    ! 
    type(TDM1), intent(in)    :: T1
    type(Hubbard), intent(in) :: Hub     

    integer             :: it, num_d
    ! Format string
    character(len=10) :: fmtstr
    
    ! Open Files
    open(12, file="GF_av_ref.txt", status="replace", form="formatted", action="write")
   
    num_d = T1%properties(1)%nclass
    !write(*,*) "Number of distance classes: ", num_d
    
    write(fmtstr,"(A1,I0,A6)") "(",num_d,"E16.8)"
    !write(*,*) "Format specifier:" , fmtstr
    
    ! Write Data to Files
    do it = 0, T1%L-1
        ! Spin-Avg GF
        write(12,fmtstr) T1%properties(1)%values(:,it,T1%avg)
    end do        
    
    ! Close Files
    close(12)
    
    ! Write parameters into another file (for sanity check when reading)
    open(12, file="ref_params.txt", status="replace", form="formatted", action="write")
        write(12,"(I0)") T1%L
        write(12,*) T1%dtau
        write(12,"(I0)") num_d
        write(12,"(I0)") Hub%nWarm
        write(12,"(I0)") Hub%nPass
        write(12,"(I0)") Hub%tausk
        write(12,"(I0)") Hub%P0%nBin
    close(12)

  end subroutine DFPT_Print_Ref_G
  
  
  
  
  
  subroutine DFPT_Read_Ref_G(T1, ref_GF)
    ! =======
    ! 
    ! Reads the averaged Green's functions (only spin-avg)
    ! from the "build_ref"-run for use as the "reference system"
    ! 
    ! =========
    ! 
    type(TDM1), intent(in)    :: T1
    real(wp), intent(inout), allocatable :: ref_GF(:,:)
    integer :: it
    ! Just for sanity check
    integer  :: L, num_d
    real(wp) :: dtau
    
    ! Open Files
    open(12, file="GF_av_ref.txt", status="old", form="formatted", action="read")

    ! Read Data
    do it = 0, T1%L-1
        ! Spin-Avg GF
        read(12,*) ref_GF(:,it)
    end do    
      
    ! Close Files
    close(12)
    
    ! Some sanity checks, read parameters from reference simulation and
    ! compare to current one
    open(12, file="ref_params.txt", status="old", form="formatted", action="read")
        read(12,*) L
        read(12,*) dtau
        read(12,*) num_d
    close(12)
    
    !write(*,*) "Checking Reference data..."
    if (L .ne. T1%L) then
        write(*,*) "Mismatch in Number of time-slices L!"
        stop
    elseif ( (abs(dtau - T1%dtau)) > 0.001 ) then
        write(*,*) "Mismatch in dtau!"
        stop
    elseif (num_d .ne. T1%properties(1)%nclass) then
        write(*,*) "Mismatch in Geometry! (Size?)"
        stop
    end if
    !write(*,*) "Reference Data probably okay!"
    
  end subroutine DFPT_Read_Ref_G
  
  
  subroutine DFPT_Read_Pert_Params(tp,mu,hsc)
    real(wp), intent(inout) :: tp, mu, hsc
    
    open(12, file="pert_params.txt", status="old", form="formatted", action="read")
        read(12,*) tp
        read(12,*) mu
        read(12,*) hsc
    close(12)

  end subroutine DFPT_Read_Pert_Params

  

  subroutine DFPT_Extend_Ref_G_ij(ref_GF,ref_GF_ij,T1)
    ! 
    ! ===========================
    ! 
    ! Extends the default reference G array, which depends on Distance class
    ! D only, to an (i,j)-matrix; Makes it easier to get g_tilde and do the FT's
    ! to (k,k').
    ! 
    ! ===========================
    ! 
    type(TDM1), intent(in)    :: T1
    ! 
    real(wp), intent(in)   , allocatable :: ref_GF(:,:)
    real(wp), intent(inout), allocatable :: ref_GF_ij(:,:,:)
    ! 
    integer :: it, ii, jj, D, N, L
    
    ! Number of sites and time slices
    L = T1%L
    N = T1%properties(1)%n
    
    do it = 0, L-1
        do ii = 1, N
            do jj = 1, N
                D = T1%properties(1)%D(ii,jj)
                ref_GF_ij(ii,jj,it) = ref_GF(D,it)
            end do
        end do
    end do
     
  
  end subroutine DFPT_Extend_Ref_G_ij
  
  
  
  
  ! =========================================================
  ! 
  ! Stuff from here is adapted from Sasha's mod_global.f90 
  ! and from the dfgreens.f 
  ! 
  ! =========================================================
  
  
  subroutine DFPT_build_dual_G0(tp,mu,hsc,N,L,k_vecs,ref_GF_kkp_tau,dual_G0,dual_G0_full,dual_G0sc,dual_G0sc_full,Nom,dtau)
    ! 
    ! Builds the "bare" dual GF from reference G from previous DQMC run
    ! and with the "perturbation matrix", for a given t' and mu with SC-fiekd hsc
    ! 
    
    ! Input
    integer,  intent(in) :: N, L, Nom
    real(wp), intent(in) :: tp, mu, hsc, dtau
    complex(wp), intent(in),    allocatable :: ref_GF_kkp_tau(:,:,:)
    real(wp),    intent(in),    allocatable :: k_vecs(:,:)
    complex(wp), intent(inout), allocatable :: dual_G0(:,:)
    complex(wp), intent(inout), allocatable :: dual_G0_full(:)
    complex(wp), intent(inout), allocatable :: dual_G0sc(:,:)
    complex(wp), intent(inout), allocatable :: dual_G0sc_full(:)
         
    ! Local variables
    integer :: kk, ll, iom, idx
    complex(wp), dimension(L)     :: ctaudata
    complex(wp), dimension(0:Nom) :: comdata
    real(wp)   , dimension(0:L-1) :: ref_GF_av
    complex(wp), dimension(N,0:Nom) :: gkw0, Gkwdf, Gkwscdf 
    real(wp), dimension(N) :: tdk, Dgapk
    complex(wp) :: a,b,c,ai,bi,ci
    
    ! Change Sign of reference GF (due to QMC notation) and transform to Matsubara Freq's
    !write(*,'(A20)') 'G_k_ref(iom)'
    do kk = 1, N
        ! Only Spin-Average of ref. GF (last index 1) (assume spin-symmetry)
        ref_GF_av(:) = ref_GF_kkp_tau(kk,kk,:)
        ! Change notation here so that: ctaudata(1) = - G_ref(tau=0)
        do ll = 0, L-1
            ctaudata(ll+1) = - ref_GF_av(ll)
        end do
        ! Matsubara
        call nfourier(ctaudata,comdata,0,L,Nom,dtau)
        
         !write(*,'(3X,A3,3X,A1,6X,A4)') 'iom', 'k', 'gkw0'
        do iom = 0, Nom
            gkw0(kk,iom) = comdata(iom)
            !write(*,'(2I5,2F15.10)') iom, kk, gkw0(kk,iom)
        end do 
        
    end do
    
    ! Now: Calculate Dual G0 in (k,iw)-space
    !write(*,'(A20)') 'G_k_dual(iom)'
    do kk = 1, N
        ! Perturbation (mu and t'-dispersion) and D-gap function
        tdk(kk) = - 4.0_wp * tp * cos(k_vecs(kk,1)) * cos(k_vecs(kk,2)) - mu
        Dgapk(kk) = 2.0_wp * hsc * (cos(k_vecs(kk,1)) - cos(k_vecs(kk,2)))
        do iom = 0, Nom
          a = tdk(kk)
          b = -a 
          c = Dgapk(kk)
          call triinv(a,b,c,ai,bi,ci)
          a = ai - gkw0(kk,iom)
          b = bi + conjg(gkw0(kk,iom)) 
          c = ci
          call triinv(a,b,c,ai,bi,ci)
          Gkwdf(kk,iom) = ai
          Gkwscdf(kk,iom) = ci 
            !write(*,'(2I5,2F15.10)') iom, kk, Gkwdf(kk,iom)
        end do
            !write(*,'(2I5,2F15.10)') 
    end do
    
    dual_G0(:,:) = Gkwdf(:,:)
    dual_G0sc(:,:) = Gkwscdf(:,:)
    
    ! 
    ! Put the dual G0 into the format G(volume) that is used later
    ! 
    do kk = 1, N
        ! Positive Frequencies
        do ll = 0, L/2 - 1
            idx = kk + ll * N
            dual_G0_full(idx) = dual_G0(kk,ll)
            dual_G0sc_full(idx) = dual_G0sc(kk,ll)
        end do
        
        ! Negative Frequencies
        do ll = L/2, L-1
            idx = kk + ll * N
            dual_G0_full(idx) = conjg(dual_G0(kk, L-1-ll ))
            dual_G0sc_full(idx) = conjg(dual_G0sc(kk, L-1-ll ))
        end do
    end do
  
  end subroutine DFPT_build_dual_G0
  
  subroutine triinv(a,b,c,ai,bi,ci)
    implicit none
    complex(wp) :: a,b,c,ai,bi,ci,det
  !    |a  c|**(-1)=  1/(ab-cc) |b   -c|
  !    |c* b|                   |-c*  a|
       det=a*b-c*conjg(c)
       ai=b/det
       bi=a/det
       ci=-c/det
  end subroutine triinv
  
  
  subroutine DFPT_Extend_Ref_G_ttp(ref_GF_ij,tm,ref_GF_full)
    ! 
    ! Extends G_ref(i,j,tau) to full G_ref(i,j,t,t') and
    ! puts it in the format G_ref(Volume,Volume) which is 
    ! required later for the 6D FFT Routine
    ! 
    type(TDM1), intent(in)    :: tm
    real(wp)  , intent(in)   , allocatable :: ref_GF_ij(:,:,:)
    real(wp)  , intent(inout), allocatable :: ref_GF_full(:,:)
    ! 
    ! Local variables
    ! 
    integer :: N, L, V
    real(wp), allocatable :: gtmp(:)
    integer :: i1, i2, iL, iL1, iL2, iks, jks
    
    N = tm%properties(1)%n
    L = tm%L
    V = N*L
    allocate(gtmp(-L+1:L-1))
    
    do i1 = 1, N
        do i2 = 1, N
            do iL = 0, L-1
                ! Take only spin-averaged Reference GF 
                gtmp(iL) = ref_GF_ij(i1,i2,iL)
            end do
            ! Reflect G (due to symmetry) for negative arguments
            do iL = 1, L-1
                gtmp(-iL) = - gtmp(L-iL)
            end do
            
            ! Put into correct slot of the big array
            do iL1 = 0, L-1
                jks = i1 + iL1 * N
                do iL2 = 0, L-1
                    iks = i2 + iL2 * N
                    ref_GF_full(iks,jks) = - gtmp(iL2 - iL1)
                
                end do
            end do
        
        
        end do
    end do 
  
  
  end subroutine DFPT_Extend_Ref_G_ttp
  
  
  subroutine DFPT_Get_Lattice_G(tp,mu,hsc,N,L,k_vecs,ref_GF_kkp_tau,Sdual,Sdualsc,Nom,dtau,lattG,lattGsc)
    ! 
    ! Obtains the lattice GF from the dual self-energy and perturbation matrix
    ! 
    
    ! Input
    integer,  intent(in) :: N, L, Nom
    real(wp), intent(in) :: tp, mu, hsc, dtau
    complex(wp), intent(in),    allocatable :: ref_GF_kkp_tau(:,:,:)
    real(wp),    intent(in),    allocatable :: k_vecs(:,:)
    complex(wp), intent(in),    allocatable :: Sdual(:)
    complex(wp), intent(inout), allocatable :: lattG(:)
    complex(wp), intent(in),    allocatable :: Sdualsc(:)
    complex(wp), intent(inout), allocatable :: lattGsc(:)
    ! Local variables
    integer :: kk, ll, iom, idx
    complex(wp), dimension(L)     :: ctaudata
    complex(wp), dimension(0:Nom) :: comdata
    real(wp)   , dimension(0:L-1) :: ref_GF_av
    complex(wp), dimension(N,0:Nom) :: gkw0
    complex(wp), dimension(N*L) :: greftmp
    real(wp), dimension(N) :: tdk, Ek
    real(wp) :: Ntot,Ntotk,omega,pi,Beta
    complex(wp) :: Xi 
    complex(wp) :: a,b,c,ai,bi,ci
    
    Xi    = dcmplx(0.0_wp,1.0_wp)
    ! 
    ! First: Transform reference GF to Matsubara Freq's 
    ! 
    
    ! Change Sign of reference GF (due to QMC notation) and transform to Matsubara Freq's
    !write(*,'(A20)') 'G_k_ref(iom)'
    do kk = 1, N
        ! Only Spin-Average of ref. GF (assume spin-symmetry)
        ref_GF_av(:) = ref_GF_kkp_tau(kk,kk,:)
        ! Change notation here so that: ctaudata(1) = -G_ref(tau=0)
        do ll = 0, L-1
              ctaudata(ll+1) = -ref_GF_av(ll)
        end do
        ! Matsubara
        call nfourier(ctaudata,comdata,0,L,Nom,dtau)
        
        !write(*,'(3X,A3,3X,A1,6X,A4)') 'iom', 'k', 'gkw0'
        do iom = 0, Nom
            gkw0(kk,iom) = comdata(iom)
            !write(*,'(2I5,2F15.10)') iom, kk, gkw0(kk,iom)
        end do
    end do
    
    
    ! 
    ! Put Gref(kk,iw) into the same shape as Sdual, i.e. G(N*L), and get
    ! the lattice GF
    ! 
    do kk = 1, N
        ! Perturbation (mu and t'-dispersion) and full hopping tk-mu assuming t=1
        tdk(kk) = - 4.0_wp * tp * cos(k_vecs(kk,1)) * cos(k_vecs(kk,2)) - mu
        Ek(kk) =-2.0_wp*(cos(k_vecs(kk,1))+cos(k_vecs(kk,2)))-4.0_wp*tp*cos(k_vecs(kk,1))*cos(k_vecs(kk,2))-mu
        ! Positive Frequencies 2x2 Nambu-Gorkiov GF (now without Hsc-field!)
        ! lattG(idx) = ( ( greftmp(idx) + Sdual(idx) )**(-1.0_wp) - tdk(kk) )**(-1.0_wp)
        do ll = 0, L/2 - 1
            idx = kk + ll * N
            greftmp(idx) = gkw0(kk,ll)
            ! 2x2 1-st inversion
            a = greftmp(idx) + Sdual(idx)
            b = -conjg(a)
            c = Sdualsc(idx) 
              call triinv(a,b,c,ai,bi,ci)
            ! 2x2 2-nd inversion without Hsc-field
            a = ai - tdk(kk)
            b = bi + tdk(kk)
            c = ci
              call triinv(a,b,c,ai,bi,ci)
            lattG(idx) = ai
            lattGsc(idx) = ci
        end do
        
        ! Negative Frequencies
        do ll = L/2, L-1
            idx = kk + ll * N
            greftmp(idx) = conjg(gkw0(kk, L-1-ll ))
            greftmp(idx) = gkw0(kk,ll)
            ! 2x2 1-st inversion
            a = greftmp(idx) + Sdual(idx)
            b = -conjg(a)
            c = Sdualsc(idx) 
              call triinv(a,b,c,ai,bi,ci)
            ! 2x2 2-nd inversion without Hsc-field
            a = ai - tdk(kk)
            b = bi + tdk(kk)
            c = ci
              call triinv(a,b,c,ai,bi,ci)
            lattG(idx) = ai
            lattGsc(idx) = ci           
        end do
    end do

    !  Symmetrize lattice GF and Sdual according to 8 elements of 2d scuare group 
    call Symmetrise_square_bz_swave(N,L,lattG) 
    call Symmetrise_square_bz_swave(N,L,Sdual) 
    !  Gor'kov SC-part has D-wave symmetry
    call Symmetrise_square_bz_Dwave(N,L,lattGsc) 
    call Symmetrise_square_bz_Dwave(N,L,Sdualsc) 

    ! Estimate Number of Fermions (occupation) via the Matsubara sum
    ! K and  Omega sum Integrate NOS Sum = Fermi Function with exp(Beta(Ek-Mu))
       Beta = Real(L) * dtau
       pi = acos(-1.0_wp)
       Ntot=0.0_wp
       do kk = 1, N
        Ntotk=0.0_wp
         do ll = 0,L/2-1
          idx = kk + ll * N
          omega=(2.0_wp*real(ll)+1.0_wp)*Pi/Beta
          Ntotk=Ntotk+Real(lattG(idx)-1.0_wp/(Xi*omega-Ek(kk)))
         end do
        Ntot=Ntot+Ntotk*2.0_wp/Beta+1.0_wp/(exp(Beta*Ek(kk))+1.0_wp)
       end do
       Ntot=Ntot/Real(N)
       print*,'N-tot per 2 spins Dop=',Ntot*2.0_wp,Ntot*2.0_wp-1.0_wp

  end subroutine DFPT_Get_Lattice_G
  

  subroutine Symmetrise_square_bz_swave(N,L,lattG)
  ! Input
    integer, intent(in) :: N, L
    complex(wp) :: lattG(:)
    
  ! Local
  integer, parameter :: nsym=8
  integer :: Nx, i, j, it, ijt, isym, idx, jdx, ijdx
  real(wp), allocatable :: kx(:), ky(:)
  complex(wp), allocatable :: G_sym(:)
  real(wp) :: two_pi, dk, xnew, ynew
  real(wp), dimension(2,2,nsym) :: symop
  real(wp), dimension(nsym) :: charac

  ! Def. Number of K-points in X(Y) directions
  Nx = int(sqrt(real(N)))
  if (mod(Nx,2) /= 0) stop 'Please use even Nx for X(Y)-directions!'

  allocate(kx(Nx), ky(Nx))
  allocate(G_sym(N*L))

  two_pi = 2.0_wp * acos(-1.0_wp)
  dk = two_pi / real(Nx,wp)

  do i=1,Nx
    kx(i) = real(i-1,wp)*dk
    ky(i) = real(i-1,wp)*dk
  end do

  ! Symmetry operators (SC-2d lattice)
  symop(:,:,1) = reshape([1.0_wp,  0.0_wp,  0.0_wp, 1.0_wp],  [2,2]) ! E
  symop(:,:,2) = reshape([0.0_wp, -1.0_wp, 1.0_wp,  0.0_wp],  [2,2]) ! C4+
  symop(:,:,3) = reshape([-1.0_wp,  0.0_wp,  0.0_wp, -1.0_wp],[2,2]) ! C2
  symop(:,:,4) = reshape([0.0_wp,  1.0_wp, -1.0_wp, 0.0_wp],  [2,2]) ! C4-
  symop(:,:,5) = reshape([-1.0_wp, 0.0_wp, 0.0_wp,  1.0_wp],  [2,2]) ! σx
  symop(:,:,6) = reshape([1.0_wp,  0.0_wp, 0.0_wp, -1.0_wp],  [2,2]) ! σy
  symop(:,:,7) = reshape([0.0_wp,  1.0_wp, 1.0_wp,  0.0_wp],  [2,2]) ! σd
  symop(:,:,8) = reshape([0.0_wp, -1.0_wp,-1.0_wp, 0.0_wp],   [2,2]) ! σd'

  charac = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]  ! s-wave
 
  ! Symmetrisation with direct index calculation
  G_sym = 0.0_wp
  do isym = 1, nsym
    do i = 1, Nx
      do j = 1, Nx
       do it = 0, L-1 
        ijt = i + (j-1) * Nx + it * N
        xnew = symop(1,1,isym)*kx(i) + symop(1,2,isym)*ky(j)
        ynew = symop(2,1,isym)*kx(i) + symop(2,2,isym)*ky(j)
        ! Wrap
        xnew = mod(xnew, two_pi)
        ynew = mod(ynew, two_pi)
        ! Nearest grid index
        idx = 1 + nint(xnew/dk)
        if (idx > Nx) idx = idx - Nx
        if (idx < 1) idx = idx + Nx
        jdx = 1 + nint(ynew/dk)
        if (jdx > Nx) jdx = jdx - Nx
        if (jdx < 1) jdx = jdx + Nx
        ijdx = idx + (jdx-1) * Nx + it * N
        G_sym(ijt) = G_sym(ijt) + charac(isym) * lattG(ijdx)
       end do
      end do
    end do
  end do

  lattG(:) = G_sym(:) / real(nsym,wp)

  deallocate(kx,ky,G_sym)
  end subroutine Symmetrise_square_bz_swave
 
 subroutine Symmetrise_square_bz_Dwave(N,L,lattG)
  ! Input
    integer, intent(in) :: N, L
    complex(wp) :: lattG(:)
    
  ! Local
  integer, parameter :: nsym=8
  integer :: Nx, i, j, it, ijt, isym, idx, jdx, ijdx
  real(wp), allocatable :: kx(:), ky(:)
  complex(wp), allocatable :: G_sym(:)
  real(wp) :: two_pi, dk, xnew, ynew
  real(wp), dimension(2,2,nsym) :: symop
  real(wp), dimension(nsym) :: charac

  ! Def. Number of K-points in X(Y) directions
  Nx = int(sqrt(real(N)))
  if (mod(Nx,2) /= 0) stop 'Please use even Nx for X(Y)-directions!'

  allocate(kx(Nx), ky(Nx))
  allocate(G_sym(N*L))

  two_pi = 2.0_wp * acos(-1.0_wp)
  dk = two_pi / real(Nx,wp)

  do i=1,Nx
    kx(i) = real(i-1,wp)*dk
    ky(i) = real(i-1,wp)*dk
  end do

  ! Symmetry operators (SC-2d lattice)
  symop(:,:,1) = reshape([1.0_wp,  0.0_wp,  0.0_wp, 1.0_wp],  [2,2]) ! E
  symop(:,:,2) = reshape([0.0_wp, -1.0_wp, 1.0_wp,  0.0_wp],  [2,2]) ! C4+
  symop(:,:,3) = reshape([-1.0_wp,  0.0_wp,  0.0_wp, -1.0_wp],[2,2]) ! C2
  symop(:,:,4) = reshape([0.0_wp,  1.0_wp, -1.0_wp, 0.0_wp],  [2,2]) ! C4-
  symop(:,:,5) = reshape([-1.0_wp, 0.0_wp, 0.0_wp,  1.0_wp],  [2,2]) ! σx
  symop(:,:,6) = reshape([1.0_wp,  0.0_wp, 0.0_wp, -1.0_wp],  [2,2]) ! σy
  symop(:,:,7) = reshape([0.0_wp,  1.0_wp, 1.0_wp,  0.0_wp],  [2,2]) ! σd
  symop(:,:,8) = reshape([0.0_wp, -1.0_wp,-1.0_wp, 0.0_wp],   [2,2]) ! σd'

  charac = [1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp,-1.0_wp,-1.0_wp]  ! D-wave (B1g)
 
  ! Symmetrisation with direct index calculation
  G_sym = 0.0_wp
  do isym = 1, nsym
    do i = 1, Nx
      do j = 1, Nx
       do it = 0, L-1 
        ijt = i + (j-1) * Nx + it * N
        xnew = symop(1,1,isym)*kx(i) + symop(1,2,isym)*ky(j)
        ynew = symop(2,1,isym)*kx(i) + symop(2,2,isym)*ky(j)
        ! Wrap
        xnew = mod(xnew, two_pi)
        ynew = mod(ynew, two_pi)
        ! Nearest grid index
        idx = 1 + nint(xnew/dk)
        if (idx > Nx) idx = idx - Nx
        if (idx < 1) idx = idx + Nx
        jdx = 1 + nint(ynew/dk)
        if (jdx > Nx) jdx = jdx - Nx
        if (jdx < 1) jdx = jdx + Nx
        ijdx = idx + (jdx-1) * Nx + it * N
        G_sym(ijt) = G_sym(ijt) + charac(isym) * lattG(ijdx)
       end do
      end do
    end do
  end do

  lattG(:) = G_sym(:) / real(nsym,wp)

  deallocate(kx,ky,G_sym)
  end subroutine Symmetrise_square_bz_Dwave 

  subroutine DFPT_Write_Lattice_G(N,L,dtau,lattG,lattGsc,Sdual,Sdualsc)
    ! 
    ! Writes the lattice Green's function  ((k,iw)-dependent) into a file
    ! 
    
    ! Input
    integer, intent(in) :: N, L
    real(wp), intent(in) :: dtau
    complex(wp), intent(in) :: lattG(:),Sdual(:)
    complex(wp), intent(in) :: lattGsc(:),Sdualsc(:)
    
    ! Local
    integer :: Nx, i, j, it, ijt
    real(wp) :: wn, pi 
    ! Format String
    character(len=10) :: fmtstr
    
    write(fmtstr,"(A1,I0,A6)") "(",9,"E16.8)"
    
    ! Open Files
    open(12, file="G_lat.txt", status="replace", form="formatted", action="write")

    ! Def. Number of K-points in X(Y) directions
    Nx = int(sqrt(real(N)))
    pi = acos(-1.0_wp)
    do i = 1, Nx
      do j = 1, Nx
        write(12,*)"#Kx Ky(2pi/N)", (i-1), (j-1)
       !write(12,*)"#Kx Ky", (i-1)/pi/real(N), (j-1)/pi/real(N)
        do it = 0, L/2-1 
          ijt = i + (j-1) * Nx + it * N
          wn = (2.0_wp*real(it)+1.0_wp)*pi/(real(L)*dtau)
        write(12,fmtstr) wn, lattG(ijt),lattGsc(ijt),Sdual(ijt),Sdualsc(ijt)
        end do
        write(12,*)
      end do  
    end do
    
    ! Close Files
    close(12)
  
  end subroutine DFPT_Write_Lattice_G

  
  subroutine nfourier(ctaudata,comdata,iflag,L,Nom,dtau)
    !  From DMFT-RMP: Fourier-Transform the Natural-Spline Interpolation
    !  of function G(tau) in order to calculate function G(omega)
    !  J. Stoer R. Bulirsch, Introduction to numerical analysis (Springer, New York, 1980)
    
    ! Input
    integer,  intent(in) :: L, Nom, iflag
    real(wp), intent(in) :: dtau
    complex(wp), intent(in)   , dimension(L)     :: ctaudata
    complex(wp), intent(inout), dimension(0:Nom) :: comdata
    
    ! Local Variables
    real(wp), dimension(L+1) :: ctaucopy
    real(wp), dimension(L+1) :: tau
    real(wp) :: cf(4,L+1)
    real(wp), dimension(L) :: a, b, c, d
    complex(wp) :: cdummy, explus, ex, xi
    real(wp) :: om
    real(wp) :: beta
    integer :: ii, jj


    beta = L * dtau
    xi = dcmplx(0.d0, 1.d0)
    ! Change notations to DMFT-RMP: G(1) = G(tau=0)
    ! 
    ! - Note: QUEST GF tau-index goes from 0:L-1, not from 1:L! Start is at tau = 0!
    ! 
    !do ii = 1, L-1
    do ii = 1, L
      !ctaucopy(ii+1) = ctaudata(ii)
      ctaucopy(ii) = ctaudata(ii)
    end do

    ! Extend to tau = beta, which is not included in QUEST GF, by using symmetry
    if (iflag .eq. 0) then
      ctaucopy(L+1) = - ONE - ctaucopy(1)
    elseif (iflag .eq. 1) then
      ctaucopy(L+1) = - ctaudata(1)
    end if

    ! Not necessary since tau = 0 already available
    !if (iflag .eq. 0) then
    !  ctaucopy(1)= - ONE - ctaudata(L)
    !else
    !  this is for Sig_Dual - antisymm????
    !  ctaucopy(1) = - ctaudata(L)
    !end if

    ! Restore G(tau=Beta)
    !if (iflag .eq. 0) then
    !  ctaucopy(L+1)= - ONE - ctaucopy(1)
    !else
    !  ctaucopy(L+1)= - ctaucopy(1)
    !end if

    ! Spline Interpolation: The Spline is given by
    ! G(tau)= a(i) + b(i) (tau-tau_i) + c(i) (...)^2 + d(i) (...)^3
    ! The following Spline subr. is from ELK-code!
    
    ! Tau
    do ii = 1, L+1
      tau(ii)=dtau*(ii-1)
    end do
    
    ! Interpolation
    call spline(L+1,tau,1,ctaucopy,cf)
    
    ! Get Coefficients
    do ii = 1, L
      a(ii) = ctaucopy(ii)
      b(ii) = cf(1,ii)
      c(ii) = cf(2,ii)
      d(ii) = cf(3,ii)
    end do

    ! The Spline multiplied by the Exponential can now be explicitly
    ! integrated. The following formulas were obtained using MATHEMATICA
    comdata = ZERO
    do ii = 0, Nom
        om=(TWO * ii + ONE)*pi/Beta
        do jj = 1, L
            cdummy = xi*om*dtau*jj
            explus = exp(cdummy)
            cdummy = xi*om*dtau*(jj-1)
            ex     = exp(cdummy)
            comdata(ii) = comdata(ii) + explus*( &
                        ( - 6.0_wp * d(jj) ) / om**4 + &
                        ( 2.0_wp * xi * c(jj) + 6.0_wp * dtau * xi * d(jj)  ) / om**3 + &
                        ( b(jj) + 2.0_wp * dtau * c(jj) + 3.0_wp * dtau**2 * d(jj) ) / om**2 + &
                        ( -xi * a(jj) - dtau * xi * b(jj) - dtau**2 * xi * c(jj) - &
                          dtau**3 * xi * d(jj)) / om )
 
            comdata(ii) = comdata(ii) + ex*( &
                          6.0_wp * d(jj) / om**4 - 2.0_wp * xi * c(jj) / om**3 &
                          -b(jj) / om**2 + xi * a(jj) / om )
        end do
    end do
    
  end subroutine nfourier
  
  
  
  subroutine invfourier(comdata,ctaudata,iflag,L,Nom,dtau)
    ! 
    ! Inverse Fourier transform:
    ! G(iw) -> G(tau)
    ! 
    ! G_t(i) = G(i*dtau) for i=0,...,L-1
    ! 
    ! G_w(n) = G(i w_n), for n=0,L/2-1
    !        - w_n = (2*n+1)pi/beta
    ! 
    ! Symmetry: G(iw_(-n)) = G(iw_(n-1))*
    
    ! Input
    integer,  intent(in) :: L, Nom, iflag
    real(wp), intent(in) :: dtau
    complex(wp), dimension(0:Nom), intent(in) :: comdata
    complex(wp), dimension(L), intent(inout)  :: ctaudata
  
    ! Local variables
    real(wp) :: ONEPI = 3.14159265359_wp
    complex(wp) :: com, exom, dum, xi
    real(wp) :: tau, etau, om, omega, muk, Beta
    integer :: i, j
    
    
    xi = dcmplx(0.d0, 1.d0)
    Beta = L * dtau
    
    ctaudata = ZERO
    do i = 1, L
        ! Here (QUEST): tau goes from 0:L-1, i.e. 0...(L-1)*dtau
        tau = (i-1) * dtau
        do j = 0, Nom
            omega = (TWO*j + ONE) * ONEPI / Beta
            om = mod(omega*tau, TWO*ONEPI)
            com = xi * om
            exom = exp(-com)
            dum = comdata(j) * exom
            ctaudata(i) = ctaudata(i) + TWO / Beta * real(dum)
        end do
    end do
    ! Special Treatment for tau = 0?
  
  end subroutine invfourier
  
  
  
  ! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU Lesser General Public
  ! License. See the file COPYING for license details.
  
  !BOP
  ! !ROUTINE: spline
  ! !INTERFACE:
  subroutine spline(n,x,ld,f,cf)
  ! !INPUT/OUTPUT PARAMETERS:
  !   n  : number of points (in,integer)
  !   x  : abscissa array (in,real(n))
  !   ld : leading dimension (in,integer)
  !   f  : input data array (in,real(ld,n))
  !   cf : cubic spline coefficients (1,2,3) and work space (4) (out,real(4,n))
  ! !DESCRIPTION:
  !   Calculates the coefficients of a cubic spline fitted to input data. In other
  !   words, given a set of data points $f_i$ defined at $x_i$, where
  !   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
  !   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
  !   is the interpolating function for $x\in[x_i,x_{i+1})$. This is done by
  !   determining the end-point coefficients $c_2^1$ and $c_2^n$ from the first
  !   and last three points, and then solving the tridiagonal system
  !   $$ d_{i-1}c_2^{i-1}+2(d_{i-1}+d_i)c_2^i+d_ic_2^{i+1}
  !    =3\left(\frac{f_{i+1}-f_i}{d_i}-\frac{f_i-f_{i-1}}{d_{i-1}}\right), $$
  !   where $d_i=x_{i+1}-x_i$, for the intermediate coefficients.
  ! 
  ! !REVISION HISTORY:
  !   Created October 2004 (JKD)
  !   Improved speed and accuracy, April 2006 (JKD)
  !   Optimisations and improved end-point coefficients, February 2008 (JKD)
  !EOP
  !BOC
  !implicit none
  ! arguments
  integer, intent(in) :: n
  real(8), intent(in) :: x(n)
  integer, intent(in) :: ld
  real(8), intent(in) :: f(ld,n)
  real(8), intent(out) :: cf(4,n)
  ! local variables
  integer i
  real(8) t1,t2,t3,t4
  if (n.le.0) then
    write(*,*)
    write(*,'("Error(spline): n <= 0 : ",I8)') n
    write(*,*)
    stop
  end if
  if (n.eq.1) then
    cf(:,1)=0.d0
    return
  end if
  if (n.eq.2) then
    cf(1,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
    cf(2:3,1)=0.d0
    cf(1,2)=cf(1,1)
    cf(2:3,2)=0.d0
    return
  end if
  cf(4,1)=1.d0/(x(2)-x(1))
  cf(1,1)=cf(4,1)*(f(1,2)-f(1,1))
  cf(4,2)=1.d0/(x(3)-x(2))
  cf(1,2)=cf(4,2)*(f(1,3)-f(1,2))
  cf(2,1)=1.d0
  ! estimate second derivative at the first point
  cf(3,1)=(cf(1,2)-cf(1,1))/(x(3)-x(1))
  ! use Gaussian elimination to solve tridiagonal system
  t1=(x(2)-x(1))*cf(4,2)
  t2=t1*cf(2,1)
  t3=1.d0/(2.d0*(t1+1.d0))
  cf(2,2)=t3
  t4=3.d0*(cf(1,2)-cf(1,1))*cf(4,2)-t2*cf(3,1)
  cf(3,2)=t4
  do i=3,n-1
    cf(4,i)=1.d0/(x(i+1)-x(i))
    cf(1,i)=cf(4,i)*(f(1,i+1)-f(1,i))
    t1=(x(i)-x(i-1))*cf(4,i)
    t2=t1*t3
    t3=1.d0/(2.d0*t1+2.d0-t2)
    cf(2,i)=t3
    t4=3.d0*(cf(1,i)-cf(1,i-1))*cf(4,i)-t2*t4
    cf(3,i)=t4
  end do
  ! estimate second derivative at the last point
  t3=(cf(1,n-1)-cf(1,n-2))/(x(n)-x(n-2))
  cf(3,n)=t3
  cf(2,n)=t3
  t3=cf(3,n-1)-t3
  cf(3,n-1)=t3
  t2=cf(2,n-1)
  cf(2,n-1)=t2*t3
  do i=n-2,2,-1
    t3=cf(3,i)-t2*t3
    cf(3,i)=t3
    t2=cf(2,i)
    cf(2,i)=t2*t3
  end do
  cf(2,1)=cf(2,1)*cf(3,1)
  do i=1,n-1
    t1=0.3333333333333333333d0*(cf(2,i+1)-cf(2,i))
    cf(3,i)=t1*cf(4,i)
    cf(1,i)=cf(1,i)-(cf(2,i)+t1)*(x(i+1)-x(i))
  end do
  ! determine end-point coefficients
  t1=x(n)-x(n-1)
  cf(1,n)=cf(1,n-1)+(2.d0*cf(2,n-1)+3.d0*cf(3,n-1)*t1)*t1
  cf(3,n)=cf(3,n-1)
  return
  end subroutine spline
  !EOC
  
end module dfpt_tools
