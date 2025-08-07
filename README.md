*********************

DF-QMC:  overcoming the "sign-problem" in a lattice Quantum Monte-Carlo 

*********************

The Dual-Fermion QMC Perturbation Theory
for 2d square Hubbard lattice

Developed by:

 - Alexander Lichtenstein (alichten@physnet.uni-hamburg.de)

 - Edin Kapetanović  (ekapetan@physnet.uni-hamburg.de)

References: S. Iskakov, M. I. Katsnelson, and A. I. Lichtenstein, npj Comput: Mater. 10, 36 (2024) https://doi.org/10.1038/s41524-024-01221-w

E. Stepanov, S. Iskakov, M. Katsnelson, and A. Lichtenstein: https://doi.org/10.48550/arXiv.2502.08635

Based on the Auxiliary Fields lattice DQMC approach by Jorge E. Hirsch, PRB 31, 4403 (1985)
 
Fortran version: QUEST code 1.4.9 developed by

 - Richard Scalettar (UC Davis) 

Compilation:

 - cd QUEST  "make"  to produce  DQMC libdqmc.a 

 - cd DFPT  "make" to produce normal and SC-field DF for 1-st order perturbation - executables are in the "bin"   

 - Compilers (both sequential and MPI) have to be specified at the top of corresponding Makefileis

 - Additional Libraries:  LAPACK and FFTW  need to be specified

The "test"-folder contain example input files

The "Hub4x4" and "Hub4x4SC" folders contain outputs for normal and SC-state 

Calculations:
 - First run.sh0 "build_ref" (sequential) DQMC for the Reference Green's function (use larger number of QMC-sweeps)
 - Then  run.sh1 DF-QMC in normal state (dfpt_mpi) or SC-state (dfsc_mpi) with maximum possible MPI processors and QMC-sweeps
 - input parameters for DFPT in "pert_params.txt" with Mu, t' and Hsd (for SC-calculatiopns)

Outputs:
 - G_lat.txt contains G(k,omega) and Sdual(k,omega) 
 - in o-out file there is estimation of the occupation number for given Mu and t'   

=================================================

Folders:

  - DFPT/SRC :  Dual Fermion Perturbation Theory developed by Alexander Lichtenstein (UHH)
	
  -	QUEST/SRC : QUEST 1.4.9 source code developed by Richard Scalettar (UC Davis) - included for complitness 
 
  -	DFPT/Hub4x4 : Example of DFPT for doped Hubbard 4x4 lattice 
 
  -	DFSC/Hub4x4 : Example of DFSC for doped  Hubbard 4x4 lattice in d-wave SC field
 
  -	DFPT/DOC_DF : DF-QMC documentation.

=================================================

                    
