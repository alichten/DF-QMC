##########################################################
#
# Makefile for the Dual Fermion Perturbation Theory Code
#
##########################################################

DFPT_DIR = $(shell pwd)

# Fortran-Compilers, sequential and MPI
#FC_SEQ = gfortran
FC_SEQ = ifort 
FC_MPI = mpifort 

ifndef FC_SEQ
  $(error "FC_SEQ" is not defined in the Makefile.)
endif

ifndef FC_MPI
  $(error "FC_MPI" is not defined in the Makefile.)
endif

# Compiler Flags
ifeq ($(FC_SEQ), gfortran)
  FC_FLAGS  = -m64 -Wall -O3 -funroll-loops -fallow-argument-mismatch
endif

ifeq ($(FC_SEQ), ifort)
  FC_FLAGS  = -m64 -warn all -O3 -unroll
endif

# Top directory of QUEST and the library
QUEST_DIR = $(shell pwd)/../QUEST/
DQMCLIB   = $(QUEST_DIR)/libdqmc.a

# BLAS/LAPACK (OpenBLAS or MKL) 
libOpenBLAS = /usr/lib64/libopenblas.a 
#LAPACKLIB = $(libOpenBLAS)
MKLPATH   = $(MKLROOT)/lib/intel64
LAPACKLIB = -L$(MKLPATH) -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

# Standard C++ library
CXXLIB = -lstdc++

# FFTW
FFTW_DIR = /home/alichten/fftw
FFTW_INC = $(FFTW_DIR)/include/
FFTW_LIB = $(FFTW_DIR)/lib/libfftw3.a
FFTW_LIB = -l fftw3 

ifndef QUEST_DIR
  $(error QUEST_DIR is not defined in the Makefile.)
endif

ifndef DQMCLIB
  $(error DQMCLIB is not defined in the Makefile.)
endif

ifndef LAPACKLIB
  $(error LAPACKLIB is not defined in the Makefile.)
endif

ifndef CXXLIB
  $(error CXXLIB is not defined in the Makefile.)
endif

ifndef FFTW_INC
  $(error FFTW_INC is not defined in the Makefile.)
endif

ifndef FFTW_LIB
  $(error FFTW_LIB is not defined in the Makefile.)
endif

# Summary of necessary libraries
LIB = $(CXXLIB) $(LAPACKLIB) $(FFTW_LIB)


# ------------------------------------------
export

all: 
	(cd SRC; $(MAKE))

.PHONY: clean

clean:
	(cd SRC; $(MAKE) clean)
	(cd bin; rm -f build_ref dfpt_mpi)
