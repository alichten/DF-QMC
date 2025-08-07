FLAGS=$(FC_FLAGS) -I$(QUEST_DIR)/SRC -I$(FFTW_INC)

.SUFFIXES : .mod .o .F90

all: $(DQMCLIB) dfpt_tools.F90 build_ref.F90
	$(FC_SEQ) $(FLAGS) -o ../bin/build_ref dfpt_tools.F90 build_ref.F90 $(DQMCLIB) $(LIB)
	$(FC_MPI) $(FLAGS) -o ../bin/dfpt_mpi dfpt_tools.F90 dfpt_mpi.F90 $(DQMCLIB) $(LIB)
	$(FC_MPI) $(FLAGS) -o ../bin/dfsc_mpi dfsc_tools.F90 dfsc_mpi.F90 $(DQMCLIB) $(LIB)



.PHONY: clean

clean:
	rm -f *.o *.mod
