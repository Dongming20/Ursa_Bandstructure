F90?=ifx
MPI?=impi

NAME_compute = compute_bandstructure

help:
	@echo "----------------------------------------------------------------------------"
	@echo "Usage: make F90=<f90> {configure,compute,all}"
	@echo ""
	@echo "Example: 'make all' will create all the Quasar executable files in ../bin"
	@echo ""
	@echo "Other possible options:"
	@echo "--- <f90> is your own Fortran90 compiler (possible choices: ifx, gfortran). Default is ifx"
	@echo ""
	@echo "REMARK: -Quasar must be linked with the FEAST eigenvalue library: www.feast-solver.org"
	@echo "          -FEAST must be compiled with the same <f90> and <mpi> options."
	@echo "        -Quasar also needs the Tetgen library (GNU-GPL License): wias-berlin.de/software/tetgen"
	@echo "          -install the executable tetgen in your PATH or copy it in the /bin Quasar directory."
	@echo "----------------------------------------------------------------------------"


compute: $(NAME_compute)

all: compute clean

############################################################
######### Fortran Compiler set-up
############################################################


# intel fortran
ifeq ($(F90),ifx)
F90FLAGS= -O3 -qopenmp 
#F90FLAGS += -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces  #-heap-arrays 1
endif

# gnu fortran
ifeq ($(F90),gfortran)
F90FLAGS= -O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none 
endif

# portland group fortran compiler
ifeq ($(F90),pgf90)
F90FLAGS= -O3 -mp
endif


##############################################################
######### include MPI directives
##############################################################


### intel mpi
ifeq ($(MPI),impi)
PF90=mpiifort -f90=$(F90)
PF90FLAGS= $(F90FLAGS) 
endif


### mpich
ifeq ($(MPI),mpich)
PF90= mpif90.mpich -fc=$(F90)
PF90FLAGS= $(F90FLAGS) 
endif

### openmpi ...requires shell environment variable "OMPI_FC=$(F90)"
ifeq ($(MPI),openmpi)
export OMPI_FC=$(F90) # for BASH shell
PF90= mpif90.openmpi
PF90FLAGS= $(F90FLAGS) 
endif

#===========================================================
# PATH To FEAST library
#===========================================================
ARCH?=x64
FEASTROOT?=$(PWD)/../..
LOCLIBS=-L$(FEASTROOT)/lib/$(ARCH)
LIB=$(LOCLIBS) -lfeast


#===========================================================
# MKL Link set-up
#===========================================================


# intel fortran
ifeq ($(F90),ifx)
LIB += -qmkl=parallel -lmkl_blacs_intelmpi_lp64  -liomp5 -lpthread -lm -ldl
endif

# gnu fortran
ifeq ($(F90),gfortran)
LIB += -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl
endif

# portland group fortran compiler (=========>not tested<==============)
ifeq ($(F90),pgf90)
FLIBS +=-Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -pgf90libs -mp -lpthread -lm -ldl
endif


###################################
#### SOURCE CODE FILES
###################################

CODE90_compute = class_linkedlist.o fem.o potential.o compute_bandstructure.o

###################################
#### OBJECTS CODE FILES
###################################

OBJS_compute =  $(CODE90_compute:.f90=p.o)



# .SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: clean all

#.f90.o:
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%p.o: %.f90
	$(PF90) $(PF90FLAGS) -c $< -o $@


#============================================================
# COMPILE and LINK
#============================================================

$(NAME_compute): $(OBJS_compute)
	$(PF90) -o $(NAME_compute) $(OBJS_compute) $(OBJS_linkedlist) $(OBJS_fem) $(OBJS_potential) $(PF90FLAGS) $(LIB)
#	mv $(NAME_compute) ../bin


#==========================================================
# Clean up directory: delete object files and modules
#==========================================================
clean:
	-@rm  *.o *.mod

