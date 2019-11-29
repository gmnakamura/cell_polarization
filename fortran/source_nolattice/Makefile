F90=gfortran-mp-9  #gfortran #ifort
FFLAGS=-O2  -fopenmp #  -fbounds-check -Wuninitialized  -fopenmp

OBJ=subroutines.o
SOURCE=subroutines.f90 

main: main.f90 $(OBJ)
	$(F90) $(FFLAGS) main.f90 $(OBJ) -o main
$(OBJ): $(SOURCE)
	$(F90) $(FFLAGS) -c $(SOURCE)
