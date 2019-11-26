F90=ifort
FFLAGS=-O2 -CB -CU #-fopenmp #-fbounds-check -finit-local-zero  -fopenmp

OBJ=subroutines.o
SOURCE=subroutines.f90 

main: main.f90 $(OBJ)
	$(F90) $(FFLAGS) main.f90 $(OBJ) -o main
$(OBJ): $(SOURCE)
	$(F90) $(FFLAGS) -c $(SOURCE)
