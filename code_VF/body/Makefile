FC = gfortran
FFLAGS = -g -Wall -fcheck=all -ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow -fbounds-check
SRC = mod_constantes.f90 mod_initialize.f90 mod_flux.f90 mod_TimeSchemes.f90 mod_ApproxSol.f90 mod_ExactSol.f90 mod_error.f90 main.f90
OBJ = ${SRC:.f90=.o}

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<


run: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean :
	@rm -f *.mod *.o run

clear :
	@rm -f *.dat