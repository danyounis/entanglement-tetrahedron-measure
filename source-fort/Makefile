FC=gfortran#<gcc/11.2.0/b1, lapack/3.9.0/b2, openblas/0.3.10/b1>

main: tetra_main tidy
mini: tetra_mini tidy
adhoc1: tetra_adhoc1 tidy

tetra_main: prec.o math.o optimize.o tetra_main.o
	$(FC) $^ -o tetra_main -fopenmp -llapack -lblas

tetra_mini: prec.o math.o optimize.o tetra_mini.o
	$(FC) $^ -o tetra_mini -fopenmp -llapack -lblas

tetra_adhoc1: prec.o math.o optimize.o tetra_adhoc1.o
	$(FC) $^ -o tetra_adhoc1 -fopenmp -llapack -lblas

prec.o: prec.F08
	$(FC) -c -DDOUBLE $<

math.o: math.f08
	$(FC) -c $< -fopenmp

optimize.o: optimize.f08
	$(FC) -c $<

tetra_main.o: tetra_main.f08
	$(FC) -c $< -fopenmp

tetra_mini.o: tetra_mini.f08
	$(FC) -c $< -fopenmp

tetra_adhoc1.o: tetra_adhoc1.f08
	$(FC) -c $< -fopenmp

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod tetra_main tetra_mini tetra_adhoc1
	@echo Done.
