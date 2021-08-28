FFC=gfortran
FFLAG=-O3 -Wall -fbounds-check -fbacktrace

SRC=main.f90
LIB=/home/yamaguchi/numerical_library/src/*.f90

main:
	${FFC} ${FFLAG} ${SRC} ${LIB}

# clean up
clean:
	rm -f *.a *.mod *~
