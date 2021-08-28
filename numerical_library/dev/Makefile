FFC=gfortran
FFLAG=-O3 -Wall -fbounds-check -fbacktrace

SRC=*.f90
OBJ=a.out
LIB=/home/yamaguchi/numerical_library/src/*.f90

main:
	${FFC} ${FFLAG} -o ${OBJ} ${SRC} ${LIB}

# clean up
clean:
	rm -f *.a *.mod *~
