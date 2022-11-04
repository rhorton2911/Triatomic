EXEC=H3Plus
OBJ=wigner.o plql.o intp.o rsg.o numbers.o basismodule.o main.o 
FC=gfortran
FLAGS= -Wall -Wno-tabs -pedantic -Werror -fimplicit-none -fopenmp
ERRFLAGS1= -g3 -Wextra -fcheck=all -fbacktrace -Waliasing -Winteger-division 
ERRFLAGS= $(ERRFLAGS1) -Wsurprising --fpe-trap=invalid,zero,overflow,underflow -fstack-check 

#Comment out when not debugging
CFLAGS=$(FLAGS) $(ERRFLAGS)

all: $(EXEC)

#Linking step: link all compiled files into an executable. 
$(EXEC): $(OBJ)
	$(FC) $(CFLAGS) $(OBJ) -o $(EXEC)

#Compiling step: call to objects in OBJ in the linking step causes these
#the be called to generate the object files.
main.o: wigner.f plql.f numbers.f90 basismodule.f90 main.f90
	$(FC) $(CFLAGS)	-c main.f90

basismodule.o: wigner.f plql.f numbers.f90 basismodule.f90 
	$(FC) $(CFLAGS) -c basismodule.f90

numbers.o: numbers.f90
	$(FC) $(CFLAGS) -c numbers.f90

wigner.o: wigner.f
	$(FC) -w --std=legacy -c wigner.f

plql.o: plql.f
	$(FC) -w -g --std=legacy -c plql.f

rsg.o: rsg.f
	$(FC) -w --std=legacy -c rsg.f

intp.o: intp.f
	$(FC) -w --std=legacy -c intp.f

clean: 
	rm -rf *.o *.mod $(EXEC)

