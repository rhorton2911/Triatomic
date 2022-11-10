EXEC=H3Plus
MCCCOBJ=numbers.o input_parser.o MPI_module.o target_data.o input_data.o \
	Legendre.o grid_radial.o form.o sturmian_class.o

        #Compiling Legendre.f90 into Legendre.o will automatically produce the .mod files for all modules within.
        #Matrix_Print.o Data_Module.o Associated_Legendre_Functions.o Special_Functions.o \  #Contained in Legendre.f90

MYOBJ=wigner.o plql.o intp.o rsg.o basismodule.o main.o 
OBJ= $(MCCCOBJ) $(MYOBJ)
FC=gfortran
FLAGS= -Wall -Wno-tabs -pedantic  -fimplicit-none -fopenmp #-Werror
ERRFLAGS1= -g -Wextra -fcheck=all -fbacktrace -Waliasing -Winteger-division -fbounds-check
ERRFLAGS= $(ERRFLAGS1) -Wsurprising --fpe-trap=invalid,zero,overflow,underflow -fstack-check 

#Comment out when not debugging
CFLAGS=$(FLAGS) $(ERRFLAGS)

all: $(EXEC)

#Linking step: link all compiled files into an executable. 
$(EXEC): $(OBJ)
	$(FC) $(CFLAGS) $(OBJ) -o $(EXEC)

#Compiling step: call to objects in OBJ in the linking step causes these
#the be called to generate the object files.
main.o: wigner.f plql.f numbers.f90 gridr.f90 sturmian.f90 input_data.f90 basismodule.f90 main.f90
	$(FC) $(CFLAGS)	-c main.f90

basismodule.o: wigner.f plql.f numbers.f90 basismodule.f90 
	$(FC) $(CFLAGS) -c basismodule.f90

numbers.o: numbers.f90
	$(FC) $(CFLAGS) -c numbers.f90

#################################MCCC Modules##########################################
sturmian_class.o: sturmian.f90 gridr.f90 form.f90
	$(FC) -w -ffree-line-length-512 -c sturmian.f90 -o sturmian_class.o

form.o: form.f90
	$(FC) -c form.f90

grid_radial.o: gridr.f90 Legendre.f90 
	$(FC) $(CFLAGS) -w -c gridr.f90 -o grid_radial.o

Legendre.o: Legendre.f90
	$(FC) $(CFLAGS) -w -c Legendre.f90

target_data.o: target_data.f90 numbers.f90
	$(FC) $(CFLAGS) -c target_data.f90

input_data.o: input_data.f90 target_data.f90 input_parser.f90 MPI_module.f90
	$(FC) -ffree-line-length-512 -c input_data.f90

MPI_module.o: MPI_module.f90
	$(FC) $(CFLAGS) -c MPI_module.f90

input_parser.o: input_parser.f90
	$(FC)  -c input_parser.f90

######################################MCCC Modules#####################################


wigner.o: wigner.f
	$(FC) $(CFLAGS) -w --std=legacy -c wigner.f

plql.o: plql.f
	$(FC) $(ERRFLAGS) -w -g --std=legacy -c plql.f

rsg.o: rsg.f
	$(FC) $(CFLAGS) -w -g  --std=legacy -c rsg.f

intp.o: intp.f
	$(FC) $(CFLAGS) -w --std=legacy -c intp.f

clean: 
	rm -rf *.o *.mod $(EXEC)

