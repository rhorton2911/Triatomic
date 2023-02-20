EXEC=H3Plus
MCCCOBJ=numbers.o input_parser.o MPI_module.o target_data.o input_data.o \
  Legendre.o grid_radial.o form.o sturmian_class.o ovlpste1me_module.o spheroidal_fn.o\
  state_class.o
STATECLASSOBJ=vnc_module.o one_electron_func_module.o state_class_module.o natorb.o  vmat_exch_module.o rearrange12.o Horbitals.o one_electron_func.o one_electron_func_group.o one_electron_func_spheroidal.o H12.o  vnc_group.o  rearrange.o  osc.o osc12.o
SCATOBJ=vnc_module.o one_electron_func_module.o state_class_module.o osc.o osc12.o vnc.o Horbitals.o rearrange12.o\
  natorb.o one_electron_func.o one_electron_func_spheroidal.o \
  vmat_exch_module.o H12.o pathpotmod.o Ps_structure.o ps_extra.o vmatrixmod.o channels.o pol_2el.o distpot.o vmat.o\
	kgrid.o ps_coefsSPD.o posvmat.o channels_info.o \
  iqpackd.o kgrid_igor.o cc.o dcoul.o spheroidal_waves.o numerov.o contwaves.o\
  onshellt_mod.o optical.o BornICS_module.o reconstruct.o vmat12.o tmatsolve.o phenom_pol.o tmatsolve_cmplx.o redistribute.o scat.o Jloop.o rearrange.o 
  #From line 5 on, dependencies of Jloop.o  (and their dependencies)


#Compiling Legendre.f90 into Legendre.o will automatically produce the .mod files for all modules within.  #Matrix_Print.o Data_Module.o Associated_Legendre_Functions.o Special_Functions.o \  #Contained in Legendre.f90
LEGACYHELPER=wigner.o plql.o intp.o rsg.o 
OBJ= $(LEGACYHELPER) $(MCCCOBJ) basismodule.o $(STATECLASSOBJ) main.o

#---------------------Choose compiler name depnding on machine ---------------------#
HOSTNAME=$(shell hostname)
$(info HOSTNAME is: $(HOSTNAME))
ifneq (,$(findstring gadi-, $(HOSTNAME)))
        #If current hostname is a gadi node, use gadi-specific parameters
	FC=gfortran
	LAPACKLIB= -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp -lm
else
	#Default, use setup for setonix
	FC=ftn 
	LAPACKLIB= -llapack
endif


#-------------------- Choose compiler flags ---------------------------#
FLAGS= -Wall -Wno-tabs -pedantic -fimplicit-none  -fopenmp -g #-Werror
ERRFLAGS1= -Wextra -fcheck=all -fbacktrace -Waliasing -Winteger-division -fbounds-check
ERRFLAGS= $(ERRFLAGS1) -Wsurprising  -fstack-check #--fpe-trap=invalid,zero,overflow
#Comment out errflags when not debugging
CFLAGS=$(FLAGS)  #$(ERRFLAGS) 
MCFLAGS = $(CFLAGS) -ffree-line-length-512

all: $(EXEC)

#Linking step: link all compiled files into an executable. 
$(EXEC): $(OBJ)
	$(FC) $(CFLAGS) $(OBJ) $(LAPACKLIB) -o $(EXEC)

#Compiling step: call to objects in OBJ in the linking step causes these
#the be called to generate the object files.
main.o: wigner.f plql.f numbers.f90 gridr.f90 sturmian.f90 input_data.f90 basismodule.f90 state_class.f90 main.f90
	$(FC) $(CFLAGS)	-c main.f90

basismodule.o: wigner.f plql.f numbers.f90 basismodule.f90 
	$(FC) $(CFLAGS) -c basismodule.f90

numbers.o: numbers.f90
	$(FC) $(CFLAGS) -c numbers.f90

#################################MCCC Modules##########################################
#MCCC code has a lot of dependencies. .o files in MCCOBJ are listed in order of depenencies
#Later .o files depend on modules defined in earlier ones
sturmian_class.o: sturmian.f90 gridr.f90 form.f90
	$(FC) -w -ffree-line-length-512 -c sturmian.f90 -o sturmian_class.o

#Contains standalone functions
form.o: form.f90
	$(FC) -c form.f90

#Contains a global variable and data type declaration
grid_radial.o: gridr.f90 Legendre.f90 
	$(FC) $(MCFLAGS) -w -c gridr.f90 -o grid_radial.o

Legendre.o: Legendre.f90
	$(FC) $(MCFLAGS) -w -c Legendre.f90

#Contains a global variable and some helper routines
target_data.o: target_data.f90 numbers.f90
	$(FC) $(MCFLAGS) -c target_data.f90

#Contains global variables and type definitions
input_data.o: input_data.f90 target_data.f90 input_parser.f90 MPI_module.f90
	$(FC) -ffree-line-length-512 -c input_data.f90

#Contains global variables
MPI_module.o: MPI_module.f90
	$(FC) $(MCFLAGS) -c MPI_module.f90

input_parser.o: input_parser.f90
	$(FC)  -c input_parser.f90

state_class.o: ovlpste1me_module.f90 spheroidal_fn.f90 state_class.f90
	$(FC) $(MCFLAGS) -c state_class.f90

#FROM HERE DOWN ARE ALL THE DEPENDENCIES FOR STATE_CLASS
#Contains global variables
ovlpste1me_module.o: ovlpste1me_module.f90 
	$(FC) $(MCFLAGS) -c ovlpste1me_module.f90 

#Contains standalone subroutines
rearrange.o: state_class_module.f90 rearrange.f90
	$(FC) $(MCFLAGS) -c rearrange.f90

rearrange12.o: rearrange12.f90
	$(FC) $(MCFLAGS) -c rearrange12.f90

#Contains global variables
state_class_module.o: state_class_module.f90
	$(FC) $(MCFLAGS) -c state_class_module.f90

#Contains global variables
one_electron_func_module.o: one_electron_func_module.f90
	$(FC) $(MCFLAGS) -c one_electron_func_module.f90

#Contains standalone subroutines
one_electron_func.o: one_electron_func.f90
	$(FC) $(MCFLAGS) -c one_electron_func.f90

#Contains standalone subroutines
one_electron_func_group.o: one_electron_func_group.f90
	$(FC) $(MCFLAGS) -c one_electron_func_group.f90

#Contains standalone subroutines
one_electron_func_spheroidal.o: one_electron_func_spheroidal.f90
	$(FC) $(MCFLAGS) -c one_electron_func_spheroidal.f90

spheroidal_fn.o: spheroidal_fn.f90
	$(FC) $(MCFLAGS) -c spheroidal_fn.f90

#Contains standalone subroutines
vnc.o: vnc.f90
	$(FC) $(MCFLAGS) -c vnc.f90

#Contains standalone subroutines
vnc_group.o: vnc_group.f90
	$(FC) $(MCFLAGS) -c vnc_group.f90

#File contains standalone functions
Horbitals.o: Horbitals.f90
	$(FC) $(MCFLAGS) -c Horbitals.f90

#File contains standalone functions
H12.o: H12.f90 basismodule.f90
	$(FC) $(MCFLAGS) -c H12.f90

#File contains standalone functions
osc.o: osc.f90
	$(FC) $(MCFLAGS) -c osc.f90

#File contains standalone functions
osc12.o: osc12.f90
	$(FC) $(MCFLAGS) -c osc12.f90

#File contains standalone functions
Jloop.o: Jloop.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c Jloop.f90

#Contains global variables
vnc_module.o: vnc_module.f90
	$(FC) $(MCFLAGS) -c vnc_module.f90

#File contains standalone functions
natorb.o: natorb.f90
	$(FC) $(MCFLAGS) -c natorb.f90

#Contains global variables
vmat_exch_module.o: vmat_exch_module.f90
	$(FC) $(MCFLAGS) -c vmat_exch_module.f90

#Contains global variables
channels.o: channels.f90
	$(FC) $(MCFLAGS) -c channels.f90

posvmat.o: posvmat.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c posvmat.f90

#Contains global variables
Ps_structure.o: Ps_structure.f90
	$(FC) -c Ps_structure.f90

#Contains global variables
vmatrixmod.o: vmatrixmod.f90
	$(FC) $(MCFLAGS) -c vmatrixmod.f90

#Contains global variables
channels_info.o: channels_info.f90
	$(FC) $(MCFLAGS) -c channels_info.f90

iqpackd.o: iqpackd.f
	$(FC) --std=legacy -c iqpackd.f

kgrid_igor.o: kgrid_igor.f
	$(FC) --std=legacy -c kgrid_igor.f

kgrid.o: kgrid.f90
	$(FC) $(MCFLAGS) -c kgrid.f90

distpot.o: distpot.f90
	$(FC) $(MCFLAGS) -c distpot.f90

#Contains MPI parallelisation, some subroutine argument are not typecast from real to complex
contwaves.o: contwaves.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c contwaves.f90

spheroidal_waves.o: spheroidal_waves.f90
	$(FC) $(MCFLAGS) -c spheroidal_waves.f90

onshellt_mod.o: onshellt_mod.f90
	$(FC) $(MCFLAGS) -c onshellt_mod.f90

optical.o: optical.f90
	$(FC) $(MCFLAGS) -c optical.f90

BornICS_module.o: BornICS_module.f90
	$(FC) $(MCFLAGS) -c BornICS_module.f90

ps_extra.o: ps_extra.f
	$(FC) --std=legacy -c ps_extra.f

pathpotmod.o: pathpotmod.f90
	$(FC)  -c pathpotmod.f90

vmat.o: vmat.f90
	$(FC) $(MCFLAGS) -c vmat.f90

pol_2el.o: pol_2el.f90
	$(FC) $(MCFLAGS) -c pol_2el.f90

ps_coefsSPD.o: ps_coefsSPD.f
	$(FC) --std=legacy -c ps_coefsSPD.f

cc.o: cc.f
	$(FC) --std=legacy -c cc.f

dcoul.o: dcoul.f
	$(FC) --std=legacy -c dcoul.f

numerov.o: numerov.f 
	$(FC) --std=legacy -c numerov.f

scat.o: scat.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c scat.f90

reconstruct.o: reconstruct.f90
	$(FC) -ffree-line-length-512  -c reconstruct.f90

vmat12.o: vmat12.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c vmat12.f90

phenom_pol.o: phenom_pol.f90
	$(FC) $(MCFLAGS) -c phenom_pol.f90

tmatsolve.o: tmatsolve.f90
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c tmatsolve.f90

tmatsolve_cmplx.o: tmatsolve_cmplx.f90
	$(FC) $(MCFLAGS) -c tmatsolve_cmplx.f90

redistribute.o: redistribute.f
	$(FC) -ffree-line-length-512 -fallow-argument-mismatch -c redistribute.f

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

