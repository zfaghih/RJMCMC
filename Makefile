##
## Receiver function inversion makefile 
## Jan Dettmer jand@uvic.ca 2012
##
## With full optimization KISPIOX
#FC      = mpif90 -fc=/opt/pgi/osx86-64/2013/bin/pgfortran -C
#FC      = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -fast -O4 -fastsse -Mstandard -Mfprelaxed -Mipa=fast,inline
FC       = mpif90 -ffree-line-length-0 -C
##############################################################################
##
##  Define binary and object paths:
##
BIN = ./bin
OBJ = ./obj
OBJS = $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o \
       $(OBJ)/pythag.o $(OBJ)/prjmh_temper_csem.o
FO = $(OBJ)/libmartin.a
##############################################################################

all: bin/prjmh_temper_csem

clean: clean.done
	rm -f *.done
	(cd ./martin/;/bin/rm -f *.o *.mod *.l core *.a;)
clean.done:
	rm -f bin/*
	rm -f *.o *.mod $(OBJ)/*

##############################################################################
#
# Targets for linking of programs:
#
# PRJMH temper Receiver:
#
$(BIN)/prjmh_temper_csem: libmartin $(OBJS) $(FO)
	$(FC)  -o $(BIN)/prjmh_temper_csem $(OBJS) $(FO)

$(OBJ)/prjmh_temper_csem.o: prjmh_temper_csem.f90 
	$(FC) -c -o $(OBJ)/prjmh_temper_csem.o prjmh_temper_csem.f90

##############################################################################
#
# Targets for separate subroutines:
#
libmartin:
	(cd ./martin/;make;)
	/bin/ln -sf ./martin/*.mod .

##############################################################################
#
# Targets to compile NR subroutines that are in separate files:
#
$(OBJ)/ludcmp.o: ludcmp.f90
	$(FC) -c -o $(OBJ)/ludcmp.o ludcmp.f90

$(OBJ)/svdcmp.o: svdcmp.f90
	$(FC) -c -o $(OBJ)/svdcmp.o svdcmp.f90

$(OBJ)/pythag.o: pythag.f90
	$(FC) -c -o $(OBJ)/pythag.o pythag.f90

##############################################################################
#
#  Targets for numerical recipes modules:
#
$(OBJ)/nrtype.o: nrtype.f90 
	$(FC) -c -o $(OBJ)/nrtype.o nrtype.f90

$(OBJ)/nr.o: nr.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nr.o nr.f90

$(OBJ)/nrutil.o: nrutil.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nrutil.o nrutil.f90

$(OBJ)/mod_sac_io.o : mod_sac_io.f90
	$(FC) -c -o $(OBJ)/mod_sac_io.o mod_sac_io.f90

##############################################################################
# ...this is the end my friend
#
# EOF
