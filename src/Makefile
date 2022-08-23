# This Makefile make with mkmkmf
# Define suffixes to look for
.SUFFIXES: .f .f90 .o .mod

SRCDIRS  =ATOM BAND COMMON COMPTON CRPA DEBYE DMDW EELS EELSMDFF ERRORMODS EXCH FF2X FMS FOVRG FULLSPECTRUM GENFMT IOMODS KSPACE LDOS MATH MKGTR MODS PAR PATH POT RDINP RHORRP RIXS SCREEN SELF SFCONV TDLDA XSPH INPGEN 
INCLUDEFLAGS=$(SRCDIRS:%=-I../%)

EXECUTABLES =atomic band compton crpa dmdw dym2feffinp eels ff2x fms fullspectrum genfmt ldos mkgtr opconsat path pot rdinp rhorrp rixs screen screenw2 sfconv xsph 
EXECDIR=../bin/Seq

$(shell [ -e Compiler.mk ] || cp Compiler.mk.default Compiler.mk)

include Compiler.mk

ifneq ($(F90),gfortran)
	FPPTASK=-DFEFF -fpp 
else
	FPPTASK=-DFEFF -cpp
endif
parFPPTASK = $(shell $(MPIF90) -show |grep -c gfortran |awk '{if($$1==0) print "-DFEFF -fpp"; else print "-DFEFF -cpp"}')

lbtool = ar ru
parmessg = "\n\n\n\tmodular src is in ../mod/MPI/\n\tPlease cd there and compile\n\twith the appropriate script."
messg = "\n\n\n\tmodular src is in ../mod/Seq/\n\tPlease cd there and compile\n\twith the appropriate script."
monmessg = "\n\n\n\tsrc is in ../mod/MONO/\n\tPlease cd there and compile\n\twith the appropriate script."
smessg1 = "\n\tBinaries are located in ../bin/"
smessg2 = "\n\tYou can use the script named\n\t\tfeff"
smessg3 = "\n\tlocated in ../bin/ to run feff"

# Source files needed by ATOM/atomic.f90.
include DEP/atomic$(DEPTYPE).mk
# Object files needed by ATOM/atomic.f90.
atomic_OBJ = ${atomicSRC:.f90=.o}
atomic_MODULES =${atomic_MODULESRC:.f90=.o}

# Source files needed by BAND/band.f90.
include DEP/band$(DEPTYPE).mk
# Object files needed by BAND/band.f90.
band_OBJ = ${bandSRC:.f90=.o}
band_MODULES =${band_MODULESRC:.f90=.o}

# Source files needed by COMPTON/compton.f90.
include DEP/compton$(DEPTYPE).mk
# Object files needed by COMPTON/compton.f90.
compton_OBJ = ${comptonSRC:.f90=.o}
compton_MODULES =${compton_MODULESRC:.f90=.o}

# Source files needed by CRPA/crpa.f90.
include DEP/crpa$(DEPTYPE).mk
# Object files needed by CRPA/crpa.f90.
crpa_OBJ = ${crpaSRC:.f90=.o}
crpa_MODULES =${crpa_MODULESRC:.f90=.o}

# Source files needed by DMDW/dmdw.f90.
include DEP/dmdw$(DEPTYPE).mk
# Object files needed by DMDW/dmdw.f90.
dmdw_OBJ = ${dmdwSRC:.f90=.o}
dmdw_MODULES =${dmdw_MODULESRC:.f90=.o}

# Source files needed by DMDW/dym2feffinp.f90.
include DEP/dym2feffinp$(DEPTYPE).mk
# Object files needed by DMDW/dym2feffinp.f90.
dym2feffinp_OBJ = ${dym2feffinpSRC:.f90=.o}
dym2feffinp_MODULES =${dym2feffinp_MODULESRC:.f90=.o}

# Source files needed by EELS/eels.f90.
include DEP/eels$(DEPTYPE).mk
# Object files needed by EELS/eels.f90.
eels_OBJ = ${eelsSRC:.f90=.o}
eels_MODULES =${eels_MODULESRC:.f90=.o}

# Source files needed by FF2X/ff2x.f90.
include DEP/ff2x$(DEPTYPE).mk
# Object files needed by FF2X/ff2x.f90.
ff2x_OBJ = ${ff2xSRC:.f90=.o}
ff2x_MODULES =${ff2x_MODULESRC:.f90=.o}

# Source files needed by FMS/fms.f90.
include DEP/fms$(DEPTYPE).mk
# Object files needed by FMS/fms.f90.
fms_OBJ = ${fmsSRC:.f90=.o}
fms_MODULES =${fms_MODULESRC:.f90=.o}

# Source files needed by FULLSPECTRUM/fullspectrum.f90.
include DEP/fullspectrum$(DEPTYPE).mk
# Object files needed by FULLSPECTRUM/fullspectrum.f90.
fullspectrum_OBJ = ${fullspectrumSRC:.f90=.o}
fullspectrum_MODULES =${fullspectrum_MODULESRC:.f90=.o}

# Source files needed by GENFMT/genfmt.f90.
include DEP/genfmt$(DEPTYPE).mk
# Object files needed by GENFMT/genfmt.f90.
genfmt_OBJ = ${genfmtSRC:.f90=.o}
genfmt_MODULES =${genfmt_MODULESRC:.f90=.o}

# Source files needed by LDOS/ldos.f90.
include DEP/ldos$(DEPTYPE).mk
# Object files needed by LDOS/ldos.f90.
ldos_OBJ = ${ldosSRC:.f90=.o}
ldos_MODULES =${ldos_MODULESRC:.f90=.o}

# Source files needed by MKGTR/mkgtr.f90.
include DEP/mkgtr$(DEPTYPE).mk
# Object files needed by MKGTR/mkgtr.f90.
mkgtr_OBJ = ${mkgtrSRC:.f90=.o}
mkgtr_MODULES =${mkgtr_MODULESRC:.f90=.o}

# Source files needed by OPCONSAT/opconsat.f90.
include DEP/opconsat$(DEPTYPE).mk
# Object files needed by OPCONSAT/opconsat.f90.
opconsat_OBJ = ${opconsatSRC:.f90=.o}
opconsat_MODULES =${opconsat_MODULESRC:.f90=.o}

# Source files needed by PATH/path.f90.
include DEP/path$(DEPTYPE).mk
# Object files needed by PATH/path.f90.
path_OBJ = ${pathSRC:.f90=.o}
path_MODULES =${path_MODULESRC:.f90=.o}

# Source files needed by POT/pot.f90.
include DEP/pot$(DEPTYPE).mk
# Object files needed by POT/pot.f90.
pot_OBJ = ${potSRC:.f90=.o}
pot_MODULES =${pot_MODULESRC:.f90=.o}

# Source files needed by RDINP/rdinp.f90.
include DEP/rdinp$(DEPTYPE).mk
# Object files needed by RDINP/rdinp.f90.
rdinp_OBJ = ${rdinpSRC:.f90=.o}
rdinp_MODULES =${rdinp_MODULESRC:.f90=.o}

# Source files needed by RHORRP/rhorrp.f90.
include DEP/rhorrp$(DEPTYPE).mk
# Object files needed by RHORRP/rhorrp.f90.
rhorrp_OBJ = ${rhorrpSRC:.f90=.o}
rhorrp_MODULES =${rhorrp_MODULESRC:.f90=.o}

# Source files needed by RIXS/rixs.f90.
include DEP/rixs$(DEPTYPE).mk
# Object files needed by RIXS/rixs.f90.
rixs_OBJ = ${rixsSRC:.f90=.o}
rixs_MODULES =${rixs_MODULESRC:.f90=.o}

# Source files needed by SCREEN/screen.f90.
include DEP/screen$(DEPTYPE).mk
# Object files needed by SCREEN/screen.f90.
screen_OBJ = ${screenSRC:.f90=.o}
screen_MODULES =${screen_MODULESRC:.f90=.o}

# Source files needed by DYNSCREEN/screenw2.f90.
include DEP/screenw2$(DEPTYPE).mk
# Object files needed by DYNSCREEN/screenw2.f90.
screenw2_OBJ = ${screenw2SRC:.f90=.o}
screenw2_MODULES =${screenw2_MODULESRC:.f90=.o}

# Source files needed by SFCONV/sfconv.f90.
include DEP/sfconv$(DEPTYPE).mk
# Object files needed by SFCONV/sfconv.f90.
sfconv_OBJ = ${sfconvSRC:.f90=.o}
sfconv_MODULES =${sfconv_MODULESRC:.f90=.o}

# Source files needed by XSPH/xsph.f90.
include DEP/xsph$(DEPTYPE).mk
# Object files needed by XSPH/xsph.f90.
xsph_OBJ = ${xsphSRC:.f90=.o}
xsph_MODULES =${xsph_MODULESRC:.f90=.o}


########################################################################################################
###                                    MAKE TARGETS                                               ###
########################################################################################################

all:
	cp PAR/sequential.src PAR/par.f90
	mkdir -p ../bin/Seq
	make $(EXECUTABLES) "F90=$(F90)" "FLAGS=$(FLAGS)" "EXECDIR=../bin/Seq/"
	. Utility/MkFeffScript seq
	cp `find ./ -name *.script` ../bin/Seq/
	echo -e $(smessg1) $(smessg2) $(smessg3)

mpi:
	cp PAR/parallel.src PAR/par.f90
	mkdir -p ../bin/MPI
	make $(EXECUTABLES) "F90=$(MPIF90)" "FLAGS=$(MPIFLAGS)" "EXECDIR=../bin/MPI/" "FPPTASK=$(parFPPTASK)"
	. Utility/MkFeffScript mpi
	cp `find ./ -name *.script` ../bin/MPI/
	echo -e $(parmessg1) $(parmessg2) $(parmessg3)

src:
	cp PAR/sequential.src PAR/par.f90
	perl Utility/make_modular_sourcefiles.pl

srcmpi:
	cp PAR/parallel.src PAR/par.f90
	perl Utility/make_modular_sourcefiles.pl mpi

deps:
	python Utility/mkdep.py

# Make clean
clean:
	rm -rf */*.o */*.mod
#####################################################################################################
###                                  END MAKE TARGETS                                             ###
#####################################################################################################
#####################################################################################################
###                                  atomic                                             ###
#####################################################################################################
atomic:$(EXECDIR)/atomic
$(EXECDIR)/atomic: $(atomic_MODULES) $(atomic_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(atomic_MODULES) $(atomic_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(atomic_MODULES) $(atomic_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END atomic                                             ###
#####################################################################################################

#####################################################################################################
###                                  band                                             ###
#####################################################################################################
band:$(EXECDIR)/band
$(EXECDIR)/band: $(band_MODULES) $(band_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(band_MODULES) $(band_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(band_MODULES) $(band_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END band                                             ###
#####################################################################################################

#####################################################################################################
###                                  compton                                             ###
#####################################################################################################
compton:$(EXECDIR)/compton
$(EXECDIR)/compton: $(compton_MODULES) $(compton_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(compton_MODULES) $(compton_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(compton_MODULES) $(compton_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END compton                                             ###
#####################################################################################################

#####################################################################################################
###                                  crpa                                             ###
#####################################################################################################
crpa:$(EXECDIR)/crpa
$(EXECDIR)/crpa: $(crpa_MODULES) $(crpa_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(crpa_MODULES) $(crpa_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(crpa_MODULES) $(crpa_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END crpa                                             ###
#####################################################################################################

#####################################################################################################
###                                  dmdw                                             ###
#####################################################################################################
dmdw:$(EXECDIR)/dmdw
$(EXECDIR)/dmdw: $(dmdw_MODULES) $(dmdw_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(dmdw_MODULES) $(dmdw_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(dmdw_MODULES) $(dmdw_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END dmdw                                             ###
#####################################################################################################

#####################################################################################################
###                                  dym2feffinp                                             ###
#####################################################################################################
dym2feffinp:$(EXECDIR)/dym2feffinp
$(EXECDIR)/dym2feffinp: $(dym2feffinp_MODULES) $(dym2feffinp_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(dym2feffinp_MODULES) $(dym2feffinp_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(dym2feffinp_MODULES) $(dym2feffinp_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END dym2feffinp                                             ###
#####################################################################################################

#####################################################################################################
###                                  eels                                             ###
#####################################################################################################
eels:$(EXECDIR)/eels
$(EXECDIR)/eels: $(eels_MODULES) $(eels_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(eels_MODULES) $(eels_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(eels_MODULES) $(eels_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END eels                                             ###
#####################################################################################################

#####################################################################################################
###                                  ff2x                                             ###
#####################################################################################################
ff2x:$(EXECDIR)/ff2x
$(EXECDIR)/ff2x: $(ff2x_MODULES) $(ff2x_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(ff2x_MODULES) $(ff2x_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(ff2x_MODULES) $(ff2x_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END ff2x                                             ###
#####################################################################################################

#####################################################################################################
###                                  fms                                             ###
#####################################################################################################
fms:$(EXECDIR)/fms
$(EXECDIR)/fms: $(fms_MODULES) $(fms_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(fms_MODULES) $(fms_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(fms_MODULES) $(fms_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END fms                                             ###
#####################################################################################################

#####################################################################################################
###                                  fullspectrum                                             ###
#####################################################################################################
fullspectrum:$(EXECDIR)/fullspectrum
$(EXECDIR)/fullspectrum: $(fullspectrum_MODULES) $(fullspectrum_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(fullspectrum_MODULES) $(fullspectrum_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(fullspectrum_MODULES) $(fullspectrum_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END fullspectrum                                             ###
#####################################################################################################

#####################################################################################################
###                                  genfmt                                             ###
#####################################################################################################
genfmt:$(EXECDIR)/genfmt
$(EXECDIR)/genfmt: $(genfmt_MODULES) $(genfmt_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(genfmt_MODULES) $(genfmt_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(genfmt_MODULES) $(genfmt_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END genfmt                                             ###
#####################################################################################################

#####################################################################################################
###                                  ldos                                             ###
#####################################################################################################
ldos:$(EXECDIR)/ldos
$(EXECDIR)/ldos: $(ldos_MODULES) $(ldos_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(ldos_MODULES) $(ldos_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(ldos_MODULES) $(ldos_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END ldos                                             ###
#####################################################################################################

#####################################################################################################
###                                  mkgtr                                             ###
#####################################################################################################
mkgtr:$(EXECDIR)/mkgtr
$(EXECDIR)/mkgtr: $(mkgtr_MODULES) $(mkgtr_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(mkgtr_MODULES) $(mkgtr_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(mkgtr_MODULES) $(mkgtr_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END mkgtr                                             ###
#####################################################################################################

#####################################################################################################
###                                  opconsat                                             ###
#####################################################################################################
opconsat:$(EXECDIR)/opconsat
$(EXECDIR)/opconsat: $(opconsat_MODULES) $(opconsat_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(opconsat_MODULES) $(opconsat_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(opconsat_MODULES) $(opconsat_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END opconsat                                             ###
#####################################################################################################

#####################################################################################################
###                                  path                                             ###
#####################################################################################################
path:$(EXECDIR)/path
$(EXECDIR)/path: $(path_MODULES) $(path_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(path_MODULES) $(path_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(path_MODULES) $(path_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END path                                             ###
#####################################################################################################

#####################################################################################################
###                                  pot                                             ###
#####################################################################################################
pot:$(EXECDIR)/pot
$(EXECDIR)/pot: $(pot_MODULES) $(pot_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(pot_MODULES) $(pot_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(pot_MODULES) $(pot_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END pot                                             ###
#####################################################################################################

#####################################################################################################
###                                  rdinp                                             ###
#####################################################################################################
rdinp:$(EXECDIR)/rdinp
$(EXECDIR)/rdinp: $(rdinp_MODULES) $(rdinp_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(rdinp_MODULES) $(rdinp_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(rdinp_MODULES) $(rdinp_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END rdinp                                             ###
#####################################################################################################

#####################################################################################################
###                                  rhorrp                                             ###
#####################################################################################################
rhorrp:$(EXECDIR)/rhorrp
$(EXECDIR)/rhorrp: $(rhorrp_MODULES) $(rhorrp_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(rhorrp_MODULES) $(rhorrp_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(rhorrp_MODULES) $(rhorrp_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END rhorrp                                             ###
#####################################################################################################

#####################################################################################################
###                                  rixs                                             ###
#####################################################################################################
rixs:$(EXECDIR)/rixs
$(EXECDIR)/rixs: $(rixs_MODULES) $(rixs_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(rixs_MODULES) $(rixs_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(rixs_MODULES) $(rixs_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END rixs                                             ###
#####################################################################################################

#####################################################################################################
###                                  screen                                             ###
#####################################################################################################
screen:$(EXECDIR)/screen
$(EXECDIR)/screen: $(screen_MODULES) $(screen_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(screen_MODULES) $(screen_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(screen_MODULES) $(screen_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END screen                                             ###
#####################################################################################################

#####################################################################################################
###                                  screenw2                                             ###
#####################################################################################################
screenw2:$(EXECDIR)/screenw2
$(EXECDIR)/screenw2: $(screenw2_MODULES) $(screenw2_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(screenw2_MODULES) $(screenw2_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(screenw2_MODULES) $(screenw2_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END screenw2                                             ###
#####################################################################################################

#####################################################################################################
###                                  sfconv                                             ###
#####################################################################################################
sfconv:$(EXECDIR)/sfconv
$(EXECDIR)/sfconv: $(sfconv_MODULES) $(sfconv_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(sfconv_MODULES) $(sfconv_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(sfconv_MODULES) $(sfconv_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END sfconv                                             ###
#####################################################################################################

#####################################################################################################
###                                  xsph                                             ###
#####################################################################################################
xsph:$(EXECDIR)/xsph
$(EXECDIR)/xsph: $(xsph_MODULES) $(xsph_OBJ) 
	@$(F90) $(FLAGS) $(INCLUDEFLAGS) -o $@ $(xsph_MODULES) $(xsph_OBJ) $(LDFLAGS)
	@echo -e $(F90) $(FLAGS) INCLUDEFLAGS -o $@ $(xsph_MODULES) $(xsph_OBJ) $(LDFLAGS)
#####################################################################################################
###                              END xsph                                             ###
#####################################################################################################

include DEP/dependencies.mk

# How to make a .o file from .f90
.f90.o:
	@echo "cd $(<D); $(F90) -c $(FLAGS) $(FPPTASK) INCLUDEFLAGS $(<F)"
	@cd $(<D); $(F90) -c $(FLAGS) $(FPPTASK) $(INCLUDEFLAGS) $(<F)

.f.o:
	@echo "cd $(<D); $(F90) -FI -c $(FLAGS) $(FPPTASK) INCLUDEFLAGS $(<F)"
	@cd $(<D); $(F90) -FI -c $(FLAGS) $(FPPTASK) $(INCLUDEFLAGS) $(<F)

.mod.o:
	

