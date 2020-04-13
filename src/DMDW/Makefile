# ------------------Macro-Defs---------------------
#F90=ifort -O3 -fpp -g -zero
#F90=pgf90 -O3 -fastsse -Mnoflushz
#F90=pgf90 -g
#F90=pgf90 -g -O0 -Mbounds -Ktrap=fp -Mchkfpstk -Mchkstk
#F90=g95 -pedantic -Wall -O3
#F90=g95 -O3 -cpp -ffree-line-length-huge -ffixed-line-length-132
F90=gfortran -g -fbacktrace -O3 -cpp -ffree-line-length-none
#F90=gfortran -g -enable-checking -Q -O0 -cpp -ffree-line-length-none
FEFFMODULE='-fpp -DFEFF'
# -------------------End-macro-Defs---------------------------

include vars.defs

all: dmdw dym2feffinp

# Here is the link step 
dmdw:m_kinds.o m_const_and_conv.o m_strings.o m_ptable.o m_math.o m_dmdw.o dmdw.o m_errorfile.o m_controls.o m_kklist.o m_strfacs.o str2dp.o str.o rdhead.o
	 $(F90) -o dmdw m_kinds.o m_const_and_conv.o m_strings.o \
	m_errorfile.o m_controls.o m_kklist.o \
	m_strfacs.o m_ptable.o m_dmdw.o m_math.o dmdw.o \
	str2dp.o str.o rdhead.o 
dym2feffinp:m_kinds.o m_const_and_conv.o  m_errorfile.o  m_controls.o m_kklist.o m_strfacs.o m_strings.o m_ptable.o m_cmdline.o m_math.o  \
	 m_dmdw.o dym2feffinp.o str2dp.o str.o rdhead.o
	 $(F90) -o dym2feffinp m_kinds.o m_const_and_conv.o  m_errorfile.o  m_controls.o m_kklist.o m_strfacs.o m_strings.o m_ptable.o \
	 m_cmdline.o m_dmdw.o m_math.o dym2feffinp.o str2dp.o str.o rdhead.o

# Here are the compile steps

m_strfacs.o:$(KSPACE)/m_strfacs.f90
	$(F90) -c $(KSPACE)/m_strfacs.f90

m_controls.o:$(KSPACE)/m_controls.f90
	$(F90) -c $(KSPACE)/m_controls.f90

m_kklist.o:$(KSPACE)/m_kklist.f90
	$(F90) -c $(KSPACE)/m_kklist.f90

m_constants.o:$(COMMON)/m_constants.f90  
	 $(F90) -c $(COMMON)/m_constants.f90 

m_errorfile.o:$(ERROR)/m_errorfile.f90
	$(F90) -c $(ERROR)/m_errorfile.f90
 
m_par.o:$(PAR)/m_par.f90  
	 $(F90) -c $(PAR)/m_par.f90 

par.o:$(PAR)/par.f90  m_par.o
	 $(F90) -c $(PAR)/par.f90 

str2dp.o:$(COMMON)/str2dp.f90
	 $(F90) -c $(COMMON)/str2dp.f90 

str.o:$(COMMON)/str.f90 str2dp.o
	 $(F90) -c $(COMMON)/str.f90 

rdhead.o:$(COMMON)/rdhead.f90  
	 $(F90) -c $(COMMON)/rdhead.f90 

m_kinds.o:$(COMMON)/m_kinds.f90  
	 $(F90) -c $(COMMON)/m_kinds.f90 

m_math.o:m_math.f90 m_kinds.o
	 $(F90) -c m_math.f90 

m_const_and_conv.o:m_const_and_conv.f90 m_kinds.o 
	 $(F90) -c ./m_const_and_conv.f90 

m_strings.o:$(INPGEN)/m_strings.f90
	 $(F90) -c $(INPGEN)/m_strings.f90 

m_ptable.o:$(INPGEN)/m_ptable.f90 m_strings.o m_kinds.o
	 $(F90) -c $(INPGEN)/m_ptable.f90 

m_cmdline.o:m_cmdline.f90
	 $(F90) -c ./m_cmdline.f90

m_dmdw.o:./m_dmdw.f90 m_kinds.o m_const_and_conv.o m_strings.o
	 $(F90) -c ./m_dmdw.f90 

dmdw.o:dmdw.f90 m_kinds.o m_const_and_conv.o m_dmdw.o m_errorfile.o
	 $(F90) -c dmdw.f90

dym2feffinp.o:dym2feffinp.f90 m_kinds.o m_const_and_conv.o m_dmdw.o m_errorfile.o 
	 $(F90) -c dym2feffinp.f90 

# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo *.exe \
	V*.f *.dbg album F.err \
	$(INST_DIR)/dmdw $(INST_DIR)/dmdw.exe dmdw \
	$(INST_DIR)/dym2feffinp $(INST_DIR)/dym2feffinp.exe dym2feffinp \
	$(INST_DIR)/ai2bdf.script $(INST_DIR)/dynG2dym.script \
	$(INST_DIR)/fchk2dym.script
  
clean_examples:
	rm -f $(TEST_DIR)/*.log $(TEST_DIR)/*.out

clean_standalone:
	rm -f -r dmdw_standalone

clean_all: clean clean_examples clean_standalone

install:
	cp dmdw dym2feffinp dynG2dym.script fchk2dym.script ai2bdf.script \
	$(INST_DIR)

examples:
	cd $(TEST_DIR); ./Run_Examples $(INST_DIR2)

STANDALONE_FILES=Makefile ai2bdf.script dmdw.f90 dym2feffinp.f90 \
fchk2dym.script dynG2dym.script m_cmdline.f90 m_const_and_conv.f90 \
m_dmdw.f90 $(INPGEN)/m_ptable.f90 $(COMMON)/m_kinds.f90 m_math.f90 \
$(ERROR)/m_errorfile.f90 $(COMMON)/m_inpmodules.f90 $(KSPACE)/m_controls.f90 \
$(KSPACE)/m_kklist.f90 $(KSPACE)/m_strfacs.f90 \
$(PAR)/m_par.f90 $(PAR)/par.f90 \
$(COMMON)/str.f90 $(COMMON)/chopen.f90 $(COMMON)/rdhead.f90 \
$(INPGEN)/m_strings.f90 $(COMMON)/str2dp.f90 $(COMMON)/m_constants.f90

standalone:
	mkdir -p dmdw_standalone/src
	mkdir -p dmdw_standalone/bin
	cp $(STANDALONE_FILES) dmdw_standalone/src
	cp standalone_vars.defs dmdw_standalone/src/vars.defs
	cp DMDW_Documentation.txt dmdw_standalone
	cp -R $(TEST_DIR) dmdw_standalone/examples
	cp -R $(TUT_DIR) dmdw_standalone
#	tar -cf dmdw_standalone.tar dmdw_standalone
#	gzip dmdw_standalone.tar
#	rm -r dmdw_standalone

