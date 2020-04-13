# FMKMF_F90 set to gfortran -O3
# FMKMF_SFTAG not set: using f90
# FMKMF_SPATH not set: using . 
# FMKMF_LINKOPTS not set: using no link options 
# Main program is pot_generator_test.f90 
# process_fsource called with arg pot_generator_test.f90 
# pot_generator_test.f90 Uses Module m_Kinds
# pot_generator_test.f90 Uses Module m_Pot_Generator
# pot_generator_test.f90 Uses Module m_PTable
# Full list of modules in pot_generator_test.f90: m_Kinds m_Pot_Generator m_PTable 
# Uses m_Kinds which is in ./m_kinds.f90
# process_fsource called with arg ./m_kinds.f90 
# Full list of modules in ./m_kinds.f90:  
# Uses m_Pot_Generator which is in ./m_pot_generator.f90
# process_fsource called with arg ./m_pot_generator.f90 
# ./m_pot_generator.f90 Uses Module m_Kinds
# ./m_pot_generator.f90 Uses Module m_Strings
# ./m_pot_generator.f90 Uses Module m_PTable
# Full list of modules in ./m_pot_generator.f90: m_Kinds m_Strings m_PTable 
# Uses m_Kinds which is in ./m_kinds.f90
# ./m_kinds.f90 already in list
# Uses m_Strings which is in ./m_strings.f90
# process_fsource called with arg ./m_strings.f90 
# Full list of modules in ./m_strings.f90:  
# Uses m_PTable which is in ./m_ptable.f90
# process_fsource called with arg ./m_ptable.f90 
# ./m_ptable.f90 Uses Module m_Kinds
# ./m_ptable.f90 Uses Module m_Strings
# Full list of modules in ./m_ptable.f90: m_Kinds m_Strings 
# Uses m_Kinds which is in ./m_kinds.f90
# ./m_kinds.f90 already in list
# Uses m_Strings which is in ./m_strings.f90
# ./m_strings.f90 already in list
# Uses m_PTable which is in ./m_ptable.f90
# ./m_ptable.f90 already in list

# ------------------Macro-Defs---------------------
F90=gfortran -O3 

# -------------------End-macro-Defs---------------------------

# Here is the link step 
pot_generator_test:../COMMON/m_kinds.o m_strings.o m_ptable.o m_pot_generator.o pot_generator_test.o 
	 $(F90) -o pot_generator_test ../COMMON/m_kinds.o m_strings.o m_ptable.o m_pot_generator.o pot_generator_test.o   

# Here are the compile steps
 
m_kinds.o:./m_kinds.f90  
	 $(F90) -c ./m_kinds.f90 

m_strings.o:./m_strings.f90  
	 $(F90) -c ./m_strings.f90 

m_ptable.o:./m_ptable.f90 m_kinds.o m_strings.o 
	 $(F90) -c ./m_ptable.f90 

m_pot_generator.o:./m_pot_generator.f90 m_kinds.o m_strings.o m_ptable.o 
	 $(F90) -c ./m_pot_generator.f90 

pot_generator_test.o:pot_generator_test.f90 m_kinds.o m_pot_generator.o m_ptable.o 
	 $(F90) -c pot_generator_test.f90 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
