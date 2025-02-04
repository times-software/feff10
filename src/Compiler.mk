##### This file contains compiler settings for the FEFF9 project
## Define below:
## - F90    :  choice of fortran90 compiler
## - FLAGS  :  compilation options
## - MPIF90 :  mpif90 compiler (optional)
## - USE_MKL:  blas/lapack libraries (optional)

########################################
#              ifort
########################################
# F90 = ifort 
# for a MacBook1,1 with a 32-bit Intel Core Duo processor ("Yonah") running Snow Leopard - generate 32bit executable:
# FLAGS = -O3  -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib/ -m32 -O2 -axsse4.2,sse3 -prec_div -132  -g 
# Regular 64-bit code e.g. for iMac, (recent-ish) MacBook Pro, etc. (10.7, 10.6) :
# Full-on debugging mode for ifort:  (note this catches different things depending on optimization settings)
# FLAGS = -O3 -check arg_temp_created -gen-interfaces -warn interfaces -g -fp-stack-check -traceback  -FR -heap-arrays -check bounds 
# FLAGS = -O0  -g -fp-stack-check -traceback -heap-arrays -check bounds 
# Production code Mac: (note this runs for Leopard and higher without performance sacrifice compared to 10.8 optimizations) :
# FLAGS =   -O3 -mmacosx-version-min=10.5 $(FCINCLUDE)
# Production code Linux:
# FLAGS = -O3

########################################
#              pgf90
########################################
# F90 = pgf90
# FLAGS = -O3

########################################
#              gfortran
########################################
 F90 =  gfortran
 #FLAGS = -O0 -ffree-line-length-none -fcheck=bounds -g  -ffpe-trap=invalid
 FLAGS = -O0 -ffree-line-length-none -fcheck=bounds -ffpe-trap=invalid -g -fallow-argument-mismatch

########################################
#              mkl
########################################
# # MAIN SWITCH (toggle mkl on/off) :
# uncomment to use MKL
# comment to use standard FEFF blas/lapack (MATH/lu.f90)
# USE_MKL = "on"

# Must edit for local mkl installation
# $MKLROOT must exist (optionally you can define it here)

## Mac MKL:
 MKLSEQUENTIALMAC = $(MKLROOT)/lib/libmkl_blas95_lp64.a $(MKLROOT)/lib/libmkl_lapack95_lp64.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_sequential.a $(MKLROOT)/lib/libmkl_core.a -lpthread -lm
 MKLTHREADEDMAC =   $(MKLROOT)/lib/libmkl_blas95_lp64.a $(MKLROOT)/lib/libmkl_lapack95_lp64.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm
## Intel MKL:
 MKLSEQUENTIALLINUX =  $(MKLROOT)/lib/em64t/libmkl_blas95_lp64.a $(MKLROOT)/lib/em64t/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread -lm
# # Windows MKL:
#    On Windows I use Visual Studio rather than this Makefile/Compiler.mk structure.
#    Just enable "Use math libraries - MKL - Sequential" in the Project Settings
#    For the time-consuming modules (pot, fms, screen, ldos, compton)

#  Choose the right OS here :
 MKL_LDFLAGS = $(MKLSEQUENTIALMAC)
# MKL_LDFLAGS = $(MKLSEQUENTIALLINUX)
 MKL_FCINCLUDE = -I$(MKLROOT)/include/em64t/lp64 -I$(MKLROOT)/include

#  If no MKL, use these defaults instead (lu.f90 provides non-optimized blas/lapack)
 FEFF_LDFLAGS = 
 FEFF_FCINCLUDE =

# #### Settings derived here :
ifdef USE_MKL
#  1/ Use the MKL settings defined above:
	LDFLAGS = $(MKL_LDFLAGS)
	FCINCLUDE = $(MKL_FCINCLUDE)
	DEPTYPE = _MKL
else
#  2/ Use the standard FEFF blas/lapack (slower but no need to install MKL):
	LDFLAGS = 
	FCINCLUDE = 
	DEPTYPE = 
endif

########################################
#              mpi
########################################
 MPIF90 = mpif90
 MPIFLAGS =  -O3 -ffree-line-length-none
# MPIFLAGS = -O3 -check arg_temp_created -gen-interfaces -warn interfaces -g -fp-stack-check -traceback  -heap-arrays -check bounds
