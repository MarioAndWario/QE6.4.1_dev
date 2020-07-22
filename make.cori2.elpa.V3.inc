# make.inc.  Generated from make.inc.in by configure.
#  1) modules/3.2.11.4                                 13) dvs/2.12_2.2.156-7.0.1.1_8.9__g5aab709e
#  2) altd/2.0                                         14) alps/6.6.58-7.0.1.1_6.4__g437d88db.ari
#  3) craype-network-aries                             15) rca/2.2.20-7.0.1.1_4.46__g8e3fb5b.ari
#  4) craype/2.6.2                                     16) atp/2.1.3
#  5) cray-libsci/19.06.1                              17) PrgEnv-intel/6.0.5
#  6) udreg/2.3.2-7.0.1.1_3.31__g8175d3d.ari           18) craype-haswell
#  7) ugni/6.0.14.0-7.0.1.1_7.33__ge78e5b0.ari         19) cray-mpich/7.7.10
#  8) pmi/5.0.14                                       20) craype-hugepages2M
#  9) dmapp/7.1.1-7.0.1.1_4.48__g38cf134.ari           21) intel/19.0.3.199
# 10) gni-headers/5.0.12.0-7.0.1.1_6.28__g3b1768f.ari  22) cray-hdf5-parallel/1.10.5.2
# 11) xpmem/2.2.20-7.0.1.1_4.10__g0475745.ari          23) cray-fftw/3.3.8.4
# 12) job/2.2.4-7.0.1.1_3.36__g36b56f4.ari             24) python/3.7-anaconda-2019.10
# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
# 	$(CPP) $(CPPFLAGS) $< -o $*.F90
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<

# .f.o and .c.o: do not modify

.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS)  -c $<



# Top QE directory, useful for locating libraries,  linking QE with plugins
# The following syntax should always point to TOPDIR:
TOPDIR = $(dir $(abspath $(filter %make.inc,$(MAKEFILE_LIST))))
# if it doesn't work, uncomment the following line (edit if needed):

# TOPDIR = /global/homes/m/mwu/softwares/QE/QE6.5

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# To use libxc (v>=3.0.1), add -D__LIBXC to DFLAGS
# See include/defs.h.README for a list of options and their meaning
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas

# MANUAL_DFLAGS  = additional precompilation option(s), if desired
#                  BEWARE: it does not work for IBM xlf! Manually edit FDFLAGS
MANUAL_DFLAGS  =
DFLAGS         =  -D__DFTI -D__EXX_ACE -D__MPI -D__OPENMP -D__SCALAPACK -D__ELPA_2017
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I$(TOPDIR)/include -I/some/other/directory/
#          the latter contains .e.g. files needed by FFT libraries
# for libxc add -I/path/to/libxc/include/

IFLAGS       = -I$(TOPDIR)/include -I$(TOPDIR)/FoX/finclude -I$(TOPDIR)/S3DE/iotk/include/ -I$(MKLROOT)/include/fftw -I$(MKLROOT)/include -I$(HDF5_DIR)/include -I$(ELPA_DIR)/include/elpa_openmp-2017.05.001.rc1/modules -I$(ELPA_DIR)/include/elpa_openmp-2017.05.001.rc1/elpa

# MOD_FLAG = flag used by f90 compiler to locate modules

MOD_FLAG      = -I

# BASEMOD_FLAGS points to directories containing basic modules,
# while BASEMODS points to the corresponding module libraries
# Each Makefile can add directories to MODFLAGS and libraries to QEMODS

BASEMOD_FLAGS= $(MOD_FLAG)$(TOPDIR)/iotk/src \
               $(MOD_FLAG)$(TOPDIR)/Modules \
               $(MOD_FLAG)$(TOPDIR)/FFTXlib \
	       $(MOD_FLAG)$(TOPDIR)/LAXlib \
	       $(MOD_FLAG)$(TOPDIR)/UtilXlib \
	       $(MOD_FLAG)$(TOPDIR)/FoX/finclude

# Compilers: fortran-90, fortran-77, C
# If a parallel compilation is desired, MPIF90 should be a fortran-90
# compiler that produces executables for parallel execution using MPI
# (such as for instance mpif90, mpf90, mpxlf90,...);
# otherwise, an ordinary fortran-90 compiler (f90, g95, xlf90, ifort,...)
# If you have a parallel machine but no suitable candidate for MPIF90,
# try to specify the directory containing "mpif.h" in IFLAGS
# and to specify the location of MPI libraries in MPI_LIBS

MPIF90         = ftn
F90            = ftn
CC             = cc
F77            = ftn

# GPU architecture (Kepler: 35, Pascal: 60, Volta: 70 )
GPU_ARCH=

# CUDA runtime (Pascal: 8.0, Volta: 9.0)
CUDA_RUNTIME=

# CUDA F90 Flags
CUDA_F90FLAGS=

# C preprocessor and preprocessing flags - for explicit preprocessing,
# if needed (see the compilation rules above)
# preprocessing flags must include DFLAGS and IFLAGS

CPP            = cpp
CPPFLAGS       = -P -traditional -Uvector $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90, F77
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = -O3 -xCORE-AVX2 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias $(DFLAGS) $(IFLAGS)
F90FLAGS       = $(FFLAGS) -nomodule -qopenmp -fpp $(FDFLAGS) $(CUDA_F90FLAGS) $(IFLAGS) $(MODFLAGS)
FFLAGS         = -O2 -xCORE-AVX2 -align array64byte -threads -heap-arrays 4096 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl -qopenmp

# compiler flags without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack

FFLAGS_NOOPT   = -O0 -assume byterecl

# compiler flag needed by some compilers when the main program is not fortran
# Currently used for Yambo

FFLAGS_NOMAIN   = -nofor_main

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty
# for libxc, set LD_LIBS=-L/path/to/libxc/lib/ -lxcf90 -lxcf03 -lxc

LD             = ftn
LDFLAGS        = -static -qopenmp
LD_LIBS        = -L$(ELPA_DIR)/lib -lelpa_openmp

# External Libraries (if any) : blas, lapack, fft, MPI

# If you have nothing better, use the local copy via "--with-netlib" :
# BLAS_LIBS = /your/path/to/espresso/LAPACK/blas.a
# BLAS_LIBS_SWITCH = internal

BLAS_LIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
BLAS_LIBS_SWITCH = external

# If you have nothing better, use the local copy via "--with-netlib" :
# LAPACK_LIBS = /your/path/to/espresso/LAPACK/lapack.a
# LAPACK_LIBS_SWITCH = internal
# For IBM machines with essl (-D__ESSL): load essl BEFORE lapack !
# remember that LAPACK_LIBS precedes BLAS_LIBS in loading order

LAPACK_LIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group
LAPACK_LIBS_SWITCH = external

SCALAPACK_LIBS = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a

# nothing needed here if the the internal copy of FFTW is compiled
# (needs -D__FFTW in DFLAGS)

FFT_LIBS =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group

# HDF5
HDF5_LIBS = -Wl,-rpath,$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib -lz -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
FOX_LIB  = -L$(TOPDIR)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common\
            -lFoX_utils -lFoX_fsys 
FOX_FLAGS = 
# For parallel execution, the correct path to MPI libraries must
# be specified in MPI_LIBS (except for IBM if you use mpxlf)

MPI_LIBS       = 

# IBM-specific: MASS libraries, if available and if -D__MASS is defined in FDFLAGS

MASS_LIBS      = 

# CUDA libraries
CUDA_LIBS=
CUDA_EXTLIBS = 

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv

AR             = ar
ARFLAGS        = ruv

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = ranlib

# all internal and external libraries - do not modify

FLIB_TARGETS   = all

LIBOBJS        = $(TOPDIR)/clib/clib.a  $(TOPDIR)/iotk/src/libiotk.a
LIBXC_LIBS     = 
QELIBS         = $(CUDA_LIBS) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FOX_LIB) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(HDF5_LIBS) $(LIBXC_LIBS) $(LD_LIBS)

# wget or curl - useful to download from network
WGET = wget -O

# Install directory - "make install" copies *.x executables there
PREFIX = /usr/local
