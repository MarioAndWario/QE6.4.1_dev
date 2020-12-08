# make.inc.  Generated from make.inc.in by configure.
  # 1) emacs/25.1             3) python/3.6             5) impi/2018.5.274
  # 2) vim/7.4                4) intel/2018.5.274.par   6) fftw/fftw
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

# TOPDIR = /admin/build/admin/rpms/frontera/BUILD/qe-6.4.1-TACC/q-e-qe-6.4.1

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# To use libxc (v>=3.0.1), add -D__LIBXC to DFLAGS
# See include/defs.h.README for a list of options and their meaning
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas

# MANUAL_DFLAGS  = additional precompilation option(s), if desired
#                  BEWARE: it does not work for IBM xlf! Manually edit FDFLAGS
MANUAL_DFLAGS  =
DFLAGS         =  -D__DFTI -D__EXX_ACE -D__MPI -D__SCALAPACK -D__ELPA_2017
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I$(TOPDIR)/include -I/some/other/directory/
#          the latter contains .e.g. files needed by FFT libraries
# for libxc add -I/path/to/libxc/include/

IFLAGS = -I$(MKLROOT)/include/fftw -I$(MKLROOT)/include -I$(TOPDIR)/include -I$(TOPDIR)/FoX/finclude -I$(TOPDIR)/S3DE/iotk/include/ -I$(ELPA_DIR)/include/2017.11.001/modules -I$(ELPA_DIR)/include/2017.11.001/elpa #-I$(HDF5_DIR)/include 

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

MPIF90         = mpiifort #mpif90
F90            = ifort
CC             = icc
F77            = ifort

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
CPPFLAGS       = -P -traditional $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90, F77
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = -O3  -xCORE-AVX512 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias $(DFLAGS) $(IFLAGS)
F90FLAGS       = $(FFLAGS) -nomodule -fpp $(FDFLAGS) $(CUDA_F90FLAGS) $(IFLAGS) $(MODFLAGS)
FFLAGS         = -O2 -align array64byte -threads -heap-arrays 4096  -xCORE-AVX512 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl

# compiler flags without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack

FFLAGS_NOOPT   = -O0 -assume byterecl -g -traceback

# compiler flag needed by some compilers when the main program is not fortran
# Currently used for Yambo

FFLAGS_NOMAIN   = -nofor_main

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty
# for libxc, set LD_LIBS=-L/path/to/libxc/lib/ -lxcf90 -lxc

LD             = mpiifort #mpif90
LDFLAGS = -static-intel -static-libgcc -static-libstdc++
LD_LIBS        = -Wl,--as-needed -liomp5 -Wl,--no-as-needed

# External Libraries (if any) : blas, lapack, fft, MPI

# If you have nothing better, use the local copy via "--with-netlib" :
# BLAS_LIBS = /your/path/to/espresso/LAPACK/blas.a
# BLAS_LIBS_SWITCH = internal

BLAS_LIBS      =  -Wl,--start-group     $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a     $(MKLROOT)/lib/intel64/libmkl_core.a     $(MKLROOT)/lib/intel64/libmkl_intel_thread.a     $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a   -Wl,--end-group
BLAS_LIBS_SWITCH = external

# If you have nothing better, use the local copy via "--with-netlib" :
# LAPACK_LIBS = /your/path/to/espresso/LAPACK/lapack.a
# LAPACK_LIBS_SWITCH = internal
# For IBM machines with essl (-D__ESSL): load essl BEFORE lapack !
# remember that LAPACK_LIBS precedes BLAS_LIBS in loading order

LAPACK_LIBS    =  -Wl,--start-group     $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a     $(MKLROOT)/lib/intel64/libmkl_core.a     $(MKLROOT)/lib/intel64/libmkl_intel_thread.a     $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a   -Wl,--end-group
LAPACK_LIBS_SWITCH = external

SCALAPACK_LIBS = $(ELPA_DIR)/lib/libelpa.a $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a

# nothing needed here if the the internal copy of FFTW is compiled
# (needs -D__DFTI -D__EXX_ACE in DFLAGS)

FFT_LIBS       =   -Wl,--start-group     $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a     $(MKLROOT)/lib/intel64/libmkl_core.a     $(MKLROOT)/lib/intel64/libmkl_intel_thread.a     $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a   -Wl,--end-group

# HDF5
HDF5_LIB = #-Wl,-rpath,$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib -lz -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
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

AR             = xiar
ARFLAGS        = ruv

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = ranlib

# all internal and external libraries - do not modify

FLIB_TARGETS   = all

LIBOBJS        = $(TOPDIR)/clib/clib.a  $(TOPDIR)/iotk/src/libiotk.a
LIBXC_LIBS     = 
QELIBS         = $(CUDA_LIBS) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FOX_LIB) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(HDF5_LIB) $(LIBXC_LIBS) $(LD_LIBS)

# wget or curl - useful to download from network
WGET = wget -O

# Install directory - "make install" copies *.x executables there
PREFIX = /usr/local

#############################################################################
# Copyright (c) 2014-2019, Intel Corporation                                #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions        #
# are met:                                                                  #
# 1. Redistributions of source code must retain the above copyright         #
#    notice, this list of conditions and the following disclaimer.          #
# 2. Redistributions in binary form must reproduce the above copyright      #
#    notice, this list of conditions and the following disclaimer in the    #
#    documentation and/or other materials provided with the distribution.   #
# 3. Neither the name of the copyright holder nor the names of its          #
#    contributors may be used to endorse or promote products derived        #
#    from this software without specific prior written permission.          #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         #
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     #
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      #
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  #
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#############################################################################
# Hans Pabst (Intel Corp.)
#############################################################################

# TBB runtime made for oldest/latest supported GCC
TBBGCC_OLD ?= 1
# ultimately disable/enable Intel TBB's malloc proxy
TBBMALLOC ?= 0

ifneq (0,$(TBBMALLOC))
  ifneq (,$(TBBROOT))
    GCC = $(notdir $(shell which gcc 2> /dev/null))
    ifneq (,$(GCC))
      GCC_VERSION_STRING = $(shell $(GCC) --version 2> /dev/null | head -n1 | sed "s/..* \([0-9][0-9]*\.[0-9][0-9]*\.*[0-9]*\)[ \S]*.*/\1/")
      GCC_VERSION_MAJOR = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f1)
      GCC_VERSION_MINOR = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f2)
      GCC_VERSION_PATCH = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f3)
      TBBLIBDIR = $(TBBROOT)/lib/intel64/gcc$(GCC_VERSION_MAJOR).$(GCC_VERSION_MINOR)
      TBBMALLOCLIB = $(wildcard $(TBBLIBDIR)/libtbbmalloc_proxy.so)
    endif
    ifeq (,$(TBBMALLOCLIB))
      ifneq (0,$(TBBGCC_OLD))
        TBBGCCDIR = $(shell ls -1 "$(TBBROOT)/lib/intel64" | tr "\n" " " | cut -d" " -f1)
      else
        TBBGCCDIR = $(shell ls -1 "$(TBBROOT)/lib/intel64" | tr "\n" " " | rev | cut -d" " -f2 | rev)
      endif
      TBBLIBDIR = $(TBBROOT)/lib/intel64/$(TBBGCCDIR)
      TBBMALLOCLIB = $(wildcard $(TBBLIBDIR)/libtbbmalloc_proxy.so)
    endif
    ifneq (,$(TBBMALLOCLIB))
      LIBS += $(TBBMALLOCLIB) $(TBBLIBDIR)/libtbbmalloc.so
      ifneq (1,$(TBBMALLOC)) # TBBMALLOC=2
        FCFLAGS += -heap-arrays
      endif
    endif
  endif
endif


#############################################################################
# Copyright (c) 2017-2019, Intel Corporation                                #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions        #
# are met:                                                                  #
# 1. Redistributions of source code must retain the above copyright         #
#    notice, this list of conditions and the following disclaimer.          #
# 2. Redistributions in binary form must reproduce the above copyright      #
#    notice, this list of conditions and the following disclaimer in the    #
#    documentation and/or other materials provided with the distribution.   #
# 3. Neither the name of the copyright holder nor the names of its          #
#    contributors may be used to endorse or promote products derived        #
#    from this software without specific prior written permission.          #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         #
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     #
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      #
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  #
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#############################################################################
# Hans Pabst (Intel Corp.)
#############################################################################

# LIBXSMM (https://github.com/hfp/libxsmm)
#
#LIBXSMMROOT = /path/to/libxsmm

ifneq (,$(LIBXSMMROOT))
ifneq (,$(wildcard $(LIBXSMMROOT)/lib/libxsmmext.a))
  WRAP ?= 1
  ifeq (2,$(WRAP))
    LDFLAGS += -Wl,--wrap=dgemm_
  else
    LDFLAGS += -Wl,--wrap=sgemm_,--wrap=dgemm_
  endif
  LDFLAGS += -Wl,--export-dynamic
  LAPACK_LIBS := $(LIBXSMMROOT)/lib/libxsmmext.a $(LIBXSMMROOT)/lib/libxsmm.a $(LAPACK_LIBS)
  ifeq (,$(OPENMP))
  ifeq (sequential,$(MKL_OMPRTL))
    LD_LIBS += -liomp5
  endif
  endif
else
  $(info LIBXSMM library not found at $(LIBXSMMROOT))
endif
endif


default: all

init_us_1.o: init_us_1.f90
	$(MPIF90) $(F90FLAGS) -O1 -c $<

new_ns.o: new_ns.f90
	$(MPIF90) $(F90FLAGS) -O1 -c $<

us_exx.o: us_exx.f90
	$(MPIF90) $(F90FLAGS) -qopenmp -qoverride_limits -c $<

wypos.o: wypos.f90
	$(MPIF90) $(F90FLAGS) -qopenmp -qoverride_limits -O0 -c $<

lr_apply_liouvillian.o: lr_apply_liouvillian.f90
	$(MPIF90) $(F90FLAGS) -qopenmp -qoverride_limits -O1 -c $<

realus.o: realus.f90
	$(MPIF90) $(F90FLAGS) -qopenmp -qoverride_limits -O1 -c $<

