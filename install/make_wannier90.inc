#=======================================
# WANNIER90 	
#=======================================

TOPDIR = /global/homes/m/mwu/softwares/QE/QE6.4.1/QE6.4.1_dev_elpa

F90=ifort
FCOPTS=$(FFLAGS) -nomodule 
LDOPTS= 

LIBS =    -lmkl_intel_lp64  -lmkl_sequential -lmkl_core 
