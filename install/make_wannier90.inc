#=======================================
# WANNIER90 	
#=======================================

TOPDIR = /home1/03355/tg826544/software/QE/QE6.4.1/QE6.4.1_dev

F90=mpif90
FCOPTS=$(FFLAGS) -nomodule 
LDOPTS= 

LIBS =    -lmkl_intel_lp64  -lmkl_sequential -lmkl_core 
