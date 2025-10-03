###############################################################################
#
# File      : common.mk
# Author    : Alex Stivala (astivala)
# Created   : July 2008
#
# $Id: Makefile 1690 2008-07-16 04:35:47Z astivala $
#
# Definitions of compilers and compile options for all Makefiles.
# Use GNU make.
#
# set MODE=DEBUG to build with debugging and profiling on. Otherwise
# default is to build with optimizations on and no debug or profile.
#
#
###############################################################################



# Includes and libraries from UMFPACK (included in SuiteSparse)
# http://www.cise.ufl.edu/research/sparse/umfpack/
# https://people.engr.tamu.edu/davis/suitesparse.html
UMFPACKINCS = -I/usr/include/suitesparse
UMFPACKLDIRS = 
UMFPACKLIBS = -lumfpack -lamd


# PARDISO libraries (Intel MKL version)
MKLROOT = /usr/local/intel/mkl/10.0.4.023
MKL = $(MKLROOT)/lib/em64t
# requires libmkl_solver_lp64 for PARDISO for MKL 10.0.4 (tango.vpac.org)
MKL_LIBS = -L/usr/local/intel/Compiler/11.0/081/bin/intel64 -L$(MKL) $(MKL)/libmkl_solver_lp64.a $(MKL)/libmkl_intel_lp64.a  -Wl,--start-group $(MKL)/libmkl_intel_thread.a $(MKL)/libmkl_core.a  -Wl,--end-group  -liomp5 -lpthread
PARDISOLIBS = $(MKL_LIBS)

# HSL 2007 MA57 libraries
# http://hsl.rl.ac.uk/hsl2007/distrib/hsl2007.html
# needs registration as a registred researcher
# also (optionally) METIS
# http://glaros.dtc.umn.edu/gkhome/views/metis
METISLIB = /home/astivala/metis-4.0/libmetis.a
HSL2007DIR = /home/astivala/Desktop/HSL2007
MA57OBJS = $(HSL2007DIR)/ma57d.o $(HSL2007DIR)/ma57dsub.o
MA57LIBS = $(MA57OBJS)  $(METISLIB)


FC         = gfortran -std=gnu -cpp
FDEBUG     = -g -fbounds-check -O0 -pg
FOPTIMIZE  = -O3 -funroll-loops
ifeq ($(MODE),DEBUG)
    FFLAGS     = $(FDEBUG) -Wall
else
    FFLAGS     = $(FOPTIMIZE) -Wall
endif
CPPFLAGS   = -D_XOPEN_SOURCE=1 -D_XOPEN_SOURCE_EXTENDED=1 $(UMFPACKINCS)

CXX        = g++

CC         = gcc -std=gnu99
CDEBUG     = -g  -O0 -pg
COPTIMIZE  = -O3 
ifeq ($(MODE),DEBUG)
    CFLAGS     = $(CDEBUG) -pedantic -Wall
else
    CFLAGS     = $(COPTIMIZE)  -pedantic -Wall
endif
#              the following warnings are not implied by -Wall
CFLAGS     += -Wextra -Wfloat-equal  \
              -Wdeclaration-after-statement -Wundef -Wshadow \
              -Wpointer-arith -Wbad-function-cast -Wcast-qual -Wcast-align\
              -Wwrite-strings -Wmissing-prototypes \
              -Wmissing-declarations -Wunreachable-code

LD         = gfortran
LDFLAGS    = $(UMFPACKLDIRS)
ifeq ($(MODE),DEBUG)
    LDFLAGS    = -g  -ffortran-bounds-check  $(UMFPACKLDIRS)
    LDFLAGS    += -pg  # for profiler gprof
endif
LDLIBS     = -lm -lblas -llapack 

