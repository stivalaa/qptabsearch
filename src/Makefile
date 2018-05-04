###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : July 2008
#
# $Id: Makefile 2094 2009-03-10 02:22:50Z astivala $
#
# Makefile for building the the C++, C,  FORTRAN-77 code for
# QP tableau searching.
#
#
# The default target "all" will compile and link the programs.
# Note there are two programs made from the same main program
# module tsrchn.f, one using sparse matrix implementation of
# QP solver and one using dense matrix implementation of
# QP solver; they just link different modules for the SOLVQP 
# subroutine.
#
# The third uses the PARDISO library included in the Intel Math Kernel
# library (MKL) - it is not included in all by default, the 'pardiso'
# target will build it also.
#
# set MODE=DEBUG to build with debugging and profiling on. Otherwise
# default is to build with optimizations on and no debug or profile.
#
# The "check" target will run ftnchek to do static checking of FORTRAN code
# (get ftnchek from http://www.dsm.fordham.edu/~ftnchek/
#  or http://netlib.org/fortran/).
#
# The "clean" target will remove object and executable files.
#
# This Makefile uses implicit rules to compile .c, .cpp and .f
# files into .o files, and some variables that may be specific to GNU make.
#
###############################################################################

include ../common.mk
-include ../local.mk

# F77 source in all programs
COMMONF77SRCS = drecip.f demvv.f dvmin.f dvmax.f rand55.f rdistm.f

# F77 source in only the numeric (Omega matrix) versions
NUMERICF77SRCS = tmatn.f tscorn.f rdtabn.f tsrchn.f

# F77 source in only the discrete (Tableaux) versions
DISCRETEF77SRCS = tmatd.f tscord.f rdtabd.f tsrchd.f

# F77 source in only the dense matrix versions
DENSEF77SRCS = dlphs1.f dlsoqp.f dlphs2.f 

# F77 source in only the sparse matrix UMFPACK versions
SPARSEF77SRCS = spsoqp.f spphs1.f spphs2.f

# F77 source in only the sparse matrix PARDISO versions
PARDISOF77SRCS = pasoqp.f paphs1.f paphs2.f

F77SRCS = $(COMMONF77SRCS) $(DENSEF77SRCS) $(SPARSEF77SRCS) \
          $(NUMERICF77SRCS) $(DISCRETEF77SRCS) $(PARDISOF77SRCS)
COMMONF77OBJS = $(COMMONF77SRCS:.f=.o)
SPARSEF77OBJS = $(SPARSEF77SRCS:.f=.o)
PARDISOF77OBJS = $(PARDISOF77SRCS:.f=.o)
DENSEF77OBJS = $(DENSEF77SRCS:.f=.o)
NUMERICF77OBJS = $(NUMERICF77SRCS:.f=.o)
DISCRETEF77OBJS = $(DISCRETEF77SRCS:.f=.o)
UMFF77OBJS = umf4_f77wrapper.o
CPPSRCS = 
CPPOBJS = $(CPPSRCS:.cpp=.o)
SPARSECSRCS =  umf4_f77wrapper.c
SPARSECOBJS = $(CSRCS:.c=.o)
SRCS = $(F77SRCS) $(CPPSRCS) $(CSRCS)
OBJS = $(COMMONF77OBJS) $(DENSEF77OBJS) $(SPARSEF77OBJS) $(CPPOBJS) \
       $(COBJS) $(UMFF77OBJS) \
       $(NUMERICF77OBJS) $(DISCRETEF77OBJS) $(PARDISOF77OBJS)
MAKEDEPEND = gcc -MM $(CPPFLAGS)
DEPENDFILE = .depend

# Program to build TAGS file for EMACS
MAKETAGS   = etags


#all: tsrchn_sparse tsrchn_dense tsrchd_sparse tsrchd_dense tsrchd_pardiso tsrchn_pardiso

all: tsrchn_sparse tsrchn_dense tsrchd_sparse tsrchd_dense

pardiso:  tsrchd_pardiso tsrchn_pardiso

tsrchn_sparse: $(COMMONF77OBJS) $(NUMERICF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tsrchn_sparse $(COMMONF77OBJS) $(NUMERICF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) \
        $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tsrchn_pardiso: $(COMMONF77OBJS) $(NUMERICF77OBJS) $(PARDISOF77OBJS)
	$(LD) -o tsrchn_pardiso $(COMMONF77OBJS) $(NUMERICF77OBJS) $(PARDISOF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(PARDISOLIBS)

tsrchn_dense: $(COMMONF77OBJS) $(NUMERICF77OBJS) $(DENSEF77OBJS) 
	$(LD) -o tsrchn_dense $(COMMONF77OBJS) $(NUMERICF77OBJS) $(DENSEF77OBJS) \
        $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tsrchd_sparse: $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tsrchd_sparse $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) \
        $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tsrchd_pardiso: $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(PARDISOF77OBJS)
	$(LD) -o tsrchd_pardiso $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(PARDISOF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(PARDISOLIBS)

tsrchd_dense: $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(DENSEF77OBJS) 
	$(LD) -o tsrchd_dense $(COMMONF77OBJS) $(DISCRETEF77OBJS) $(DENSEF77OBJS) \
        $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

check: $(F77SRCS)
	ftnchek -nonovice $+

depend: $(CSRCS) $(CPPSRCS) $(HDRS)
	$(MAKEDEPEND) $(CSRCS) $(CPPSRCS) $(HDRS) > $(DEPENDFILE)

# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) tsrchn_sparse tsrchn_dense
	$(RM) $(OBJS)

include $(DEPENDFILE)


