###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : July 2008
#
# $Id: Makefile 3490 2010-03-17 22:04:12Z alexs $
#
# Makefile for building the test programs for FORTRAN implementation of
# QP tableau searching.
#
# The default "all" target is just the "tests" target.
# The "tests" target will make test harness programs and run regression tests
# The "pardiso_tests" target builds and runs the test programs for the
# implementation using the PARDISO solver.
#
# The "clean" target will remove object and executable files.
#
#
# This Makefile uses implicit rules to compile .c, .cpp and .f
# files into .o files, and some variables that may be specific to GNU make.
#
###############################################################################

include ../common.mk
-include ../local.mk


# location of the objects to test ie to link test programs with
OBJDIR = ../src
COMMONF77OBJS = $(OBJDIR)/drecip.o $(OBJDIR)/demvv.o $(OBJDIR)/dvmin.o \
                $(OBJDIR)/dvmax.o \
                $(OBJDIR)/tmatn.o $(OBJDIR)/tscorn.o $(OBJDIR)/rand55.o \
			    $(OBJDIR)/rdtabn.o
DENSEF77OBJS = $(OBJDIR)/dlphs1.o $(OBJDIR)/dlsoqp.o $(OBJDIR)/dlphs2.o 
SPARSEF77OBJS = $(OBJDIR)/spsoqp.o $(OBJDIR)/spphs1.o $(OBJDIR)/spphs2.o

UMFF77OBJS = $(OBJDIR)/umf4_f77wrapper.o

PARDISOF77OBJS = $(OBJDIR)/pasoqp.o $(OBJDIR)/paphs1.o $(OBJDIR)/paphs2.o

MA57F77OBJS = $(OBJDIR)/masoqp.o $(OBJDIR)/maphs1.o $(OBJDIR)/maphs2.o

TESTPROGS = testp1 testp2 dltest tmtst1 tmtst2 tmtst3 tmtst4 tmtst5 tmtst6 \
            testumfpack umftst tstsp1 tstsp2 sptest \
			tmtss1 tmtss2 tmtss3 tmtss4 tmtss5 tmtss6 
TESTOBJS = $(TESTPROGS:=.o)

PARDISO_TESTPROGS =  testpardiso pasptest patstsp1 patstsp2
PARDISO_TESTOBJS = $(PARDISO_TESTPROGS:=.o)

MA57_TESTPROGS = masptest matstsp1 matstsp2
MA57_TESTOBJS = $(MA57_TESTPROGS:=.o)

#all: tests pardiso_tests

all: tests ma57_tests

tests: $(TESTPROGS)
	@for i in $+ ; do \
		echo $$i ; \
		./regress.sh $$i $$i.correct.out ; \
	done

pardiso_tests: $(PARDISO_TESTPROGS)
	@for i in $+ ; do \
		echo $$i ; \
		./regress.sh $$i $$i.correct.out ; \
	done

ma57_tests: $(MA57_TESTPROGS)
	@for i in $+ ; do \
		echo $$i ; \
		./regress.sh $$i $$i.correct.out ; \
	done

umftst: umftst.o $(COMMONF77OBJS) $(UMFF77OBJS) $(SPARSEF77OBJS)
	$(LD) -o umftst  umftst.o $(COMMONF77OBJS) $(UMFF77OBJS)  $(SPARSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)


testumfpack: testumfpack.o $(COMMONF77OBJS) $(UMFF77OBJS) $(SPARSEF77OBJS)
	$(LD) -nofor-main -o testumfpack  testumfpack.o $(COMMONF77OBJS) $(UMFF77OBJS)  $(SPARSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)
#	$(LD) -Mnomain  -o testumfpack  testumfpack.o $(COMMONF77OBJS) $(UMFF77OBJS)  $(SPARSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)

testpardiso: testpardiso.o 
	$(LD) -o testpardiso  testpardiso.o $(PARDISOOBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)


sptest: sptest.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o sptest sptest.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)

tstsp1: tstsp1.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tstsp1  tstsp1.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 

tstsp2: tstsp2.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tstsp2  tstsp2.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 

tmtst6: tmtst6.o $(COMMONF77OBJS)  $(DENSEF77OBJS)
	$(LD) -o tmtst6  tmtst6.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tmtst5: tmtst5.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o tmtst5  tmtst5.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tmtst4: tmtst4.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o tmtst4  tmtst4.o $(COMMONF77OBJS) $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tmtst3: tmtst3.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o tmtst3  tmtst3.o $(COMMONF77OBJS) $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tmtst2: tmtst2.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o tmtst2  tmtst2.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

tmtst1: tmtst1.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o tmtst1  tmtst1.o $(COMMONF77OBJS) $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

dltest: dltest.o $(COMMONF77OBJS) $(DENSEF77OBJS)
	$(LD) -o dltest dltest.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

testp1: testp1.o $(COMMONF77OBJS)  $(DENSEF77OBJS)
	$(LD) -o testp1  testp1.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

testp2: testp2.o $(COMMONF77OBJS)  $(DENSEF77OBJS)
	$(LD) -o testp2  testp2.o $(COMMONF77OBJS)  $(DENSEF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)


tmtss6: tmtss6.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss6  tmtss6.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 

tmtss5: tmtss5.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss5  tmtss5.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)

tmtss4: tmtss4.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss4  tmtss4.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 

tmtss3: tmtss3.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss3  tmtss3.o $(COMMONF77OBJS) $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS)

tmtss2: tmtss2.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss2  tmtss2.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 

tmtss1: tmtss1.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS)
	$(LD) -o tmtss1  tmtss1.o $(COMMONF77OBJS)  $(SPARSEF77OBJS) $(UMFF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(UMFPACKLDIRS) $(UMFPACKLIBS) $(LDLIBS) 




pasptest: sptest.o $(COMMONF77OBJS) $(PARDISOF77OBJS)
	$(LD) -o pasptest sptest.o $(COMMONF77OBJS) $(PARDISOF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(PARDISOLIBS)

patstsp1: patstsp1.o $(COMMONF77OBJS) $(PARDISOF77OBJS)
	$(LD) -o patstsp1  patstsp1.o $(COMMONF77OBJS) $(PARDISOF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(PARDISOLIBS)

patstsp2: patstsp2.o $(COMMONF77OBJS) $(PARDISOF77OBJS)
	$(LD) -o patstsp2  patstsp2.o $(COMMONF77OBJS) $(PARDISOF77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(PARDISOLIBS)




masptest: sptest.o $(COMMONF77OBJS) $(MA57F77OBJS)
	$(LD) -o masptest sptest.o $(COMMONF77OBJS) $(MA57F77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(MA57LIBS)

matstsp1: matstsp1.o $(COMMONF77OBJS) $(MA57F77OBJS)
	$(LD) -o matstsp1  matstsp1.o $(COMMONF77OBJS) $(MA57F77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(MA57LIBS)

matstsp2: matstsp2.o $(COMMONF77OBJS) $(MA57F77OBJS)
	$(LD) -o matstsp2  matstsp2.o $(COMMONF77OBJS) $(MA57F77OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS) $(MA57LIBS)



clean:
	$(RM) $(TESTPROGS) $(TESTOBJS) $(PARDISO_TESTPROGS) $(PARDISO_TESTOBJS) $(MA57_TESTOBJS)



