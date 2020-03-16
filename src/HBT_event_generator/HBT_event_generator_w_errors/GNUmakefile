# ============================================================================
#  Makefile HBTeg                             Chris Plumberg, October 31, 2018
# ============================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
#running with icpc
#CFLAGS= $(shell ../../crank/configure_compiler_C++_Flags.sh)
#running with g++
CFLAGS= -O3 -fopenmp -std=c++11 -g
#CFLAGS += -Wall
#CFLAGS += -lgsl -lgslcblas -g -Wall -Warray-bounds -Wchar-subscripts -Wcomment \
#						-Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized \
#						-Wmain -Wnarrowing -Wnonnull -Wparentheses \
#						-Wreorder -Wreturn-type -Wsequence-point -Wtrigraphs \
#						-Wunused-function -Wunused-but-set-variable -Wunused-variable \
#						-Wwrite-strings
#CFLAGS += -Wsign-compare -Wpointer-sign
#CFLAGS += -lgsl -lgslcblas -g
###CFLAGS += -lgsl -lgslcblas -g -profile-functions -profile-loops=all -profile-loops-report=2

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
LDFLAGS			+=		-L/home/blixen/plumberg/lib -lshifter
SYSTEMFILES     =       $(SRCGNU)
SRCDIR		=	src

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	run_HBT_event_generator.e
endif

ifeq "$(CHECKS)" ""
CHECKS		=	run_checks.e
endif

SRC			=	main.cpp \
				$(SRCDIR)/Arsenal.cpp \
				$(SRCDIR)/ParameterReader.cpp \
				$(SRCDIR)/gauss_quadrature.cpp \
				$(SRCDIR)/NumeratorAndDenominator.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_rotInv.cpp \
				$(SRCDIR)/NumeratorAndDenominator_wBA.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_wBA.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_rotInv_wBA.cpp \
				$(SRCDIR)/HBT_event_generator.cpp

CHECKSRC	=	run_checks.cpp \
				$(SRCDIR)/Arsenal.cpp \
				$(SRCDIR)/ParameterReader.cpp \
				$(SRCDIR)/gauss_quadrature.cpp \
				$(SRCDIR)/NumeratorAndDenominator.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_rotInv.cpp \
				$(SRCDIR)/NumeratorAndDenominator_wBA.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_wBA.cpp \
				$(SRCDIR)/NumeratorAndDenominator1D_rotInv_wBA.cpp \
				$(SRCDIR)/HBT_event_generator.cpp

INC			= 	main.h \
				$(SRCDIR)/EventRecord.h \
				$(SRCDIR)/ParticleRecord.h \
				$(SRCDIR)/Stopwatch.h \
				$(SRCDIR)/Arsenal.h \
				$(SRCDIR)/ParameterReader.h \
				$(SRCDIR)/gauss_quadrature.h \
				$(SRCDIR)/ensemble.h \
				$(SRCDIR)/random_events.h \
				$(SRCDIR)/estimate_error.h \
				$(SRCDIR)/HBT_event_generator.h

CHECKINC 	= 	run_checks.h \
				$(SRCDIR)/EventRecord.h \
				$(SRCDIR)/ParticleRecord.h \
				$(SRCDIR)/Stopwatch.h \
				$(SRCDIR)/Arsenal.h \
				$(SRCDIR)/ParameterReader.h \
				$(SRCDIR)/gauss_quadrature.h \
				$(SRCDIR)/ensemble.h \
				$(SRCDIR)/random_events.h \
				$(SRCDIR)/estimate_error.h \
				$(SRCDIR)/HBT_event_generator.h

# -------------------------------------------------

OBJDIR			=	.
#SRCFILES 		= 	$(SRC) $(INC) GNUmakefile
#CHECKSRCFILES 	= 	$(CHECKSRC) $(CHECKINC) GNUmakefile
OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
					$(basename $(SRC))))
CHECKOBJS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
					$(basename $(CHECKSRC))))
TARGET			=	$(MAIN)
INSTPATH		=	..
#INSTPATH		=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install target checks

all:		mkobjdir $(TARGET) $(CHECKS)

target:		mkobjdir $(TARGET)

checks:		mkobjdir $(CHECKS)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

$(CHECKS):	$(CHECKOBJS)	
		$(CC) $(CHECKOBJS) -o $(CHECKS) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm -f $(OBJECTS)
		-rm -f $(CHECKOBJS)

distclean:	
		-rm -f $(TARGET)
		-rm -f $(CHECKS)
		-rm -f $(OBJECTS)
		-rm -f $(CHECKOBJS)
#		-rm $(MAIN)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: $(SRCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: $(SRCDIR)/ParameterReader.h $(SRCDIR)/Arsenal.h
$(SRCDIR)/HBT_event_generator.cpp: $(SRCDIR)/HBT_event_generator.h $(SRCDIR)/Arsenal.h \
									$(SRCDIR)/Stopwatch.h $(SRCDIR)/EventRecord.h \
									$(SRCDIR)/ParticleRecord.h $(SRCDIR)/ensemble.h
./main.cpp: $(SRCDIR)/Stopwatch.h $(SRCDIR)/ParameterReader.h main.h $(SRCDIR)/EventRecord.h \
			$(SRCDIR)/ParticleRecord.h $(SRCDIR)/ensemble.h
./run_checks.cpp: $(SRCDIR)/Stopwatch.h $(SRCDIR)/ParameterReader.h run_checks.h $(SRCDIR)/EventRecord.h \
			$(SRCDIR)/ParticleRecord.h $(SRCDIR)/ensemble.h

