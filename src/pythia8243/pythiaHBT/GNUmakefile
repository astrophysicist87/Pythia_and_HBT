# GNUmakefile for pythiaHBT
# -------------------------

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl -O3

CXX := g++

RM				=	rm -f
O               =   .o
SRCDIR			=	src
INCDIR			=	include
OBJDIR			=	obj

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	main_BEeffects_OpenMP
endif

ifeq "$(MAIN2)" ""
MAIN2		=	main_BEeffects_HepMC
endif

MAINSRC		=	main_BEeffects_OpenMP.cpp
MAIN2SRC	=	main_BEeffects_HepMC.cpp

SRC			=	$(SRCDIR)/Arsenal.cpp \
				$(SRCDIR)/ParameterReader.cpp

INC			= 	$(INCDIR)/Arsenal.h \
				$(INCDIR)/ParameterReader.h

# -------------------------------------------------

#OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
#					$(basename $(SRC))))
OBJECTS				=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
						$(notdir $(basename $(SRC)))))

TARGET			=	$(MAIN)
TARGET2			=	$(MAIN2)
INSTPATH		=	..

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

#%.cpp:
#	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install target target2

all:		mkobjdir $(TARGET) $(TARGET2)

target:		mkobjdir $(TARGET)

target2:	mkobjdir $(TARGET2)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

#$(TARGET):	$(OBJECTS)
#		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
$(TARGET): $(MAINSRC) $(OBJECTS) $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< $(OBJECTS) -o $@ -std=c++11 -fopenmp -g $(CXX_COMMON)
#	@echo 'Inside new GNUmakefile'

$(TARGET2): $(MAIN2SRC) $(OBJECTS) $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< $(OBJECTS) -o $@ -std=c++11 -fopenmp -g \
				-I../../HepMC/build/include  -L../../HepMC/build/lib \
				-Wl,-rpath,../../HepMC/build/lib -lHepMC $(CXX_COMMON)

clean:		
		-rm -f $(OBJECTS)

distclean:	
		-rm -f $(TARGET)
		-rm -f $(TARGET2)
		-rm -f $(OBJECTS)

install:	$(TARGET) $(TARGET2)
		cp $(TARGET) $(INSTPATH)
		cp $(TARGET2) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: 				$(INCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: 		$(INCDIR)/ParameterReader.h \
									$(INCDIR)/Arsenal.h
main_BEeffects_OpenMP.cpp: 			$(INCDIR)/ParameterReader.h \
									$(INCDIR)/Arsenal.h
main_BEeffects_HepMC.cpp: 			$(INCDIR)/ParameterReader.h \
									$(INCDIR)/Arsenal.h
