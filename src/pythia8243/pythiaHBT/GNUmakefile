# GNUmakefile for pythiaHBT
# -------------------------

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl 

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

MAINSRC		=	main_BEeffects_OpenMP.cpp

SRC			=	$(SRCDIR)/Arsenal.cpp \
				$(SRCDIR)/ParameterReader.cpp

INC			= 	$(INCDIR)/Arsenal.h \
				$(INCDIR)/ParameterReader.h

# -------------------------------------------------

#OBJDIR			=	.
#SRCFILES 		= 	$(SRC) $(INC) GNUmakefile
OBJECTS			=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
					$(basename $(SRC))))
TARGET			=	$(MAIN)
INSTPATH		=	..

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

#%.cpp:
#	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install target

all:		mkobjdir $(TARGET)

target:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

#$(TARGET):	$(OBJECTS)
#		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)
$(TARGET): $(MAINSRC) $(OBJECTS) $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< $(OBJECTS) -o $@ -std=c++11 -fopenmp -g $(CXX_COMMON)
#	@echo 'Inside new GNUmakefile'

clean:		
		-rm -f $(OBJECTS)

distclean:	
		-rm -f $(TARGET)
		-rm -f $(OBJECTS)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/Arsenal.cpp: 				$(INCDIR)/Arsenal.h
$(SRCDIR)/ParameterReader.cpp: 		$(INCDIR)/ParameterReader.h \
									$(INCDIR)/Arsenal.h
main_BEeffects_OpenMP.cpp: 			$(INCDIR)/ParameterReader.h \
									$(INCDIR)/Arsenal.h
