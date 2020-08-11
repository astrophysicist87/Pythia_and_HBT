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
CFLAGS= -O3 -fopenmp -std=c++0x -lgsl -lgslcblas -g

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)
SRCDIR		=	.

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	run.e
endif

SRC		=	test.cpp \
			fit_v1.cpp \
			fit_v2.cpp

INC		= 	fit_v1.h \
			fit_v2.h

# -------------------------------------------------

OBJDIR		=	.
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	..


# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi


# -------------------------------------------------


.PHONY:		all mkobjdir clean distclean install


all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm $(OBJECTS)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

fit_v1.cpp: fit_v1.h
fit_v2.cpp: fit_v2.h
test.cpp: fit_v1.h fit_v2.h
