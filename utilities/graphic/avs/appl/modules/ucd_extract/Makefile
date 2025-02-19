# This module was obtained from the International AVS Center
# For more information, please send email to avsemail@ncsc.org

# Makefile for     : ucd_minmax
# Module Author    : Ian Curington
# Module Version # : 1.000
# Platform         : SUN 4/490
# Makefile Date    : 06/20/92

#AVS_PATH=/usr/avs
INC_FILE=$(AVS_PATH)/include/Makeinclude
MLIB=-lm
include $(INC_FILE)

AVS_LIBS=$(AVS_PATH)/lib
BASELIBS=-lgeom -lutil -lm
FBASELIBS=-lgeom -lutil
FLOWLIBS=-L$(AVS_LIBS) -lflow_c $(BASELIBS) $(LASTLIBS)
CSIMLIBS=-L$(AVS_LIBS) -lsim_c $(BASELIBS) $(LASTLIBS)
MODLIBS=-L$(AVS_LIBS) -lmdata -lmfilt -lmmapp -lmrend -lrf
AVS_INC=-I. -I$(AVS_PATH)/include
CFLAGS=$(ACFLAGS) $(AVS_INC)

# In order to use C++, you will need to modify the CPP and CPPLIBS definitions
# to match your local system conventions.

CPP_MAIN_LIB=-lC
CPPFLOWLIBS=$(FLOWLIBS) -lc $(CPP_MAIN_LIB) -lc
CPPSIMLIBS=$(CSIMLIBS) -lc $(CPP_MAIN_LIB) -lc
CPPHOME=/usr
CPPLIBS=-L$(CPPHOME)/lib
CPP=$(CPPHOME)/bin/CC
CPPFLAGS=$(CFLAGS) $(CPPLIBS)



DEFINITIONS=
OBJECTS=ucd_extract.o
LIBRARIES=

ucd_extract:	$(OBJECTS)
	cc $(CPPFLAGS) -o ucd_extract $(OBJECTS) $(FLOWLIBS) $(LIBRARIES)

.o:
	cc $(CPPFLAGS) $(DEFINITIONS) -c $<

.f.o:
	$(F77) $(F77FLAGS)  -c $<

clean:
	rm *.o
