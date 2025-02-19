#
# Makefile for create_mpeg
#

#AVS_PATH=/usr/avs
INC_FILE=$(AVS_PATH)/include/Makeinclude
include $(INC_FILE)

AVS_LIBS = $(AVS_PATH)/lib
INC = $(AVS_PATH)/include
BASELIBS=-lgeom -lutil -lm $(LASTLIBS)
CFLAGS = -I$(INC) $(AOPTCFLAGS) $(LOCAL_CFLAGS) -DAVS
CFLOWLIBS=-L$(AVS_LIBS) -lflow_c -L/usr/lib/X11R5 $(BASELIBS)
CSIMLIBS=-L$(AVS_LIBS) -lsim_c -L/usr/lib/X11R5 $(BASELIBS)
CC=cc

MODULES=create_mpeg

create_mpeg:	create_mpeg.o
	$(CC) $(CFLAGS) -o create_mpeg create_mpeg.o $(CFLOWLIBS)

create_mpeg.o:	create_mpeg.c
	$(CC) $(CFLAGS) -c create_mpeg.c


