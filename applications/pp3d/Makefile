#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR= -I.

SRC= avs3d.f bmul.f bndry.f coeff.f dfkt.f error.f fpost.f \
  gupwd.f indat3d.f init1.f mgrout.f nsdef.f optcnl.f parq3d.f pp3d.f \
  projma.f prostp.f rdparm.f supwdg.f trsort.f util.f \
  xmrout.f zvalue1.f 


EXEC=pp3d-$(ID)

FEATLIB= feat3d feat2d sysutils 

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src

test: $(EXEC) run
