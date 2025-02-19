#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR=-I.

SRC= avs3d.f bmul.f bndry.f cc3d.f coeff.f dfkt.f error.f \
     fpost.f gupwd.f indat3d.f init1.f matmul.f mgrout.f mgstp.f \
     nsdef.f optcnl.f parq3d.f rdparm.f supwdg.f util.f vanca.f \
     xmrout.f zvalue1.f

EXEC=cc3d-$(ID)

FEATLIB= feat3d feat2d sysutils 

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src

test: $(EXEC) run
