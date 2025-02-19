#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR=-I.

SRC= avs2d.f  dfkt.f   indat2d.f  optcnl.f  projma.f  trsort.f \
     bmul.f   error.f  init1.f    prostp.f  util.f \
     bndry.f  fpost.f  mgrout.f   parq2d.f  rdparm.f  xmrout.f \
     coeff.f  gupwd.f  nsdef.f    pp2d.f    supwdg.f  zvalue1.f

EXEC= pp2d-$(ID)

FEATLIB= feat2d sysutils

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src

test: $(EXEC) run
