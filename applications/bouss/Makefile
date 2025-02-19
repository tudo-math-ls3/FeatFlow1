#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR=-I. -Iinclude

SRC= avs2d.f  dfkt.f   indat2d.f  optcnl.f  projma.f  trsort.f \
     bmul.f   error.f  init1.f    prostp.f  util.f \
     bndry.f  fpost.f  mgrout.f   parq2d.f  rdparm.f  xmrout.f \
     coeff.f  gupwd.f  nsdef.f    bouss.f   supwdg.f  zvalue1.f \
     hwahl.f  mgcort.f oralB.f    tcalc.f   tdef.f    tdfkt.f \
     trand.f

EXEC= bouss-$(ID)

FEATLIB= feat2d sysutils

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src
vpath %.inc include

test: $(EXEC) run
