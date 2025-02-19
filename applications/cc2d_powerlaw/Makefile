#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR=-I. -Iinclude

SRC= avs2d.f  coeff.f  gupwd.f    mgrout.f  util.f \
     bmul.f   dfkt.f   indat2d.f  mgstp.f   parq2d.f  vanca.f \
     bndry.f  error.f  init1.f    nsdef.f   rdparm.f  xmrout.f \
     cc2d.f   fpost.f  matmul.f   optcnl.f  supwdg.f  zvalue1.f \
     oseen2.f RestProlMatModif.f  stabil.f  anisoref.f

EXEC=cc2d_powerlaw-$(ID)

FEATLIB= umfpack2 feat2d sysutils
FEATLIB+= $(if $(BLASLIB),,lapack)

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src
vpath %.inc include

test: $(EXEC) run
