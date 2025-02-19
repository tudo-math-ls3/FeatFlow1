#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

INCDIR=-I.

SRC= parq2d.f rdparm.f tr2to3.f util.f

EXEC=tr2to3-$(ID)

FEATLIB= feat2d sysutils

include $(FEATFLOW)/Rules_apps.mk
-include Deps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src

test: $(EXEC) run




