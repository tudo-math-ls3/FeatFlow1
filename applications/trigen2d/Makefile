#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

INCDIR=-I.

SRC= avstr.f bdrchk.f cref1.f cref2.f cref3.f cref4.f cref5.f cref6.f \
     owsc1.f parq2d.f rdparm.f trigen2d.f xcor.f

EXEC=trigen2d-$(ID)

FEATLIB= feat2d sysutils

include $(FEATFLOW)/Rules_apps.mk
-include Deps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src

test: $(EXEC) run




