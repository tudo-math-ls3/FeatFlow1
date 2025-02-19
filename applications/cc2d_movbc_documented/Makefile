#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

INCDIR=-I. -I./include

SRC1= cc2d.f dfkt.f mgstp.f nonstead.f nscalc.f nsdef.f optcnl.f \
      vanca.f zvalue1.f 
      
SRC2= bmul.f bndry.f coeff.f em30.f em31.f exxx.f gupwd.f ltrafo.f \
      matrest.f parq2domega.f parq2dwrap.f qtrafo.f supwdg.f ttrafo.f\
      xmrout.f
     
SRC3= arcalc.f bdrchk.f invert.f matmul.f util.f

SRC4= m011.f mgrout.f mgrout10.f mgrout30.f mgrout31.f

SRC5= avs2d.f bdintegral.f draglift.f draglift2.f error.f fpost.f\
      gmv2d.f
      
SRC6= alcgen.f datutil.f done.f init1.f rdini.f rdparm.f 
      
SRCUSER=fictbdry.f indat2d.f parq2d.f
      
SRC= $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5) $(SRC6) $(SRCUSER)

EXEC=cc2d-$(ID)

FEATLIB= feat2d sysutils

include $(FEATFLOW)/Rules_apps.mk

# The user defineable files and the include files are in the main
# directory. The source files for the black-box solver must be
# found in the "src"-subdirectory:

vpath %.f src src/discret src/init src/misc src/multigrid src/postproc
vpath %.inc include 

test: $(EXEC) run
