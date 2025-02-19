#!/usr/bin/env gmake

FEATFLOW=../..
include $(FEATFLOW)/Globals.mk

all: bench
bench: sysinfo cc2d cc3d pp2d pp3d finish
	@(date >>$(BENCHLOG))

.NOTPARALLEL: cc2d pp2d cc3d pp3d finish

cc2d pp2d cc3d pp3d:
#
# We process one test with Streamline Diffusion and one test
# with Upwind.
# Some of the old source code files had errors that were fixed
# over the time.
# To be compatible to the old timing results, we replace
# some of the source code files at compile time 'old' versions.
# These 'old versions' can be found in the 'altsrc'
# subdirectory. The MAKE utility will automatically take these
# instead of the files in the application directory because
# of the VPATH setting. The ADDITIONALSRC variable contains
# filenames of 'old' files from the FEAT2D/3D libraries that are not
# present in the SRC= specifier in the makefiles of the applications.
# ADDITIONALSRC is added to SRC and therefore replaces 'new' FEAT
# library routines at compile time.
#
	@($(MAKE) -f ../$@/Makefile ADDITIONALSRC=sveb.f VPATH=./altsrc/$@:../$@/src clean)
	@echo | tee -a $(BENCHLOG)
	@echo "Processing benchmark:" $@ "with streamline diffusion ("$@"-SDG)" | tee -a $(BENCHLOG)
	@(cp \#data/$@.dat_sdg \#data/$@.dat)
	-@($(MAKE) -f ../$@/Makefile ADDITIONALSRC=sveb.f VPATH=./altsrc/$@:../$@/src run)
	@(cat \#data/$@.sdg >>$(BENCHLOG))
	@echo | tee -a $(BENCHLOG)
	@echo "Processing benchmark:" $@ "with upwinding ("$@"-UPW)" | tee -a $(BENCHLOG)
	@(cp \#data/$@.dat_upw \#data/$@.dat)
	-@($(MAKE) -f ../$@/Makefile ADDITIONALSRC=sveb.f VPATH=./altsrc/$@:../$@/src run)
	@(cat \#data/$@.upw >>$(BENCHLOG))

finish:
	@echo '============================================================'
	@echo "Benchmarks done."
	@echo 
	@echo "The benchmark protocol was saved in the file" 
	@echo 
	@echo '          '$(BENCHLOGFILENAME)
	@echo 
	@echo 'in the FeatFlow installation directory.'
	@echo 
	@echo 'Please, email the protocol file to featflow@featflow.de'
	@echo 'if you like your information to be published on our homepage.'
	@echo 

sysinfo: id
	@(echo 'Featflow ' $(FFVER) ' benchmark' >$(BENCHLOG))
	@(date >>$(BENCHLOG))
	@(uname -a   >>$(BENCHLOG))
	@($(MAKE) id >>$(BENCHLOG))
	-@($(FEATFLOW)/bin/info_cpu $(ID) >>$(BENCHLOG))
	-@($(FEATFLOW)/bin/info_f77 $(ID) $(FC) >>$(BENCHLOG))

clean clean_exec purge purge_all: 
	@($(MAKE) -f ../cc2d/Makefile VPATH=./altsrc/cc2d:../cc2d/src $@)
	@($(MAKE) -f ../cc3d/Makefile VPATH=./altsrc/cc3d:../cc3d/src $@)
	@($(MAKE) -f ../pp2d/Makefile VPATH=./altsrc/pp2d:../pp2d/src $@)
	@($(MAKE) -f ../pp3d/Makefile VPATH=./altsrc/pp3d:../pp3d/src $@)

id: .id

