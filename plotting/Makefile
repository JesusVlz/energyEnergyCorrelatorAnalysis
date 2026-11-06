PROGRAM       = DrawingClasses
version       = EEC
CXX           = g++
ROOT_CONFIG   ?= root-config
# Use the dictionary generator from the same ROOT installation as its headers/libs.
ROOTCLING     = $(shell $(ROOT_CONFIG) --bindir)/rootcling
CXXFLAGS      = -O -Wall -fPIC -D$(version) # Add -g for debugger
CXXFLAGS     += $(shell $(ROOT_CONFIG) --cflags) $(INCS)
LDFLAGS       = $(shell $(ROOT_CONFIG) --libs)

HDRSDICT = SplitCanvas.h
HDRS += $(HDRSDICT) AlgorithmLibrary.h EECHistogramManager.h EECCard.h EECDrawer.h EECComparingDrawer.h EECBackgroundScale.h EECSignalToBackgroundUnfoldingScale.h EECUnfoldConfiguration.h SystematicUncertaintyOrganizer.h HybridModelHistogramManager.h HolguinHistogramManager.h CoLBTHistogramManager.h JewelHistogramManager.h

OBJS = $(HDRS:.h=.o) nanoDict.o
DEPS = $(OBJS:.o=.d)
ROOT_BUILD_CONFIG = .root-build-config

.PHONY: all $(PROGRAM) clean cl FORCE
.DELETE_ON_ERROR:

all: $(PROGRAM).so

# Keep the historical explicit target while tracking the file actually produced.
$(PROGRAM): $(PROGRAM).so

$(PROGRAM).so: $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(CXX) -shared $(CXXFLAGS) $(OBJS) $(LDFLAGS) -lEG -o $@
	@echo "done"

# Recheck the environment each time, but only change the stamp when it differs.
# cmsenv alone does not invalidate files made by a different ROOT/compiler.
$(ROOT_BUILD_CONFIG): FORCE
	@set -e; { \
	  printf '%s\n' "ROOT_CONFIG=$$(command -v $(ROOT_CONFIG))"; \
	  $(ROOT_CONFIG) --version --prefix --bindir; \
	  printf '%s\n' "ROOTCLING=$(ROOTCLING)" "CXX=$(CXX)" \
	    "CXX_PATH=$$(command -v $(firstword $(CXX)))" \
	    "CPPFLAGS=$(CPPFLAGS)" "CXXFLAGS=$(CXXFLAGS)" "LDFLAGS=$(LDFLAGS)"; \
	  $(CXX) --version; \
	} > $@.tmp
	@if ! cmp -s $@.tmp $@; then mv -f $@.tmp $@; else rm -f $@.tmp; fi

FORCE:

$(OBJS): $(ROOT_BUILD_CONFIG) Makefile

%.o: %.cxx %.h
	$(COMPILE.cc) -MMD -MP $(OUTPUT_OPTION) $<

nanoDict.o: nanoDict.cc nanoDict_rdict.pcm
	$(COMPILE.cc) -MMD -MP $(OUTPUT_OPTION) $<

# Both files come from one rootcling invocation (GNU make >= 4.3).
nanoDict.cc nanoDict_rdict.pcm &: $(HDRSDICT) LinkDef.h $(ROOT_BUILD_CONFIG) Makefile
	@echo "Generating dictionary with $(ROOTCLING) ..."
	$(ROOTCLING) -f nanoDict.cc -c -D$(version) $(CPPFLAGS) $(INCS) $(HDRSDICT) LinkDef.h

clean:
	rm -f $(OBJS) $(DEPS) $(PROGRAM).o $(PROGRAM).so \
	  nanoDict.cc nanoDict_rdict.pcm $(ROOT_BUILD_CONFIG) $(ROOT_BUILD_CONFIG).tmp

# Recursive makes order clean before the build even with make -j.
cl:
	$(MAKE) clean
	$(MAKE) all

-include $(DEPS)
