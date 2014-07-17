ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = $(ROOTCFLAGS) -g -Wall -fPIC

LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

GLIBS         = $(filter-out -lz, $(ROOTGLIBS))
GLIBS         += -lMinuit -lRooFitCore -lRooFit -lRooStats

EXTHEADERS    = -I../util
EXTLIBS       = lib/libutil.so

# -- Default rules
$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^ 

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $< 

ANA = plotHpt.o # anaH.o 

DICTFILES = ${ANA:.o=Dict.o}
DICTHEADERS = ${ANA:.o=Dict.h}


# ================================================================================
all: dir links $(addprefix obj/,$(ANA) $(DICTFILES)) 
# --------------------------------------------------
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA) $(DICTFILES)) $(EXTLIBS) -o lib/libh0.so 


# -- create directories if not yet existing
dir:
	mkdir -p obj bin lib

# -- create links if not yet existing
links:
	cd lib && ln -f -s ../../util/lib/libutil.so && cd - 

# ======================================================================
# -- Executables
# ======================================================================

clean:
	rm -f $(addprefix obj/,$(ANA) $(DICTFILES)) 
	rm -f $(DICTHEADERS) 
	rm -f bin/runH
	rm -f lib/*
