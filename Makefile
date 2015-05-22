# prerequisites: 
# DELPHES, an environment variable pointing to a DELPHES-3.1.2 release build
# ../util, containing utilities

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -- simple non-optimized compilation
CXXFLAGS      = -g -O0 -Wall -fPIC -pipe -Wuninitialized -O 
LD            = $(CXX)
LDFLAGS       = -g 
SOFLAGS       = -g -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) 
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) -lTMVA -lRooFitCore -lRooFit -lRooStats

EXTHEADERS    = -I../util -I$(DELPHES)
LIBPATH       = $(shell pwd)/lib

READER = anaH.o
ANA = plotHpt.o RooOneSidedProfileLL.o hstat.o
EXE =  runH.o runPlot.o

DICTFILES = ${ANA:.o=Dict.o}
DICTHEADERS = ${ANA:.o=Dict.h}

# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^ 

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $< 

.PHONY: prep all clean vars

# ================================================================================
all: vars prep lib bin
# -----------------------------------------------------------------------

# -- library
lib: $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) $(LIBPATH)/libDelphes.so $(LIBPATH)/libutil.so -o lib/libh0.so 

# -- binaries
bin: lib/libh0.so obj/runH.o  obj/runPlot.o
	$(LD) $(LDFLAGS) -o bin/runH $(GLIBS) obj/runH.o $(LIBPATH)/libDelphes.so $(LIBPATH)/libutil.so $(LIBPATH)/libh0.so 
	$(LD) $(LDFLAGS) -o bin/runPlot $(GLIBS) obj/runPlot.o  $(LIBPATH)/libDelphes.so $(LIBPATH)/libutil.so lib/libh0.so


# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../../util/lib/libutil.so && cd - 
	cd lib && ln -f -s $(DELPHES)/libDelphes.so && cd - 

clean:
	rm -f $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	rm -f $(DICTHEADERS) 
	rm -f bin/runH bin/runPlot
	rm -f lib/*

# -- ensure that the environment variable DELPHES is set
vars:
ifndef DELPHES
    $(error DELPHES is undefined, please set it to your local Delphes installation)
endif
