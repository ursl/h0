# required ENV variables: 
# DELPHES

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#CXXFLAGS      = -g -O3 -Wall -fPIC -pipe -Wuninitialized
CXXFLAGS      = -g -O0 -Wall -fPIC -pipe -Wuninitialized -O 
LD            = $(CXX)
LDFLAGS       = -g 
SOFLAGS       = -g -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) 
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) -lTMVA -lRooFitCore -lRooFit -lRooStats

EXTHEADERS    = -I../util
LIBPATH       = `pwd`/lib

READER = anaH.o ExRootTreeReader.o runH.o
ANA = plotHpt.o 

DICTFILES = ${ANA:.o=Dict.o}
DICTHEADERS = ${ANA:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^ 

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $< 


# ================================================================================
all: prep lib bin
# -----------------------------------------------------------------------

# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../../util/lib/libutil.so && cd - 
	cd delphes && ln -f -s  $(DELPHES)/classes && cd - 
	cd lib && ln -f -s $(DELPHES)/libDelphes.so && cd - 

# -- library
lib: $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) $(EXTLIBS) -o lib/libh0.so 

# -- binaries
bin: lib/libh0.so obj/runH.o 
	$(LD) $(LDFLAGS) -o bin/runH $(GLIBS) obj/runH.o -L$(LIBPATH) -lh0 -lDelphes -lutil

# -- other stuff
obj/ExRootTreeReader.o: 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c delphes/ExRootTreeReader.cc -o $@


clean:
	rm -f $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	rm -f $(DICTHEADERS) 
	rm -f delphes/classes
	rm -f bin/runH
	rm -f lib/*
