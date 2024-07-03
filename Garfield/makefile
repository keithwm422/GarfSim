OBJDIR  = $(GARFIELD_HOME)/build/Object
SRCDIR  = $(GARFIELD_HOME)/Source
INCDIR  = $(GARFIELD_HOME)/install/include/Garfield
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR  = $(GARFIELD_HOME)/install/lib64
DCsimINCDIR  = $(DCsim_HOME)/Include
DCsimLIBDIR  = $(DCsim_HOME)/Library

# Compiler (and Linker)
CXX = $(shell root-config --cxx)

# Compiler (and Linker) flags
CXXFLAGS  = $(shell root-config --cflags)
CXXFLAGS += -I$(INCDIR) -I$(DCsimINCDIR) -I$(HEEDDIR)
# Code Generation and Optimization flags
CXXFLAGS += -fno-common -fPIC -O3
# Warning flags
CXXFLAGS += -W -Wall -Wextra -Wno-long-long
# Debug flags
CXXFLAGS += -g

# Libraries (note: the order matters)
LDLIBS  = -L$(LIBDIR) -Wl,--whole-archive -lGarfield  -L$(DCsimLIBDIR) -l DCsim -Wl,--no-whole-archive
# Fortran libraries
LDLIBS += -lgfortran
# ROOT libraries
LDLIBS += $(shell root-config --glibs)
LDLIBS +=  -lGeom

garfroot: garfroot.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile: MakeGasFile.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_2: MakeGasFile_2.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_3: MakeGasFile_3.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_flight: MakeGasFile_flight.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
jet: jet.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet: smalljet.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
gasgainstudy: gas_gain_study.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_fieldplots: smalljet_fieldplots.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
printgas: printgas.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
full_tpc_sim:	full_tpc_sim.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
filterjob:  filterjob.cc
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
