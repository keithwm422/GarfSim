OBJDIR  = $(GARFIELD_HOME)/build/Object
SRCDIR  = $(GARFIELD_HOME)/Source
INCDIR  = $(GARFIELD_HOME)/install/include
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
CXXFLAGS += -fno-common -fPIC -O3 -fopenmp
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
Csv_lorentz_vs_temp: Csv_lorentz_vs_temp.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_driftv_vs_bfield: Csv_driftv_vs_bfield.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_drift_vs_pressure: Csv_drift_vs_pressure.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_drift_vs_pressure_zoomed: Csv_drift_vs_pressure_zoomed.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_driftv_vs_efield: Csv_driftv_vs_efield.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_drift_vs_Temp: Csv_drift_vs_Temp.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
Csv_drift_vs_gascomp: Csv_drift_vs_gascomp.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
LoadGasFile: LoadGasFile.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile: MakeGasFile.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_2: MakeGasFile_2.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_3: MakeGasFile_3.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_flight: MakeGasFile_flight.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFile_multiB: MakeGasFile_multiB.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFineGridTemperature: MakeGasFineGridTemperature.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
MakeGasFineGridTemperature_Bfieldoff: MakeGasFineGridTemperature_Bfieldoff.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
jet: jet.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet: smalljet.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_drifttime: smalljet_drifttime.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_check_num_e: smalljet_check_num_e.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_drifttime_electron_only: smalljet_drifttime_electron_only.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_drifttime_electron_vs_coordinates: smalljet_drifttime_electron_vs_coordinates.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_drifttime_electron_vs_coordinates_finegrid: smalljet_drifttime_electron_vs_coordinates_finegrid.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
isochron_generator: isochron_generator.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
IsoHEATB: IsoHEATB.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
IsoMRE: IsoMRE.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
field_map_gen: field_map_gen.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
eField_plotter: eField_plotter.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
eField_plotter_v2: eField_plotter_v2.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
VoltageMap_ICRC2025: VoltageMap_ICRC2025.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
eField_plotter_forkenichi: eField_plotter_forkenichi.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
eField_plotter_BESS: eField_plotter_BESS.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
driftField_plotter_forkenichi: driftField_plotter_forkenichi.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
track_resolution_calculator: track_resolution_calculator.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
track_energy_deltaElectrons: track_energy_deltaElectrons.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
track_plotter: track_plotter.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
track_numelectrons: track_numelectrons.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
driftv_vsEField_plotter: driftv_vsEField_plotter.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
isochron_gen_p: isochron_gen_p.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
GasFileGenerator: GasFileGenerator.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
thread_TGraph2D_tester: thread_TGraph2D_tester.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
isochron_gen_p_full: isochron_gen_p_full.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
ViewIsochroneDCT: ViewIsochroneDCT.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
ViewBFieldExample: ViewBFieldExample.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
HijackedIsochroneDCT: HijackedIsochroneDCT.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
HijackedIsochroneDCT_allwires: HijackedIsochroneDCT_allwires.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
TScaleIso_run: TScaleIso_run.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
HijackedIsochroneDCT_topwires: HijackedIsochroneDCT_topwires.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
HijackedIsochroneDCT_testBmap: HijackedIsochroneDCT_testBmap.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
avalanch_gain_calc: avalanch_gain_calc.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
gain_calculator: gain_calculator.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
avalanch_gain_calc_v2: avalanch_gain_calc_v2.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
avalanch_drifttime_vs_temp: avalanch_drifttime_vs_temp.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
avalanch_drifttime_vs_temp_v2: avalanch_drifttime_vs_temp_v2.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
basic_drift_time_calculator: basic_drift_time_calculator.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
basic_drift_time_calculator_MC: basic_drift_time_calculator_MC.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
plotDriftVelocitye: plotDriftVelocitye.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
gasgainstudy: gas_gain_study.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
smalljet_fieldplots: smalljet_fieldplots.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
printgas: printgas.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
full_tpc_sim:	full_tpc_sim.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
make_potential_map_original: make_potential_map_original.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
make_potential_map: make_potential_map.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
make_potential_map_wires_for_cathode: make_potential_map_wires_for_cathode.C
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
filterjob:  filterjob.cc
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDLIBS)
