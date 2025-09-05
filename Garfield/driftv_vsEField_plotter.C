#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1.h>
#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "SolidBox.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "ViewCell.hh"
#include "ViewField.hh"
#include "ViewMedium.hh"
#include "TrackSimple.hh"
#include "ViewDrift.hh"
#include "TrackHeed.hh"
#include "DriftLineRKF.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "ViewSignal.hh"
#include "Random.hh"
#include <TROOT.h>
#include <TRint.h>
#include "DCsim.hh"
#include <fstream>
#include "Plotting.hh"
#include "TMath.h"
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <omp.h>
#include <chrono>

using namespace Garfield;

int main(int argc, char * argv[]) {
  auto start = std::chrono::high_resolution_clock::now();

  bool realtimeplots = true;
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  std::stringstream wid;
  // TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
  const double zpos = invals[0]; // HELIX zposition slice to plot Efield in millimeters
  std::cout << "zpos slicing is [mm] : " << zpos << std::endl;
  std::stringstream efieldfilename;
  efieldfilename << "efieldprofile_zpos_" << zpos << ".root";      
  Garfield::plottingEngine.SetDefaultStyle();
  MediumMagboltz * gas = new MediumMagboltz();
  // Setup the gas
  const double pressure = 760.; //Torr
  const double temperature = 293.15; //K
  // Set the temperature [K] and pressure [Torr]
  gas->LoadGasFile("Flight2024_Bon1T_P_755.865_T_293.15_10Ar_90CO2_multiE.gas");
  // lets just print out the drift velocity to a file?
  ViewMedium mediumView;
  mediumView.SetMedium(gas);
  mediumView.EnableExport("driftv_vs_e.txt");
  mediumView.PlotElectronVelocity('e');
  std::stringstream fieldfilename;
  fieldfilename << "fieldprofile_v_vs_E.root";      
  TPad * mypad=mediumView.GetCanvas();
  mypad->SaveAs(fieldfilename.str().c_str());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  app->Run(kTRUE);
  return 0;
}
