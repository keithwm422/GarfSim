#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include "TMath.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  // TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
//  const double pressure = 1 * AtmosphericPressure; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr) 
    MediumMagboltz master, sparse;
    master.LoadGasFile("Master_Dense_Drift.gas");
    sparse.LoadGasFile("Comparison_Sparse.gas");

    double E_nominal = 740.0; // Your operating field
    double vx1, vy1, vz1, vx2, vy2, vz2;
    double dummy = 0;

    // Get velocities (vz is usually the drift direction)
    master.ElectronVelocity(E_nominal, 0, 0, 0, 0, 0, vx1, vy1, vz1);
    sparse.ElectronVelocity(E_nominal, 0, 0, 0, 0, 0, vx2, vy2, vz2);

    double diff = TMath::Abs(vx1 - vx2); // in cm/us
    double drift_dist = 1.0; // assume a 1 cm drift cell
    double time_error = drift_dist * (1/vx1 - 1/vx2); // in microseconds
    double spatial_error_microns = TMath::Abs(time_error * vx1) * 10000;

    printf("Master Velocity: %.6f cm/ns\n", vx1);
    printf("Sparse Velocity: %.6f cm/ns\n", vx2);
    printf("Expected Systematic Shift: %.6f microns\n", spatial_error_microns);
}
