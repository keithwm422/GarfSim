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
  const double pressure = invals[0]*51.7149; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr) 
  std::cout << "pressure in torr: " << pressure << std::endl;
  const double temperature = 273.15 + invals[1];
  std::cout << "temperature in kelvin: " << temperature << std::endl;

  std::stringstream outfilename;
  outfilename << "DriftStudy_Boff_P_" << pressure <<"_T_" << temperature << ".gas";
  //outfilename << "Flight2024_Bon_P_" << pressure <<"_T_" << temperature << "_10Ar_90CO2_multiE.gas";

  MediumMagboltz *gasMaster = new MediumMagboltz();
  gasMaster->SetComposition("CO2", 90., "AR", 10.);
  gasMaster->SetTemperature(temperature);
  gasMaster->SetPressure(pressure); // Use your calculated pressure
  // DENSE GRID: Focus on the 100 - 5000 V/cm range
  // This ensures the spline has enough "anchors" to be perfectly smooth.
  // from IsoHEATB code
  /*
  const int nFields = 5;
  const double E_not = 984.25;
  const double emin = E_not-E_not;
  const double emax = E_not+E_not;
  // Flag to request logarithmic spacing.
  const bool useLog = false;
  const double bmin=0;
  const double bmax=2.5; // do we need magnetic field on?
  const int nBFields=4;
  const double amin=0;
  const double amax=TMath::Pi()/2.0; // do we need magnetic field on?
  const int nAFields=4;

  gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,amin,amax,nAFields);
  // Turn on penning transfer?
  gas->EnablePenningTransfer();
  gas->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gas->GetNumberOfLevels();
  const int ncoll = 5;
  gas->GenerateGasTable(ncoll);
  */ 


  const int nFields = 30; 
  const double emin = 100.0;
  const double emax = 5000.0;
  const bool useLogE = false; // Linear is better for narrow, high-precision ranges
  const double bmin=0;
  const double bmax=0; // do we need magnetic field on?
  const int nBFields=1;
  //gasMaster->SetFieldGrid(emin, emax, nFields, useLogE);
  gasMaster->SetFieldGrid(emin, emax, nFields, useLogE, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 
  // High statistics for the "Truth" file
  const int ncoll = 15;
  gasMaster->EnablePenningTransfer();
  gasMaster->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gasMaster->GetNumberOfLevels();
  //const int ncoll = 10;
  gasMaster->GenerateGasTable(ncoll);
  gasMaster->WriteGasFile("Master_Dense_Drift_isoheatb.gas");

  MediumMagboltz *gasSparse = new MediumMagboltz();
  gasSparse->SetComposition("CO2", 90., "AR", 10.);
  gasSparse->SetTemperature(temperature);
  gasSparse->SetPressure(pressure);

  // SPARSE GRID: Mimics the "5 points up to 200k" logic
  const int nFieldsSparse = 5; 
  const double eminSparse = 0;
  const double emaxSparse = 2.0*984.25; 
  const bool useLogESparse = false;
  gasSparse->SetFieldGrid(eminSparse, emaxSparse, nFieldsSparse, useLogESparse);
  gasSparse->SetFieldGrid(eminSparse, emaxSparse, nFieldsSparse, useLogESparse, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 

  // Use lower ncoll to see how statistical noise affects the curve
  const int ncollSparse = 5; 
  gasSparse->EnablePenningTransfer();
  gasSparse->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gasSparse->GetNumberOfLevels();
  gasSparse->GenerateGasTable(ncollSparse);
  gasSparse->WriteGasFile("Comparison_Sparse_isoheatb.gas");


    //master.LoadGasFile("Master_Dense_Drift.gas");
    //sparse.LoadGasFile("Comparison_Sparse.gas");

    double E_nominal = 740.0; // Your operating field
    double vx1, vy1, vz1, vx2, vy2, vz2;
    double dummy = 0;

    // Get velocities (vz is usually the drift direction)
    gasMaster->ElectronVelocity(E_nominal, 0, 0, 0, 0, 0, vx1, vy1, vz1);
    gasSparse->ElectronVelocity(E_nominal, 0, 0, 0, 0, 0, vx2, vy2, vz2);

    double diff = TMath::Abs(vx1 - vx2); // in cm/us
    double drift_dist = 1.0; // assume a 1 cm drift cell
    double time_error = drift_dist * (1/vx1 - 1/vx2); // in microseconds
    double spatial_error_microns = TMath::Abs(time_error * vx1) * 10000;

    printf("Master Velocity: %.6f cm/us\n", vx1);
    printf("Sparse Velocity: %.6f cm/us\n", vx2);
    printf("Expected Systematic Shift: %.6f microns\n", spatial_error_microns);

}
