#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include "TMath.h"

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

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
  outfilename << "Flight2024_P_" << pressure <<" _T_" << temperature << "_.gas";      

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("CO2", 90.,"AR", 10.);
 
  // Set the field range to be covered by the gas table. 
  const int nFields = 20;
  const double emin =    100.;
  const double emax = 100000.;
  // Flag to request logarithmic spacing.
  const bool useLog = true;
  const double bmin=0;
  const double bmax=0;
  const int nBFields=1;
  gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 

  const int ncoll = 10;
  // Switch on debugging to print the Magboltz output.
  gas->EnableDebugging();
  // Run Magboltz to generate the gas table.
  gas->GenerateGasTable(ncoll);
  gas->DisableDebugging();
  // Save the table. 
  gas->WriteGasFile(outfilename.str());

  // app.Run(kTRUE);

}
