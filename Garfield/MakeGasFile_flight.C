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
  //outfilename << "Flight2024_Boff_P_" << pressure <<"_T_" << temperature << "_.gas";
  outfilename << "Flight2024_midway_Boff_P_" << pressure <<"_T_" << temperature << "_multiE_90CO2_10Ar_01232024_v2.gas";

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("CO2", 90.,"AR", 10.);

  // Set the field range to be covered by the gas table. 
  //const int nFields = 1;
  //const double E_not = 984.25;
  //const double emin = E_not;
  //const double emax = E_not;
  //for efield drift study
  const int nFields = 11;
  //make E_not the midpt
  const double E_not = 984.25;
  const double emin = E_not-E_not;
  const double emax = E_not+E_not;

  // Flag to request logarithmic spacing.
  const bool useLog = false;
  const double bmin=0;
  const double bmax=0; // do we need magnetic field on?
  const int nBFields=1;
  gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 

  // Turn on penning transfer?
  gas->EnablePenningTransfer();
  gas->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gas->GetNumberOfLevels();
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
