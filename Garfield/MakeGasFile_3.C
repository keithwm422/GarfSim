#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include "TMath.h"

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  // TApplication app("app", &argc, argv);
 
  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 273.15;
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("CO2", 95.,"AR", 5.);
 
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
  gas->WriteGasFile("keith_co2_95_AR_5_T273.gas");

  // app.Run(kTRUE);

}
