#include <iostream>
#include <sstream>
#include <iomanip> // Required for setprecision
#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include "TMath.h"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/FundamentalConstants.hh"
#include <fstream>

using namespace Garfield;

int main(int argc, char * argv[]) {

   TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
  std::ofstream outFile("driftve.csv");
// 2. Check if the file opened successfully
    if (outFile.is_open()) {
        // 3. Write data to the file just like you use std::cout
        outFile << "B,theta,E,vE,vexB,vB,lorentz" << std::endl;

        
        std::cout << "Successfully wrote to the file." << std::endl;
    } else {
        // Handle potential errors (e.g., permissions or disk full)
        std::cerr << "Error: Could not open the file for writing." << std::endl;
        return 1;
    }
  double myvelocity;
  MediumMagboltz * gas = new MediumMagboltz();
  const double pressure = 14.616*51.7149; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr)
  gas->SetTemperature(299.0); // from CLI
  gas->SetPressure(pressure);
  gas->SetComposition("CO2", 90.,"AR", 10.);
  const int nFields = 5;
  const double E_not = 720.0;
  //const double E_not = 984.25;
  const double emin = E_not-E_not;
  const double emax = E_not+E_not;
  // Flag to request logarithmic spacing.
  const bool useLog = false;
  const double bmin=0;
  const double bmax=2.5; // do we need magnetic field on?
  const int nBFields=4;
  // angles
  const double amin=(5.0/6.0)*TMath::Pi()/2.0;
  const double amax=TMath::Pi()/2.0; // do we need magnetic field on?
  const int naFields=3;

  gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,amin,amax,naFields);
  // Turn on penning transfer?
  gas->EnablePenningTransfer();
  gas->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gas->GetNumberOfLevels();
  const int ncoll = 5;
  gas->GenerateGasTable(ncoll);
  // lets just print out the drift velocity to a file?
  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  //gas->PrintGas();
  std::vector<double> efields;
  std::vector<double> bfields;
  std::vector<double> angles;
  gas->GetFieldGrid(efields, bfields, angles);
  const auto nE = efields.size();
  const auto nB = bfields.size();
  const auto nA = angles.size();
      std::cout << "B,theta,E,vE,vexB,vB,lorentz\n";
  for (size_t j = 0; j < nB; ++j) {
    for (size_t k = 0; k < nA; ++k) {
      for (size_t i = 0; i < nE; ++i) {
        double ve = 0.;
        double vexB = 0.;
        double vB = 0.;
        gas->GetElectronVelocityE(i, j, k, ve);
        gas->GetElectronVelocityExB(i,j,k, vexB);
        gas->GetElectronVelocityB(i,j,k, vB);
        // Convert from cm/ns to cm/us.
        ve *= 1.e3;
        vexB *= 1.e3;
        vB *= 1.e3;
        double lorentzangle = 0.;
        gas->GetElectronLorentzAngle(i, j, k, lorentzangle);
        //alpha = exp(alpha);
        std::printf("%10.4f,%10.4f,%10.4f,%10.10f,%10.10f,%10.10f,%10.10f\n", bfields[j], angles[k] * RadToDegree, efields[i], ve, vexB, vB, lorentzangle);
        outFile << std::fixed << std::setprecision(8) << bfields[j] << "," << angles[k] * RadToDegree << "," << efields[i] << "," << ve << "," << vexB << "," << vB << "," << lorentzangle << std::endl;
      }
    }
  }
  //ViewMedium view;
  //view.SetMedium(gas);
  //view.PlotElectronVelocity('e'); 
        // 4. Close the file to free up resources
        outFile.close();
  
   //app.Run(kTRUE);

}
