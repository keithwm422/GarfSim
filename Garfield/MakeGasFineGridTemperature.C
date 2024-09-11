#include <iostream>
#include <sstream>
#include <fstream>

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
  const double temperature_lower = 273.15 + invals[1];
  const double temperature_upper = 273.15 + invals[2];
  // make 100 data points?
  const double temp_step=(temperature_upper-temperature_lower)/100.0;

  double mylorentz[101] =   {0};
  double     myvel[101] =   {0};
  double    mytemp[101] =   {0};
  std::ofstream outfile;
  outfile.open("lorentz_vs_temp.csv");
  outfile << "T,la" << std::endl;
  std::ofstream outfile2;
  outfile2.open("driftv_vs_temp_fine.csv");
  outfile2 << "T,v" << std::endl;


  double temperature_i=temperature_lower;
  int gasfile_i=0;
  while(gasfile_i<101){
    // Setup the gas.
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetTemperature(temperature_i);
    mytemp[gasfile_i]=temperature_i;
    gas->SetPressure(pressure);
    gas->SetComposition("CO2", 90.,"AR", 10.);
    // Set the field range to be covered by the gas table. 
    const int nFields = 1;
    const double E_not = 984.25;
    const double emin = E_not;
    const double emax = E_not;
    const bool useLog = false;
    const double bmin=1;
    const double bmax=1; // do we need magnetic field on?
    const int nBFields=1;
    gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 
    const int ncoll = 10;
    // Switch on debugging to print the Magboltz output.
    gas->EnableDebugging();
    // Run Magboltz to generate the gas table.
    gas->GenerateGasTable(ncoll);
    gas->DisableDebugging();
    // now parse the gas for the thing we want and put into csv in Temp,la
    //gas->PrintGas();
    std::vector<double> efields;
    std::vector<double> bfields;
    std::vector<double> angles;
    gas->GetFieldGrid(efields, bfields, angles);
    const auto nE = efields.size();
    const auto nB = bfields.size();
    const auto nA = angles.size();
    for (size_t j = 0; j < nB; ++j) {
      for (size_t k = 0; k < nA; ++k) {
        std::cout << "B = " << bfields[j] << " T, theta = "
                  << angles[k] * RadToDegree << " degree\n";
        std::cout << "   E [V/cm]     vE [cm/us]    alpha [1/cm]\n";
        for (size_t i = 0; i < nE; ++i) {
          double ve = 0.;
          gas->GetElectronVelocityE(i, j, k, ve);
          // Convert from cm/ns to mm/us.
          ve *= 1.e4;
          myvel[gasfile_i]=ve;
          double alpha = 0.;
          gas->GetElectronTownsend(i, j, k, alpha);
          alpha = exp(alpha);
          double lorentzangle = 0.;
          gas->GetElectronLorentzAngle(i, j, k, lorentzangle);
          mylorentz[gasfile_i]=lorentzangle;
          std::printf("%10.3f    %10.3f    %10.3f\n", efields[i], ve, lorentzangle);
        }
      }
    }
    gasfile_i++;
    temperature_i+=temp_step;
    delete gas;
  }
  for(int j=0;j<101;j++){
    std::printf("%10.3f,%10.3f\n", mytemp[j], mylorentz[j]);
    outfile << mytemp[j] << "," << mylorentz[j] << std::endl;
    outfile2 << mytemp[j] << "," << myvel[j] << std::endl;
  }
}
