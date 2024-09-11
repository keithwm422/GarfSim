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
  double myvelocity;
  MediumMagboltz * gas = new MediumMagboltz();
//  gas->LoadGasFile("Flight2024_Bgrid_P_755.865_T_293.15_.gas");
  gas->LoadGasFile("Flight2024_BvsTstudy_P_755.865_T_283.15_.gas");
  gas->PrintGas();
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
      std::cout << "   E [V/cm]     vE [cm/us]    lorentz []\n";
      for (size_t i = 0; i < nE; ++i) {
        double ve = 0.;
        gas->GetElectronVelocityE(i, j, k, ve);
        // Convert from cm/ns to cm/us.
        ve *= 1.e3;
        double lorentzangle = 0.;
        gas->GetElectronLorentzAngle(i, j, k, lorentzangle);
        //alpha = exp(alpha);
        std::printf("%10.3f    %10.3f    %10.3f\n", efields[i], ve, lorentzangle);
      }
    }
  } 
  
  // app.Run(kTRUE);

}
