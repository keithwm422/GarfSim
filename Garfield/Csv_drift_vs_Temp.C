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

  double myvelocity[5] =   {0};
  double    mytemps[5] =   {0};

  MediumMagboltz * gas = new MediumMagboltz();
  std::ofstream outfile;
  outfile.open("driftv_vs_Temperature.csv");
  outfile << "T,v" << std::endl;
  int gasfile_i=0;
  while(gasfile_i<5){
    switch (gasfile_i)
    {
    case 0:
      gas->LoadGasFile("Flight2024_P_755.038_T_283.15_.gas");
      mytemps[gasfile_i]=283.15;
      break;
    case 1:
      gas->LoadGasFile("Flight2024_P_755.038_T_288.15_.gas");
      mytemps[gasfile_i]=288.15;
      break;
    case 2:
      gas->LoadGasFile("Flight2024_P_755.038_T_293.15_.gas");
      mytemps[gasfile_i]=293.15;
      break;
    case 3:
      gas->LoadGasFile("Flight2024_P_755.038_T_298.15_.gas");
      mytemps[gasfile_i]=298.15;
      break;
    case 4:
      gas->LoadGasFile("Flight2024_P_755.038_T_303.15_.gas");
      mytemps[gasfile_i]=303.15;
      break;
    }
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
          double alpha = 0.;
          myvelocity[gasfile_i]=ve;
          gas->GetElectronTownsend(i, j, k, alpha);
          alpha = exp(alpha);
          std::printf("%10.3f    %10.3f    %10.3f\n", efields[i], ve, alpha);
        }
      }
    }
    gasfile_i++;
  }
  for(int j=0;j<5;j++){
    std::printf("%10.3f,%10.3f\n", mytemps[j], myvelocity[j]);
    outfile << mytemps[j] << "," << myvelocity[j] << std::endl;
  }

  // app.Run(kTRUE);

}
