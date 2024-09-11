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
/* files used
Gascomp_90_10_P_620.579_T_293.15_.gas
Gascomp_90_10_P_646.436_T_293.15_.gas
Gascomp_90_10_P_672.294_T_293.15_.gas
Gascomp_90_10_P_698.151_T_293.15_.gas
Gascomp_90_10_P_724.009_T_293.15_.gas
Gascomp_90_10_P_749.866_T_293.15_.gas
Gascomp_90_10_P_775.724_T_293.15_.gas
Gascomp_90_10_P_801.581_T_293.15_.gas
Gascomp_90_10_P_827.438_T_293.15_.gas
Gascomp_90_10_P_879.153_T_293.15_.gas
Gascomp_90_10_P_905.011_T_293.15_.gas
Gascomp_90_10_P_930.868_T_293.15_.gas
*/
int main(int argc, char * argv[]) {

  double myvelocity[13] =   {0};
  double mypressure[13] =   {0};

  MediumMagboltz * gas = new MediumMagboltz();
  std::ofstream outfile;
  outfile.open("driftv_vs_pressure.csv");
  outfile << "P,v" << std::endl;
  int gasfile_i=0;
  while(gasfile_i<5){
    switch (gasfile_i)
    {
    case 0:
      gas->LoadGasFile("Gascomp_90_10_P_724.009_T_293.15_.gas"); // 14.0
      mypressure[gasfile_i]=724.009;
      break;
    case 1:
      gas->LoadGasFile("Gascomp_90_10_P_736.937_T_293.15_.gas"); // 14.25
      mypressure[gasfile_i]=736.937;
      break;
    case 2:
      gas->LoadGasFile("Gascomp_90_10_P_749.866_T_293.15_.gas"); // 14.5
      mypressure[gasfile_i]=749.866;
      break;
    case 3:
      gas->LoadGasFile("Gascomp_90_10_P_762.795_T_293.15_.gas"); // 14.75
      mypressure[gasfile_i]=762.795;
      break;
    case 4:
      gas->LoadGasFile("Gascomp_90_10_P_775.724_T_293.15_.gas"); // 15
      mypressure[gasfile_i]=775.724;
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
    std::printf("%10.3f,%10.3f\n", mypressure[j], myvelocity[j]);
    outfile << mypressure[j] << "," << myvelocity[j] << std::endl;
  }

  // app.Run(kTRUE);

}
