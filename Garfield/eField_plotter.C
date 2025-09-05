#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1.h>
#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "SolidBox.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "ViewCell.hh"
#include "ViewField.hh"
#include "TrackSimple.hh"
#include "ViewDrift.hh"
#include "TrackHeed.hh"
#include "DriftLineRKF.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "ViewSignal.hh"
#include "Random.hh"
#include <TROOT.h>
#include <TRint.h>
#include "DCsim.hh"
#include <fstream>
#include "Plotting.hh"
#include "TMath.h"
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <omp.h>
#include <chrono>

using namespace Garfield;

double getMean(const std::vector<double> &input){
  double sum = std::accumulate(input.begin(), input.end(), 0.0);
  double mean = sum / input.size();
  return mean;
}

double getSigma(const std::vector<double> &input){
  //std::cout << " size is " << static_cast<double>(input.size()) << std::endl;
  //double sum = std::accumulate(input.begin(), input.end(), 0.0);
  //double mean = sum / input.size();
  double mean = getMean(input);
  std::vector<double> diff(input.size());

  std::transform(input.begin(), input.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / input.size());
  return stdev;
}

int main(int argc, char * argv[]) {
  auto start = std::chrono::high_resolution_clock::now();

  bool realtimeplots = true;
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  std::stringstream wid;



  int pid = getpid();
  timeval t;
  gettimeofday(&t, NULL);
  int seed = pid*t.tv_usec;
  std::cout << "Random Seed: " << seed << std::endl;
  randomEngine.Seed(seed);

  for (int iplane=0;iplane<1;iplane++){
    for (int iw =0;iw<7;iw++){
      int iadd = iplane*7 + iw;
      wid.str("");
      wid << "a_" <<iplane << "_" << iw ;
      std::cout << iadd << " " << wid.str() << std::endl;
      std::string str(wid.str());
      const char * name = str.c_str();
    }
  }
  
  Garfield::plottingEngine.SetDefaultStyle();
  MediumMagboltz * gas = new MediumMagboltz();

  // Setup the gas
  const double pressure = 760.; //Torr
  const double temperature = 293.15; //K
 
  // Set the temperature [K] and pressure [Torr]
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("co2", 85, "ar", 15);

//  gas->LoadGasFile("co2_90_AR_10_T273.gas");
//  gas->LoadGasFile("keith_co2_85_AR_15_T273.gas");
  gas->LoadGasFile("Flight2024_Boff_P_755.865_T_293.15_multiE_90CO2_10Ar.gas");

  // lets just print out the drift velocity to a file?

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
//  cmp->SetMagneticField(0.,0.,0.0);
  cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,4.7,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
  const double vPotential= -2700;
  const double rPotential= 175e-4;
  
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;
  //const double cathodesep =0.2;
  //const double cathodesep =0.1;
  const double fieldsep = 0.2;
   
  for (int iplane=0;iplane<2;iplane++){
    for(int iw=0;iw<17;iw++){
//    for(int iw=0;iw<34;iw++){
//    for(int iw=0;iw<68;iw++){
      float y = 3.2-iw*cathodesep;
      float x = 7.62-15.24*iplane;
      cmp->AddWire(x,y,2 * rCathode, vCathode, "c");
      //std::cout << " wire " << x << " " << y << " " << vCathode << " " << "c" << std::endl;
    }
  }
  
  int nwire = 7;
  for (int iplane=0;iplane<1;iplane++){
    for(int iw=0;iw<nwire;iw++){
      wid << "a_" << iplane << "_" << iw;      
      float y = 2.4-iw*anodesep;
      float sign = 1.0;
      if(iw%2==0) sign = -1.0;
      float x = sign*300e-4;
      cmp->AddWire(x,y,2 * rAnode, vAnode, "a");
      //std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
    }
    for(int iw=0;iw<nwire+1;iw++){
      float y = 2.8-iw*potentialsep;
      float x = 0;
      cmp->AddWire(x,y,2 * rPotential, vPotential, "p");
      std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
    }
    cmp->AddWire(0,3.2,2 * rPotential, vPotential, "p2");
    cmp->AddWire(0,-3.2,2 * rPotential,vPotential, "p2");
  }

  //cmp->AddPlaneX(-8,-10000,"cP1"); 
  //cmp->AddPlaneX(8,-10000,"cP2");
  
  float v;
  float tweak;
  nwire=(int)(16.0/fieldsep) -1;
  for (int iplane=0;iplane<2;iplane++){
    for(int iw=0;iw<nwire;iw++){
      float y = 3.6-7.2*iplane;
      float x = 8-(iw+1)*fieldsep;
      tweak = -0.15;
      if(x>=0){ 
	      v= (1+tweak*(8-x)/8 )*( vCathode +  (vPotential-vCathode)*(8-x)/8.)  ;
      }
      else{
	      v=(1+tweak*(8+x)/8)*(vPotential + (vCathode-vPotential)*(-x)/8.);
	    }
      std::cout << " wire " << x << " " << y << " " << v << " " << "T" << std::endl;
      cmp->AddWire(x,y,2*rCathode, v, "T");
    }
  }

  Sensor * sensor = new Sensor;
  //ViewSignal * vs1 = new ViewSignal;

  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,20000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");
  //vs1->SetSensor(sensor);
  //double ex,ey,ez;
  //int stat_efield;
  //sensor->ElectricField(-1.0,-1.0,0,ex,ey,ez,gas,stat_efield); ///const double x, const double y, const double z,
                           //double &ex, double &ey, double &ez, double &v,
                           //Medium *&medium, int &status


  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetComponent(cmp);
  //  fieldView.SetElectricFieldRange(-800., 800.0);
  //constexpr bool plotProfile = true;
  if (true) {
    fieldView.PlotProfile(-7.5, 0.2, 0.0,7.5 , 0.2, 0., "ex"); // this line is for a field-shaping/potential wire region near (0,0.4,0). 
    //fieldView.PlotProfile(-7.5, 0.0, 0.0,7.5 , 0.0, 0., "ex"); // this line is for a sense wire region near (0,0,0). 

  } else {
    fieldView.SetArea(-7.5, -2.4, 7.5, 2.4);
    //fieldView.SetArea(-1.5, -0.4, 1.5, 0.4);
    fieldView.PlotContour("ex");
    ViewCell cellView;
    //cellView.SetCanvas(&canvas);
    //cellView.SetComponent(&cmp);
    //cellView.SetArea(1.1 * xMin, 1.1 * yMin, -1., 1.1 * xMax, 1.1 * yMax, 1.);
    //cellView.Plot2d();
  }
  std::cout << "Efield found: \n";
   auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  app->Run(kTRUE);
  return 0;
}
