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



using namespace Garfield;

double getSigma(const std::vector<double> &input){
  std::cout << " size is " << static_cast<double>(input.size()) << std::endl;
  double sum = std::accumulate(input.begin(), input.end(), 0.0);
  double mean = sum / input.size();

  std::vector<double> diff(input.size());

  std::transform(input.begin(), input.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / input.size());
  return stdev;
}

int main(int argc, char * argv[]) {

  char * simoutFile;
  double trackx,trackang;
  int ntrack;
  char * cdum;
  cdum = getenv("DCSimNtrack");
  ntrack  = atoi(cdum);
  std::cout << " # tracks " << ntrack << std::endl;
  cdum= getenv("DCSimtrackx");
  trackx = atof(cdum);
  std::cout << " Track Start X " << trackx << std::endl;
  cdum = getenv("DCSimtrackang");
  trackang = atof(cdum);
  std::cout << " Track Angle " << trackang << std::endl;
  simoutFile = getenv("DCSimOutFile");
  bool realtimeplots = true;
  int maxclustersize = 10000;
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  TH1D *SensewireSig[180];
  TFile * Outfile = new TFile(simoutFile,"recreate");

  TH1F ncluster("ncluster","ncluster",100,0.,2000.);
  TH1F nclusterused("nclusterused","nclusterused",100,0.,2000.);

  TH1F clustersize("clustersize","clustersize",1000,0.,100.);
  TH1F clustersizeused("clustersizeused","clustersizeused",1000,0.,100.);

  std::stringstream wid;
  double x_i=7.60249, y_i=0, z_i=0.2, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0
  std::ofstream outputfilecathode("cathode_wires_notracks.txt");
  std::ofstream outputfilesense("sense_wires_notracks.txt");
  std::ofstream outputfilegetelectron("electronstartpoints_notracks.txt");
  std::ofstream outputfilegetelectronendpoint("electronendpoints_notracks.txt");
  std::ofstream outputfiledrifttimes("electrondrifttimes_298_notracks_BON_zi02.txt");
  std::ofstream outputfile("electron_average_drifttimes_298_notracks_BON_zi02.txt");
  outputfile << "num_electrons" << "," << "running_average" << "," << "new_average" << "," << "running_sum" << "," << "curr_sig" << "," << "ne" << "," << "ni" << std::endl;
  std::vector<double> electrons_drifted;
  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);
  outputfiledrifttimes << std::fixed << std::showpoint;
  outputfiledrifttimes << std::setprecision(10);
  //name for drifttimes vs drift distance


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
      SensewireSig[iadd] = new TH1D(name,name,5000,0,10000);
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
  gas->LoadGasFile("Flight2024_BvsTstudy_P_755.865_T_298.15_.gas");

  // lets just print out the drift velocity to a file?

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
//  cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,4.7,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
  const double vPotential= -250;
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
      outputfilecathode << " wire " << x << " " << y << " " << vCathode << " " << "c" << std::endl;
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
      outputfilesense << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
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
  ViewSignal * vs1 = new ViewSignal;

  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,10000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");
  vs1->SetSensor(sensor);
  
  ViewDrift * vd = new ViewDrift();
  TCanvas* canvas3 = new TCanvas("hye");
  vs1->SetCanvas(canvas3);
    TCanvas* canvas1 = new TCanvas();
    TCanvas* canvas2 = new TCanvas();
    TCanvas* canvas4 = new TCanvas();
  if(realtimeplots){
    ViewCell * view = new ViewCell();
    ViewField * viewfield = new ViewField();
    view->SetComponent(cmp);
    view->DisableWireMarkers();
    viewfield->SetComponent(cmp);
    viewfield->SetSensor(sensor);
    //    viewfield->PlotSurface("e");
    ViewDrift * vd = new ViewDrift();
    TCanvas* canvas1 = new TCanvas();
    TCanvas* canvas2 = new TCanvas();
    TCanvas* canvas4 = new TCanvas();
    vd->SetCanvas(canvas2);
    viewfield->SetElectricFieldRange(0.,3e5);
    viewfield->SetCanvas(canvas4);
    viewfield->PlotContour();  
    view->SetCanvas(canvas1);
    view->Plot2d();
    canvas1->Update();
    canvas1->Print("View2D.pdf","pdf");
  }
 
  AvalancheMC * driftline = new AvalancheMC();
  //  driftline->EnableDebugging();  
  driftline->SetDistanceSteps(0.001);
  //driftline->EnableMagneticField();
  driftline->EnableDiffusion();
  driftline->SetSensor(sensor);
  //  driftline->EnablePlotting(vd);
  //  driftline->EnableSignalCalculation();
  unsigned int ne=0, ni=0;
    std::cout << __LINE__ << std::endl;

    //driftline->AvalancheElectron(-300e-4+rAnode+100e-4,2.4,0,0);
    //std::cout << __LINE__ << std::endl;
    //driftline->GetAvalancheSize(ne, ni);
    //std::cout << "avalanche # electrons= " << ne << " # ions= " << ni << std::endl;

  for (int iplane=0;iplane<1;iplane++){
    for (int iw =0;iw<7;iw++){
      int iadd = iplane*7 + iw;
      SensewireSig[iadd]->Reset();
    }
  }
  double r=0.01;
  // really the x_i here is 7.62 but also subtract off the diameter of the cathode wire (2*rCathode)
  int num_electrons=0;
  double min_variation=0.001;
  double running_sum=0;
  double running_average=0;
  double new_average=0;
  double curr_sig=0;
  bool keep_running=true;
  while(keep_running){
    double xendpoint = 0, yendpoint = 0, zendpoint=0;
    double xendpoint2 = 0, yendpoint2 = 0, zendpoint2=0;
    double tendpoint = 0, tendpoint2 = 0;
	  //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
    int i=0;
    int stat=0;
	  //outputfilegetelectron << "startpoint" << i  << " of " << ncl << "  electrons " << x << " " << y << " " << z << " " << t << " " << e << " " << dx << " " << dy << " " << dz << std::endl;
    //std::cout << __LINE__ << std::endl;
	  driftline->DriftElectron(x_i,y_i,z_i,0);
    //std::cout << __LINE__ << std::endl;
	  //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
	  driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint, xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
    curr_sig=tendpoint2-tendpoint;
    electrons_drifted.push_back(curr_sig);
    num_electrons++;
    running_sum+=curr_sig;
    new_average=running_sum/num_electrons;
    if(num_electrons>300 && TMath::Abs((new_average-running_average)/running_average)<=min_variation) keep_running=false;
    else{
      // compute the average some more
      outputfile << num_electrons << "," << running_average << "," << new_average << "," << running_sum << "," << curr_sig << "," << ne << "," << ni << std::endl;
      //std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " new average is " << new_average << std::endl;    
      running_average=new_average;
      outputfiledrifttimes << i << "," << tendpoint2-tendpoint << std::endl;
    }
  }
  outputfilecathode.close();
  outputfilesense.close();
  outputfilegetelectron.close();
  outputfilegetelectronendpoint.close();
  outputfiledrifttimes.close();
  std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " RMS is : " << getSigma(electrons_drifted) << std::endl;
  outputfile.close();

  app->Run(kTRUE);
}
