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



using namespace Garfield;

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
  std::ofstream outputfilecathode("cathode_wires.txt");
  std::ofstream outputfilesense("sense_wires.txt");
  std::ofstream outputfilegetelectron("electronstartpoints.txt");
  std::ofstream outputfilegetelectronendpoint("electronendpoints.txt");
  std::ofstream outputfiledrifttimes("electrondrifttimes_288_762_x4c.txt");
  //name for drifttimes vs drift distance

  outputfilecathode << "wire cathode information" << std::endl;
  outputfilesense << "wire sense information" << std::endl;
  outputfiledrifttimes << "e,t" << std::endl;

  //outputfile << "num_electrons" << "," << "running_average" << "," << "new_average" << "," << "running_sum" << "," << "curr_sig" << "," << "ne" << "," << "ni" << std::endl;

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
  gas->LoadGasFile("Flight2024_P_755.038_T_288.15_.gas");

  // lets just print out the drift velocity to a file?

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
//  cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  //SolidBox * enclosure = new SolidBox(0,0,0,10,4.7,11);
  SolidBox * enclosure = new SolidBox(0,0,0,25,25,25);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  const double vCathode = -7500;
  const double rCathode = 175e-4;
  const double vAnode = 0;
  const double rAnode = 20e-4;
  const double vPotential= -2700;
  const double rPotential = 175e-4;
  const double anode_cathode_sep=7.62;
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;

  for (int iplane=0;iplane<4;iplane++){
    for(int iw=0;iw<121;iw++){
      float y = 24-iw*cathodesep;
      float x = 24-16.0*iplane;
      cmp->AddWire(x,y,2 * rCathode, vCathode, "c");
std::cout << " wire " << x << " " << y << " " << vCathode << " " << "c" << std::endl;
    }
  }
 
  for (int iplane=0;iplane<3;iplane++){
    for(int iw=0;iw<60;iw++){
      wid << "a_" << iplane << "_" << iw;
      
      float y = 23.6-iw*anodesep;
      float sign = 1.0;
      if(iw%2==0) sign = -1.0;
      float x = 16-iplane*16+sign*300e-4;
      cmp->AddWire(x,y,2 * rAnode, vAnode, "a");
      std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;

    }
    for(int iw=0;iw<61;iw++){
      float y = 24-iw*potentialsep;
      float x = 16-iplane*16;
      cmp->AddWire(x,y,2 * rPotential, vPotential, "p");
std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
    }
  }

  cmp->AddPlaneX(-25,-7500,"cP"); 
  cmp->AddPlaneX(25,-7500,"cP");
  
  float v;
  for (int iplane=0;iplane<2;iplane++){
    for(int iw=0;iw<241;iw++){
      float y = 24.5-49*iplane;
      float x = 24-iw*cathodesep/2;
      if(x>16)
	{ 
	  v=vCathode + (vPotential-vCathode)*(24-x)/8.;
	}
      else if(x<=16 && x>8) 
	{
	  v=vPotential + (vCathode-vPotential)*(16-x)/8.;
	}
      else if(x<=8 && x>0)
	{ 
	  v=vCathode + (vPotential-vCathode)*(8-x)/8.;
	}
      else if(x<=0 && x>-8) 
	{
	  v=vPotential + (vCathode-vPotential)*(-x)/8.;
	}
      else if(x<=-8 && x>-16)
	{ 
	  v=vCathode + (vPotential-vCathode)*(-8-x)/8.;
	}
      else if(x<=-16 && x>-24) 
	{
	  v=vPotential + (vCathode-vPotential)*(-16-x)/8.;
	}
      std::cout << " wire " << x << " " << y << " " << v << " " << "T" << std::endl;
      cmp->AddWire(x,y,2 * rCathode, v, "T");
    }
  }

  Sensor * sensor = new Sensor;
  ViewSignal * vs1 = new ViewSignal;
  std::cout << __LINE__ << std::endl;

  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,5000);
  std::cout << __LINE__ << std::endl;
  
  cmp->AddReadout("a"); 
  sensor->AddElectrode(cmp,"a");
  vs1->SetSensor(sensor);
  std::cout << __LINE__ << std::endl;
  
  ViewDrift * vd = new ViewDrift();
  TCanvas* canvas3 = new TCanvas("hye");
  vs1->SetCanvas(canvas3);
  TCanvas* canvas1 = new TCanvas();
  TCanvas* canvas2 = new TCanvas();
  TCanvas* canvas4 = new TCanvas();
  std::cout << __LINE__ << std::endl;
  if(realtimeplots){
    ViewCell * view = new ViewCell();
    ViewField * viewfield = new ViewField();
    view->SetComponent(cmp);
    view->DisableWireMarkers();
    viewfield->SetComponent(cmp);
    viewfield->SetSensor(sensor);
    //viewfield->PlotSurface("e");
    ViewDrift * vd = new ViewDrift();
    vd->SetCanvas(canvas2);
    //viewfield->SetElectricFieldRange(0.,3e5);
    viewfield->SetCanvas(canvas4);
    //viewfield->PlotContour("v");  
    viewfield->Plot("v","colz");  
    view->SetCanvas(canvas4);
    view->Plot2d();
    canvas1->Update();
    canvas1->Print("Canvas1.pdf","pdf");
    canvas2->Update();
    canvas2->Print("Canvas2.pdf","pdf");
    canvas3->Update();
    canvas3->Print("Canvas3.pdf","pdf");
    canvas4->Update();
    canvas4->Print("Canvas4.pdf","pdf");

  }
  outputfilecathode.close();
  outputfilesense.close();
  outputfilegetelectron.close();
  outputfilegetelectronendpoint.close();
  outputfiledrifttimes.close();
  //cmp->PrintCell();
  app->Run(kTRUE);
}
