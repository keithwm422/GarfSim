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
#include <TGraph.h>
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

  bool realtimeplots = false;
  bool try_getting_sig_from_sensor=true;

  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  TFile * Outfile = new TFile(simoutFile,"recreate");

  std::stringstream wid;
  
  std::ofstream outputfile("gasgain_study_rev1.txt");
  outputfile << "num_electrons" << "," << "running_average" << "," << "new_average" << "," << "running_sum" << "," << "curr_sig" << "," << "ne" << "," << "ni" << std::endl;

  TGraph * mygraph;
  double average_sig[3]={0};
  double gas_comp[3]={15,10,5};
  int pid = getpid();
  timeval t;
  gettimeofday(&t, NULL);
  int seed = pid*t.tv_usec;
  std::cout << "Random Seed: " << seed << std::endl;
  randomEngine.Seed(seed);
  Garfield::plottingEngine.SetDefaultStyle();

  // Setup the gas- all mixes same
  const double pressure = 760.; //Torr
  const double temperature = 293.15; //K
  char * IonData = getenv("GARFIELD_IONDATA") ;

  // pointers to instances for calculating the avalanches in
  MediumMagboltz         * gas[3];
  ComponentAnalyticField * cmp[3]; 
  GeometrySimple         * geo[3];
  SolidBox               * enclosure[3];

  // parameters about the enclosure and wires
  const double vCathode       = -7500;
  const double rCathode       = 175e-4;
  const double vAnode         = 0;
  const double rAnode         = 20e-4;
  const double vPotential     = -2000;
  const double rPotential     = 175e-4;
  const double anodesep       = 0.8;
  const double potentialsep   = 0.8;
  const double cathodesep     = 0.4;
  const double fieldsep       = 0.2;
  int nwire                   = 7;
  const int nTwires           = (int)(16.0/fieldsep) -1;

  //  sensors for recording signals
  Sensor      * sensor[3];
  ViewSignal  * vs1[3];

  // canvas and histograms for all plots
  TCanvas     * canvas3[3];
  TH1D        * wire_sig[nwire];

  // Avalanch electrons at same location in each iteration
  AvalancheMC * driftline[3];

  for(int i=0;i<3;i++){
    gas[i]          = new MediumMagboltz();
    cmp[i]          = new ComponentAnalyticField();
    geo[i]          = new GeometrySimple();
    enclosure[i]    = new SolidBox(0,0,0,10,4.7,11);
    sensor[i]       = new Sensor;
    vs1[i]          = new ViewSignal;
    canvas3[i]      = new TCanvas();
    driftline[i]    = new AvalancheMC();
    // Set the temperature [K] and pressure [Torr]
    gas[i]->SetTemperature(temperature);
    gas[i]->SetPressure(pressure);
    //gas->SetComposition("co2", 85, "ar", 15);
    //gas->LoadGasFile("co2_90_AR_10_T273.gas");
    //gas->LoadGasFile("keith_co2_85_AR_15_T273.gas");
    if(i==0)      gas[i]->LoadGasFile("keith_co2_85_AR_15_T273.gas");
    else if(i==1) gas[i]->LoadGasFile("keith_co2_90_AR_10_T273.gas");
    else if(i==2) gas[i]->LoadGasFile("keith_co2_95_AR_5_T273.gas");
    gas[i]->LoadIonMobility(IonData);
    gas[i]->PrintGas();
    // now cmp, basically all the same except the gas added is different in each
    cmp[i]->SetMagneticField(0.,0.,1.0);
    geo[i]->AddSolid(enclosure[i], gas[i]);
    cmp[i]->SetGeometry(geo[i]);
    for (int iplane=0;iplane<2;iplane++){
      for(int iw=0;iw<17;iw++){
        float y = 3.2-iw*cathodesep;
        float x = 8-16.0*iplane;
        cmp[i]->AddWire(x,y,2 * rCathode, vCathode, "c");
        std::cout << " wire " << x << " " << y << " " << vCathode << " " << "c" << std::endl;
      }
    }
    for (int iplane=0;iplane<1;iplane++){
      for(int iw=0;iw<nwire;iw++){
        wid << "a_" << iplane << "_" << iw;      
        float y    = 2.4-iw*anodesep;
        float sign = 1.0;
        if(iw%2==0) sign = -1.0;
        float x    = sign*300e-4;
        cmp[i]->AddWire(x,y,2 * rAnode, vAnode, "a");
        std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
      }
      for(int iw=0;iw<nwire+1;iw++){
        float y = 2.8-iw*potentialsep;
        float x = 0;
        cmp[i]->AddWire(x,y,2 * rPotential, vPotential, "p");
        std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
      }
      cmp[i]->AddWire(0,3.2,2 * rPotential, vPotential, "p2");
      cmp[i]->AddWire(0,-3.2,2 * rPotential,vPotential, "p2");
    }
    //cmp->AddPlaneX(-8,-10000,"cP1"); 
    //cmp->AddPlaneX(8,-10000,"cP2");
    float v;
    float tweak;
    for (int iplane=0;iplane<2;iplane++){
      for(int iw=0;iw<nTwires;iw++){
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
        cmp[i]->AddWire(x,y,2*rCathode, v, "T");
      }
    }
    sensor[i]->AddComponent(cmp[i]);
    sensor[i]->SetTimeWindow(0,2,5000);
    cmp[i]->AddReadout("a"); 
    sensor[i]->AddElectrode(cmp[i],"a");
    vs1[i]->SetSensor(sensor[i]);
    vs1[i]->SetCanvas(canvas3[i]);
    /*if(realtimeplots){
      ViewCell  * view = new ViewCell();
      ViewField * viewfield = new ViewField();
      view->SetComponent(cmp[i]);
      view->DisableWireMarkers();
      viewfield->SetComponent(cmp[i]);
      viewfield->SetSensor(sensor[i]);
      //    viewfield->PlotSurface("e");
      ViewDrift * vd      = new ViewDrift();
      TCanvas   * canvas1 = new TCanvas();
      TCanvas   * canvas2 = new TCanvas();
      TCanvas   * canvas4 = new TCanvas();
      vd->SetCanvas(canvas2);
      viewfield->SetElectricFieldRange(0.,3e5);
      viewfield->SetCanvas(canvas4);
      viewfield->PlotContour();  
      view->SetCanvas(canvas1);
      view->Plot2d();
      canvas1->Update();
      canvas1->Print("View2D.pdf","pdf");
    }*/
    //  driftline[i]->EnableDebugging();  
    driftline[i]->SetDistanceSteps(0.001);
    //driftline[i]->EnableMagneticField();
    driftline[i]->EnableDiffusion();
    driftline[i]->SetSensor(sensor[i]);
    //  driftline[i]->EnablePlotting(vd);
    //  driftline[i]->EnableSignalCalculation();
    unsigned int ne=0, ni=0;
    int num_electrons=0;
    double min_variation=0.00001;
    double running_sum=0;
    double running_average=0;
    double new_average=0;
    double curr_sig=0;
    bool keep_running=true;
    while(keep_running){
      driftline[i]->AvalancheElectron(-300e-4+rAnode+100e-4,2.4,0,0);
      driftline[i]->GetAvalancheSize(ne, ni);
      curr_sig = sensor[i]->GetSignal("a",0,0);
      sensor[i]->ClearSignal();
      num_electrons++;
      running_sum+=curr_sig;
      new_average=running_sum/num_electrons;
      if(num_electrons>10 && TMath::Abs((new_average-running_average)/running_average)<=min_variation) keep_running=false;
      else{
        // compute the average some more
        outputfile << num_electrons << "," << running_average << "," << new_average << "," << running_sum << "," << curr_sig << "," << ne << "," << ni << std::endl;
        //std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " new average is " << new_average << std::endl;    

        running_average=new_average;
      }
    }
    std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << std::endl;
    average_sig[i]=running_average;
    if(realtimeplots){
      vs1[i]->PlotSignal("a");
      //   vd->Plot(true, false);
      //    canvas1->Update();
      //  canvas1->Print("newtrack.gif");
    }
    //wire_sig[i]->Add(vs1[i]->GetHistogram());
  }
  average_sig[2]=average_sig[2]/average_sig[1];
  average_sig[0]=average_sig[0]/average_sig[1];
  average_sig[1]=average_sig[1]/average_sig[1];
  mygraph = new TGraph(3,gas_comp,average_sig);
  mygraph->Write();
  Outfile->Close();
  outputfile.close();
  app->Run(kTRUE);
}
