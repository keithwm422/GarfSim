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
  gas->LoadGasFile("keith_co2_90_AR_10_T273.gas");
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
  const double rCathode= 175e-4;
  const double vAnode= 0;
  const double rAnode= 20e-4;
  const double vPotential= -2450;
  const double rPotential= 175e-4;
  
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;
  const double fieldsep = 0.2;
   
  for (int iplane=0;iplane<2;iplane++){
    for(int iw=0;iw<17;iw++){
      float y = 3.2-iw*cathodesep;
      float x = 8-16.0*iplane;
      cmp->AddWire(x,y,2 * rCathode, vCathode, "c");
      std::cout << " wire " << x << " " << y << " " << vCathode << " " << "c" << std::endl;
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
      std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
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
  sensor->SetTimeWindow(0,2,5000);
  
  cmp->AddReadout("a"); 
  sensor->AddElectrode(cmp,"a");
  vs1->SetSensor(sensor);
  
  //ViewDrift * vd = new ViewDrift();
  TCanvas* canvas3 = new TCanvas();
  vs1->SetCanvas(canvas3);

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

    driftline->AvalancheElectron(-300e-4+rAnode+100e-4,2.4,0,0);
    std::cout << __LINE__ << std::endl;
    driftline->GetAvalancheSize(ne, ni);
    std::cout << "avalanche # electrons= " << ne << " # ions= " << ni << std::endl;

/*  AvalancheMC * driftline_i = new AvalancheMC();
  driftline_i->SetDistanceSteps(0.001);
  //driftline_i->EnableMagneticField();
  driftline_i->SetSensor(sensor);
  driftline_i->EnableSignalCalculation();

  TrackHeed * track = new TrackHeed();
  track->EnableElectricField();
  track->EnableMagneticField();
  track->EnableDeltaElectronTransport();
    
  for (int itrack=0;itrack<ntrack;itrack++){
    int clustercount=0;
    int usedclustercount=0;
    //  float trackenergy = 1e9; 
    float trackenergy = 8*1e9;     // for Be
    float trackangle = trackang;
    float trackstartX = trackx;
    float trackstartY = 3.5;
    track->SetSensor(sensor);
    // track->SetParticle("p");
    track->SetParticleUser(8.4375e9,4);    // for Be
    track->SetKineticEnergy(trackenergy);
    //    track->EnablePlotting(vd);
    std::cout << " NEW TRACK STARTED " << itrack << " x: " << trackstartX << std::endl;

    double dirx = trackangle/TMath::Sqrt(1+trackangle*trackangle);
    double diry = -1.0/TMath::Sqrt(1+trackangle*trackangle);

    track->NewTrack(trackstartX,trackstartY,0,0,dirx,diry,0);
    
    for (int iplane=0;iplane<1;iplane++){
      for (int iw =0;iw<7;iw++){
        int iadd = iplane*7 + iw;
        SensewireSig[iadd]->Reset();
      }
    }
    double xcl, ycl, zcl, tcl, ecl, extra;
    int ncl;
    double r=0.01;
    double x=0, y=0, z=0, t=0, e=0, dx=0, dy=0, dz=0;
    double xendpoint = 0, yendpoint = 0, zendpoint=0;
    double xendpoint2 = 0, yendpoint2 = 0, zendpoint2=0;
    double tendpoint = 0, tendpoint2 = 0;
    int i=0;
    int stat;
    std::cout << "rough # clusters: " << track->GetClusterDensity()*5.6 << std::endl;
    while (track->GetCluster(xcl, ycl, zcl, tcl, ncl, ecl, extra)){
      std::cout << " new cluster: # used " << usedclustercount << " Total: " << clustercount << " " << xcl << " " << ycl << std::endl;
      clustercount++;
      clustersize.Fill(ncl);
      if(ncl<maxclustersize){
        usedclustercount++;
	      clustersizeused.Fill(ncl);
	      for(i = 0; i < ncl; i++){
	        track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
	        //std::cout << " ****  " << i  << " of " << ncl << "  electrons " << x << " " << y << " " << z << " " << t << " " << e << " " << dx << " " << dy << " " << dz << std::endl;
          //std::cout << __LINE__ << std::endl;
	        driftline->DriftElectron(x,y,z,0);
          //std::cout << __LINE__ << std::endl;
	        //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
	        driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint, xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
          //std::cout << " ion start point "  << xendpoint2 << " " << yendpoint2 << " " << tendpoint2 << std::endl;
	        double angle = RndmGaussian(0,1.4);
	        driftline_i->DriftIon(xendpoint2 + r*sin(angle), yendpoint2 + r*cos(angle), zendpoint2,tendpoint2); 
	        int iplane, iw;
	        iplane=0;
	        iw = (int)((yendpoint2+2.5)/anodesep);  
      	  if(iw>=0 && iw <7){
      	    int iadd = iw;
      	    std::cout << " signal on plane " << iplane << " wire " << iw <<  " " << iadd  << std::endl;
      	    vs1->PlotSignal("a");
	          if(realtimeplots){
	          //vd->Plot();
	          //canvas2->Update();
	          //canvas1->Update();
	          }
	          if(iadd>=0 && iadd<7){
	            SensewireSig[iadd]->Add(vs1->GetHistogram());
	          }
	        }
	        //else std::cout << " invalid wire " << iw << std::endl;
	        sensor->ClearSignal();
	      }
	      if(clustercount>40000) break;
      }
    }
    for(int iw=0;iw<7;iw++){
      SensewireSig[iw]->Write();
    }
    ncluster.Fill(clustercount);
    nclusterused.Fill(usedclustercount);
  }*/
  ncluster.Write();
  clustersize.Write();
  nclusterused.Write();
  clustersizeused.Write();
  TH1D * wire_sig[2]; 
  wire_sig[0] = new TH1D("wire_sig0","wire_sig0",5000,0,10000);
  wire_sig[1] = new TH1D("wire_sig1","wire_sig1",5000,0,10000);

  //vs1->GetHistogram()->Write();
  if(realtimeplots){
    vs1->PlotSignal("a");
    //   vd->Plot(true, false);
    //    canvas1->Update();
    //  canvas1->Print("newtrack.gif");
  }
  wire_sig[0]->Add(vs1->GetHistogram());
  wire_sig[0]->Write();
  Outfile->Close();
  app->Run(kTRUE);
}
