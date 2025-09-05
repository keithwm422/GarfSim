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
#include "ViewMedium.hh"
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


  // TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
  const double zpos = invals[0]; // HELIX zposition slice to plot Efield in millimeters
  std::cout << "zpos slicing is [mm] : " << zpos << std::endl;

  TFile * OutFile = new TFile("trackres_test1.root","recreate");
  TH1D ncluster("ncluster","ncluster",100,0.,2000.);
  TH1D nclusterused("nclusterused","nclusterused",100,0.,2000.);
  TH1D clustersize("clustersize","clustersize",1000,0.,100.);
  TH1D clustersizeused("clustersizeused","clustersizeused",1000,0.,100.);
  TH1D clusterenergy("clusterEnergy","ClusterEnergy",2000,-0.5,1999.5);
  TH1D ElectronKE("ElectronKE","ElectronKE",2000,-0.5,1999.5);
  TH1D ClusterElectronsSize("ClusterElectronsSize","ClusterElectronsSize",1000,-0.5,999.5);
  TH1D NumElectronsPerCluster("NumElectronsPerCluster","NumElectronsPerCluster",1000,-0.5,999.5);
  TH1D NumPhotonsPerCluster("NumPhotonsPerCluster","NumPhotonsPerCluster",1000,-0.5,999.5);
  TH1D NumIonsPerCluster("NumIonsPerCluster","NumIonsPerCluster",1000,-0.5,999.5);

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
  std::stringstream efieldfilename;
  efieldfilename << "efieldprofile_zpos_" << zpos << ".root";      

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
  //gas->LoadGasFile("Flight2024_Bon1T_P_755.865_T_293.15_10Ar_90CO2_multiE.gas");

  // lets just print out the drift velocity to a file?
  //ViewMedium mediumView;
  //mediumView.SetMedium(gas);
  //mediumView.PlotElectronVelocity('e');
  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
  const double vPotential= -2700;
  const double rPotential= 175e-4;
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;
  const double fieldsep = 0.2;
  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,20,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  // full column goes
 /*HLX_Geometry referenced: center of wires is (y,z) in millimeters
   +0.3,+284
   -0.3,+276
   ...
   +0.3,-276
   -0.3,-284
  */
 int nwire = 72;
 float last_y_anode=0;
 float last_y_potential=0;
 float starting_y=28.4;
 for (int iplane=0;iplane<1;iplane++){
   for(int iw=0;iw<nwire;iw++){
     wid << "a_" << iplane << "_" << iw;      
     float y = starting_y-iw*anodesep;
     float sign = 1.0;
     if(iw%2==0) sign = -1.0;
     float x = sign*300e-4;
     //float x = 0;
     cmp->AddWire(x,y,2 * rAnode, vAnode, "a");
     std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
     last_y_anode=y;
   }
   for(int iw=0;iw<nwire;iw++){
     float y = (starting_y-(potentialsep/2.0))-iw*potentialsep;
     float x = 0;
     cmp->AddWire(x,y,2 * rPotential, vPotential, "p");
     std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
     last_y_potential=y;
   }
   // add in 3 potential wires at the bottom and top now separated by 4mm from eachother and by 4mm from the bottom or top wire
   cmp->AddWire(0,starting_y+(potentialsep/2.0),2 * rPotential, vPotential, "pT");
   cmp->AddWire(0,starting_y+(2.0*potentialsep/2.0),2 * rPotential, vAnode, "pT");
   cmp->AddWire(0,starting_y+(3.0*potentialsep/2.0),2 * rPotential, vPotential, "pT");
   std::cout << " wire " << 0 << " " << starting_y+(potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
   std::cout << " wire " << 0 << " " << starting_y+(2.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
   std::cout << " wire " << 0 << " " << starting_y+(3.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
   cmp->AddWire(0,last_y_potential-(potentialsep/2.0),2 * rPotential, vAnode, "pT");
   cmp->AddWire(0,last_y_potential-(2.0*potentialsep/2.0),2 * rPotential, vPotential, "pT");
   cmp->AddWire(0,last_y_potential-(3.0*potentialsep/2.0),2 * rPotential, vAnode, "pT");
   std::cout << " wire " << 0 << " " << last_y_potential-(potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
   std::cout << " wire " << 0 << " " << last_y_potential-(2.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
   std::cout << " wire " << 0 << " " << last_y_potential-(3.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
 }
 cmp->AddPlaneX(-7.62,-7500,"cP1"); 
 cmp->AddPlaneX(7.62,-7500,"cP2");

 float v=0;
 float tweak=0;
 float drift_dist_max=7.62*1.0;
 int nstrips_per_PCB=40; // from field-shaping PCBs
 int nstrips=2*nstrips_per_PCB; // 2 pcbs for a drift cell (1 on left and 1 on right)
 int count_me=0;
 double strip_sep=drift_dist_max/(double) (nstrips_per_PCB+2);  // plus 2 for the separation from the -7500V and GND that the PCBs are offset from the planes by
 double size_tweak=0.95;
 double strip_size=size_tweak*strip_sep;
 float starting_y_strips=starting_y+(4.0*potentialsep/2.0);
 float last_y_strips=last_y_potential-(4.0*potentialsep/2.0);
 float starting_x_strips=drift_dist_max-(1.5*strip_sep); // one pad away is where the first division should start
 float voltage_step=(vAnode-vCathode)/(nstrips_per_PCB+1);
 // 2 planes because top and bottom of chamber
 // starting_x
 // PCB thickness 1/16 inch
 double pcb_thiqness=0.15875; // in cm
 int num_wires_for_rep=4;

 for (int iplane=0;iplane<2;iplane++){
   for(int ipad=0;ipad<nstrips_per_PCB;ipad++){

     float y = starting_y_strips; // top or btm
     if(iplane >0) y=last_y_strips;
     float sign = 0.0;
     if(ipad%2==0) sign = -1.0;
     y = y+(sign*pcb_thiqness);
     float x = starting_x_strips-((double)(ipad)*strip_sep); // left or right
     float left_x = (-1.0*starting_x_strips)+((double)(ipad)*strip_sep); // left or right
     //v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vPotential-vCathode)*(drift_dist_max-x)/drift_dist_max)  ;
     //v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vAnode-vCathode)*(drift_dist_max-x)/drift_dist_max);
     v=vCathode +  voltage_step*(ipad+1);
     //std::cout << " wire " << x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
     //std::cout << " wire " << left_x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
     count_me++;
     //cmp->AddWire(x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
     //cmp->AddWire(left_x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
     // add many wires per pad strip?
     float pad_r=strip_size/(2.0*(float)(num_wires_for_rep));
     for(int wire_in_pad=0;wire_in_pad<num_wires_for_rep;wire_in_pad++){
       float x_delt=x-((float)(num_wires_for_rep-1)*pad_r)+((float)(wire_in_pad)*2.0*pad_r);
       float leftx_delt=left_x-((float)(num_wires_for_rep-1)*pad_r)+((float)(wire_in_pad)*2.0*pad_r);
       cmp->AddWire(x_delt,y,pad_r,v,"T");
       cmp->AddWire(leftx_delt,y,pad_r,v,"T");
       std::cout << " wire " << x_delt << " " << pad_r << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
     }
   }
 }
 std::cout << " strip_sep is: "  << strip_sep  << std::endl;
 std::cout << " estimated diameter should be "  <<  2.0*drift_dist_max/(double) (2.0*(nstrips_per_PCB-1)) << std::endl;
 std::cout << " diameter of cathode is: "  << 2.0*rCathode << std::endl;
 std::cout << " strip_size is: "  << strip_size  << std::endl;

 Sensor * sensor = new Sensor;
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


  //TCanvas canvas("c", "", 600, 600);
  TrackHeed track(sensor);
  track.EnableElectricField();
  // turn on delta electron transport
  track.EnableDeltaElectronTransport();
  //track.DisableDeltaElectronTransport();
  //track.EnableCoulombScattering();
  //track.SetParticle("muon");
  track.SetParticleUser(8.4375e9,4);    // for Be
  track.SetKineticEnergy(8*1e9);
  //track.SetEnergy(1700.e9); // in what units?
  DriftLineRKF drift(sensor);
  // Create a canvas.
  TCanvas * cD = new TCanvas("cD", "", 600, 600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  drift.EnablePlotting(&driftView);
  track.EnablePlotting(&driftView);
  ViewCell cellView;
  cellView.SetCanvas(cD);
  cellView.SetComponent(cmp);

  //const double rTrack = ;
  
  double xcl, ycl, zcl, tcl, ecl, extra;
  int nel,nph,nio;
  int ncl;
  double r=0.01;
  double x=0, y=0, z=0, t0=0, e=0, dx=0, dy=0, dz=0;
  const double x0 = 7.5;
  const double y0 = 20;
  double totalcluster_energy=0;
  int num_deltas=0;
  int num_electrons=0;
  track.NewTrack(x0, y0, 0, 0, -0.1, -20, 0); //const double x0, const double y0, const double z0, t0, dx0, dy0, dz0
  // Loop over the clusters along the track.
  std::cout << "clusters along track have size: "  << track.GetClusters().size() << std::endl;
  std::cout << "First cluster x,y,z,t is : "  << track.GetClusters().front().x << "," << track.GetClusters().front().y << "," << track.GetClusters().front().z << "," << track.GetClusters().front().t << std::endl;
  std::cout << "Last cluster x,y,z,t is : "  << track.GetClusters().back().x << "," << track.GetClusters().back().y << "," << track.GetClusters().back().z << "," << track.GetClusters().back().t << std::endl;

  for (const auto& cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    clusterenergy.Fill(cluster.energy);
    totalcluster_energy+=cluster.energy;
    //ClusterElectronsSize.Fill(cluster.electrons.size());
    NumElectronsPerCluster.Fill(cluster.electrons.size());
    num_electrons+=cluster.electrons.size();
    NumPhotonsPerCluster.Fill(cluster.photons.size());
    NumIonsPerCluster.Fill(cluster.ions.size());    
    //if(cluster.y<13 && cluster.y>-13){
      //std::cout << "clusters in center x,y,z,t, E: "  << cluster.x << "," << cluster.y << "," << cluster.z << "," << cluster.t << "," << cluster.energy << std::endl;
      //while loop on the clusters themselves
      //while (track.GetCluster(xcl, ycl, zcl, tcl, ncl, ecl, extra)){
      //while (track.GetCluster(xcl, ycl, zcl, tcl, nel,nph,nio, ecl, extra)){
        //for(int i = 0; i < cluster.electrons.size(); i++){
	      //  track.GetElectron(i,x,y,z,t0,e,dx,dy,dz);
        for(const auto& electron : cluster.electrons){
          if(electron.e>0){ // than this should be a delta electron...
            //std::cout << "electron: "  << electron.e << std::endl;
            ElectronKE.Fill(electron.e);
            num_deltas++;
          }
        }
      //for (const auto& electron : cluster.electrons) {
        //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
      //  ElectronKE.Fill(electron.e);
        //std::cout << "electron cluster size is: "  << cluster.electrons.size() << std::endl;
        //drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);
      //}  
      //}
    //}
  }
    cellView.Plot2d();
  constexpr bool twod = true;
  constexpr bool drawaxis = true;
  driftView.Plot(twod, drawaxis, true);
  std::cout << "Efield found: \n";
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  std::cout << "Average Electron KE: " << ElectronKE.GetMean() << " eV\n";
  std::cout << "Total cluster energy: " << totalcluster_energy << " eV\n";
  std::cout << "Total delta electrons in the track: " << num_deltas << "\n";
  std::cout << "Total electrons in the track: " << num_electrons << "\n";

  OutFile->cd();
  clusterenergy.Write();
  ElectronKE.Write();
  NumElectronsPerCluster.Write();
  NumPhotonsPerCluster.Write();
  NumIonsPerCluster.Write();
  ClusterElectronsSize.Write();
  OutFile->Close();
  app->Run(kTRUE);
  return 0;
}
