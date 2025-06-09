#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <ROOT/TThreadedObject.hxx>
#include <TThread.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
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
#include "RandomEngine.hh"
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
#include <random>
#include <thread>

using namespace Garfield;


int whichWireID_int(double ypos){
  return std::round((double)(((71)/2.0)-ypos/8.0));
}

double whichYPos(int thiswireID){
  return (double)((8.0 *((71)/2.0 - (double)(thiswireID))));
}

int main(int argc, char * argv[]) {
  ROOT::EnableThreadSafety();
  // TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
  const int wireIDOfInterest = invals[0];
  const double temperature   = invals[1];
  const double BFieldValue   = invals[2];
  std::stringstream ingasfilename;
  //ingasfilename << "/home/kmcbride/garfield/keiths_code/GarfSim/Garfield/FlightGasFiles/BOFF/Flight2024_Boff_P_755.865_T_" << temperature << ".15_multiE_90CO2_10Ar_01122024.gas";
  ingasfilename << "/home/kmcbride/garfield/keiths_code/GarfSim/Garfield/FlightGasFiles/BON/Flight2024_Bon_1Tesla_P_755.865_T_" << temperature << ".15_multiE_90CO2_10Ar_01122024.gas";  // /home/kmcbride/garfield/keiths_code/GarfSim/Garfield/FlightGasFiles/BON/Flight2024_Bon_1Tesla_P_755.865_T_279.15_multiE_90CO2_10Ar_01122024.gas
  MediumMagboltz * gas = new MediumMagboltz();
  gas->LoadGasFile(ingasfilename.str());
  // lets just print out the drift velocity to a file?

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


  auto start = std::chrono::high_resolution_clock::now();
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);

  std::stringstream wid;

  // now apply the y coordinates (which is HELIX z coordinates) correctly based on which wireID
  // maximum lorentz angle will be at ~2Tesla, so deflection of up to 18mm? if a sense wire is located at 4mm then need to go up and down by 20 mm to be safe
  // first find where the sense wire is
  const double YPosOfInterest = whichYPos(wireIDOfInterest)/10.0; // now in cm D'oh
  double zone_y_min=YPosOfInterest-2.0; // lowest for top
  double zone_y_max=YPosOfInterest+2.0; // highest for top
  if(zone_y_min < -30.0) zone_y_min=-29.6;
  if(zone_y_max > 29.8) zone_y_max=29.6;
  // if top 2 or bottom two wires, we need to expand the z-axis a little more
  if(wireIDOfInterest==0 || wireIDOfInterest==1) zone_y_min=zone_y_min-1.0; // need 1 more cm in vertical to accnt for Efield distortion
  if(wireIDOfInterest==70 || wireIDOfInterest==71) zone_y_max=zone_y_max+1.0; // need 1 more cm in vertical to accnt for Efield distortion

  std::cout << "zone_y_min and max are: " << zone_y_min << ", " << zone_y_max << std::endl;
  // for the drift time vs temp calculation of the DCT
  //const double max_y_i=28.8-rCathode, max_x_i=max_x_in, stepy=0.2,stepx=0.25, min_y_i=-28.8+rCathode;
  const double max_y_i=zone_y_max, max_x_i=7.6, stepy=0.1,stepx=0.2, min_y_i=zone_y_min;
  const double x_i=-7.6, y_i=zone_y_min, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)
  
  // correct the stepping and grid if we have negative numbers inputted for x?
  
  const double realstepx = stepx;
  // find the size of the vector we will be storing
  int size_of_vec_in_y=(int)((max_y_i-y_i)/stepy)+1;
  // for parallel, find the number of iterations we will be using to meet the x grid
  int size_of_x_grid=(int)((max_x_i-x_i)/realstepx)+1;
  std::cout << "size of x_grid will be: " << size_of_x_grid << std::endl;
  std::cout << "size of vectors will be: " << size_of_vec_in_y << std::endl;
  std::cout << "max_y_i minus y_i will be: " << max_y_i-y_i << std::endl;
  std::cout << "max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
  std::cout << "casted max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
   //<< "," << new_average << "," << getSigma(electrons_drifted) << "," << y_i << "," << z_i << "," << x_i << std::endl;
  std::cout << "wire id and yposition found are: " << wireIDOfInterest << " , " << YPosOfInterest << std::endl;
  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);
  int max_electrons=100;
  std::stringstream outrootfilename;
  outrootfilename << "BON_iso_persensewire/isochron_" << BFieldValue << "_" << wireIDOfInterest << "_" << temperature  << "_.root";
  TFile * Outfile = new TFile(outrootfilename.str().c_str(),"recreate");

  int pid = getpid();
  timeval t;
  gettimeofday(&t, NULL);
  int seed = pid*t.tv_usec;
  std::cout << "Random Seed: " << seed << std::endl;
  gRandom->TRandom::GetSeed();

  for (int iplane=0;iplane<1;iplane++){
    for (int iw =0;iw<7;iw++){
      int iadd = iplane*7 + iw;
      wid.str("");
      wid << "a_" <<iplane << "_" << iw ;
      //std::cout << iadd << " " << wid.str() << std::endl;
      std::string str(wid.str());
      const char * name = str.c_str();
    }
  }
  

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
//  cmp->SetMagneticField(0.,0.,0.0);
  cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,31,11);
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
      float sign = -1.0;
      if(iw%2==0) sign = 1.0;
      float x = sign*300e-4;
      //float x = 0;
      cmp->AddWire(x,y,2 * rAnode, vAnode, "a");
      //std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
      last_y_anode=y;
    }
    for(int iw=0;iw<nwire;iw++){
      float y = (starting_y-(potentialsep/2.0))-iw*potentialsep;
      float x = 0;
      cmp->AddWire(x,y,2 * rPotential, vPotential, "p");
      //std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
      last_y_potential=y;
    }
    // add in 3 potential wires at the bottom and top now separated by 4mm from eachother and by 4mm from the bottom or top wire
    cmp->AddWire(0,starting_y+(potentialsep/2.0),2 * rPotential, vPotential, "pT");
    cmp->AddWire(0,starting_y+(2.0*potentialsep/2.0),2 * rPotential, vAnode, "pT");
    cmp->AddWire(0,starting_y+(3.0*potentialsep/2.0),2 * rPotential, vPotential, "pT");
    //std::cout << " wire " << 0 << " " << starting_y+(potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << starting_y+(2.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << starting_y+(3.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
    cmp->AddWire(0,last_y_potential-(potentialsep/2.0),2 * rPotential, vAnode, "pT");
    cmp->AddWire(0,last_y_potential-(2.0*potentialsep/2.0),2 * rPotential, vPotential, "pT");
    cmp->AddWire(0,last_y_potential-(3.0*potentialsep/2.0),2 * rPotential, vAnode, "pT");
    //std::cout << " wire " << 0 << " " << last_y_potential-(potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << last_y_potential-(2.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << last_y_potential-(3.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
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
        //std::cout << " wire " << x_delt << " " << pad_r << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
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
  unsigned int ne=0, ni=0;
  //make some histograms by default for debugging at first
  Outfile->cd();
  TH1D * numE_thiswire_h = new TH1D("numEthiswire","numEthiswire",200,-0.5,199.5);
  TH1D * numdiff_h = new TH1D("numdiff","numdiff",100,-5,495);
  TH1D * xendpoint2_h = new TH1D("xendpoint2","xendpoint2",1540,-7.7,7.7);
  TH1D * yendpoint2_h = new TH1D("yendpoint2","yendpoint2",6000,-30,30);
  TH1D * zendpoint2_h = new TH1D("zendpoint2","zendpoint2",6000,-30,30);
  TH1D * xendpoint_h = new TH1D("xendpoint","xendpoint",1540,-7.7,7.7);
  TH1D * yendpoint_h = new TH1D("yendpoint","yendpoint",6000,-30,30);
  TH1D * zendpoint_h = new TH1D("zendpoint","zendpoint",6000,-30,30);
  TH1D * xdiff_h = new TH1D("xdiff","xdiff",1540,-7.7,7.7);
  TH1D * ydiff_h = new TH1D("ydiff","ydiff",6000,-30,30);
  TH1D * zdiff_h = new TH1D("zdiff","zdiff",6000,-30,30);

  TH1D * xendpoint2All_h = new TH1D("xendpoint2All","xendpoint2All",1540,-7.7,7.7);
  TH1D * yendpoint2All_h = new TH1D("yendpoint2All","yendpoint2All",6000,-30,30);
  TH1D * zendpoint2All_h = new TH1D("zendpoint2All","zendpoint2All",6000,-30,30);
  TH1D * xendpointAll_h = new TH1D("xendpointAll","xendpointAll",1540,-7.7,7.7);
  TH1D * yendpointAll_h = new TH1D("yendpointAll","yendpointAll",6000,-30,30);
  TH1D * zendpointAll_h = new TH1D("zendpointAll","zendpointAll",6000,-30,30);
  TH1D * xdiffAll_h = new TH1D("xdiffAll","xdiffAll",1540,-7.7,7.7);
  TH1D * ydiffAll_h = new TH1D("ydiffAll","ydiffAll",6000,-30,30);
  TH1D * zdiffAll_h = new TH1D("zdiffAll","zdiffAll",6000,-30,30);

  TH1D * xendpoint2Up_h = new TH1D("xendpoint2Up","xendpoint2Up",1540,-7.7,7.7);
  TH1D * yendpoint2Up_h = new TH1D("yendpoint2Up","yendpoint2Up",6000,-30,30);
  TH1D * zendpoint2Up_h = new TH1D("zendpoint2Up","zendpoint2Up",6000,-30,30);
  TH1D * xendpointUp_h = new TH1D("xendpointUp","xendpointUp",1540,-7.7,7.7);
  TH1D * yendpointUp_h = new TH1D("yendpointUp","yendpointUp",6000,-30,30);
  TH1D * zendpointUp_h = new TH1D("zendpointUp","zendpointUp",6000,-30,30);
  TH1D * xdiffUp_h = new TH1D("xdiffUp","xdiffUp",1540,-7.7,7.7);
  TH1D * ydiffUp_h = new TH1D("ydiffUp","ydiffUp",6000,-30,30);
  TH1D * zdiffUp_h = new TH1D("zdiffUp","zdiffUp",6000,-30,30);

  TH1D * xendpoint2Low_h = new TH1D("xendpoint2Low","xendpoint2Low",1540,-7.7,7.7);
  TH1D * yendpoint2Low_h = new TH1D("yendpoint2Low","yendpoint2Low",6000,-30,30);
  TH1D * zendpoint2Low_h = new TH1D("zendpoint2Low","zendpoint2Low",6000,-30,30);
  TH1D * xendpointLow_h = new TH1D("xendpointLow","xendpointLow",1540,-7.7,7.7);
  TH1D * yendpointLow_h = new TH1D("yendpointLow","yendpointLow",6000,-30,30);
  TH1D * zendpointLow_h = new TH1D("zendpointLow","zendpointLow",6000,-30,30);
  TH1D * xdiffLow_h = new TH1D("xdiffLow","xdiffLow",1540,-7.7,7.7);
  TH1D * ydiffLow_h = new TH1D("ydiffLow","ydiffLow",6000,-30,30);
  TH1D * zdiffLow_h = new TH1D("zdiffLow","zdiffLow",6000,-30,30);
  // want the bin centers to be such that steps in x and y (garfield coordinates) are placing drift locations simulated at the bin centers of x and y, so we need to shift the bins by half a step size in each direction
  // when constructing this histogram, careful with bin numbers and edges
  int nbinsx=std::round((max_x_i-x_i)/(stepx/2.0))+1;
  double x_low_edge=x_i-(stepx/4.0);
  double x_high_edge=max_x_i+(stepx/4.0);
  int nbinsy=std::round((max_y_i-y_i)/(stepy/2.0))+1;
  double y_low_edge=y_i-(stepy/4.0);
  double y_high_edge=max_y_i+(stepy/4.0);
  std::cout << "nbinsx, x_low_edge,x_high_edge     " << nbinsx << "," << x_low_edge << ", " << x_high_edge << std::endl;
  std::cout << "nbinsy, y_low_edge,y_high_edge     " << nbinsy << "," << y_low_edge << ", " << y_high_edge << std::endl;
  std::stringstream h3dnames;
  h3dnames.str( std::string());
  h3dnames.clear();
  // this is the wire of interest
  h3dnames << "sensewire_" << wireIDOfInterest;
  TH3D * iso_per_sensewire = new TH3D(h3dnames.str().c_str(),h3dnames.str().c_str(),1500,-0.5,14999.5,nbinsy,y_low_edge,y_high_edge,nbinsx,x_low_edge,x_high_edge);
  // now get the wire above it
  double UpYPosOfInterest = 1000.0; // now in cm D'oh
  if(wireIDOfInterest-1<0){
    UpYPosOfInterest=32.0;
  }
  else{
    UpYPosOfInterest = whichYPos(wireIDOfInterest-1)/10.0; // now in cm D'oh
  }
  //const double UpYPosOfInterest = whichYPos(wireIDOfInterest-1)/10.0; // now in cm D'oh
  const double zone_y_minUp=UpYPosOfInterest-2.0; // lowest for top
  const double zone_y_maxUp=UpYPosOfInterest+2.0; // highest for top
  int nbinsyUp=std::round((zone_y_maxUp-zone_y_minUp)/(stepy/2.0))+1;
  double y_low_edgeUp=zone_y_minUp-(stepy/4.0);
  double y_high_edgeUp=zone_y_maxUp+(stepy/4.0);
  h3dnames.str( std::string());
  h3dnames.clear();
  h3dnames << "sensewire_" << wireIDOfInterest-1; // subtract one from wireID means we get the wire above it (and if wireid was 0, we get -1 which means we fill this histogram with information about non-sense wire endpoints? interesting
  TH3D * iso_per_sensewireUp = new TH3D(h3dnames.str().c_str(),h3dnames.str().c_str(),1500,-0.5,14999.5,nbinsyUp,y_low_edgeUp,y_high_edgeUp,nbinsx,x_low_edge,x_high_edge);
  // now get the wire below it
  double LowYPosOfInterest = 1000.0; // now in cm D'oh
  if(wireIDOfInterest+1>71){
    LowYPosOfInterest=-32.0;
  }
  else{
    LowYPosOfInterest = whichYPos(wireIDOfInterest+1)/10.0; // now in cm D'oh
  }
  const double zone_y_minLow=LowYPosOfInterest-2.0; // lowest for top
  const double zone_y_maxLow=LowYPosOfInterest+2.0; // highest for top
  int nbinsyLow=std::round((zone_y_maxLow-zone_y_minLow)/(stepy/2.0))+1;
  double y_low_edgeLow=zone_y_minLow-(stepy/4.0);
  double y_high_edgeLow=zone_y_maxLow+(stepy/4.0);
  h3dnames.str( std::string());
  h3dnames.clear();
  h3dnames << "sensewire_" << wireIDOfInterest+1; // add one to wireID means we get the wire below it (and if wireid was >71, we get 72 which means we fill this histogram with information about non-sense wire endpoints? interesting
  TH3D * iso_per_sensewireLow = new TH3D(h3dnames.str().c_str(),h3dnames.str().c_str(),1500,-0.5,14999.5,nbinsyLow,y_low_edgeLow,y_high_edgeLow,nbinsx,x_low_edge,x_high_edge);
  // really the x_i here is 7.62 but also subtract off the diameter of the cathode wire (2*rCathode)
  //while(x_i<=max_x_i){
  #pragma omp parallel for collapse(2)
  //for(double x_delt=x_i; x_delt<=max_x_i;x_delt+=stepx){
  for(int ix=0;ix<size_of_x_grid;ix++){
    for(int iy=0;iy<size_of_vec_in_y;iy++){
      //double x_delt=x_i+((double)(ix)*stepx);
      double x_delt=x_i+((double)(ix)*realstepx);
      double y_delt=y_i+((double)(iy)*stepy);
      if(TMath::Abs(x_delt)<1.5*300e-4) continue;
      AvalancheMC * driftline = new AvalancheMC();
      //  driftline->EnableDebugging();
      //driftline->SetTimeSteps(step_size/1000.0);
      driftline->SetDistanceSteps(40.0/1000.0);  // 400um step size for outer region
      if(TMath::Abs(x_delt)<1.1){
        driftline->SetDistanceSteps(1.0/1000.0); // 10um for inner region
        std::cout << " step size will be 1 instead of 40" << std::endl;
      }
      //driftline->SetDistanceSteps(0.001);
      //driftline->EnableMagneticField();
      driftline->EnableDiffusion();
      driftline->SetSensor(sensor);
      //  driftline->EnablePlotting(vd);
      //  driftline->EnableSignalCalculation();
      unsigned int ne=0, ni=0;
      double num_electrons_y;
      double average_y, stddev_y, y_y, z_y, x_y, avg_x, avg_y, avg_z, std_x, std_y, std_z, e_a, i_a;
      int iter_y=0;
      int num_electrons=0;
      int num_electrons_thiswire=0;
      int num_electrons_thiswireUp=0;
      int num_electrons_thiswireLow=0;
      double min_variation=0.1;
      double running_sum=0;
      double running_average=0;
      double new_average=0;
      double curr_sig=0;
      bool keep_running=true;
      bool is_pos_borked=false;
      std::stringstream histonames;
      //histonames.str( std::string());
      //histonames.clear();
      //TH1D * drifttime_h = new TH1D(histonames.str().c_str(),histonames.str().c_str(),1500,-0.5,14999.5);
      //histonames.str( std::string());
      //histonames.clear();
      //histonames << "thiswires_num_electrons_" << y_delt*10.0 << "_y_" << x_delt*10.0;
      //histonames << "status_e_z_" << y_delt*10.0 << "_y_" << x_delt*10.0;
      //TH1D * stat_h = new TH1D(histonames.str().c_str(),histonames.str().c_str(),200,-0.5,199.5); // fill it with the difference between electrons that are drifted and don't end up on this wire, vs those that do end up on this wire (so always positive)
      //TH1D * stat_h = new TH1D(histonames.str().c_str(),histonames.str().c_str(),20,-19.5,0.5);
      //histonames.str( std::string());
      //histonames.clear();
      //histonames << "iso_per_wire_" << y_delt*10.0 << "_" << x_delt*10.0;
      //TH3D * iso_per_sensewire = new TH3D(histonames.str().c_str(),histonames.str().c_str(),1540,-7.7,7.7,6000,-30,30,1500,-0.5,14999.5);
      while(keep_running && !is_pos_borked){
        double xendpoint = 0, yendpoint = 0, zendpoint=0;
        double xendpoint2 = 0, yendpoint2 = 0, zendpoint2=0;
        double tendpoint = 0, tendpoint2 = 0;
        double energy0=0,energy1=0;
	      //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
        int stat=0;
	bool did_it_drift=driftline->DriftElectron(x_delt,y_delt,z_i,0);
        if(!did_it_drift) is_pos_borked=true;
	      //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
	      driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint,xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
        //std::cout << "stat is " << stat << std::endl;
        //if(ne_av>1) std::cout << "ne_av,ni_av is: " << ne_av << "," << ni_av << std::endl;
        //if(stat!=0) keep_running=false;
        curr_sig=tendpoint2-tendpoint;
        if(stat==-5){
          int which_wireID=whichWireID_int(yendpoint2*10.0);
          xendpointAll_h->Fill(xendpoint);
          xendpoint2All_h->Fill(xendpoint2);
          xdiffAll_h->Fill(xendpoint2-xendpoint);
          yendpointAll_h->Fill(yendpoint);
          yendpoint2All_h->Fill(yendpoint2);
          ydiffAll_h->Fill(yendpoint2-yendpoint);
          zendpointAll_h->Fill(zendpoint);
          zendpoint2All_h->Fill(zendpoint2);
          zdiffAll_h->Fill(zendpoint2-zendpoint);
          if(which_wireID==wireIDOfInterest){
            //iso_per_sensewire->Fill(xendpoint,yendpoint,curr_sig); // order is x, y, time
            iso_per_sensewire->Fill(curr_sig,yendpoint,xendpoint); // order is NOW time, y, x
            xendpoint_h->Fill(xendpoint);
            xendpoint2_h->Fill(xendpoint2);
            xdiff_h->Fill(xendpoint2-xendpoint);
            yendpoint_h->Fill(yendpoint);
            yendpoint2_h->Fill(yendpoint2);
            ydiff_h->Fill(yendpoint2-yendpoint);
            zendpoint_h->Fill(zendpoint);
            zendpoint2_h->Fill(zendpoint2);
            zdiff_h->Fill(zendpoint2-zendpoint);
            num_electrons_thiswire++;
            running_sum+=curr_sig;
            new_average=running_sum/num_electrons_thiswire;
          }
          else if(which_wireID==wireIDOfInterest-1){
            iso_per_sensewireUp->Fill(curr_sig,yendpoint,xendpoint); // order is NOW time, y, x
            xendpointUp_h->Fill(xendpoint);
            xendpoint2Up_h->Fill(xendpoint2);
            xdiffUp_h->Fill(xendpoint2-xendpoint);
            yendpointUp_h->Fill(yendpoint);
            yendpoint2Up_h->Fill(yendpoint2);
            ydiffUp_h->Fill(yendpoint2-yendpoint);
            zendpointUp_h->Fill(zendpoint);
            zendpoint2Up_h->Fill(zendpoint2);
            zdiffUp_h->Fill(zendpoint2-zendpoint);
            num_electrons_thiswireUp++;
          }
          else if(which_wireID==wireIDOfInterest+1){
            iso_per_sensewireLow->Fill(curr_sig,yendpoint,xendpoint); // order is NOW time, y, x
            xendpointLow_h->Fill(xendpoint);
            xendpoint2Low_h->Fill(xendpoint2);
            xdiffLow_h->Fill(xendpoint2-xendpoint);
            yendpointLow_h->Fill(yendpoint);
            yendpoint2Low_h->Fill(yendpoint2);
            ydiffLow_h->Fill(yendpoint2-yendpoint);
            zendpointLow_h->Fill(zendpoint);
            zendpoint2Low_h->Fill(zendpoint2);
            zdiffLow_h->Fill(zendpoint2-zendpoint);
            num_electrons_thiswireLow++;
          }
          //drifttime_h->Fill(curr_sig);
          //stat_h->Fill(stat);
          //iso_per_sensewire->Fill(xendpoint,yendpoint,curr_sig);
          num_electrons++;
        }
        //}
        if(num_electrons%30==0){
          std::cout << "num_electrons: " << num_electrons << " simulated for x=" << x_delt << " y=" << y_delt << std::endl;
        }
        if(num_electrons_thiswire%30==0){
          std::cout << "num_electrons_thiswire: " << num_electrons_thiswire << " simulated for x=" << x_delt << " y=" << y_delt << std::endl;
        }
        if(num_electrons_thiswireUp%30==0){
          std::cout << "num_electrons_thiswireUp: " << num_electrons_thiswireUp << " simulated for x=" << x_delt << " y=" << y_delt << std::endl;
        }
        if(num_electrons_thiswireLow%30==0){
          std::cout << "num_electrons_thiswireLow: " << num_electrons_thiswireLow << " simulated for x=" << x_delt << " y=" << y_delt << std::endl;
        }
        if(num_electrons_thiswire>max_electrons){
          if(TMath::Abs((new_average-running_average)/running_average)<=min_variation){
            keep_running=false; // successful checkout of this wireIDs real drift position acceptance
          }
          else is_pos_borked=true;
        }
        if(num_electrons_thiswireUp>max_electrons) keep_running=false;
        if(num_electrons_thiswireLow>max_electrons) keep_running=false;
        if(num_electrons>max_electrons && num_electrons_thiswire<10 && num_electrons_thiswireUp<10 && num_electrons_thiswireLow<10) is_pos_borked=true;
        else{
          // compute the average some more
          running_average=new_average;
        }
      }// end of while loop? now we write to file?

      // write these out to vectors to spit out to file at the end
      //outputfile << num_electrons << "," << new_average << "," << getSigma(electrons_drifted) << "," << y_delt << "," << z_i << "," << x_delt << std::endl;
      if(is_pos_borked){
        std::cout << "x,y borked is " << x_delt << "," << y_delt << std::endl;
        //std::cout << "ypos borked is " << iter_y << std::endl;
        iter_y++;
      }
      else{
        num_electrons_y=num_electrons;
        average_y=new_average;
        y_y=y_delt;
        z_y=z_i;
        x_y=x_delt;
        iter_y++;
      }
      std::cout << "x = " << ix << ", y= " << iy << ", threadId = "<< omp_get_thread_num() << std::endl; //i, j, omp_get_thread_num());
      //y_i=min_y_i;
      Outfile->cd();
      //int checking_writing = drifttime_h->Write();
      //if(checking_writing<=0){
      //  std::cout << "FAILED TO WRITE, error: " << checking_writing << ", and means and RMS of the histo are " << drifttime_h->GetMean() << " , " << drifttime_h->GetRMS() << " \n";
      //}
      //int checking_writing_2 = iso_per_sensewire->Write();
      //if(checking_writing_2<=0){
      //  std::cout << "FAILED TO WRITE, error: " << checking_writing_2 << " \n";
      //}
      //stat_h->Write();
      numdiff_h->Fill(num_electrons-num_electrons_thiswire);
      numE_thiswire_h->Fill(num_electrons_thiswire);
      delete driftline;
      //delete drifttime_h;
      //delete stat_h;
      //delete iso_per_sensewire;
      //delete aval;
    }// end of for loop y. move file writing inside
  }// end of for loop x
  //std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " RMS is : " << getSigma(electrons_drifted) << std::endl;
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  std::cout << "Closing file and finished run" << std::endl;
  Outfile->cd();
  /*for(int hi=0;hi<72;hi++){
    //iso_per_sensewire[hi]->Write();
    auto newh=iso_per_sensewire[hi]->Merge();
    newh->Write();
  }*/
    //auto newh=iso_per_sensewire_h.Merge();
    //newh->Write();
  numE_thiswire_h->Write();
  numdiff_h->Write();
  xendpoint2_h->Write();
  yendpoint2_h->Write();
  zendpoint2_h->Write();
  xendpoint_h->Write();
  yendpoint_h->Write();
  zendpoint_h->Write();
  xdiff_h->Write();
  ydiff_h->Write();
  zdiff_h->Write();

  xendpoint2Up_h->Write();
  yendpoint2Up_h->Write();
  zendpoint2Up_h->Write();
  xendpointUp_h->Write();
  yendpointUp_h->Write();
  zendpointUp_h->Write();
  xdiffUp_h->Write();
  ydiffUp_h->Write();
  zdiffUp_h->Write();

  xendpoint2Low_h->Write();
  yendpoint2Low_h->Write();
  zendpoint2Low_h->Write();
  xendpointLow_h->Write();
  yendpointLow_h->Write();
  zendpointLow_h->Write();
  xdiffLow_h->Write();
  ydiffLow_h->Write();
  zdiffLow_h->Write();

  xendpoint2All_h->Write();
  yendpoint2All_h->Write();
  zendpoint2All_h->Write();
  xendpointAll_h->Write();
  yendpointAll_h->Write();
  zendpointAll_h->Write();
  xdiffAll_h->Write();
  ydiffAll_h->Write();
  zdiffAll_h->Write();
  iso_per_sensewire->Write();
  iso_per_sensewireUp->Write();
  iso_per_sensewireLow->Write();
  Outfile->Close();
  return 0;
}
