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
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackSimple.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Random.hh"
#include <TROOT.h>
#include <TRint.h>
#include "DCsim.hh"
#include <fstream>
#include "TMath.h"
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include "Garfield/ViewIsochrons.hh"
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
  TApplication app("app", &argc, argv);


    std::stringstream wid;
  MediumMagboltz* gas = new MediumMagboltz();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();

  //cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,31.,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
//  const double rAnode= 175e-4;
  const double vPotential= -2700;
//  const double rPotential= 175e-4;
  const double rPotential= 175e-4;
  
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;
  //const double cathodesep =0.2;
  //const double cathodesep =0.1;
  const double fieldsep = 0.2;

  /* 
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
*/
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
  double size_tweak=0.9;
  double strip_size=size_tweak*strip_sep;
  float starting_y_strips=starting_y+(4.0*potentialsep/2.0);
  float last_y_strips=last_y_potential-(4.0*potentialsep/2.0);
  float starting_x_strips=drift_dist_max-(1.5*strip_sep); // one pad away is where the first division should start
  float voltage_step=(vAnode-vCathode)/(nstrips_per_PCB+1);
  // 2 planes because top and bottom of chamber
  // starting_x
  for (int iplane=0;iplane<2;iplane++){
    for(int ipad=0;ipad<nstrips_per_PCB;ipad++){
      float y = starting_y_strips; // top or btm
      if(iplane >0) y=last_y_strips;
      float x = starting_x_strips-((double)(ipad)*strip_sep); // left or right
      float left_x = (-1.0*starting_x_strips)+((double)(ipad)*strip_sep); // left or right
	    //v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vPotential-vCathode)*(drift_dist_max-x)/drift_dist_max)  ;
	    //v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vAnode-vCathode)*(drift_dist_max-x)/drift_dist_max);
      v=vCathode +  voltage_step*(ipad+1);
      std::cout << " wire " << x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
      std::cout << " wire " << left_x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;

      count_me++;
      cmp->AddWire(x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
      cmp->AddWire(left_x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"


    }
    /*for(int iw=0;iw<nstrips+1;iw++){
      float y = starting_y_strips-(2.0*last_y_strips)*iplane;
      float x = drift_dist_max-(iw)*strip_sep;
      tweak = 0;
      if(x>=0){ 
	      //v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vPotential-vCathode)*(drift_dist_max-x)/drift_dist_max)  ;
	      v= (1+tweak*(drift_dist_max-x)/drift_dist_max )*( vCathode +  (vAnode-vCathode)*(drift_dist_max-x)/drift_dist_max);
      }
      else{
	      //v=(1+tweak*(drift_dist_max+x)/drift_dist_max)*(vPotential + (vCathode-vPotential)*(-x)/drift_dist_max);
	      v=(1+tweak*(drift_dist_max+x)/drift_dist_max)*(vAnode + (vCathode-vAnode)*(-x)/drift_dist_max);
	    }
      std::cout << " wire " << x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
      count_me++;
      cmp->AddWire(x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
    }*/
  }
  
  std::cout << " strip_sep is: "  << strip_sep  << std::endl;
  std::cout << " estimated diameter should be "  <<  2.0*drift_dist_max/(double) (2.0*(nstrips_per_PCB-1)) << std::endl;
  std::cout << " diameter of cathode is: "  << 2.0*rCathode << std::endl;
  std::cout << " strip_size is: "  << strip_size  << std::endl;

  // Load the field map.
  ComponentGrid * cmpB = new ComponentGrid();
  //cmp.SetCylindricalCoordinates(); // /c/Users/keith/HELIX/thermal_chamber_temps
  //cmp.LoadMagneticField("solenoid.txt", "XZ"); // come up with a file that has x,y,z in cm and bx,by,bz in Tesla
  cmpB->SetGeometry(geo);
  cmpB->LoadMagneticField("/home/kmcbride/master/08212025/helix-tools/00build/HEATModel_xslice_-200_column_0.csv", "XYZ"); // come up with a file that has x,y,z in cm and bx,by,bz in Tesla

  Sensor * sensor = new Sensor;
  //ViewSignal * vs1 = new ViewSignal;

  sensor->AddComponent(cmp);
  sensor->AddComponent(cmpB);
  sensor->SetTimeWindow(0,2,20000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");

  ViewField view;
  view.SetSensor(sensor);
  view.SetMagneticFieldRange(-2, 2);
  view.PlotProfile(-7, 0, 0., 7, 0., 0., "bz"); 
  app.Run(true);
}
