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
#include <TGraph.h>
#include <TH1.h>
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentGrid.hh"
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

using namespace Garfield;
 
// removed
//#include "Garfield/Plotting.hh"

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
  const double thistilt = invals[0]; // HELIX zposition slice to plot Efield in millimeters
  std::cout << "Tilt will be : " << thistilt << "um" <<  std::endl;
  const double thisOldWay = invals[1]; // HELIX zposition slice to plot Efield in millimeters
  bool input_oldWay = (thisOldWay >= 0) ? false : true; // if negative, we will do oldway
  std::cout << "Old way : " << input_oldWay <<  std::endl;
  const int real_in_factor = invals[2]; // HELIX zposition slice to plot Efield in millimeters
  int in_factor = 1;
  if(real_in_factor >1) in_factor = real_in_factor;
  std::cout << "Grid Factor multiplication " << in_factor <<  std::endl;


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
      std::cout << iadd << " " << wid.str() << std::endl;
      std::string str(wid.str());
      const char * name = str.c_str();
    }
  }

  //Garfield::plottingEngine.SetDefaultStyle();
  MediumMagboltz * gas = new MediumMagboltz();

  // Setup the gas
  const double pressure = 760.; //Torr
  const double temperature = 293.15; //K
 
  // Set the temperature [K] and pressure [Torr]
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("co2", 85, "ar", 15);

  gas->LoadGasFile("/home/kmcbride/garfield/keiths_code/GarfSim/Garfield/FlightGasFiles/BOFF/Flight2024_Boff_P_755.865_T_299.15_multiE_90CO2_10Ar_01122024.gas");
  //gas->LoadGasFile("Flight2024_Boff_P_755.865_T_299.15_logE.gas");
  //gas->LoadGasFile("/home/kmcbride/garfield/keiths_code/GarfSim/Garfield/FlightGasFiles/BON/Flight2024_Bon_1Tesla_P_755.865_T_299.15_multiE_90CO2_10Ar_01122024.gas");

  // lets just print out the drift velocity to a file?

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
   //   SolidBox * enclosure = new SolidBox(0,0,0,10,31.,11);
  SolidBox * enclosure = new SolidBox(0,0,0,20,31.0,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);
  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
  //const double rAnode= 40e-4;
//  const double rAnode= 175e-4;
  const double vPotential= -2700;
//  const double rPotential= 175e-4;
  //const double rPotential= 350e-4;
  const double rPotential= 250e-4;
  
  const double anodesep = 0.8;
  const double potentialsep = 0.8;
  const double cathodesep =0.4;
  //const double cathodesep =0.2;
  //const double cathodesep =0.1;
  const double fieldsep = 0.2;
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
      //wid << "a_" << iplane << "_" << iw;      
      float y = starting_y-iw*anodesep;
      float sign = -1.0;
      if(iw%2==0) sign = 1.0;
      float x = sign*300e-4;
      //float x = 0;
      cmp->AddWire(x,y,rAnode, vAnode, "a");
      std::cout << " wire " << x << " " << y << " " << vAnode << " " << "a" << std::endl;
      last_y_anode=y;
    }
    for(int iw=0;iw<nwire;iw++){
      float y = (starting_y-(potentialsep/2.0))-iw*potentialsep;
      float x = 0;
      //float x = 100e-4;
      cmp->AddWire(x,y, rPotential, vPotential, "p");
      std::cout << " wire " << x << " " << y << " " << vPotential << " " << "p" << std::endl;
      last_y_potential=y;
    }
    // add in 3 potential wires at the bottom and top now separated by 4mm from eachother and by 4mm from the bottom or top wire
    cmp->AddWire(0,starting_y+(potentialsep/2.0),  rPotential, vPotential, "pT");
    cmp->AddWire(0,starting_y+(2.0*potentialsep/2.0),  rPotential, vAnode, "pT");
    cmp->AddWire(0,starting_y+(3.0*potentialsep/2.0),  rPotential, vPotential, "pT");
    //std::cout << " wire " << 0 << " " << starting_y+(potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << starting_y+(2.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << starting_y+(3.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;

    cmp->AddWire(0,last_y_potential-(potentialsep/2.0),  rPotential, vAnode, "pT");
    cmp->AddWire(0,last_y_potential-(2.0*potentialsep/2.0),  rPotential, vPotential, "pT");
    cmp->AddWire(0,last_y_potential-(3.0*potentialsep/2.0),  rPotential, vAnode, "pT");
    //std::cout << " wire " << 0 << " " << last_y_potential-(potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << last_y_potential-(2.0*potentialsep/2.0) << " " << vPotential << " " << "pT" << std::endl;
    //std::cout << " wire " << 0 << " " << last_y_potential-(3.0*potentialsep/2.0) << " " << vAnode << " " << "pT" << std::endl;

  }
  double the_final_distance=7.62;
  bool oldWay = true;
  bool gaussian=false;
  if(oldWay){
    cmp->AddPlaneX(-1.0 * the_final_distance,vCathode,"cP1"); 
    cmp->AddPlaneX(the_final_distance,vCathode,"cP2");
  }
  else{
    // need to do cmp->AddWire() like above
    // but now for 2 mil diameter (50um) on 3 mil centers/pitch (75um).
    double cathode_start  = starting_y+(4.0*potentialsep/2.0); // we will need to address this probably
    double cathode_mesh_sep = 0.0075;  // 75 microns which is 0.075 mm 
    int nmesh = 8000; // crazy but maybe correct...
    double cathode_mesh_diameter = 0.005; // 50 microns which is 0.05 mm 
    for(int iw=0;iw<nmesh+1;iw++){
      float y = cathode_start-(iw*cathode_mesh_sep);
      float xLeft = -1.0 * the_final_distance; // this is the leftmost edge of the chamber, we will approximate the cathode plane with a bunch of wires, we can adjust the diameter and spacing to get a good approximation to a plane
      float xRight = the_final_distance; // this is the rightmost edge of the chamber, we will approximate the cathode plane with a bunch of wires, we can adjust the diameter and spacing to get a good approximation to a plane
      // turning on a deflection?
      // Parameters for deflection
      if(gaussian){
        double peakDeflection = xLeft * 0.03; // 1% deflection
        double yMid = 0.0;            // Center of the range
        double sigma = cathode_start / 2.0;           // Spread (adjust to taste)
        // Apply smooth deflection
        double exponent = -std::pow(y - yMid, 2) / (2.0 * std::pow(sigma, 2));
        xLeft += peakDeflection * std::exp(exponent);
      }
      //else if (!gaussian){ // linear tilt? symmetric about 0
      //  double tilt_amnt=thistilt/10000.0; // in cm from microns
      //  xLeft += ((tilt_amnt)/(cathode_start)*y);
      //  xRight += ((tilt_amnt)/(cathode_start)*y);
      //}
      cmp->AddWire(xLeft,y,cathode_mesh_diameter, vCathode, "cP1");
      cmp->AddWire(xRight,y,cathode_mesh_diameter, vCathode, "cP2");
      //std::cout << " wire " << xLeft << " " << y << " " << vCathode << " " << "cP1" << std::endl;
      //std::cout << " wire " << xRight << " " << y << " " << vCathode << " " << "cP2" << std::endl;
    }
  }



  float v=0;
  float tweak=0;
  float drift_dist_max=the_final_distance*1.0; // adjusted
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
      //std::cout << " wire " << x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;
      //std::cout << " wire " << left_x << " " << y << " " << v << " " << "T" << " count: " << count_me  << "out of " << nstrips << std::endl;

      count_me++;
      cmp->AddWire(x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
      cmp->AddWire(left_x,y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
      // lets try adding in the other side (subtract off two of the drift_dist_max from the left_x and add on two drift_dist_max to the x)
      cmp->AddWire(x+(2.0*drift_dist_max),y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"
      cmp->AddWire(left_x-(2.0*drift_dist_max),y,strip_size, v, "T"); // arguments are "xloc, yloc, diameter, voltage, label"


    }
  }
  
  std::cout << " strip_sep is: "  << strip_sep  << std::endl;
  std::cout << " estimated diameter should be "  <<  2.0*drift_dist_max/(double) (2.0*(nstrips_per_PCB-1)) << std::endl;
  std::cout << " diameter of cathode is: "  << 2.0*rCathode << std::endl;
  std::cout << " strip_size is: "  << strip_size  << std::endl;

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
  // plotting helix E_y (E_x here) versus helix z (which is y here) at different Helix y values (x positions here) to see if the plane compared to the mesh matters

  ComponentGrid * cmpE = new ComponentGrid();
  cmpE->SetGeometry(geo);
  std::stringstream elecfilename;
  //elecfilename << "generated_electric_field_grid.txt";
  elecfilename << "generated_electric_field_grid_" << in_factor << ".txt"; // generated_electric_field_grid_1.txt
  cmpE->LoadElectricField(elecfilename.str(), "XYZ",false,false);

  std::stringstream TF1name; // for accessing FullBField files
  TFile * thisguy; // output file
  TF1name.str("");
  // BFieldIsoMaps_wireID_1.root
  TF1name << "drift_compareGrid_andAnalytic_" << in_factor << ".root";

  std::cout << "writing to file: " << TF1name.str().c_str() << std::endl;
  thisguy = TFile::Open(TF1name.str().c_str(),"RECREATE");
  if(!thisguy || !thisguy->IsOpen() || thisguy->IsZombie()){
    return -1;
  }
  thisguy->cd();
  double y_value_to_plot = 6.0;
  while(y_value_to_plot>4.0){
    double x_min = 0.0;
    double x_max = 7.6;
    double x_step=0.0025; // 25um steps? // this will be 10k points?
    TGraph * thisgEy = new TGraph();
    TGraph * thisgVy = new TGraph();
    std::stringstream thisname;
    thisname << "Ey_vs_y_" << std::fixed << std::setprecision(0) << y_value_to_plot*10.0;
    thisgEy->SetName(thisname.str().c_str());
    thisgEy->SetTitle("Diff Drift Electric Field vs y position; y(mm);Ey(V/cm)");
    thisname.str("");
    thisname << "Vy_vs_y_" << std::fixed << std::setprecision(0) << y_value_to_plot*10.0;
    thisgVy->SetName(thisname.str().c_str());
    thisgVy->SetTitle("Diff Drift Velocity vs y position; y(mm);percent diff in v");
    while(x_min<x_max+x_step/2.0){
      Medium * medium = nullptr;
      double ex = 0., ey = 0., ez = 0.;
      int status;
      sensor->ElectricField(x_min,y_value_to_plot,0,ex,ey,ez,medium,status);
      double ex_g = 0., ey_g = 0., ez_g = 0.;
      int status_g;
      cmpE->ElectricField(x_min,y_value_to_plot,0,ex_g,ey_g,ez_g,medium,status_g);
      if(ex<1.0 && ex>-1.0){
        std::cout << "whyyy" << std::endl;
      }
      thisgEy->AddPoint(x_min*10.0,ex-ex_g);
      if(status!=0 || x_min==0.03){
        std::cout << "oof " << status << "at x step " << x_min << std::endl;
              x_min+=x_step;
        continue;
      }
      // get drift velocity
      double vx = 0., vy = 0., vz = 0.;
      double vx_g = 0., vy_g = 0., vz_g = 0.;
      gas->ElectronVelocity(ex, ey, ez, 0, 0, 0, vx, vy, vz);
      gas->ElectronVelocity(ex_g, ey_g, ez_g, 0, 0, 0, vx_g, vy_g, vz_g);
      
      thisgVy->AddPoint(x_min*10.0,TMath::Abs(((vx*1000.0) - (vx_g*1000.0))/(vx*1000.0)*100.0));
      x_min+=x_step;
    }
    thisgEy->Write();
    thisgVy->Write();
    delete thisgEy;
    delete thisgVy;
    y_value_to_plot-=0.8;
  }
  std::cout << "Efields found: \n";
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  thisguy->Close();
  //app->Run(kTRUE);
  return 0;
}
