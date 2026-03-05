#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TLegend.h>
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
#include "Garfield/ViewMedium.hh"
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
#include <cstdlib> // For system()

using namespace Garfield;

#define numhists 20


void extractFileFromTarGz(const std::string& archivePath, const std::string& fileNameToExtract, const std::string& outputPath) {
    // Construct the command to extract the specific file
    std::string command = "tar -xzvf " + archivePath + " " + " -C " + outputPath + " " + fileNameToExtract; //tar -xzvf HEATModelGarfieldFiles.tar.gz -C /home/kmcbride/HEATModel_files_for$    // Execute the command
    int result = std::system(command.c_str());

    if (result == 0) {
        std::cout << "Successfully extracted '" << fileNameToExtract << "' to '" << outputPath << "'." << std::endl;
    } else {
        std::cerr << "Error extracting file: " << result << std::endl;
    }
}

double getMean(const std::vector<int> &input){
  double sum = std::accumulate(input.begin(), input.end(), 0.0);
  double mean = sum / input.size();
  return mean;
}

double getSigma(const std::vector<int> &input){
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
  //TRint* app = new TRint("Garfield", &argc, argv, 0, 0);

  std::stringstream wid;
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
 // const int wireIDOfInterest = invals[0];
  const double temperature   = invals[0];
  const double BFieldValue   = invals[1];
  const double input_tstep   = invals[2];
  const int which_column     = invals[3];
  const double x_slice     = invals[4];
  std::cout << "Will simulate gas gain" << std::endl;
  //
  int pid = getpid();
  timeval t;
  gettimeofday(&t, NULL);
  int seed = pid*t.tv_usec;
  std::cout << "Random Seed: " << seed << std::endl;
  gRandom->TRandom::GetSeed();

// gasfile stuff
//std::cout << "RSS before gasfile allocation: " << getCurrentRSS() / 1024 << " KB" << std::endl;
  MediumMagboltz * gas = new MediumMagboltz();
  // Setup the gas
  //const double temperature=299.15;
  const double pressure = 14.616*51.7149; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr) 
  gas->SetTemperature(temperature); // from CLI
  gas->SetPressure(pressure);
  gas->SetComposition("CO2", 90.,"AR", 10.);
  const int nFields = 5;
  const double E_not = 984.25;
  const double emin = E_not-E_not;
  const double emax = E_not+E_not;
  // Flag to request logarithmic spacing.
  const bool useLog = false;
  const double bmin=0;
  const double bmax=2; // do we need magnetic field on?
  const int nBFields=4;
  gas->SetFieldGrid(emin, emax, nFields, useLog, bmin,bmax,nBFields,TMath::Pi()/2.0,TMath::Pi()/2.0,1); 
  // Turn on penning transfer?
  gas->EnablePenningTransfer();
  gas->SetMaxElectronEnergy(200);
  std::cout << "number of levels: " << gas->GetNumberOfLevels();
  const int ncoll = 5;
  gas->GenerateGasTable(ncoll);
  // lets just print out the drift velocity to a file?
  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();
  //std::cout << "RSS after gasfile allocation: " << getCurrentRSS() / 1024 << " KB" << std::endl;
  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  //cmp->SetMagneticField(0.,0.,BFieldValue);
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
      if(iw==36) cmp->AddWire(x,y,2 * rAnode, vAnode, "asig");
      else cmp->AddWire(x,y,2 * rAnode, vAnode, "a");
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
  }
  
  std::cout << " strip_sep is: "  << strip_sep  << std::endl;
  std::cout << " estimated diameter should be "  <<  2.0*drift_dist_max/(double) (2.0*(nstrips_per_PCB-1)) << std::endl;
  std::cout << " diameter of cathode is: "  << 2.0*rCathode << std::endl;
  std::cout << " strip_size is: "  << strip_size  << std::endl;

// geometry detector wires

  std::cout << " strip_sep is: "  << strip_sep  << std::endl;
  std::cout << " estimated diameter should be "  <<  2.0*drift_dist_max/(double) (2.0*(nstrips_per_PCB-1)) << std::endl;
  std::cout << " diameter of cathode is: "  << 2.0*rCathode << std::endl;
  std::cout << " strip_size is: "  << strip_size  << std::endl;

  // need to tar extract the file so we dont need all of them unloaded:
  std::stringstream magfilename;
  magfilename << "HEATModelForGarfield/HEATModel_xslice_";
  magfilename << std::fixed << std::setprecision(0) << x_slice << "_column_" << which_column << ".csv";

  std::string archive = "HEATModelGarfieldFiles.tar.gz";
  //std::string fileToExtract = "HEATModelForGarfield/HEATModel_xslice_0_column_1.csv"; // Path within the archive // looks like HEATModelForGarfield/HEATModel_xslice_0_column_0.csv
  std::string outputDir = "/home/kmcbride/garfield/isoHEATB_codes/GarfSim/Garfield/HEATModel_files_for_garfield";
  extractFileFromTarGz(archive, magfilename.str(), outputDir);
  // Now you can read the extracted file:
  std::string extractedFilePath = outputDir + "/" + magfilename.str(); // Adjust if file path within archive differs from output path

  // Load the field map.
  ComponentGrid * cmpB = new ComponentGrid();
  cmpB->SetGeometry(geo);
  //cmpB->LoadMagneticField("garfield_HEAT_example_v2.csv", "XYZ"); // come up with a file that has x,y,z in cm and bx,by,bz in Tesla
  if(which_column==0 || which_column==1 || which_column==-1){
    //cmpB->LoadMagneticField(magfilename.str() , "XYZ");
    cmpB->LoadMagneticField(extractedFilePath , "XYZ");
  }
  //cmpB->LoadMagneticField("garfield_HEAT_example_v2.csv", "XYZ");
  //else if(which_column==-1) cmpB->LoadMagneticField("garfield_left_column.csv", "XYZ");
  //else if(which_column==1) cmpB->LoadMagneticField("garfield_right_column.csv", "XYZ");
  //cmpB->LoadMagneticField("/home/kmcbride/master/08212025/helix-tools/00build/HEATModel_xslice_-200_column_0.csv" , "XYZ");
  //else if(which_column==3) cmpB->LoadMagneticField("/home/kmcbride/master/08212025/helix-tools/00build/HEATModel_xslice_200_column_0.csv" , "XYZ");
  //else if(which_column==4) cmpB->LoadMagneticField("/home/kmcbride/master/08212025/helix-tools/00build/HEATModel_xslice_-100_column_0.csv" , "XYZ");
  else cmpB->LoadMagneticField("garfield_HEAT_example_v2.csv", "XYZ"); // come up with a file that has x,y,z in cm and bx,by,bz in Tesla



  Sensor * sensor = new Sensor;
  sensor->AddComponent(cmpB);
  sensor->AddComponent(cmp);


  sensor->SetTimeWindow(0,2,20000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("asig");
  sensor->AddElectrode(cmp,"asig");
  sensor->EnableComponent(0, true);
  sensor->EnableComponent(1, true);
  AvalancheMicroscopic * aval = new AvalancheMicroscopic();
  aval->EnableSignalCalculation(false);
  //aval->EnableMagneticField();
  //driftline->EnableDiffusion();
  //driftline->EnableAttachment();
  aval->SetSensor(sensor);
  std::cout << "Avalanch size limit is: " << (int) (aval->GetAvalancheSizeLimit()) << std::endl;
  aval->EnableAvalancheSizeLimit(100000);
  std::cout << "Avalanch size limit after setting is: " << (int) (aval->GetAvalancheSizeLimit()) << std::endl;
  //bool did_it_drift=driftline->DriftElectron(x_delt,y_delt,z_i,0);
  std::vector<int> electrons_drifted_aval;
  std::vector<int> ions_drifted_aval;
  for(int i=0;i<100;i++){
    double xendpoint, yendpoint, zendpoint, tendpoint, energy0,xendpoint2, yendpoint2, zendpoint2, tendpoint2, energy1;
    int stat;
    int ne_av,ni_av;
    bool did_it_drift=aval->AvalancheElectron(1,1,0,0,0,0,0,0);
  
    //std::cout << __LINE__ << std::endl;
    //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
    aval->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint, energy0,xendpoint2, yendpoint2, zendpoint2, tendpoint2, energy1, stat);
    aval->GetAvalancheSize(ne_av,ni_av);
    electrons_drifted_aval.push_back(ne_av);
    ions_drifted_aval.push_back(ni_av);
    //std::cout << "Avalanch : " << ne_av << "," << ni_av << std::endl;

  }
    std::cout << "Avalanch : " << getMean(electrons_drifted_aval) << "," << getMean(ions_drifted_aval) << std::endl;

      TCanvas * cD = new TCanvas("cD", "", 600, 600);

  //aval->EnablePathLengthComputation();
  //aval->EnableAvalancheSizeLimit(1000);
  //  driftline->EnablePlotting(vd);
  //  driftline->EnableSignalCalculation();
 return 0;
}
