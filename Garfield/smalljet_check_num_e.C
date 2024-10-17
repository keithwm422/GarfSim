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
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  TH1D *SensewireSig[180];

  TH1F ncluster("ncluster","ncluster",100,0.,2000.);
  TH1F nclusterused("nclusterused","nclusterused",100,0.,2000.);

  TH1F clustersize("clustersize","clustersize",1000,0.,100.);
  TH1F clustersizeused("clustersizeused","clustersizeused",1000,0.,100.);

  std::stringstream wid;
  double max_y_i=0.5, max_z_i=0.5, stepy=0.05,stepz=0.05, min_y_i=-0.5;
  double x_i=7.60249, y_i=min_y_i, z_i=-0.5, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0
  // find the size of the vector we will be storing
  int size_of_vec_in_y=(int)((max_y_i-y_i)/stepy);
  std::cout << "size of vectors will be: " << size_of_vec_in_y << std::endl;
  std::cout << "max_y_i minus y_i will be: " << max_y_i-y_i << std::endl;
  std::cout << "max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
  std::cout << "casted max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
   //<< "," << new_average << "," << getSigma(electrons_drifted) << "," << y_i << "," << z_i << "," << x_i << std::endl;

  std::ofstream outputfile("electron_drifttimes_coordinates_BOFF_with_distances_small_1000e_270P.txt");
  outputfile << "num_electrons,average,stddev,y,z,x,dx,dy,dz,sigx,sigy,sigz" << std::endl;
  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);
  outputfile << std::fixed << std::showpoint;
  outputfile << std::setprecision(10);
  //name for drifttimes vs drift distance
  std::ofstream statfile("electron_status_values_270P.txt");
  std::ofstream numfile("electron_num_endpoints.txt");


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
  gas->LoadGasFile("Flight2024_Boff_P_755.865_T_298.15_.gas");

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
  const double vPotential= -270;
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
  ViewSignal * vs1 = new ViewSignal;

  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,10000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");
  vs1->SetSensor(sensor);
  
 
  AvalancheMC * driftline = new AvalancheMC();
  //  driftline->EnableDebugging();  
  driftline->SetDistanceSteps(0.001);
  //driftline->EnableMagneticField();
  driftline->EnableDiffusion();
  driftline->SetSensor(sensor);
  //  driftline->EnablePlotting(vd);
  //  driftline->EnableSignalCalculation();
  unsigned int ne=0, ni=0;

  for (int iplane=0;iplane<1;iplane++){
    for (int iw =0;iw<7;iw++){
      int iadd = iplane*7 + iw;
      SensewireSig[iadd]->Reset();
    }
  }
  // really the x_i here is 7.62 but also subtract off the diameter of the cathode wire (2*rCathode)
  while(z_i<=max_z_i){
    std::vector<double> num_electrons_y(size_of_vec_in_y);
    std::vector<double> average_y(size_of_vec_in_y);
    std::vector<double> stddev_y(size_of_vec_in_y);
    std::vector<double> y_y(size_of_vec_in_y);
    std::vector<double> z_y(size_of_vec_in_y);
    std::vector<double> x_y(size_of_vec_in_y);
    std::vector<double> avg_x(size_of_vec_in_y);
    std::vector<double> avg_y(size_of_vec_in_y);
    std::vector<double> avg_z(size_of_vec_in_y);
    std::vector<double> std_x(size_of_vec_in_y);
    std::vector<double> std_y(size_of_vec_in_y);
    std::vector<double> std_z(size_of_vec_in_y);

    int iter_y=0;
    while(y_i<=max_y_i){
      int num_electrons=0;
      double min_variation=0.1;
      double running_sum=0;
      double running_average=0;
      double new_average=0;
      double curr_sig=0;
      double curr_x=0;
      double curr_y=0;
      double curr_z=0;
      bool keep_running=true;
      std::vector<double> electrons_drifted_t;
      std::vector<double> electrons_drifted_x;
      std::vector<double> electrons_drifted_y;
      std::vector<double> electrons_drifted_z;
      while(keep_running){
        double xendpoint = 0, yendpoint = 0, zendpoint=0;
        double xendpoint2 = 0, yendpoint2 = 0, zendpoint2=0;
        double tendpoint = 0, tendpoint2 = 0;
	      //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
        int stat=0;
	      //outputfilegetelectron << "startpoint" << i  << " of " << ncl << "  electrons " << x << " " << y << " " << z << " " << t << " " << e << " " << dx << " " << dy << " " << dz << std::endl;
        //std::cout << __LINE__ << std::endl;
	      driftline->DriftElectron(x_i,y_i,z_i,0);
        //std::cout << __LINE__ << std::endl;
	      int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
          numfile << num_electrons << "," << nelectronpoints << "," << y_i << "," << z_i << "," << curr_sig << "," << curr_x << "," << curr_y << "," << curr_z << std::endl;
	      driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint, xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
        curr_sig=tendpoint2-tendpoint;
        curr_x=xendpoint2-xendpoint;
        curr_y=yendpoint2-yendpoint;
        curr_z=zendpoint2-zendpoint;
        electrons_drifted_t.push_back(curr_sig);
        electrons_drifted_x.push_back(curr_x);
        electrons_drifted_y.push_back(curr_y);
        electrons_drifted_z.push_back(curr_z);
        if(stat<0){
          statfile << num_electrons << "," << stat << "," << y_i << "," << z_i << "," << curr_sig << "," << curr_x << "," << curr_y << "," << curr_z << std::endl;
        }
        num_electrons++;
        running_sum+=curr_sig;
        new_average=running_sum/num_electrons;
        if(num_electrons>1000 && TMath::Abs((new_average-running_average)/running_average)<=min_variation) keep_running=false;
        else{
          // compute the average some more
          running_average=new_average;
        }
      }
      // write these out to vectors to spit out to file at the end
      //outputfile << num_electrons << "," << new_average << "," << getSigma(electrons_drifted) << "," << y_i << "," << z_i << "," << x_i << std::endl;
      num_electrons_y[iter_y]=num_electrons;
      average_y[iter_y]=new_average;
      stddev_y[iter_y]=getSigma(electrons_drifted_t);
      y_y[iter_y]=y_i;
      z_y[iter_y]=z_i;
      x_y[iter_y]=x_i;
      avg_x[iter_y]=getMean(electrons_drifted_x);
      avg_y[iter_y]=getMean(electrons_drifted_y);
      avg_z[iter_y]=getMean(electrons_drifted_z);
      std_x[iter_y]=getSigma(electrons_drifted_x);
      std_y[iter_y]=getSigma(electrons_drifted_y);
      std_z[iter_y]=getSigma(electrons_drifted_z);
      y_i+=stepy;
      iter_y++;
    }
    // dump to file
    auto itA = num_electrons_y.begin();
    auto itB = average_y.begin();
    auto itC = stddev_y.begin();
    auto itD = y_y.begin();
    auto itE = z_y.begin();
    auto itF = x_y.begin();
    auto itG = avg_x.begin();
    auto itH = avg_y.begin();
    auto itI = avg_z.begin();
    auto itJ = std_x.begin();
    auto itK = std_y.begin();
    auto itL = std_z.begin();

    int myval=0;
    while(myval<iter_y){
      outputfile << num_electrons_y[myval] << "," << average_y[myval] << "," << stddev_y[myval] << "," 
      << y_y[myval] << "," << z_y[myval] << "," << x_y[myval] << "," << avg_x[myval] << "," << avg_y[myval] << "," << avg_z[myval]
      << "," << std_x[myval] << "," << std_y[myval] << "," << std_z[myval]
      << std::endl;
      myval++;
    }
    /*for (auto& [a,b,c,d,e,f] : zip(num_electrons_y, average_y,stddev_y,y_y,z_y,x_y)) {
      outputfile << a << "," << b << "," << c << "," << d << "," << e << "," << f << std::endl;

    }*/
    z_i+=stepz;
    y_i=min_y_i;
  }
  //std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " RMS is : " << getSigma(electrons_drifted) << std::endl;
  
  outputfile.close();
  statfile.close();
  numfile.close();
  app->Run(kTRUE);
  return 0;
}
