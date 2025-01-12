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
  //keeping x constant (y in HELIX)
  //double max_y_i=0.25, max_z_i=0.25, stepy=0.01,stepz=0.01, min_y_i=-0.25;
  //double x_i=7.60249, y_i=min_y_i, z_i=-0.25, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0
  // keeping z constant (x in HELIX, along wire)

  //double max_y_i=1.0, max_x_i=0.03+0.5, stepy=0.1,stepx=0.005, min_y_i=-1.0;

  //this was for isochron_generator_v2_BON
  //double max_y_i=0.3, max_x_i=0.03+0.5, stepy=0.1,stepx=0.005, min_y_i=-0.3;
  //double x_i=0.03+0.05, y_i=min_y_i, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)

  //double x_i=0.03, y_i=min_y_i, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)

  //double max_y_i=0.5, max_x_i=7.5, stepy=0.1,stepx=0.25, min_y_i=-0.5;
  //double x_i=0.5, y_i=min_y_i, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)

  // for the full view of the isochrons
  const double max_y_i=2.45, max_x_i=7.5, stepy=0.05,stepx=0.1, min_y_i=-2.4;
  const double x_i=-7.5, y_i=min_y_i, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)

  // find the size of the vector we will be storing
  int size_of_vec_in_y=(int)((max_y_i-y_i)/stepy)+1;
  // for parallel, find the number of iterations we will be using to meet the x grid
  int size_of_x_grid=(int)((max_x_i-x_i)/stepx)+1;
  std::cout << "size of vectors will be: " << size_of_vec_in_y << std::endl;
  std::cout << "max_y_i minus y_i will be: " << max_y_i-y_i << std::endl;
  std::cout << "max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
  std::cout << "casted max_y_i minus y_i/stepy will be: " << (max_y_i-y_i)/stepy << std::endl;
   //<< "," << new_average << "," << getSigma(electrons_drifted) << "," << y_i << "," << z_i << "," << x_i << std::endl;

  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);
  std::ofstream outputfile("isochron_BOFF_T293_p_full_compare_nopenning.txt",std::ios_base::app);
  outputfile << "num_electrons,average,stddev,y,z,x,dx,dy,dz,sigx,sigy,sigz" << std::endl;
  outputfile << std::fixed << std::showpoint;
  outputfile << std::setprecision(10);
  outputfile.close();
  //name for drifttimes vs drift distance
  //std::ofstream statfile("iso_status_checke_p.txt");


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

  // Gas not used for microscopic tracking
//  const double pressure = 755.865; //Torr
//  const double temperature = 303.15; //K
 
  // Set the temperature [K] and pressure [Torr]
//  gas->SetTemperature(temperature);
//  gas->SetPressure(pressure);
//  gas->SetComposition("co2", 90, "ar", 10);

//  gas->LoadGasFile("co2_90_AR_10_T273.gas");
//  gas->LoadGasFile("keith_co2_85_AR_15_T273.gas");
//  gas->LoadGasFile("Flight2024_Bon_P_755.865_T_303.15_.gas");
//  gas->LoadGasFile("Flight2024_Bon_P_755.865_T_303.15_multiE.gas");
//  gas->LoadGasFile("Flight2024_Boff_P_755.865_T_303.15_multiE_maxE2000.gas"); // works well! provides output max drift time same as data!
  gas->LoadGasFile("Flight2024_Boff_P_755.865_T_293.15_multiE_90CO2_10Ar_nopenning.gas");
  // lets just print out the drift velocity to a file?

  std::cout << "LOADED GAS FILE" << std::endl;

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  cmp->SetMagneticField(0.,0.,0.0);
  //cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,4.7,11);
  geo->AddSolid(enclosure, gas);
  cmp->SetGeometry(geo);

  const double vCathode= -7500;
  const double rCathode= 175e-4; // is this in centimeters? seems so
  const double vAnode= 0;
  const double rAnode= 20e-4;
  const double vPotential= -2700;
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
  //ViewSignal * vs1 = new ViewSignal;

  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,20000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");
  //vs1->SetSensor(sensor);
  
 
  //AvalancheMC * driftline = new AvalancheMC();
  //  driftline->EnableDebugging();  
  //driftline->SetDistanceSteps(0.001);
  //driftline->EnableMagneticField();
  //driftline->EnableDiffusion();
  //driftline->SetSensor(sensor);
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
  //while(x_i<=max_x_i){
  omp_lock_t writelock;
  omp_init_lock(&writelock);
  #pragma omp parallel for
  //for(double x_delt=x_i; x_delt<=max_x_i;x_delt+=stepx){
  for(int ix=0;ix<size_of_x_grid;ix++){
    double x_delt=x_i+((double)(ix)*stepx);
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
    AvalancheMC * driftline = new AvalancheMC();
    //  driftline->EnableDebugging();  
    driftline->SetDistanceSteps(0.001);
    //driftline->EnableMagneticField();
    driftline->EnableDiffusion();
    driftline->SetSensor(sensor);
    //  driftline->EnablePlotting(vd);
    //  driftline->EnableSignalCalculation();
    unsigned int ne=0, ni=0;

    //while(y_i<=max_y_i){
    for(int iy=0;iy<size_of_vec_in_y;iy++){
      double y_delt=y_i+((double)(iy)*stepy);
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
      bool is_pos_borked=false;
      std::vector<double> electrons_drifted_t;
      std::vector<double> electrons_drifted_x;
      std::vector<double> electrons_drifted_y;
      std::vector<double> electrons_drifted_z;
      while(keep_running && !is_pos_borked){
        double xendpoint = 0, yendpoint = 0, zendpoint=0;
        double xendpoint2 = 0, yendpoint2 = 0, zendpoint2=0;
        double tendpoint = 0, tendpoint2 = 0;
	      //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
        int stat=0;
	      //outputfilegetelectron << "startpoint" << i  << " of " << ncl << "  electrons " << x << " " << y << " " << z << " " << t << " " << e << " " << dx << " " << dy << " " << dz << std::endl;
        //std::cout << __LINE__ << std::endl;
	bool did_it_drift=driftline->DriftElectron(x_delt,y_delt,z_i,0);
        if(!did_it_drift) is_pos_borked=true;
        //std::cout << __LINE__ << std::endl;
	      //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
	      driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint, xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
        //std::cout << "stat is: " << stat << ", x_delt is " << x_delt << ", y_delt is " << y_delt << std::endl;
        //std::cout << "num_electrons is: " << num_electrons << std::endl;
        //if(stat!=0) keep_running=false;
        curr_sig=tendpoint2-tendpoint;
        curr_x=xendpoint2-xendpoint;
        curr_y=yendpoint2-yendpoint;
        curr_z=zendpoint2-zendpoint;
        electrons_drifted_t.push_back(curr_sig);
        electrons_drifted_x.push_back(curr_x);
        electrons_drifted_y.push_back(curr_y);
        electrons_drifted_z.push_back(curr_z);
        //if(stat<0){
        //  statfile << num_electrons << "," << stat << "," << y_delt << "," << z_i << "," << curr_sig << "," << curr_x << "," << curr_y << "," << curr_z << std::endl;
        //  is_pos_borked=true;
        //}
        num_electrons++;
        running_sum+=curr_sig;
        new_average=running_sum/num_electrons;
        if(num_electrons>1000){
          if(TMath::Abs((new_average-running_average)/running_average)<=min_variation) keep_running=false;
          else is_pos_borked=true;
        }
        //if(num_electrons>1000 && TMath::Abs((new_average-running_average)/running_average)<=min_variation) keep_running=false;
        //if(num_electrons>2) keep_running=false;
        else{
          // compute the average some more
          running_average=new_average;
        }
      }
      // write these out to vectors to spit out to file at the end
      //outputfile << num_electrons << "," << new_average << "," << getSigma(electrons_drifted) << "," << y_delt << "," << z_i << "," << x_delt << std::endl;
      if(is_pos_borked){
        std::cout << "x,y borked is " << x_delt << "," << y_delt << std::endl;
        //std::cout << "ypos borked is " << iter_y << std::endl;
        num_electrons_y[iter_y]=-666.0;
        average_y[iter_y]=-666.0;
        stddev_y[iter_y]=-666.0;
        y_y[iter_y]=y_delt;
        z_y[iter_y]=z_i;
        x_y[iter_y]=x_delt;
        avg_x[iter_y]=-666.0;
        avg_y[iter_y]=-666.0;
        avg_z[iter_y]=-666.0;
        std_x[iter_y]=-666.0;
        std_y[iter_y]=-666.0;
        std_z[iter_y]=-666.0;
        //y_delt+=stepy;
        iter_y++;
      }
      else{
        num_electrons_y[iter_y]=num_electrons;
        average_y[iter_y]=new_average;
        stddev_y[iter_y]=getSigma(electrons_drifted_t);
        y_y[iter_y]=y_delt;
        z_y[iter_y]=z_i;
        x_y[iter_y]=x_delt;
        avg_x[iter_y]=getMean(electrons_drifted_x);
        avg_y[iter_y]=getMean(electrons_drifted_y);
        avg_z[iter_y]=getMean(electrons_drifted_z);
        std_x[iter_y]=getSigma(electrons_drifted_x);
        std_y[iter_y]=getSigma(electrons_drifted_y);
        std_z[iter_y]=getSigma(electrons_drifted_z);
        //std::cout << "iter_y" << iter_y << std::endl;
        //y_i+=stepy;
        iter_y++;
      }
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
    //if(is_pos_borked) iter_y=0;
    omp_set_lock(&writelock);
    outputfile.open("isochron_BON_T303_p_full_compare.txt",std::ios_base::app);
    while(myval<iter_y){
      std::cout << "ix written is " << ix << " out of: " << size_of_x_grid << std::endl;
      //std::cout << "iter_y in writing output is " << iter_y << std::endl;
      //std::cout << "myval in writing value is " << myval << std::endl;
      outputfile << num_electrons_y[myval] << "," << average_y[myval] << "," << stddev_y[myval] << "," 
      << y_y[myval] << "," << z_y[myval] << "," << x_y[myval] << "," << avg_x[myval] << "," << avg_y[myval] << "," << avg_z[myval]
      << "," << std_x[myval] << "," << std_y[myval] << "," << std_z[myval]
      << std::endl;
      myval++;
    }
    /*for (auto& [a,b,c,d,e,f] : zip(num_electrons_y, average_y,stddev_y,y_y,z_y,x_y)) {
      outputfile << a << "," << b << "," << c << "," << d << "," << e << "," << f << std::endl;

    }*/
    //x_i+=stepx;
    outputfile.close();
    omp_unset_lock(&writelock);
    //y_i=min_y_i;
    delete driftline;
  }
  //std::cout << "# Avalanched electrons: " << num_electrons << " ave sig is: " << running_average << " RMS is : " << getSigma(electrons_drifted) << std::endl;
 auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  std::cout << "Closing file and finished run" << std::endl;
  omp_destroy_lock(&writelock);
  //outputfile.close();
  //statfile.close();
  app->Run(kTRUE);
  return 0;
}
