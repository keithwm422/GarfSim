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

  // TApplication app("app", &argc, argv);
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
//  const double pressure = 1 * AtmosphericPressure; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr)
  const double pressure = invals[0]*51.7149; // in torr- we were at 14.6 psi (1 psi = 51.7149 torr)
  std::cout << "pressure in torr: " << pressure << std::endl;
  const double temperature = invals[1];
  std::cout << "temperature in kelvin: " << temperature << std::endl;
  std::stringstream ingasfilename;
  //outfilename << "Flight2024_Boff_P_" << pressure <<"_T_" << temperature << "_.gas";
  ingasfilename << "FlightGasFiles/BOFF/Flight2024_Boff_P_755.865_T_" << temperature << ".15_multiE_90CO2_10Ar_01122024.gas"; // path FlightGasFiles/BOFF/Flight2024_Boff_P_755.865_T_279.15_multiE_90CO2_10Ar_01122024.gas
  Garfield::plottingEngine.SetDefaultStyle();
  MediumMagboltz * gas = new MediumMagboltz();
//  gas->LoadGasFile("co2_90_AR_10_T273.gas");
//  gas->LoadGasFile("keith_co2_85_AR_15_T273.gas");
//  gas->LoadGasFile("Flight2024_Bon_P_755.865_T_298.15_.gas");
  gas->LoadGasFile(ingasfilename.str());

  // lets just print out the drift velocity to a file?

  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();
  //return 0;


  auto start = std::chrono::high_resolution_clock::now();
  bool realtimeplots = true;
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);

  std::stringstream wid;

  // for the drift time vs temp calculation of the DCT
  const double max_y_i=30.0, max_x_i=7.619, stepy=0.1,stepx=0.25, min_y_i=-30.0;
  const double x_i=-7.619, y_i=min_y_i, z_i=0.0, t_i=0, e_i=0, dx_i=0, dy_i=0, dz_i=0; // -300e-4+rAnode+100e-4 , 2.4 ,0 ,0 (don't forget the stagger in x! add extra 0.03)

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
  std::stringstream outputfilename;
  //outfilename << "Flight2024_Boff_P_" << pressure <<"_T_" << temperature << "_.gas";
  outputfilename << "isochronGenerator_p_fullEfield_T" << temperature << ".txt";
  std::ofstream outputfile(outputfilename.str().c_str(),std::ios_base::app);

  outputfile << "num_electrons,average,stddev,y,z,x,dx,dy,dz,sigx,sigy,sigz" << std::endl;
  outputfile << std::fixed << std::showpoint;
  outputfile << std::setprecision(10);
  outputfile.close();
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
  

  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  cmp->SetMagneticField(0.,0.,0.0);
//  cmp->SetMagneticField(0.,0.,1.0);

  GeometrySimple * geo = new GeometrySimple();
 
  SolidBox * enclosure = new SolidBox(0,0,0,10,31,11);
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

  Sensor * sensor = new Sensor;
  sensor->AddComponent(cmp);
  sensor->SetTimeWindow(0,2,20000); // might need to change this, its start, step size, number of steps
  cmp->AddReadout("a");
  sensor->AddElectrode(cmp,"a");
  unsigned int ne=0, ni=0;
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
    std::vector<double> e_a(size_of_vec_in_y);
    std::vector<double> i_a(size_of_vec_in_y);
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

    for(int iy=0;iy<size_of_vec_in_y;iy++){
      double y_delt=y_i+((double)(iy)*stepy);
      int num_electrons=0;
      int ne_av=0, ni_av=0;
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
        double energy0=0,energy1=0;
	      //track->GetElectron(i,x,y,z,t,e,dx,dy,dz);
        int stat=0;
	      //outputfilegetelectron << "startpoint" << i  << " of " << ncl << "  electrons " << x << " " << y << " " << z << " " << t << " " << e << " " << dx << " " << dy << " " << dz << std::endl;
        //std::cout << __LINE__ << std::endl;
	bool did_it_drift=driftline->DriftElectron(x_delt,y_delt,z_i,0);
        if(!did_it_drift) is_pos_borked=true;
        //std::cout << __LINE__ << std::endl;
	      //int nelectronpoints = driftline->GetNumberOfElectronEndpoints();
	      driftline->GetElectronEndpoint(0, xendpoint, yendpoint, zendpoint, tendpoint,xendpoint2, yendpoint2, zendpoint2, tendpoint2, stat);
        //if(ne_av>1) std::cout << "ne_av,ni_av is: " << ne_av << "," << ni_av << std::endl;
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
    outputfile.open(outputfilename.str().c_str(),std::ios_base::app);
    while(myval<iter_y){
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
    //delete aval;
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
