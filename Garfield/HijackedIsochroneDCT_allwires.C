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
#include <TParameter.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph2D.h>
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
#include <vector>
#include <numeric>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include "ViewIsochrons.hh"
using namespace Garfield;




int whichWireID_int(double ypos){
  return std::round((double)(((71)/2.0)-ypos/8.0));
}

double whichYPos(int thiswireID){
  return (double)((8.0 *((71)/2.0 - (double)(thiswireID))));
}

double Interpolate(const std::vector<double>& y,
                   const std::vector<double>& x, const double xx) {

  const double tol = 1.e-6 * fabs(x.back() - x.front());
  if (xx < x.front()) return y.front();
  const auto it1 = std::upper_bound(x.cbegin(), x.cend(), xx);
  if (it1 == x.cend()) return y.back();
  const auto it0 = std::prev(it1);
  const double dx = (*it1 - *it0);
  if (dx < tol) return y[it0 - x.cbegin()];
  const double f = (xx - *it0) / dx;
  return y[it0 - x.cbegin()] * (1. - f) + f * y[it1 - x.cbegin()];
}

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
  double invals[10]={0};
  for(int i = 1; i < argc; i++){
    invals[i-1] = atof(argv[i]);
    std::cout << invals[i-1] << std::endl;
  }
  const double zpos = invals[0]; // HELIX zposition slice to plot Efield in millimeters
  std::cout << "zpos slicing is [mm] : " << zpos << std::endl;
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
  gas->LoadGasFile("Flight2024_Boff_P_755.865_T_293.15_multiE_90CO2_10Ar.gas");
  char * IonData = getenv("GARFIELD_IONDATA") ;
  gas->LoadIonMobility(IonData);
  gas->PrintGas();
  ComponentAnalyticField * cmp = new ComponentAnalyticField();
  cmp->SetMagneticField(0.,0.,0.0);
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


 // We are essentially copying all of the computeDriftLines function inside of the ViewIsochrone Class
  // this function looks like : ComputeDriftLines(tstep, points, driftLines, startPoints, endPoints, statusCodes, rev);
    //-----------------------------------------------------------------------
  //   DRFEQT - The main routine (DRFEQT) accumulates equal drift time data
  //   DRFEQP   which is plotted as a set of contours in the entry DRFEQP.
  //-----------------------------------------------------------------------
  // ComputeDriftLines requires a few input parameters we expose here:
  const double tstep = 10.0; // 10 ns time steps when interpolating the resulting drift times
  const bool rev = true; // this appears to be needed to get the drift times to be constructed properly (output becomes how long it takes for electron to be at endpoint)
  // the input points which were constructed from a for loop
  std::vector<std::array<double, 3> > points;
  unsigned int nPoints = 600;
  const double ymin=-30.0;
  const double ymax=30.0;
  //unsigned int nPoints = 20;
  //const double ymin=-1.0;
  //const double ymax=1.0;
  double start_drift =-7.62;
  for (unsigned int i = 0; i < nPoints; ++i) {
    const double x0 = start_drift;
    const double y0 = ymin + i * (ymax - ymin) / nPoints;
    std::array<double, 3> p0 = {x0, y0, 0.};
    points.push_back(std::move(p0));
  }
  start_drift =7.62;  // other side of the wires
  for (unsigned int i = 0; i < nPoints; ++i) {
    const double x0 = start_drift;
    const double y0 = ymin + i * (ymax - ymin) / nPoints;
    std::array<double, 3> p0 = {x0, y0, 0.};
    points.push_back(std::move(p0));
  }

  // passed as empty to ComputerDriftLines
  std::vector<std::vector<std::array<double, 4> > > driftLines;
  std::vector<std::array<double, 3> > startPoints;  
  std::vector<std::array<double, 3> > endPoints;  
  std::vector<int> statusCodes;


  // steps when running through ComputeDriftLines
  DriftLineRKF drift; 
  drift.SetSensor(sensor);
  drift.SetMaximumStepSize(); // need to think about this potentially
  drift.EnableSignalCalculation(false);
// everything else is one big for loop around the starting points
  for (const auto& point : points) {
    // always use DriftElectron 
    drift.DriftElectron(point[0], point[1], point[2], 0.);
    const unsigned int nu = drift.GetNumberOfDriftLinePoints();
    // Check that the drift line has enough points.
    if (nu < 3) continue;
    int status = 0;
    double xf = 0., yf = 0., zf = 0., tf = 0.;
    drift.GetEndPoint(xf, yf, zf, tf, status);
    // Find the number of points to be stored.
    const unsigned int nSteps = static_cast<unsigned int>(tf / tstep);
    if (nSteps == 0) continue;
    std::vector<double> xu(nu, 0.);
    std::vector<double> yu(nu, 0.);
    std::vector<double> zu(nu, 0.);
    std::vector<double> tu(nu, 0.);
    for (unsigned int i = 0; i < nu; ++i) {
      drift.GetDriftLinePoint(i, xu[i], yu[i], zu[i], tu[i]);
    }
    if (rev) {
      for (auto& t : tu) t = tf - t;
      std::reverse(std::begin(xu), std::end(xu)); 
      std::reverse(std::begin(yu), std::end(yu)); 
      std::reverse(std::begin(zu), std::end(zu)); 
      std::reverse(std::begin(tu), std::end(tu)); 
    }
    std::vector<std::array<double, 4> > tab;
    // Interpolate at regular time intervals.
    for (unsigned int i = 0; i < nSteps; ++i) {
      const double t = (i + 1) * tstep;
      // tab.push_back(PLACO3(Interpolate(xu, tu, t),
      //                      Interpolate(yu, tu, t),
      //                      Interpolate(zu, tu, t)));
      std::array<double, 4> step = {Interpolate(xu, tu, t),
                                    Interpolate(yu, tu, t),
                                    Interpolate(zu, tu, t),t};
      tab.push_back(step);
    }
    driftLines.push_back(std::move(tab));
    std::array<double, 3> start = {xu[0], yu[0], zu[0]};
    std::array<double, 3> end = {xu[nu - 1], yu[nu - 1], zu[nu - 1]};
    startPoints.push_back(std::move(start));
    endPoints.push_back(std::move(end));
    // Store the drift line return code.
    if (rev) {
      statusCodes.push_back(status);
    } else {
      statusCodes.push_back(0);
    }
  }
  //
  // now what to do with these vectors?
  // eh figure that out afterwards...
  std::cout << "Isochrones completed: \n";
  std::cout << "Some Stats: " << std::endl;
  std::cout << "Size Statuses: " << statusCodes.size() << std::endl;
  std::cout << "startpoints: "<< startPoints.size() << std::endl;
  std::cout << "endpoints: "<< endPoints.size() << std::endl;
  std::cout << "DriftLines: "<< driftLines.size() << std::endl;
  //for (const auto& driftLine : driftLines) {
  //  std::cout << "driftline has size : " << driftLine.size() << std::endl;
  //}
  // need a TGraph2D for each of the sense wires available?
  TGraph2D * senses[72];
  TGraph2D * invsenses[72];
  //TH2F * sensesH[72];
  TH2F * invsensesH[72];
  std::stringstream graphNames;
  std::stringstream TitleNames;

  const double time_step_histo=tstep;
  int numtimebins=std::round(((double)(15000.0)-(double)(-15000.0))/(time_step_histo/1.0))+1;
  double lowTimeEdge=-15000.0-(time_step_histo/2.0);
  double highTimeEdge=15000.0+(time_step_histo/2.0);
  const double y_range=2.0;
  const double stepy=0.1;
  for(int h=0;h<72;h++){
    const double YPosOfInterest = whichYPos(h)/10.0; // now in cm D'oh
    double zone_y_min=YPosOfInterest-y_range; // lowest for top
    double zone_y_max=YPosOfInterest+y_range; // highest for top
    std::cout << "zone_y_min and max are: " << zone_y_min << ", " << zone_y_max << std::endl;
    int nbinsy=std::round((zone_y_max-zone_y_min)/(stepy/1.0))+1;
    double y_low_edge=zone_y_min-(stepy/2.0);
    double y_high_edge=zone_y_max+(stepy/2.0);
    if(zone_y_min < -30.0) zone_y_min=-29.6;
    if(zone_y_max > 29.8) zone_y_max=29.6;
    graphNames.str("");
    graphNames << "iso_inverse_" << h;
    TitleNames.str("");
    TitleNames << "Inverse sense wire " << h << ";y(cm);z(cm);time(ns)";
    invsensesH[h]= new TH2F(graphNames.str().c_str(),TitleNames.str().c_str(),numtimebins,lowTimeEdge,highTimeEdge,nbinsy,y_low_edge,y_high_edge);
  }
  for(int s=0;s<72;s++){
    senses[s] = new TGraph2D();
    graphNames.str("");
    graphNames << "Graphsensewire_" << s;
    senses[s]->SetName(graphNames.str().c_str());
    graphNames.str("");
    graphNames << "Sense Wire " << s << "; y(cm); z(cm)";
    senses[s]->SetTitle(graphNames.str().c_str());
    invsenses[s] = new TGraph2D();
    graphNames.str("");
    graphNames << "GraphInvertedsensewire_" << s;
    invsenses[s]->SetName(graphNames.str().c_str());
    graphNames.str("");
    graphNames << "Inverse Isochrone Sense Wire " << s << "; y(cm); z(cm)";
    invsenses[s]->SetTitle(graphNames.str().c_str());
  }

  TGraph * testISO = new TGraph();
  testISO->SetName("Iso5000");
  testISO->SetTitle("Isochrone example; y(cm);z(cm)");

  for (int iter=0; iter<statusCodes.size();iter++) {
    if(statusCodes[iter]>0 && statusCodes[iter]<=72){
      int thisWireIndex=statusCodes[iter]-1;
      for(const auto& drift: driftLines[iter]){
        senses[thisWireIndex]->AddPoint(drift[0],drift[1],drift[3]);
        invsenses[thisWireIndex]->AddPoint(drift[3],drift[1],drift[0]);
        if(drift[0]>0) invsensesH[thisWireIndex]->SetBinContent(invsensesH[thisWireIndex]->GetXaxis()->FindBin(drift[3]),invsensesH[thisWireIndex]->GetYaxis()->FindBin(drift[1]),drift[0]);
        else if(drift[0]<0) invsensesH[thisWireIndex]->SetBinContent(invsensesH[thisWireIndex]->GetXaxis()->FindBin(-1.0*drift[3]),invsensesH[thisWireIndex]->GetYaxis()->FindBin(drift[1]),drift[0]);
        if(drift[3]<5009 && drift[3] >4991) testISO->AddPoint(drift[0],drift[1]);
      }
    }
  }
  // these comments just show how to access the times and shit
  //const auto drift_ex=driftLines[0];
  //for(const auto& drift_arr: drift_ex){
  //  std::cout << drift_arr[0] << "," << drift_arr[1] << "," << drift_arr[2] << "," << drift_arr[3] << std::endl;
  //}

  std::cout << "And now min and max status codes are: ";
    // Use std::minmax_element to get iterators to min and max elements
  auto minmax_it = std::minmax_element(statusCodes.begin(), statusCodes.end());
  // Dereference the iterators to get the actual values
  int min_value = *minmax_it.first;
  int max_value = *minmax_it.second;
  double step = 0.5; // for ints..
  int nbins=std::round(((double)(max_value)-(double)(min_value))/(step/1.0))+1;
  double low_edge=min_value-(step/2.0);
  double high_edge=max_value+(step/2.0);
  TH1F * status_H = new TH1F("statusCodes","statusCodes", nbins,low_edge,high_edge);
  // do this per histogram to investigate profiles  
  double Startstepsizes[3] ={0.1,0.1,0.1}; // for histos..
  int StartnbinsXYZ[3]={0};
  double Startlow_edgesXYZ[3]={0};
  double Starthigh_edgesXYZ[3]={0};
  // Initialize min_val and max_val with the first element of the specified dimension
  for(int dimension_index=0;dimension_index<3;dimension_index++){
    double minXYZ = startPoints[0][dimension_index];
    double maxXYZ = startPoints[0][dimension_index];
    // Iterate through the rest of the vector to find the true min and max
    for (size_t i = 1; i < startPoints.size(); ++i) {
        minXYZ = std::min(minXYZ, startPoints[i][dimension_index]);
        maxXYZ = std::max(maxXYZ, startPoints[i][dimension_index]);
    }
    // now fill my arrays for determining appropriate bins
    StartnbinsXYZ[dimension_index]=std::round(((double)(maxXYZ)-(double)(minXYZ))/(Startstepsizes[dimension_index]/1.0))+1;
    Startlow_edgesXYZ[dimension_index]=minXYZ-(Startstepsizes[dimension_index]/2.0);
    Starthigh_edgesXYZ[dimension_index]=maxXYZ+(Startstepsizes[dimension_index]/2.0);
  }
  TH1F * startx_H = new TH1F("startX","startX", StartnbinsXYZ[0],Startlow_edgesXYZ[0],Starthigh_edgesXYZ[0]);
  TH1F * starty_H = new TH1F("startY","startY", StartnbinsXYZ[1],Startlow_edgesXYZ[1],Starthigh_edgesXYZ[1]);
  TH1F * startz_H = new TH1F("startZ","startZ", StartnbinsXYZ[2],Startlow_edgesXYZ[2],Starthigh_edgesXYZ[2]);
  double Endstepsizes[3] ={0.1,0.1,0.1}; // for histos..
  int EndnbinsXYZ[3]={0};
  double Endlow_edgesXYZ[3]={0};
  double Endhigh_edgesXYZ[3]={0};
  // Initialize min_val and max_val with the first element of the specified dimension
  for(int dimension_index=0;dimension_index<3;dimension_index++){
    double minXYZ = endPoints[0][dimension_index];
    double maxXYZ = endPoints[0][dimension_index];
    // Iterate through the rest of the vector to find the true min and max
    for (size_t i = 1; i < endPoints.size(); ++i) {
        minXYZ = std::min(minXYZ, endPoints[i][dimension_index]);
        maxXYZ = std::max(maxXYZ, endPoints[i][dimension_index]);
    }
    // now fill my arrays for determining appropriate bins
    EndnbinsXYZ[dimension_index]=std::round(((double)(maxXYZ)-(double)(minXYZ))/(Endstepsizes[dimension_index]/1.0))+1;
    Endlow_edgesXYZ[dimension_index]=minXYZ-(Endstepsizes[dimension_index]/2.0);
    Endhigh_edgesXYZ[dimension_index]=maxXYZ+(Endstepsizes[dimension_index]/2.0);
  }
  TH1F * endx_H = new TH1F("endX","endX", EndnbinsXYZ[0],Endlow_edgesXYZ[0],Endhigh_edgesXYZ[0]);
  TH1F * endy_H = new TH1F("endY","endY", EndnbinsXYZ[1],Endlow_edgesXYZ[1],Endhigh_edgesXYZ[1]);
  TH1F * endz_H = new TH1F("endZ","endZ", EndnbinsXYZ[2],Endlow_edgesXYZ[2],Endhigh_edgesXYZ[2]);


  
  for (int thisstatus : statusCodes) {
      status_H->Fill(thisstatus);
  }
      // Loop through the vector and fill my histos (index is x y or z)
  for (const auto& arr : startPoints) {
      startx_H->Fill(arr[0]); 
      starty_H->Fill(arr[1]); 
      startz_H->Fill(arr[2]); 
  }
  for (const auto& arr : endPoints) {
      endx_H->Fill(arr[0]); 
      endy_H->Fill(arr[1]); 
      endz_H->Fill(arr[2]); 
  }

  // trying to understand driftline stuff


  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds\n";
  std::stringstream outrootfilename;
  outrootfilename << "DriftLineRKFTest_allwires.root";
  TFile * Outfile = new TFile(outrootfilename.str().c_str(),"recreate");
  Outfile->cd();
  status_H->Write();
  startx_H->Write();
  starty_H->Write();
  startz_H->Write();
  endx_H->Write();
  endy_H->Write();
  endz_H->Write();
  for(int s =0;s<72;s++){
    senses[s]->Write();
    //invsenses[s]->Write();
    invsensesH[s]->Write();
  }
  //testISO->Write();
  TParameter bfield("Bfield", 0.0);
  TParameter temperature("Temperature", 293.0);
  TParameter versionflag("Isochrone_version", 5);
  bfield.Write();
  temperature.Write();
  versionflag.Write();
  Outfile->Close();
  app->Run(kTRUE);
  return 0;
}
