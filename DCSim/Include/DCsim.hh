#ifndef G_DCsim
#define G_DCsim

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TRandom3.h>
#include <TGeoManager.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

class DCsim  {

 public:
  // Constructor
  DCsim(TRandom3 r);
  // Destructor
  virtual ~DCsim();

  void Copy(DCsim * insim);
  void MakeTree(char * file);
  void SetNwires(int iw);
  void SetWirePos(int iw, float x, float y);
  void Digitize(int iw, TH1D * wavehist);
  void Digitize(TH1D* rawWaveform, TH1D* digitizedWaveform); 
  void InitDigitization();
  void CreateNoise(TH1D* NoiseWave);

  void LoadAmpParameters(double r_w,
			 double vrf,
			 double ven,
			 double vth,
			 double vin,
			 double vdet,
			 double f_z,
			 double f_p,
			 double f0 ) {R_w=r_w; Vrf=vrf; Ven=ven; Vth=vth; Vin=vin; Vdet=vdet; F_z=f_z; F_p=f_p; F_0=f0;};
			 
  double AmpVoltageNoise(double freq);
  void Filter(TH1D* rawWave, TH1D* filteredWave, TH1D* Diff1Sig, TH1D* Int2Sig, double Int1timeconstant, double Dif1timeconstant, double Int2timeconstant);
  void AddNoise(TH1D* InSig, TH1D* NoiseSig);
  void ApplyGainandNoise(TH1D* InSig, TH1D* NoiseSig, TH1D* OutSig);
  double GetWeightingIntegralI1(TH1D* Int2Sig, double Tf);
  double GetWeightingIntegralI2(TH1D* Int2Sig, double Tf);
  double GetBallisticDeficit(TH1D* PreampOut,TH1D* ShaperOut);

  int FirstThreshold(int iw);
  int FirstThreshold(int iw, TH1D* waveform);
  
  int FirstThresholdPreDig(TH1D* Waveform);
  double GetPreSampleRMS(TH1D* Waveform, int threshbin, int nsamples);
  double GetLeadingEdgeSlope(TH1D * Waveform, double rms, int threshbin, int nsamples);

  int GetSample(int iw, int isample);
  double GetSample(int isample, TH1D* digitizedWaveform);
  int GetnSamples(){return nedgesamples;};
  void FitLeadingEdges();
  void FitLeadingEdges(int iw, TH1D* digitizedWaveform);
  void FillTree(double tx, double ty, double tang);
  void WriteTree();
  TFile * DCsimOutFile(){return outfile;};
 private:
  std::string m_className;

  // Options
  bool m_debug;
  std::string m_label;
  TRandom3 rloc;
  float digitPeriod;
  float histbin;
  int waveformLength;
  int nwires;
  double gasGain;
  double transimpedanceGain;
  int fullscaleADC;
  double fullscaleVoltage;
  double ADCthresh;
  int nedgesamples;
  double noiseRMS;
  int firstedgesample;

  double R_w;
  double Vrf;
  double Ven;
  double Vth;
  double Vin;
  double Vdet;
  double F_z;
  double F_p;
  double F_0;


  TFile * outfile;
  TTree * DCtree;
  double int1;
  double int2;
  double dif1;
  float trackstartx;
  float trackstarty;
  float trackangle;
  float tthresh[7];
  float tedgefit[7];
 

  struct wireDataPacket {
    int threshbin;
    std::vector<int> Sample;
  };

  struct wireStruct {
    double x, y;
  };

  struct DCreco{
    double wireT[7];
    

  };

  struct wireSig{
    std::vector<int>Sample;
  };

  std::vector<wireStruct> wireInfo;
  std::vector<wireSig*> DCDigitization;

public: 
  wireDataPacket DCDataPacket[7]; 
  
};

#endif
