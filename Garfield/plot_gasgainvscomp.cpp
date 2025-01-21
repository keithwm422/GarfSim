#include "attrib.C"
Int_t* GetColors( Int_t type, Int_t ndata );
void plot_gasgainvscomp( ){
  TString fname;
  TCanvas *cb = new TCanvas( "cb", "Control", 0, 830, 400, 600 );
  TFile * _file;
  TH1F * h[3];
  fname="DCsimX1p0A0p0.root";
  _file = TFile::Open( fname );
  if(!_file || !_file->IsOpen() || _file->IsZombie()){
    return;
  }
  Double_t maxes[3];
  Double_t gases[3]={15,10,5};
  for(int i=0;i<3;i++){
    if(i==0) h[i]      = (TH1F*)_file->Get("wire_sig0");
    else if(i==1) h[i] = (TH1F*)_file->Get("wire_sig1");
    else if(i==2) h[i] = (TH1F*)_file->Get("wire_sig2");
    maxes[i]=h[i]->GetMinimum();
  }
  maxes[0]=maxes[0]/maxes[1];
  maxes[2]=maxes[2]/maxes[1];
  maxes[1]=maxes[1]/maxes[1];
  TGraph *all_gains = new TGraph(3,gases,maxes);
  Attrib( 20, 24, 24, 1.2, 1.7 );
  //Int_t *col = GetColors( kBird, 3 );
  //h[4]->SetOptStat(0);
  //h[4]->SetLineColor(6);
  // AttribH( h[4], kBlue+3,0.4,0.6);
  TH1F *frame = gPad->DrawFrame(0, 0, 30, 2 );
  frame->GetXaxis()->SetTitle( "Composition (% Argon)" );
  frame->GetYaxis()->SetTitle( "Gas Gain (rel to 10%)" );

  AttribG( all_gains, 1, 5, 1,3,4);
   //TF1 *f3 = gre3->GetFunction("pol3");
  all_gains->Draw("p");
  cb->SaveAs( "gasgain_vs_argon.pdf" );
}

Int_t* GetColors( Int_t type, Int_t ndata )
{
  gStyle->SetPalette( type );
  Int_t ncol = TColor::GetNumberOfColors();
  Int_t FI = TColor::GetColorPalette( 0 );
  Int_t *col = new Int_t[ndata];
  for( Int_t i = 0; i < ndata; i++ )
    col[i] = FI+ncol/ndata*i;

  return col;
}
