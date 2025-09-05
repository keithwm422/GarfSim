void grabIsochrone(){
  std::stringstream TF1name; // for accessing files
  std::stringstream histoname; // for writing histogram and TGraph names
  std::stringstream isonames; // for accessing histo names inside of the input file(s)
  TFile * outputfile; // output file
  TFile * isoInputfile; // input file
  TGraph * isoOfInterest = new TGraph();
  isoOfInterest->SetName("isoRotated");
  isoOfInterest->SetTitle("Isochrone 1532 Rotated; z(cm); y(cm)");
  // save in file.
  TF1name.str("");
  TF1name << "grabbedIsos.root";
  cout << "writing to file: " << TF1name.str().c_str() << endl;
  outputfile = TFile::Open(TF1name.str().c_str(),"RECREATE");
  if(!outputfile || !outputfile->IsOpen() || outputfile->IsZombie()){
    return;
  } 
  TF1name.str("");
  TF1name << "thisisochrone_testDCT.root";
  cout << "opening file: " << TF1name.str().c_str() << endl;
  isoInputfile = TFile::Open(TF1name.str().c_str());
  if(!isoInputfile || !isoInputfile->IsOpen() || isoInputfile->IsZombie()){
    return;
  }
  gROOT->cd();
  TCanvas * thisc = (TCanvas *)isoInputfile->Get("c2");
  TList * listc2 = thisc->GetListOfPrimitives();
  //listc2->Print();
  int i=0;
  histoname.str("");
  for (TObject* obj : *listc2) {
    isoInputfile->cd();
    if ( obj->IsA()->InheritsFrom( TGraph::Class() ) ) {
      histoname.str("");
      histoname << "thisIso" << i;
      TGraph * mygraph = (TGraph*) obj;
      std::cout << mygraph->GetN() << std::endl;
      int numentries=mygraph->GetN();
      if(i==1532){
        for(int n=0;numentries;n++){
          double y;
          double z;
          mygraph->GetPoint(n,y,z);
          outputfile->cd();
          isoOfInterest->AddPoint(z,y);
          //std::cout << n << std::endl;
        }
      }
      mygraph->SetName(histoname.str().c_str());
      histoname.str("");
      histoname << "Isochrone " << i << "; y(cm); z(cm)";
      mygraph->SetTitle(histoname.str().c_str());
      outputfile->cd();
      mygraph->Write();
      std::cout << i << std::endl;
      i++;
    }
  }
  outputfile->cd();
  isoOfInterest->Write();
  delete isoOfInterest;
  outputfile->Close();
  isoInputfile->cd();
  isoInputfile->Close();
}