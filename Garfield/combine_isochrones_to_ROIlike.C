/// \file
/// \ingroup tutorial_io
/// \notebook -nodraw
/// Macro to add histogram files
/// This macro is kept for didactical purposes only: use instead the executable $ROOTSYS/bin/hadd !
/// 
/// This macro will add histograms from a list of root files and write them
/// to a target root file. The target file is newly created and must not be
/// identical to one of the source files.
/// This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
/// which had a problem with directories more than one level deep.
/// The macro from Sven has been enhanced by Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
/// to automatically add Trees (via a chain of trees).
///
/// \macro_code
///
/// \author Sven A. Schmidt, sven.schmidt@cern.ch, 13.2.2001


#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <cstring>

namespace fs = std::filesystem;

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );

std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}


void combine_isochrones_to_ROIlike() {

  TString dir = "iso_foranalysis";  
  //TFile* fme=new TFile("flight_isochrones_output.root","recreate");
  //--- filenames are unique so we can use a set
  set<fs::path> sorted_by_name;

  for (auto &entry : fs::directory_iterator(dir.Data()))
    sorted_by_name.insert(entry.path());

  //--- print the files sorted by filename
  FileList = new TList();
  for (auto &filename : sorted_by_name) {
    TString fname = filename.c_str();
    cout << fname << endl; 
    if( fname.Index("isochron_distr_") == -1 ) continue; //isochron_distr_2.5_2.9_40_.root

    cout << fname << endl; 
    FileList->Add( TFile::Open( fname ) );
    TFile *first_source = (TFile*)FileList->First();
  }
  Target = TFile::Open( "combined_iso_ROI.root", "RECREATE" );
  MergeRootfile( Target, FileList );
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  cout << "path is " << path << endl;

  // new TGraphto store things
  auto glow = new TGraph2D();
  glow->SetTitle("myisochrone299");
  glow->SetName("myisochrone299");

  //TH2D * isolow = new TH2D("iso293","iso293",31,-77.5,77.5,289,-289,289);
  TH1D * roidistr[72];
  TH1D * Normalroidistr[72];
  TH1D * NormalLroidistr[72];
  TH1D * NormalRroidistr[72];

  TH1D * Lroidistr[72];
  TH1D * Rroidistr[72];
  ostringstream hID;
  ostringstream hTitle;
  for(int i=0;i<72;i++){  
    hID.str( std::string());
    hID.clear();
    hID << "Wire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << " Wire " << i << " ROIlike; time[ns]; counts";
    roidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);
    hID.str( std::string());
    hID.clear();
    hID << "LeftWire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << "LeftWire " << i << " ROIlike; time[ns]; counts";
    Lroidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);
    hID.str( std::string());
    hID.clear();
    hID << "RightWire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << "RightWire " << i << " ROIlike; time[ns]; counts";
    Rroidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);
    hID.str( std::string());
    hID.clear();
    hID << "NormalWire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << "NormalWire " << i << " ROIlike; time[ns]; counts";
    Normalroidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);
    hID.str( std::string());
    hID.clear();
    hID << "NormalLeftWire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << "NormalLeftWire " << i << " ROIlike; time[ns]; counts";
    NormalLroidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);
    hID.str( std::string());
    hID.clear();
    hID << "NormalRightWire_" << i << "_timedistr";
    hTitle.str( std::string());
    hTitle.clear(); 
    hTitle << "NormalRightWire " << i << " ROIlike; time[ns]; counts";
    NormalRroidistr[i] = new TH1D(hID.str().c_str(),hTitle.str().c_str(),1500,-0.5,14999.5);

  }
  TH2D * isolow = new TH2D("iso299_inverted","iso299_inverted",6000,-14997.5,14997.5,578,-289,289);

  for(int bint=1;bint<6001;bint++){
    for(int binz=1;binz<578;binz++){
      isolow->SetBinContent(bint,binz,-1000);
    }
  }
  int iter_i=0;
  // set up the TList ->TFile iterator?
  TIter next(sourcelist);
  while (TObject *first_source_obj = next()){
  TFile *first_source = (TFile*) first_source_obj;
  TDirectory *current_sourcedir = gDirectory; // this is needed once we open stuff?
  // loop over all keys in this directory
  std::string delimiter="_";

  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  first_source->cd( path ); // change to that file?
  const char * filename_to_split=first_source->GetName();
  std::string str_to_split(filename_to_split);
  std::vector<std::string> splitted_filename= split(str_to_split, delimiter); // ["scott", "tiger", "mushroom"]
  double temperature = std::stod(splitted_filename[6]);
  std::cout << "temperature is : " << temperature << std::endl;


  char dritftime_find='D';
  // goes through all the TH1Ds in the file open currently
  for (TObject* keyAsObj : *gDirectory->GetListOfKeys()){
    auto key = dynamic_cast<TKey*>(keyAsObj);
    if (strchr(key->GetName(), dritftime_find) != nullptr){
      std::cout << "Key name: " << key->GetName() << " Type: " << key->GetClassName() << std::endl;
      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TH1D::Class() ) ) {
        TH1D * thish= (TH1D*)obj;
        const char * idk=thish->GetName();
        //cout << "opening " << idk << endl;
        std::string str(idk);
        // split on underscore _
        std::vector<std::string> splitted= split(str, delimiter); // ["scott", "tiger", "mushroom"]
        //cout << "Result z,y: " << splitted[2] << ", " << splitted[4] << endl;
        if ( thish ) {
          glow->SetPoint(iter_i,std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
          if(std::stod(splitted[4]) <0) isolow->SetBinContent(isolow->GetXaxis()->FindBin(-1.0*thish->GetMean()),isolow->GetYaxis()->FindBin(std::stod(splitted[2])),std::stod(splitted[4]));
          else isolow->SetBinContent(isolow->GetXaxis()->FindBin(thish->GetMean()),isolow->GetYaxis()->FindBin(std::stod(splitted[2])),std::stod(splitted[4]));  
          iter_i++;
          // find which wire and add its contents to existing hist if so +-3.9
          double zlocation = std::stod(splitted[2]); // if 287.9 to 280.1, id=0, if 279.9 to 270.1 id=1
          double which_id=(zlocation-284)/-8.0;
          int this_id=(int)which_id;
          if(this_id<72){
            roidistr[this_id]->Add(thish);
            double deltaID=which_id-(double)this_id;
            if(TMath::Abs(deltaID)<0.05) Normalroidistr[this_id]->Add(thish);
            // also check if its left or right
            double ylocation=std::stod(splitted[4]);
            if(ylocation<0) Lroidistr[this_id]->Add(thish);
            else Rroidistr[this_id]->Add(thish);
            if(TMath::Abs(deltaID)<0.05 && ylocation < 0) NormalLroidistr[this_id]->Add(thish);
            else if(TMath::Abs(deltaID)<0.05 && ylocation > 0) NormalRroidistr[this_id]->Add(thish);

          }
        }
      }
    }
  }

  }
  //TFile *first_source = (TFile*)sourcelist->First(); 
   /*
   while ( (key = (TKey*)nextkey())) {
      // read object from first source file
      first_source->cd( path ); // change to that file?
      TObject *obj = key->ReadObj(); // nextkey is a TIter for this file name to readobject from the file?
      if ( obj->IsA()->InheritsFrom( TH1D::Class() ) ) {
        TH1D * thish= (TH1D*)obj;
        const char * idk=first_source->GetName();
        cout << "opening " << idk << endl;
        std::string str(idk);
        // split on underscore _
        std::vector<std::string> splitted= split(str, delimiter); // ["scott", "tiger", "mushroom"]
        cout << "Result: " << splitted[4] << endl;
        if ( thish ) {
          cout << thish->GetName() << endl;
          //target->cd();
          //myf->Write(str.substr(lastSlash + 2, dot - lastSlash - 2).c_str());
          //obj->Write(key->GetName());
          //try to write the object to our new
          //g->SetPoint(iter_i,thisx,thisy,thish->GetMean());
          cout << "written" << endl;
        }
         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
        TFile *nextsource = (TFile*)sourcelist->After( first_source );
        while ( nextsource ) {
            //cout << "opening " << nextsource->GetName() << endl;
            const char * idk=nextsource->GetName();
            cout << "opening " << idk << endl;
            std::string str(idk);
            // split on underscore _
            std::vector<std::string> splitted= split(str, delimiter); // ["scott", "tiger", "mushroom"]
            cout << "Result: " << splitted[4] << endl;
            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(thish->GetName());
            if (key2) {
                TObject *obj2 = key2->ReadObj();
                if ( obj->IsA()->InheritsFrom( TH1D::Class() ) ) {
                  TH1D *thish2 = (TH1D*)obj2;
                  cout << thish2->GetName() << endl;
                  // clone it>
                  //target->cd();
                  //myf2->Write(str.substr(lastSlash + 2, dot - lastSlash - 2).c_str());
                //obj->Write(key->GetName());
                 //try to write the object to our new
                 cout << "written" << endl;
               }
            }
            nextsource = (TFile*)sourcelist->After( nextsource );
        }
      }
      else {
         // object is of no type that we know or can handle
         cout << "Unknown object type, name: "
         << obj->GetName() << " title: " << obj->GetTitle() << endl;
      }

   }

   */
   // save modifications to target file
  target->cd();
  glow->Write();

  isolow->Write();
  for(int i=0;i<72;i++){
    roidistr[i]->Write();
    Lroidistr[i]->Write();
    Rroidistr[i]->Write();
    Normalroidistr[i]->Write();
    NormalLroidistr[i]->Write();
    NormalRroidistr[i]->Write();

  }

  target->SaveSelf(kTRUE);
}