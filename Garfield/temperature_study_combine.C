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


void temperature_study_combine() {

  TString dir = "temperature_isochrone_study";  
  //TFile* fme=new TFile("flight_isochrones_output.root","recreate");
  //--- filenames are unique so we can use a set
  set<fs::path> sorted_by_name;

  for (auto &entry : fs::directory_iterator(dir.Data()))
    sorted_by_name.insert(entry.path());

  //--- print the files sorted by filename
  FileList = new TList();
  for (auto &filename : sorted_by_name) {
    TString fname = filename.c_str();
    //cout << fname << endl; 
    if( fname.Index("isochron_distr_") == -1 ) continue; //isochron_distr_2.5_2.9_40_.root

    //cout << fname << endl; 
    FileList->Add( TFile::Open( fname ) );
    TFile *first_source = (TFile*)FileList->First();
  }
  Target = TFile::Open( "temperature_study_isochrones_output.root", "RECREATE" );
  MergeRootfile( Target, FileList );
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  cout << "path is " << path << endl;

  // new TGraphto store things
  auto glow = new TGraph2D();
  TH3D * statush= new TH3D("status","status",301,-301,301,31,-77.5,77.5,20,-20.5,-0.5);
  TH1D * statush2= new TH1D("status2","status2",200,-20.5,-0.5);
  TH2D * isolow = new TH2D("iso283","iso283",31,-77.5,77.5,289,-289,289);
  TH2D * isohigh = new TH2D("iso303","iso303",31,-77.5,77.5,289,-289,289);

  glow->SetTitle("myisochrone283");
  glow->SetName("myisochrone283");
  auto ghigh = new TGraph2D();
  ghigh->SetTitle("myisochrone303");
  ghigh->SetName("myisochrone303");
  std::string delimiter="_";
  int iter_i=0;

  // set up the TList ->TFile iterator?
  TIter next(sourcelist);
  while (TObject *first_source_obj = next()){
    TFile *first_source = (TFile*) first_source_obj;
    TDirectory *current_sourcedir = gDirectory; // this is needed once we open stuff?
    // loop over all keys in this directory
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key, *oldkey=0;
    first_source->cd( path ); // change to that file?
    //std::cout << "first_source name: " << first_source->GetName() << " Title: " << first_source->GetTitle() << std::endl;
    const char * filename_to_split=first_source->GetName();
    //cout << "opening " << idk << endl;
    std::string str_to_split(filename_to_split);
    // split on underscore _
    std::vector<std::string> splitted_filename= split(str_to_split, delimiter); // ["scott", "tiger", "mushroom"]
    double temperature = std::stod(splitted_filename[7]);
    //std::cout << "temperature is : " << temperature << std::endl;

    char dritftime_find='D';
    // goes through all the TH1Ds in the file open currently
    for (TObject* keyAsObj : *gDirectory->GetListOfKeys()){
      auto key = dynamic_cast<TKey*>(keyAsObj);
      if (strchr(key->GetName(), dritftime_find) != nullptr){
        //std::cout << "Key name: " << key->GetName() << " Type: " << key->GetClassName() << std::endl;
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
            //cout << thish->GetName() << endl;
            //target->cd();
            //myf->Write(str.substr(lastSlash + 2, dot - lastSlash - 2).c_str());
            //obj->Write(key->GetName());
            //try to write the object to our new
            if(temperature>290){
              ghigh->SetPoint(iter_i,std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
              isohigh->Fill(std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
            }
            else{
              glow->SetPoint(iter_i,std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
              isolow->Fill(std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
            }
            //glow->SetPoint(iter_i,std::stod(splitted[4]),std::stod(splitted[2]),thish->GetMean());
            //if(std::stod(splitted[2]) <-10 && std::stod(splitted[2]) >-60 && std::stod(splitted[4]) >6 && std::stod(splitted[4])<19){
            //  cout << "Mean is " << thish->GetMean() << " for z,y: " << std::stod(splitted[2]) << ", " << std::stod(splitted[4]) << " for point " << iter_i << endl;
            //  Double_t x,y,z;
            //  glow->GetPoint(iter_i,x,y,z);
            //  cout << "Graph2D now says " << x << "," << y << "," << z << endl;
            //}
            //cout << "written" << endl;
            iter_i++;
          }
        }
      }
      else{
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
            if(thish->GetMean()!=-5.0){
              std::cout << "first_source name: " << first_source->GetName() << " Title: " << first_source->GetTitle() << std::endl;
              std::cout << "maybe we have one to look at! " << thish->GetMean() << " , " << std::stod(splitted[3]) << " , " << std::stod(splitted[5]) << std::endl;
              statush->Fill(std::stod(splitted[3]),std::stod(splitted[5]),thish->GetMean());
              statush2->Fill(thish->GetMean());
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
  ghigh->Write();
  statush->Write();
  statush2->Write();
  isolow->Write();
  isohigh->Write();
  target->SaveSelf(kTRUE);
  std::cout << "number of histograms: " << iter_i << std::endl;
}