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

namespace fs = std::filesystem;

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );

void hadd_Efield_profiles_v3() {

  TString dir = "efield_profs_v3";  
  TFile* fme=new TFile("newfile.root","recreate");
  //--- filenames are unique so we can use a set
  set<fs::path> sorted_by_name;

  for (auto &entry : fs::directory_iterator(dir.Data()))
    sorted_by_name.insert(entry.path());

  //--- print the files sorted by filename
  FileList = new TList();
  for (auto &filename : sorted_by_name) {

    TString fname = filename.c_str();
    cout << fname << endl; 
    if( fname.Index("efieldprofile_zpos_") == -1 ) continue;

    cout << fname << endl; 
    FileList->Add( TFile::Open( fname ) );
    TFile *first_source = (TFile*)FileList->First();
  }

  fme->Close();
  Target = TFile::Open( "profiles_merged_v6.root", "RECREATE" );
  MergeRootfile( Target, FileList );
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {
     cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );
   cout << "path is " << path << endl;
   TFile *first_source = (TFile*)sourcelist->First();
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   // loop over all keys in this directory
   TChain *globChain = 0;
   std::string prename="zis_";
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {
      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TCanvas::Class() ) ) {
        TCanvas *c1 = (TCanvas*)obj;
        cout << c1->GetTitle() << endl;
        TObject * thisf=c1->FindObject("fProfile_0");
        TF1 * myf = (TF1*) thisf;
        const char * idk=first_source->GetName();
        cout << "opening " << idk << endl;
        std::string str(idk);
        // Find the position of the last '/' character
        size_t lastSlash = str.find_last_of('s');
        size_t dot = str.find('.', lastSlash + 1);
        if (lastSlash != string::npos && dot != string::npos) {
          string result = str.substr(lastSlash + 2, dot - lastSlash - 2);
          cout << "Result: " << result << endl;
        } else {
          cout << "Invalid format." << endl;
        }

        if ( myf ) {
                target->cd();
                myf->Write(str.substr(lastSlash + 2, dot - lastSlash - 2).c_str());
                //obj->Write(key->GetName());
                 //try to write the object to our new
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
            // Find the position of the last '/' character
            size_t lastSlash = str.find_last_of('s');
            size_t dot = str.find('.', lastSlash + 1);
            if (lastSlash != string::npos && dot != string::npos) {
              string result = str.substr(lastSlash + 2, dot - lastSlash - 2);
              cout << "Result: " << result << endl;
            } else {
              cout << "Invalid format." << endl;
            }
            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(c1->GetName());
            if (key2) {
                TObject *obj2 = key2->ReadObj();
                if ( obj->IsA()->InheritsFrom( TCanvas::Class() ) ) {
                  TCanvas *c2 = (TCanvas*)obj2;
                  cout << c2->GetTitle() << endl;
                  TObject * thisf2=c2->FindObject("fProfile_0");
                  TF1 * myf2 = (TF1*) thisf2;
                  cout << myf2->GetName() << endl;
                  // clone it>
                  target->cd();
                  myf2->Write(str.substr(lastSlash + 2, dot - lastSlash - 2).c_str());
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
   // save modifications to target file
   target->SaveSelf(kTRUE);
}
