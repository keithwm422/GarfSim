#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <algorithm>


void TransformGraph(TGraph* g) {

  double slope = 0.008855515147092161;
  double intercept = -0.40797048541793746;
  // at each E field value, do v = E*slope+intercept

    if (!g) return;

    int n = g->GetN();
    double x, y;

    for (int i = 0; i < n; i++) {
        // 1. Get current coordinates
        g->GetPoint(i, x, y);
        
        // 2. Apply linear transformation: y' = (y * slope) + intercept
        double newY = (y * slope) + intercept;
        
        // 3. Update the point (X remains unchanged)
        g->SetPoint(i, x, newY);
    }
}

void TransformAndSubtract(TGraph* target, TGraph* reference) {
    if (!target || !reference) {
        printf("Error: Null pointer passed to TransformAndSubtract\n");
        return;
    }

    int n = target->GetN();
    if (n != reference->GetN()) {
        printf("Warning: Graphs have different numbers of points (%d vs %d)\n", n, reference->GetN());
        // We'll proceed using the smaller count to avoid crashes
        n = std::min(n, reference->GetN());
    }

    double x_t, y_t, x_r, y_r;

    for (int i = 0; i < n; i++) {
        target->GetPoint(i, x_t, y_t);
        reference->GetPoint(i, x_r, y_r);

        // 1. Transform the target value
        // 2. Subtract the reference value
        // Formula: y' = (y_target * slope + intercept) - y_reference
        double newY = y_t - y_r;

        target->SetPoint(i, x_t, newY);
    }
}


void Plot_compareWires(){
  //BFieldCoeff_wireID_29_10.0_.root
  std::stringstream TF1name; // for accessing iso_inverse files
  std::stringstream histoname; // for writing histogram and TGraph names
  TFile * _filem[4]; // input file we iterate over
  double tilt=0.0;
  int tilt_iter=0;
  TGraph * theseOG[4][7]; // 0um, 100, 200, 400
  while(tilt_iter<2){
    TF1name.str("");
    if(tilt_iter==0)TF1name << "Tilt_0um_electrode.root"; //Tilt_0um.root
    if(tilt_iter==1)TF1name << "Tilt_0um_electrode_oldwires.root"; //Tilt_0um.root
    cout << "opening file: " << TF1name.str().c_str() << endl;
    _filem[tilt_iter] = TFile::Open(TF1name.str().c_str());
    if(!_filem[tilt_iter] || !_filem[tilt_iter]->IsOpen() || _filem[tilt_iter]->IsZombie()){
      return;
    }
    // Ey_vs_z_-6
    gROOT->cd();
    double y_dist = 7.0;
    int iter=0;
    while (y_dist>0){
      histoname.str("");
      histoname << "Ey_vs_z_" << std::fixed << std::setprecision(0) << y_dist;
      theseOG[tilt_iter][iter] = (TGraph *)_filem[tilt_iter]->Get(histoname.str().c_str()); // iso_inverse_62
      TransformGraph(theseOG[tilt_iter][iter]);
      iter++;
      y_dist -=1.0;
    }
   // now we load in the other guys
     if(tilt_iter==0) tilt=0;
     else if(tilt_iter==1) tilt=100.0;
     else if(tilt_iter==2) tilt=400.0;
     tilt_iter++;
  }

  // 4. Set "Color Schemes" (Line and Marker Styles)
  TCanvas *c1 = new TCanvas("c1", "MultiGraph Example", 800, 600);

    // Graph 1: Red
    theseOG[0][0]->SetLineColor(kRed);
    theseOG[0][0]->SetMarkerColor(kRed);
    theseOG[0][0]->SetMarkerStyle(20); // Circle

    // Graph 2: Blue
    theseOG[1][0]->SetLineColor(kBlue);
    theseOG[1][0]->SetMarkerColor(kBlue);
    theseOG[1][0]->SetMarkerStyle(21); // Square

    // Graph 3: Green (using +2 to make it darker/readable)
    /*theseOG[2][0]->SetLineColor(kGreen + 2);
    theseOG[2][0]->SetMarkerColor(kGreen + 2);
    theseOG[2][0]->SetMarkerStyle(22); // Triangle
    */

    // 5. Create TMultiGraph and add graphs
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(theseOG[0][0], "lp"); // "lp" means draw line and points
    mg->Add(theseOG[1][0], "lp");
    //mg->Add(theseOG[2][0], "lp");

    // Set Global Titles (must be done before or during Draw)
    mg->SetTitle("Comparison of Datasets;z(cm);v(mm/us)");

    // 6. Draw everything
    mg->Draw("A"); // "A" draws the axes first

    // 7. Add a Legend
    TLegend *leg = new TLegend(0.15, 0.7, 0.35, 0.85); // (x1, y1, x2, y2) in NDC coordinates
    leg->AddEntry(theseOG[0][0], "Elec-new", "lp");
    leg->AddEntry(theseOG[1][0], "Elec-old", "lp");
    //leg->AddEntry(theseOG[2][0], "Elec-old", "lp");
    leg->Draw();

    // 8. Update and Save
    c1->SetGrid();
    c1->Update();
    c1->SaveAs("MeshCompareElectrode_oldwires.pdf");
    _filem[0]->Close();
    _filem[1]->Close();
  
}

