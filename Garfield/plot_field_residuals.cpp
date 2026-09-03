#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <cmath>

void plot_field_residuals() {
    // 1. Open the two ROOT files

    TFile *file1 = TFile::Open("Tilt_0um_newdiameters.root", "READ");
    TFile *file2 = TFile::Open("Tilt_0um_electrode.root", "READ");

    if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie()) {
        std::cerr << "Error: One or both ROOT files could not be opened!" << std::endl;
        return;
    }

    // 2. Create canvas with 1x2 split
    TCanvas *c1 = new TCanvas("c1", "E-Field Residuals", 1200, 600);
    c1->Divide(2, 1);

    // 3. Loop over drift distances from +7 to -7 (skipping 0)
    for (int d = 7; d >= -7; --d) {
        if (d == 0) continue;

        TString ey_name = TString::Format("Ey_vs_z_%d", d);
        TString ez_name = TString::Format("Ez_vs_z_%d", d);

        TGraph *g_ey_wire = (TGraph*)file1->Get(ey_name);
        TGraph *g_ez_wire = (TGraph*)file1->Get(ez_name);

        TGraph *g_ey_plane = (TGraph*)file2->Get(ey_name);
        TGraph *g_ez_plane = (TGraph*)file2->Get(ez_name);

        if (!g_ey_wire || !g_ey_plane || !g_ez_wire || !g_ez_plane) {
            std::cerr << "Warning: Missing graphs for drift distance " << d << " cm. Skipping..." << std::endl;
            continue;
        }

        // --- Calculate Ey Residuals ---
        int n_points_ey = g_ey_wire->GetN();
        TGraph *g_ey_diff = new TGraph(n_points_ey);
        
        for (int i = 0; i < n_points_ey; ++i) {
            double z, e_wire;
            g_ey_wire->GetPoint(i, z, e_wire);
            
            // Interpolate the plane graph value at this exact z position
            double e_plane = g_ey_plane->Eval(z); 
            double diff = e_wire - e_plane;
            
            g_ey_diff->SetPoint(i, z, diff);
        }

        // --- Calculate Ez Residuals ---
        int n_points_ez = g_ez_wire->GetN();
        TGraph *g_ez_diff = new TGraph(n_points_ez);
        
        for (int i = 0; i < n_points_ez; ++i) {
            double z, e_wire;
            g_ez_wire->GetPoint(i, z, e_wire);
            
            // Interpolate the plane graph value at this exact z position
            double e_plane = g_ez_plane->Eval(z); 
            double diff = e_wire - e_plane;
            
            g_ez_diff->SetPoint(i, z, diff);
        }

        // --- Plot Ey Residuals (Left Subplot) ---
        c1->cd(1);
        g_ey_diff->SetMarkerStyle(20);
        g_ey_diff->SetMarkerColor(kBlack);
        g_ey_diff->SetLineColor(kBlack);
        
        TString title_ey = TString::Format("#DeltaE_{y} vs z (Drift = %d cm);z [cm];#DeltaE_{y} (Wire - Plane) [V/cm]", d);
        g_ey_diff->SetTitle(title_ey);
        g_ey_diff->Draw("APL"); // Draws axes, points, and connecting lines to trace patterns

        // --- Plot Ez Residuals (Right Subplot) ---
        c1->cd(2);
        g_ez_diff->SetMarkerStyle(20);
        g_ez_diff->SetMarkerColor(kRed);
        g_ez_diff->SetLineColor(kRed);
        
        TString title_ez = TString::Format("#DeltaE_{z} vs z (Drift = %d cm);z [cm];#DeltaE_{z} (Wire - Plane) [V/cm]", d);
        g_ez_diff->SetTitle(title_ez);
        g_ez_diff->Draw("APL");

        // Force canvas render update
        c1->Update();

        // 4. Save frame as PNG
        TString sign = (d > 0) ? "plus" : "minus";
        TString out_name = TString::Format("diff_frame_%s_%02d.png", sign.Data(), std::abs(d));
        c1->SaveAs(out_name);

        // Clean up current frame's newly allocated TGraphs to avoid RAM leaks
        delete g_ey_diff;
        delete g_ez_diff;
    }

    // Clean up open files
    file1->Close();
    file2->Close();
}
