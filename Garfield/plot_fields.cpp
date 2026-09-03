#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

void plot_fields() {
    // 1. Open the two ROOT files (Update these paths to your actual file names)
    TFile *file1 = TFile::Open("Tilt_0um_newdiameters.root", "READ");
    TFile *file2 = TFile::Open("Tilt_0um_electrode.root", "READ");
 // Tilt_0um_electrode.root           Tilt_0um_newdiameters.root
    if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie()) {
        std::cerr << "Error: One or both ROOT files could not be opened!" << std::endl;
        return;
    }

    // 2. Create a canvas with a 1x2 split subplot layout
    TCanvas *c1 = new TCanvas("c1", "E-Field Comparison", 1200, 600);
    c1->Divide(2, 1);

    // 3. Loop over drift distances from +7 down to -7
    for (int d = 7; d >= -7; --d) {
        if (d == 0) continue; // Skip 0 as requested

        // Format object names based on the current drift distance
        TString ey_name = TString::Format("Ey_vs_z_%d", d);
        TString ez_name = TString::Format("Ez_vs_z_%d", d);

        // Fetch TGraphs from File 1
        TGraph *g_ey1 = (TGraph*)file1->Get(ey_name);
        TGraph *g_ez1 = (TGraph*)file1->Get(ez_name);

        // Fetch TGraphs from File 2
        TGraph *g_ey2 = (TGraph*)file2->Get(ey_name);
        TGraph *g_ez2 = (TGraph*)file2->Get(ez_name);

        // Quick sanity check to ensure the graphs exist for this iteration
        if (!g_ey1 || !g_ey2 || !g_ez1 || !g_ez2) {
            std::cerr << "Warning: Missing graphs for drift distance " << d << " cm. Skipping..." << std::endl;
            continue;
        }

        // --- Style File 1 Graphs (Black Markers) ---
        g_ey1->SetMarkerStyle(20); g_ey1->SetMarkerColor(kBlack); g_ey1->SetLineColor(kBlack);
        g_ez1->SetMarkerStyle(20); g_ez1->SetMarkerColor(kBlack); g_ez1->SetLineColor(kBlack);

        // --- Style File 2 Graphs (Red Markers) ---
        g_ey2->SetMarkerStyle(24); g_ey2->SetMarkerColor(kRed);   g_ey2->SetLineColor(kRed);
        g_ez2->SetMarkerStyle(24); g_ez2->SetMarkerColor(kRed);   g_ez2->SetLineColor(kRed);

        // 4. Setup Left Subplot (Ey vs z)
        c1->cd(1);
        TMultiGraph *mg_ey = new TMultiGraph();
        TString title_ey = TString::Format("E_{y} vs z (Drift Distance = %d cm);z [cm];E_{y} [V/cm]", d);
        mg_ey->SetTitle(title_ey);
        
        mg_ey->Add(g_ey1, "AP"); // "P" for markers, "L" for lines if you prefer
        mg_ey->Add(g_ey2, "P");
        mg_ey->Draw("A");

        // Add Legend to Left Subplot
        TLegend *leg_ey = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg_ey->AddEntry(g_ey1, "Wire Model", "p");
        leg_ey->AddEntry(g_ey2, "Plane Model", "p");
        leg_ey->Draw();

        // 5. Setup Right Subplot (Ez vs z)
        c1->cd(2);
        TMultiGraph *mg_ez = new TMultiGraph();
        TString title_ez = TString::Format("E_{z} vs z (Drift Distance = %d cm);z [cm];E_{z} [V/cm]", d);
        mg_ez->SetTitle(title_ez);
        
        mg_ez->Add(g_ez1, "AP");
        mg_ez->Add(g_ez2, "P");
        mg_ez->Draw("A");

        // Add Legend to Right Subplot
        TLegend *leg_ez = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg_ez->AddEntry(g_ez1, "Wire Model", "p");
        leg_ez->AddEntry(g_ez2, "Plane Model", "p");
        leg_ez->Draw();

        // Update the canvas to force rendering
        c1->Update();

        // 6. Save frame as PNG. Pad the number with a sign to maintain smooth bash sorting
        TString sign = (d > 0) ? "plus" : "minus";
        TString out_name = TString::Format("frame_%s_%02d.png", sign.Data(), std::abs(d));
        c1->SaveAs(out_name);

        // Clean up heap-allocated objects for this iteration to avoid memory leaks
        delete mg_ey;
        delete mg_ez;
        delete leg_ey;
        delete leg_ez;
    }

    // Clean up files
    file1->Close();
    file2->Close();
}
