#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <vector>

void plot_diffVelocities_Grid() {
    // List of index numbers corresponding to your file suffix
    std::vector<int> indices = {1, 2, 4,8}; 
    
    // Distinct color palette for the curves
    std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2};

    auto *mg = new TMultiGraph();
    auto *legend = new TLegend(0.65, 0.65, 0.88, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);

    for (size_t i = 0; i < indices.size(); ++i) {
        int idx = indices[i];
        TString fileName = Form("drift_compareGrid_andAnalytic_%d.root", idx);
        
        TFile *file = TFile::Open(fileName, "READ");
        if (!file || file->IsZombie()) {
            printf("Error: Could not open file %s\n", fileName.Data());
            continue;
        }

        auto *graph = (TGraph*)file->Get("Vy_vs_y_44");
        if (!graph) {
            printf("Error: TGraph 'Vy_vs_y_44' not found in %s\n", fileName.Data());
            file->Close();
            continue;
        }

        // Apply visual styling (cycles through colors)
        int color = colors[i % colors.size()];
        graph->SetLineColor(color);
        graph->SetLineWidth(2);
        graph->SetMarkerColor(color);
        graph->SetMarkerStyle(20 + (i % 5));

        mg->Add(graph, "LP"); // Line + Points
        legend->AddEntry(graph, Form("%d", idx * 100), "lp");
    }

    // Create canvas and apply configuration
    auto *canvas = new TCanvas("c1", "Vy vs y Multigraph", 800, 600);
    canvas->SetGrid();
    canvas->SetLogy();

    mg->Draw("A"); // "A" draws axes around the multigraph
    mg->SetTitle("percent change in Vy vs y Comparison;y (mm); % change in V_{y}");
    
    legend->Draw();
    canvas->Update();
}
