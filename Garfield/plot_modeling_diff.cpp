#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <cmath>

double FieldToVelocity(double E) {
    double slope = 0.008855515147092161;
    double intercept = -0.40797048541793746;
    double v_mm_per_us = (E * slope) + intercept;
    return v_mm_per_us * 0.1; // Convert mm/us to cm/us
}

// Struct to pass back both the true integrated time and the final landed Z position
struct DriftResult {
    double total_time;
    double final_z;
};

// Computes 2D trajectory tracking from start_d down to 0 along the y-axis
DriftResult GetIntegratedDriftProperties(TFile* file_plane, int start_d, double initial_z) {
    double total_time = 0.0;
    double current_z = initial_z;
    
    int sign = (start_d > 0) ? 1 : -1;
    int steps = std::abs(start_d);
    double dy = 1.0; // 1.0 cm slices

    // Step from the starting distance down to the mesh at 0
    for (int i = steps; i > 0; --i) {
        int current_slice = sign * i;
        
        TString ey_name = TString::Format("Ey_vs_z_%d", current_slice);
        TString ez_name = TString::Format("Ez_vs_z_%d", current_slice);
        
        TGraph* g_ey = (TGraph*)file_plane->Get(ey_name);
        TGraph* g_ez = (TGraph*)file_plane->Get(ez_name);
        
        if (!g_ey || !g_ez) continue;

        // Evaluate the local fields at our current, dynamically shifting z-position
        double E_y = g_ey->Eval(current_z);
        double E_z = g_ez->Eval(current_z);

        double v_y = FieldToVelocity(E_y);
        double v_z = FieldToVelocity(E_z);

        if (std::abs(v_y) > 0) {
            // 1. dt is determined completely by the Ey component moving along dy
            double dt = dy / std::abs(v_y); 
            total_time += dt;

            // 2. The amount of z displacement caused by the transverse field component
            double dz = v_z * dt;
            current_z += dz; 
        }
    }

    DriftResult result = {total_time, current_z};
    return result;
}

void plot_modeling_diff() {

    TFile *file1 = TFile::Open("Tilt_0um_newdiameters.root", "READ");
    TFile *file2 = TFile::Open("Tilt_0um_electrode.root", "READ");

    if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "2D Trajectory Tracker Error", 1200, 600);
    c1->Divide(2, 1);

    for (int d = 7; d >= -7; --d) {
        if (d == 0) continue;

        TString ey_name = TString::Format("Ey_vs_z_%d", d);
        TString ez_name = TString::Format("Ez_vs_z_%d", d);

        TGraph *g_ey_wire = (TGraph*)((TGraph*)file1->Get(ey_name))->Clone();
        TGraph *g_ez_wire = (TGraph*)((TGraph*)file1->Get(ez_name))->Clone();
        TGraph *g_ey_plane = (TGraph*)((TGraph*)file2->Get(ey_name))->Clone();
        TGraph *g_ez_plane = (TGraph*)((TGraph*)file2->Get(ez_name))->Clone();

        if (!g_ey_wire || !g_ey_plane || !g_ez_wire || !g_ez_plane) continue;

        int n_points = g_ey_wire->GetN();
        TGraph *g_dy_res = new TGraph(n_points);
        TGraph *g_dz_res = new TGraph(n_points);

        // For every initial z coordinate along the wire mesh grid
        for (int i = 0; i < n_points; ++i) {
            double z_initial, E_wire_y, E_wire_z;
            g_ey_wire->GetPoint(i, z_initial, E_wire_y);
            g_ez_wire->GetPoint(i, z_initial, E_wire_z);

            // Calculate the actual drift time and path deformation using the plane file reference
            DriftResult plane_drift = GetIntegratedDriftProperties(file2, d, z_initial);

            // Convert raw field locally to velocity for the wire model at this slice
            double v_wire_y = FieldToVelocity(E_wire_y);
            double v_wire_z = FieldToVelocity(E_wire_z);


           // Compare to plane model values at the final landed position
           double E_plane_y = g_ey_plane->Eval(plane_drift.final_z); // This is an E-field!
           double E_plane_z = g_ez_plane->Eval(plane_drift.final_z); // This is an E-field!

           // Convert those interpolated fields to velocities
           double v_plane_y = FieldToVelocity(E_plane_y);
           double v_plane_z = FieldToVelocity(E_plane_z);

            // Spatial tracking residual calculations
            double delta_y = (v_wire_y - v_plane_y) * plane_drift.total_time;
            double delta_z = (v_wire_z - v_plane_z) * plane_drift.total_time;

            g_dy_res->SetPoint(i, z_initial, delta_y);
            g_dz_res->SetPoint(i, z_initial, delta_z);
        }

        // --- Plotting ---
        c1->cd(1);
        g_dy_res->SetMarkerStyle(20); g_dy_res->SetMarkerColor(kBlack); g_dy_res->SetLineColor(kBlack);
        g_dy_res->SetTitle(TString::Format("#Delta y Spatial Error (Drift = %d cm);Initial z [cm];#Delta y [cm]", d));
        g_dy_res->Draw("APL");

        c1->cd(2);
        g_dz_res->SetMarkerStyle(20); g_dz_res->SetMarkerColor(kRed); g_dz_res->SetLineColor(kRed);
        g_dz_res->SetTitle(TString::Format("#Delta z Spatial Error (Drift = %d cm);Initial z [cm];#Delta z [cm]", d));
        g_dz_res->Draw("APL");

        c1->Update();

        TString sign = (d > 0) ? "plus" : "minus";
        c1->SaveAs(TString::Format("trajectory_error_%s_%02d.png", sign.Data(), std::abs(d)));

        delete g_ey_wire; delete g_ey_plane; delete g_ez_wire; delete g_ez_plane;
        delete g_dy_res; delete g_dz_res;
    }

    file1->Close(); file2->Close();
}
