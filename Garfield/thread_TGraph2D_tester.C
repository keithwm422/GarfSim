#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <random>

#include "TFile.h"
#include "TGraph2D.h"
#include "TROOT.h"

// A mutex to protect access to the shared vector of TGraph2Ds.
std::mutex g_mutex;
std::vector<TGraph2D*> g_graphs;

// The function to be executed by each thread.
void process_data(int thread_id, const std::vector<std::vector<double>>& data) {
    // Create a thread-local TGraph2D.
    TGraph2D* graph = new TGraph2D();

    // Add points from the data vector.
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i].size() == 3) {
            graph->SetPoint(i, data[i][0], data[i][1], data[i][2]);
        }
    }

    std::cout << "Thread " << thread_id << " created TGraph2D with " << graph->GetN() << " points." << std::endl;

    // Safely add the thread-local graph to the shared vector.
    // Lock the mutex before modifying the shared resource.
    {
        std::lock_guard<std::mutex> lock(g_mutex);
        g_graphs.push_back(graph);
    }
}

int main() {
    // Enable ROOT's thread safety.
    ROOT::EnableThreadSafety();

    const int num_threads = 4;
    const int points_per_thread = 100;
    std::vector<std::thread> threads;
    std::vector<std::vector<std::vector<double>>> all_data(num_threads);

    // Seed a random number generator for example data.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 10.0);

    // Prepare example data for each thread.
    for (int i = 0; i < num_threads; ++i) {
        for (int j = 0; j < points_per_thread; ++j) {
            all_data[i].push_back({dis(gen), dis(gen), dis(gen)});
        }
    }

    // Start the worker threads.
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(process_data, i, std::ref(all_data[i]));
    }

    // Wait for all threads to finish.
    for (auto& t : threads) {
        t.join();
    }

    std::cout << "\nAll threads have completed. Merging and writing to file." << std::endl;

    // --- Recombination and file writing section (single-threaded) ---

    // Open a ROOT file for writing from the main thread.
    TFile* file = new TFile("graphs.root", "RECREATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Write all collected TGraph2Ds to the file.
    {
        std::lock_guard<std::mutex> lock(g_mutex); // Ensure no threads are still accessing the vector.
        for (size_t i = 0; i < g_graphs.size(); ++i) {
            g_graphs[i]->SetName(TString::Format("graph_%zu", i));
            g_graphs[i]->Write();
        }
    }

    file->Close();
    std::cout << "Successfully wrote all TGraph2Ds to graphs.root." << std::endl;

    // Clean up allocated TGraph2D objects.
    for (TGraph2D* graph : g_graphs) {
        delete graph;
    }

    return 0;
}
