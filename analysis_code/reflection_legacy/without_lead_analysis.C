#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TCanvas.h>

/**
 * Convert charge in pC to Number of Photo-Electrons (NPE) using raw Gain.
 * Formula: NPE = Q [pC] / (Gain * 1.602e-7)
 */
double charge_to_NPE(double charge_pC, double gain) {
    if (gain <= 0) return 0;
    const double e_charge_pC = 1.60217663e-7; // Charge of 1 PE in pC for Gain=1
    return charge_pC / (gain * e_charge_pC);
}

/**
 * Event selection logic based on NPE.
 * Returns true if both channels are above their respective NPE thresholds.
 */
bool is_selected_npe(double npe0, double npe1, double thr0_npe, double thr1_npe) {
    return (npe0 > thr0_npe && npe1 > thr1_npe);
}

/**
 * Main analysis function.
 */
void analysis() {
    const char* input_file = "output.root";
    
    TFile *f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open " << input_file << std::endl;
        return;
    }

    // --- Hardcoded Gain Values ---
    const double gain[2] = {1.0E+7, 1.0E+7}; 

    // Initialize TTreeReader
    TTreeReader reader("T_Charge", f);
    
    TTreeReaderValue<unsigned int> eventNum(reader, "eventNum");
    TTreeReaderArray<double> charge_pC(reader, "charge_pC");

    // --- Thresholds in NPE ---
    double threshold_npe_ch0 = 1.0;
    double threshold_npe_ch1 = 1.0;

    TH1D* hist_npe[2];
    for (int i = 0; i < 2; ++i) {
        hist_npe[i] = new TH1D(Form("h_npe_ch%d", i), Form("NPE Ch%d (Gain=%.1E);NPE;Counts", i, gain[i]), 400, -10, 90);
    }

    TH1D* hist_total = new TH1D("Total NPE", "Total NPE", 400, -10, 90);

    std::cout << ">>> Starting Analysis Loop (NPE-based selection applied) <<<" << std::endl;
    
    long long entryCount = 0;
    long long selectedCount = 0;

    while (reader.Next()) {
        entryCount++;

        // 1. pC -> NPE 변환 선행
        double npe[2];
        for (int i = 0; i < 2; ++i) {
            npe[i] = charge_to_NPE(charge_pC[i], gain[i]);
        }

        // 2. NPE 기반 이벤트 선택
        if (!is_selected_npe(npe[0], npe[1], threshold_npe_ch0, threshold_npe_ch1)) {
            continue;
        }
        
        selectedCount++;

        // 3. 히스토그램 채우기
        hist_npe[0]->Fill(npe[0]);
        hist_npe[1]->Fill(npe[1]);

        hist_total->Fill(npe[0] + npe[1]);

        if (selectedCount <= 10) {
            std::cout << "Event " << *eventNum << ": Ch0 NPE = " << npe[0] 
                      << ", Ch1 NPE = " << npe[1] << " (Selected by NPE > 1.0)" << std::endl;
        }
    }

    std::cout << ">>> Analysis Complete <<<" << std::endl;
    std::cout << "Total entries: " << entryCount << std::endl;
    std::cout << "Selected entries: " << selectedCount << std::endl;

    // Visualization
    TCanvas* c1 = new TCanvas("c1", "NPE Distribution", 800, 600);
    c1->SetLogy();
    hist_npe[0]->Draw();
    hist_npe[1]->SetLineColor(kRed);
    hist_npe[1]->Draw("same");

    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hist_npe[0], "Channel 0", "l");
    leg->AddEntry(hist_npe[1], "Channel 1", "l");
    leg->Draw();

    TCanvas* c2 = new TCanvas("c2", "NPE Distribution", 800, 600);
    c2->SetLogy();
    hist_total->Draw();

    // Save histograms
    TFile *f_out = new TFile("output_analysis.root", "RECREATE");
    for (int i = 0; i < 2; ++i) {
        hist_npe[i]->Write();
    }
    std::cout << "Results saved to output_analysis.root" << std::endl;
}
