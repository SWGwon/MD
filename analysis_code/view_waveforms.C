#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <iostream>

void view_waveforms(std::string fileName) {
    TFile *f = TFile::Open(fileName.c_str());
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");

    int bins = 100;
    
    UInt_t ADC[8][bins];
    T->SetBranchAddress("ADC", ADC);

    std::cout << "========================================" << std::endl;
    std::cout << "  Waveform Viewer (Event-by-Event)      " << std::endl;
    std::cout << "  - Close the window to see NEXT event  " << std::endl;
    std::cout << "  - Press Ctrl+C in terminal to EXIT    " << std::endl;
    std::cout << "========================================" << std::endl;

    for (Long64_t i = 0; i < T->GetEntries(); ++i) {
        T->GetEntry(i);
        
        // Use a fixed name for the canvas to track it
        TCanvas *c1 = new TCanvas("c_view", Form("Waveform Viewer - Event %lld", i), 900, 500);
        
        TH1D *hWave = new TH1D("hWave_tmp", Form("Event %lld, Channel 0;Sample;ADC Value", i), bins, 0, bins);
        for (int s = 0; s < bins; ++s) {
            hWave->SetBinContent(s + 1, ADC[0][s]);
        }
        
        hWave->SetLineColor(kBlue + 1);
        hWave->SetLineWidth(2);
        hWave->Draw("HIST");
        c1->Update();

        // Loop until the canvas "c_view" is closed by the user
        while (gROOT->GetListOfCanvases()->FindObject("c_view")) {
            gSystem->ProcessEvents(); // Keep GUI responsive
            gSystem->Sleep(20);       // Reduce CPU usage
        }

        // Clean up before next event
        delete hWave;
        std::cout << "Proceeding to next event..." << std::endl;
    }
}
