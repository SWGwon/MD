#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLeaf.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>

void calculate_charge(std::string fileName, int numChannels = -1) {
    // --- User Configuration ---
    double dynamicRange = 2.0;      // Volts
    int    numBits      = 14;       // ADC Resolution
    double resistance   = 50.0;     // Ohms
    double samplingTime = 2.0e-9;   // Seconds (2 ns)
    UInt_t maxADC       = (1 << numBits) - 1;

    // Open Input File
    TFile *fIn = TFile::Open(fileName.c_str());
    if (!fIn || fIn->IsZombie()) return;
    TTree *TIn = (TTree*)fIn->Get("T");

    // Determine bins and channels dynamically from the ADC branch title
    TLeaf *lADC = TIn->GetLeaf("ADC");
    if (!lADC) {
        std::cerr << "Error: Branch 'ADC' not found in " << fileName << std::endl;
        fIn->Close();
        return;
    }

    TString title = lADC->GetTitle(); // e.g., "ADC[8][100]"
    int fileChannels = 1;
    int bins = 1;

    if (title.Contains("[") && title.Contains("]")) {
        sscanf(title.Data(), "ADC[%d][%d]", &fileChannels, &bins);
    } else {
        fileChannels = 1; 
        bins = lADC->GetLen();
    }

    if (numChannels <= 0 || numChannels > fileChannels) {
        numChannels = fileChannels;
    }

    std::cout << "File has " << fileChannels << " channels with " << bins << " bins each." << std::endl;
    std::cout << "Processing " << numChannels << " channels..." << std::endl;

    int totalBinsInFile = lADC->GetLen();
    std::vector<UInt_t> adcData(totalBinsInFile);
    TIn->SetBranchAddress("ADC", adcData.data());

    UInt_t eventNumIn;
    if (TIn->GetBranch("EventNumber")) {
        TIn->SetBranchAddress("EventNumber", &eventNumIn);
    } else {
        eventNumIn = 0;
    }

    // Create Output File and Tree
    TString outFileName = "output_charge_" + TString(fileName);
    TFile *fOut = new TFile(outFileName, "RECREATE");
    TTree *TOut = new TTree("T_Charge", "Comprehensive Pulse Analysis");
    
    // Output Variables
    std::vector<double> charge_pC(numChannels);
    std::vector<double> amplitude_mV(numChannels);
    std::vector<double> peakTime_ns(numChannels);
    std::vector<double> riseTime_ns(numChannels);
    std::vector<double> fallTime_ns(numChannels);
    std::vector<double> fwhm_ns(numChannels);
    std::vector<double> baseline_ADC(numChannels);
    std::vector<double> baselineRMS_ADC(numChannels);
    std::vector<int>    isSaturated(numChannels);
    UInt_t eventNumOut;

    TOut->Branch("eventNum",       &eventNumOut,          "eventNum/i");
    TOut->Branch("charge_pC",      charge_pC.data(),      Form("charge_pC[%d]/D", numChannels));
    TOut->Branch("amplitude_mV",   amplitude_mV.data(),   Form("amplitude_mV[%d]/D", numChannels));
    TOut->Branch("peakTime_ns",    peakTime_ns.data(),    Form("peakTime_ns[%d]/D", numChannels));
    TOut->Branch("riseTime_ns",    riseTime_ns.data(),    Form("riseTime_ns[%d]/D", numChannels));
    TOut->Branch("fallTime_ns",    fallTime_ns.data(),    Form("fallTime_ns[%d]/D", numChannels));
    TOut->Branch("fwhm_ns",        fwhm_ns.data(),        Form("fwhm_ns[%d]/D", numChannels));
    TOut->Branch("baseline_ADC",   baseline_ADC.data(),   Form("baseline_ADC[%d]/D", numChannels));
    TOut->Branch("baselineRMS_ADC",baselineRMS_ADC.data(), Form("baselineRMS_ADC[%d]/D", numChannels));
    TOut->Branch("isSaturated",    isSaturated.data(),    Form("isSaturated[%d]/I", numChannels));

    // Conversion factors
    double adcToVolts = dynamicRange / (double)maxADC;
    double unitChargeFactor = (adcToVolts * samplingTime / resistance) * 1e12; // pC
    double adcToMillivolts = adcToVolts * 1000.0;
    double samplingTime_ns = samplingTime * 1e9;

    // Average Waveform Histograms
    std::vector<TH1D*> hAvgWave(numChannels);
    std::vector<TH1D*> hCharge(numChannels);
    for (int ch = 0; ch < numChannels; ++ch) {
        hAvgWave[ch] = new TH1D(Form("hAvgWave_Ch%d", ch), 
                               Form("Average Waveform Ch %d;Time [ns];Amplitude [ADC]", ch), 
                               bins, 0, bins * samplingTime_ns);
        hCharge[ch] = new TH1D(Form("hCharge_Ch%d", ch), 
                               Form("Charge Spectrum Ch %d;Charge [pC];Counts", ch), 
                               600, -1, 5);
        hAvgWave[ch]->SetDirectory(0);
        hCharge[ch]->SetDirectory(0);
    }

    Long64_t nEntries = TIn->GetEntries();
    std::vector<UInt_t> samples(bins);

    std::cout << "Starting processing " << nEntries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        TIn->GetEntry(i);
        eventNumOut = eventNumIn;

        for (int ch = 0; ch < numChannels; ++ch) {
            int offset = ch * bins;
            isSaturated[ch] = 0;

            // 1. Baseline calculation (Median)
            for (int s = 0; s < bins; ++s) {
                samples[s] = adcData[offset + s];
                if (samples[s] >= maxADC - 1) isSaturated[ch] = 1;
            }
            std::vector<UInt_t> sortedSamples = samples;
            std::nth_element(sortedSamples.begin(), sortedSamples.begin() + bins / 2, sortedSamples.end());
            double baseline = (double)sortedSamples[bins / 2];
            baseline_ADC[ch] = baseline;

            // RMS calculation
            double sumSq = 0;
            for (int s = 0; s < bins; ++s) {
                sumSq += std::pow((double)samples[s] - baseline, 2);
            }
            baselineRMS_ADC[ch] = std::sqrt(sumSq / bins);

            // 2. Find Peak Position & Amplitude
            double maxAmpADC = -1.0e9;
            int peakIdx = -1;
            for (int s = 0; s < bins; ++s) {
                double pulseVal = baseline - (double)samples[s];
                if (pulseVal > maxAmpADC) {
                    maxAmpADC = pulseVal;
                    peakIdx = s;
                }
                hAvgWave[ch]->Fill(s * samplingTime_ns, pulseVal); // Accumulate for average
            }

            amplitude_mV[ch] = maxAmpADC * adcToMillivolts;
            peakTime_ns[ch]  = peakIdx * samplingTime_ns;

            // 3. Pulse Shape Analysis (Rise/Fall Time, FWHM)
            double t10 = -1, t90 = -1, t50_1 = -1, t50_2 = -1, t90_fall = -1, t10_fall = -1;
            for (int s = 0; s < bins; ++s) {
                double val = baseline - (double)samples[s];
                double time = s * samplingTime_ns;
                
                if (t10 < 0 && val >= 0.1 * maxAmpADC) t10 = time;
                if (t50_1 < 0 && val >= 0.5 * maxAmpADC) t50_1 = time;
                if (t90 < 0 && val >= 0.9 * maxAmpADC) t90 = time;
                
                if (s > peakIdx) {
                    if (t90_fall < 0 && val <= 0.9 * maxAmpADC) t90_fall = time;
                    if (t50_2 < 0 && val <= 0.5 * maxAmpADC) t50_2 = time;
                    if (t10_fall < 0 && val <= 0.1 * maxAmpADC) t10_fall = time;
                }
            }
            riseTime_ns[ch] = (t90 > 0 && t10 > 0) ? (t90 - t10) : 0;
            fallTime_ns[ch] = (t10_fall > 0 && t90_fall > 0) ? (t10_fall - t90_fall) : 0;
            fwhm_ns[ch]     = (t50_2 > 0 && t50_1 > 0) ? (t50_2 - t50_1) : 0;

            // 4. Integrate Charge
            double totalCharge = 0;
            for (int s = 0; s < bins; ++s) {
                double temp_amp = baseline - (double)samples[s];
                //if (temp_amp > 0)
                    totalCharge += temp_amp * unitChargeFactor;
            }

            charge_pC[ch] = totalCharge;
            hCharge[ch]->Fill(totalCharge);
        }
        TOut->Fill();

        if (i > 0 && i % 100000 == 0) std::cout << "Processing event " << i << " / " << nEntries << "..." << std::endl;
    }

    // Scale Average Waveforms
    for (int ch = 0; ch < numChannels; ++ch) {
        hAvgWave[ch]->Scale(1.0 / nEntries);
    }

    // --- Visualization ---
    TCanvas *c1 = new TCanvas("c1", "Charge Spectra", 1200, 800);
    if (numChannels <= 1) c1->Divide(1, 1);
    else if (numChannels <= 2) c1->Divide(2, 1);
    else if (numChannels <= 4) c1->Divide(2, 2);
    else if (numChannels <= 6) c1->Divide(3, 2);
    else c1->Divide(4, 2);

    for (int ch = 0; ch < numChannels; ++ch) {
        c1->cd(ch + 1);
        gPad->SetLogy();
        hCharge[ch]->SetLineColor(kBlue + (ch % 4));
        hCharge[ch]->Draw();
    }

    TString plotName = TString(fileName);
    plotName.ReplaceAll(".root", Form("_ch%d_charge.png", numChannels));
    c1->SaveAs(plotName);

    TCanvas *c2 = new TCanvas("c2", "Average Waveforms", 1200, 800);
    if (numChannels <= 1) c2->Divide(1, 1);
    else if (numChannels <= 2) c2->Divide(2, 1);
    else if (numChannels <= 4) c2->Divide(2, 2);
    else if (numChannels <= 6) c2->Divide(3, 2);
    else c2->Divide(4, 2);

    for (int ch = 0; ch < numChannels; ++ch) {
        c2->cd(ch + 1);
        hAvgWave[ch]->SetLineColor(kRed + (ch % 4));
        hAvgWave[ch]->Draw("HIST");
    }
    TString avgPlotName = TString(fileName);
    avgPlotName.ReplaceAll(".root", Form("_ch%d_avg_wave.png", numChannels));
    c2->SaveAs(avgPlotName);
    
    // Save data and close
    fOut->cd();
    TOut->Write();
    for (int ch = 0; ch < numChannels; ++ch) {
        hCharge[ch]->Write();
        hAvgWave[ch]->Write();
    }
    fOut->Close();
    fIn->Close();

    std::cout << "Analysis Complete!" << std::endl;
    std::cout << "1. Data saved to: " << outFileName << std::endl;
    std::cout << "2. Charge plots saved to: " << plotName << std::endl;
    std::cout << "3. Avg waveform plots saved to: " << avgPlotName << std::endl;
}
