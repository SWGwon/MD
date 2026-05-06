#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TLeaf.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>

// 펄스 분석을 위한 구조체
struct PulseAnalysis {
    double charge_pC;
    double amplitude_mV;
    double peakTime_ns;
    double riseTime_ns;
    double fallTime_ns;
    double fwhm_ns;
    double baseline_ADC;
    double baselineRMS_ADC;
    int isSaturated;
    std::vector<double> waveform; // For average waveform calculation
};

void production(std::string fileName, int numChannels = -1) {
    // 멀티스레딩 활성화
    ROOT::EnableImplicitMT();

    // --- User Configuration ---
    double dynamicRange = 2.0;      // Volts
    int    numBits      = 14;       // ADC Resolution
    double resistance   = 50.0;     // Ohms
    double samplingTime = 2.0e-9;   // Seconds (2 ns)
    UInt_t maxADC       = (1 << numBits) - 1;

    // Open Input File to detect structure
    TFile *fIn = TFile::Open(fileName.c_str());
    if (!fIn || fIn->IsZombie()) return;
    TTree *TIn = (TTree*)fIn->Get("T");
    TLeaf *lADC = TIn->GetLeaf("ADC");
    if (!lADC) {
        std::cerr << "Error: Branch 'ADC' not found in " << fileName << std::endl;
        fIn->Close();
        return;
    }

    TString title = lADC->GetTitle(); 
    int fileChannels = 1, bins = 1;
    if (title.Contains("[") && title.Contains("]")) {
        sscanf(title.Data(), "ADC[%d][%d]", &fileChannels, &bins);
    } else {
        fileChannels = 1; 
        bins = lADC->GetLen();
    }
    if (numChannels <= 0 || numChannels > fileChannels) numChannels = fileChannels;
    fIn->Close(); // RDataFrame will reopen it

    std::cout << "File has " << fileChannels << " channels with " << bins << " bins each." << std::endl;
    std::cout << "Processing " << numChannels << " channels using multi-threading..." << std::endl;

    // Conversion factors
    double adcToVolts = dynamicRange / (double)maxADC;
    double unitChargeFactor = (adcToVolts * samplingTime / resistance) * 1e12; // pC
    double adcToMillivolts = adcToVolts * 1000.0;
    double samplingTime_ns = samplingTime * 1e9;

    // RDataFrame Setup
    ROOT::RDataFrame df("T", fileName);

    // 1. Define analysis logic
    auto processEvent = [=](const ROOT::RVec<UInt_t>& adcData) {
        std::vector<PulseAnalysis> results(numChannels);
        for (int ch = 0; ch < numChannels; ++ch) {
            int offset = ch * bins;
            auto& res = results[ch];
            res.isSaturated = 0;
            res.waveform.resize(bins);

            // Copy and sort for median baseline
            std::vector<UInt_t> sortedSamples(bins);
            for(int s=0; s<bins; ++s) {
                sortedSamples[s] = adcData[offset + s];
                if (sortedSamples[s] >= maxADC - 1) res.isSaturated = 1;
            }
            std::vector<UInt_t> copyForMedian = sortedSamples;
            std::nth_element(copyForMedian.begin(), copyForMedian.begin() + bins / 2, copyForMedian.end());
            double baseline = (double)copyForMedian[bins / 2];
            res.baseline_ADC = baseline;

            // RMS calculation
            double sumSq = 0;
            for (int s = 0; s < bins; ++s) {
                sumSq += std::pow((double)sortedSamples[s] - baseline, 2);
            }
            res.baselineRMS_ADC = std::sqrt(sumSq / bins);

            // Peak Position & Amplitude
            double maxAmpADC = -1.0e9;
            int peakIdx = -1;
            for (int s = 0; s < bins; ++s) {
                double pulseVal = baseline - (double)sortedSamples[s];
                res.waveform[s] = pulseVal;
                if (pulseVal > maxAmpADC) {
                    maxAmpADC = pulseVal;
                    peakIdx = s;
                }
            }
            res.amplitude_mV = maxAmpADC * adcToMillivolts;
            res.peakTime_ns  = peakIdx * samplingTime_ns;

            // Pulse Shape Analysis
            double t10 = -1, t90 = -1, t50_1 = -1, t50_2 = -1, t90_fall = -1, t10_fall = -1;
            for (int s = 0; s < bins; ++s) {
                double val = res.waveform[s];
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
            res.riseTime_ns = (t90 > 0 && t10 > 0) ? (t90 - t10) : 0;
            res.fallTime_ns = (t10_fall > 0 && t90_fall > 0) ? (t10_fall - t90_fall) : 0;
            res.fwhm_ns     = (t50_2 > 0 && t50_1 > 0) ? (t50_2 - t50_1) : 0;

            // Integrate Charge
            double totalCharge = 0;
            for (int s = 0; s < bins; ++s) {
                totalCharge += res.waveform[s] * unitChargeFactor;
            }
            res.charge_pC = totalCharge;
        }
        return results;
    };

    // 2. Apply analysis and extract branches
    auto df_analyzed = df.Define("results", processEvent, {"ADC"});

    // Define individual result columns
    auto df_final = df_analyzed;
    std::vector<std::string> colNames = {
        "charge_pC", "amplitude_mV", "peakTime_ns", "riseTime_ns", 
        "fallTime_ns", "fwhm_ns", "baseline_ADC", "baselineRMS_ADC", "isSaturated"
    };

    // Lambda helpers to extract vector elements from vector of PulseAnalysis
    #define EXTRACT_COL(NAME, TYPE) \
        df_final = df_final.Define(#NAME, [](const std::vector<PulseAnalysis>& res) { \
            std::vector<TYPE> v; v.reserve(res.size()); \
            for(auto& r : res) v.push_back(r.NAME); \
            return v; \
        }, {"results"});

    EXTRACT_COL(charge_pC, double)
    EXTRACT_COL(amplitude_mV, double)
    EXTRACT_COL(peakTime_ns, double)
    EXTRACT_COL(riseTime_ns, double)
    EXTRACT_COL(fallTime_ns, double)
    EXTRACT_COL(fwhm_ns, double)
    EXTRACT_COL(baseline_ADC, double)
    EXTRACT_COL(baselineRMS_ADC, double)
    EXTRACT_COL(isSaturated, int)

    // Add eventNum if exists
    if (df.GetColumnNames().end() != std::find(df.GetColumnNames().begin(), df.GetColumnNames().end(), "EventNumber")) {
        df_final = df_final.Define("eventNum", "EventNumber");
    } else {
        df_final = df_final.Define("eventNum", [](){ static std::atomic<UInt_t> i(0); return i++; });
    }

    // 3. Book histograms and calculate average waveforms
    // Charge histograms
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hChargePtrs;
    for (int ch = 0; ch < numChannels; ++ch) {
        hChargePtrs.push_back(df_final.Define(Form("charge_ch%d", ch), [ch](const std::vector<double>& v){ return v[ch]; }, {"charge_pC"})
                                      .Histo1D({Form("hCharge_Ch%d", ch), Form("Charge Spectrum Ch %d;Charge [pC];Counts", ch), 600, -1, 5}, Form("charge_ch%d", ch)));
    }

    // Average waveforms: Aggregate RVecs element-wise
    using RVecD = ROOT::VecOps::RVec<double>;
    auto df_wave = df_final.Define("waveforms", [](const std::vector<PulseAnalysis>& res) {
        std::vector<RVecD> waves;
        for(auto& r : res) waves.push_back(RVecD(r.waveform));
        return waves;
    }, {"results"});

    auto avgWavePtr = df_wave.Aggregate(
        [numChannels, bins](std::vector<RVecD>& acc, const std::vector<RVecD>& waves){
            for(int i=0; i<numChannels; ++i) acc[i] += waves[i];
            return acc;
        },
        [numChannels](std::vector<RVecD>& acc1, const std::vector<RVecD>& acc2){
            for(int i=0; i<numChannels; ++i) acc1[i] += acc2[i];
            return acc1;
        },
        "waveforms",
        std::vector<RVecD>(numChannels, RVecD(bins, 0.0))
    );

    auto count = df_final.Count();

    // 4. Execution & Output
    TString input = fileName;
    TString baseName = gSystem->BaseName(input); 
    TString outFileName = "production_output_" + baseName;
    
    std::cout << "Running analysis loop..." << std::endl;
    df_final.Snapshot("T_Charge", outFileName, {"eventNum", "charge_pC", "amplitude_mV", "peakTime_ns", "riseTime_ns", "fallTime_ns", "fwhm_ns", "baseline_ADC", "baselineRMS_ADC", "isSaturated"});

    // After Snapshot, results are ready
    long long nEntries = *count;
    std::cout << "Analysis Complete! Processed " << nEntries << " entries." << std::endl;

    // --- Visualization & Scaling ---
    TFile *fOut = TFile::Open(outFileName, "UPDATE");
    TCanvas *c1 = new TCanvas("c1", "Charge Spectra", 1200, 800);
    int cols = (numChannels <= 1) ? 1 : (numChannels <= 2) ? 2 : (numChannels <= 4) ? 2 : (numChannels <= 6) ? 3 : 4;
    int rows = (numChannels <= 1) ? 1 : (numChannels <= 4) ? 2 : 2;
    c1->Divide(cols, rows);

    for (int ch = 0; ch < numChannels; ++ch) {
        c1->cd(ch + 1);
        gPad->SetLogy();
        hChargePtrs[ch]->SetLineColor(kBlue + (ch % 4));
        hChargePtrs[ch]->Draw();
        hChargePtrs[ch]->Write();
    }
    TString plotName = TString(baseName).ReplaceAll(".root", Form("_ch%d_charge.png", numChannels));
    c1->SaveAs(plotName);

    TCanvas *c2 = new TCanvas("c2", "Average Waveforms", 1200, 800);
    c2->Divide(cols, rows);
    auto allAvgWaves = *avgWavePtr;
    for (int ch = 0; ch < numChannels; ++ch) {
        c2->cd(ch + 1);
        auto& avgWaveVec = allAvgWaves[ch];
        TH1D *hAvg = new TH1D(Form("hAvgWave_Ch%d", ch), Form("Average Waveform Ch %d;Time [ns];Amplitude [ADC]", ch), bins, 0, bins * samplingTime_ns);
        for(int s=0; s<bins; ++s) hAvg->SetBinContent(s+1, avgWaveVec[s] / nEntries);
        hAvg->SetLineColor(kRed + (ch % 4));
        hAvg->Draw("HIST");
        hAvg->Write();
    }
    TString avgPlotName = TString(baseName).ReplaceAll(".root", Form("_ch%d_avg_wave.png", numChannels));
    c2->SaveAs(avgPlotName);

    //fOut->Close();
    std::cout << "Data saved to: " << outFileName << std::endl;
}
