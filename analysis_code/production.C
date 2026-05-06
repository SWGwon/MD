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

struct EventResult {
    ROOT::VecOps::RVec<double> charge;
    ROOT::VecOps::RVec<double> amp;
    ROOT::VecOps::RVec<double> pTime;
    ROOT::VecOps::RVec<double> rTime;
    ROOT::VecOps::RVec<double> fTime;
    ROOT::VecOps::RVec<double> fwhm;
    ROOT::VecOps::RVec<double> bl;
    ROOT::VecOps::RVec<double> blRMS;
    ROOT::VecOps::RVec<int> sat;
    std::vector<ROOT::VecOps::RVec<double>> waves;
    UInt_t evNum;
};

void production(std::string fileName, int numChannels = -1) {
    ROOT::EnableImplicitMT();
    double dynamicRange = 2.0, resistance = 50.0, samplingTime = 2.0e-9;
    int numBits = 14; UInt_t maxADC = (1 << numBits) - 1;

    TFile *fCheck = TFile::Open(fileName.c_str());
    TTree *TIn = (TTree*)fCheck->Get("T");
    TLeaf *lADC = TIn->GetLeaf("ADC");
    TString title = lADC->GetTitle(); 
    int fCh = 1, bins = 1;
    if (title.Contains("[") && title.Contains("]")) sscanf(title.Data(), "ADC[%d][%d]", &fCh, &bins);
    else { fCh = 1; bins = lADC->GetLen(); }
    if (numChannels <= 0 || numChannels > fCh) numChannels = fCh;
    Long64_t nEntries = TIn->GetEntries();
    bool hasEv = (TIn->GetBranch("EventNumber") != nullptr);
    fCheck->Close();

    double factor = (dynamicRange/maxADC * samplingTime / resistance) * 1e12;
    double toMV = (dynamicRange/maxADC) * 1000.0;
    double ts_ns = samplingTime * 1e9;

    ROOT::RDataFrame df("T", fileName);
    auto analyze = [=](const ROOT::VecOps::RVec<UInt_t>& adc, UInt_t ev) {
        EventResult res; res.evNum = ev;
        res.charge.resize(numChannels); res.amp.resize(numChannels); res.pTime.resize(numChannels);
        res.rTime.resize(numChannels); res.fTime.resize(numChannels); res.fwhm.resize(numChannels);
        res.bl.resize(numChannels); res.blRMS.resize(numChannels); res.sat.resize(numChannels);
        res.waves.resize(numChannels);
        for (int ch = 0; ch < numChannels; ++ch) {
            int off = ch * bins; res.waves[ch].resize(bins);
            std::vector<UInt_t> s(bins);
            for(int i=0; i<bins; ++i) {
                s[i] = adc[off + i];
                if (s[i] >= maxADC - 1) res.sat[ch] = 1;
            }
            std::vector<UInt_t> srt = s;
            std::nth_element(srt.begin(), srt.begin() + bins/2, srt.end());
            double bl = (double)srt[bins/2]; res.bl[ch] = bl;
            double sq = 0, mA = -1e9; int pI = -1;
            for (int i = 0; i < bins; ++i) {
                double v = (double)s[i], pV = bl - v;
                sq += std::pow(v - bl, 2); res.waves[ch][i] = pV;
                if (pV > mA) { mA = pV; pI = i; }
            }
            res.blRMS[ch] = std::sqrt(sq / bins); res.amp[ch] = mA * toMV; res.pTime[ch] = pI * ts_ns;
            double t10=-1, t50_1=-1, t90=-1, t90_f=-1, t50_2=-1, t10_f=-1;
            for (int i = 0; i < bins; ++i) {
                double v = res.waves[ch][i], t = i * ts_ns;
                if (t10 < 0 && v >= 0.1 * mA) t10 = t;
                if (t50_1 < 0 && v >= 0.5 * mA) t50_1 = t;
                if (t90 < 0 && v >= 0.9 * mA) t90 = t;
                if (i > pI) {
                    if (t90_f < 0 && v <= 0.9 * mA) t90_f = t;
                    if (t50_2 < 0 && v <= 0.5 * mA) t50_2 = t;
                    if (t10_f < 0 && v <= 0.1 * mA) t10_f = t;
                }
            }
            res.rTime[ch] = (t90 > 0 && t10 > 0) ? (t90 - t10) : 0;
            res.fTime[ch] = (t10_f > 0 && t90_f > 0) ? (t10_f - t90_f) : 0;
            res.fwhm[ch]  = (t50_2 > 0 && t50_1 > 0) ? (t50_2 - t50_1) : 0;
            res.charge[ch] = std::accumulate(res.waves[ch].begin(), res.waves[ch].end(), 0.0) * factor;
        }
        return res;
    };

    auto df_1 = hasEv ? df : df.Define("EventNumber", [](){ static std::atomic<UInt_t> i(0); return i++; });
    auto df_res = df_1.Define("r", analyze, {"ADC", "EventNumber"});

    using RVecD = ROOT::VecOps::RVec<double>;
    auto avgWPtr = df_res.Aggregate(
        [numChannels, bins](std::vector<RVecD>& acc, const EventResult& r){
            for(int i=0; i<numChannels; ++i) acc[i] += r.waves[i]; return acc;
        },
        [numChannels](std::vector<RVecD>& a1, const std::vector<RVecD>& a2){
            for(int i=0; i<numChannels; ++i) a1[i] += a2[i]; return a1;
        }, "r", std::vector<RVecD>(numChannels, RVecD(bins, 0.0)));

    std::vector<ROOT::RDF::RResultPtr<TH1D>> hPtrs;
    for (int ch = 0; ch < numChannels; ++ch) {
        hPtrs.push_back(df_res.Define(Form("c%d", ch), [ch](const EventResult& r){ return r.charge[ch]; }, {"r"})
                             .Histo1D({Form("hCharge_Ch%d", ch), "Charge;pC", 600, -1, 5}, Form("c%d", ch)));
    }

    // Extract columns for TTree (memory efficient Take)
    auto df_cols = df_res.Define("charge_v", [](const EventResult& r){ return r.charge; }, {"r"})
                         .Define("amp_v",    [](const EventResult& r){ return r.amp; }, {"r"})
                         .Define("pTime_v",  [](const EventResult& r){ return r.pTime; }, {"r"})
                         .Define("rTime_v",  [](const EventResult& r){ return r.rTime; }, {"r"})
                         .Define("fTime_v",  [](const EventResult& r){ return r.fTime; }, {"r"})
                         .Define("fwhm_v",   [](const EventResult& r){ return r.fwhm; }, {"r"})
                         .Define("bl_v",     [](const EventResult& r){ return r.bl; }, {"r"})
                         .Define("blRMS_v",  [](const EventResult& r){ return r.blRMS; }, {"r"})
                         .Define("sat_v",    [](const EventResult& r){ return r.sat; }, {"r"})
                         .Define("ev_v",     [](const EventResult& r){ return r.evNum; }, {"r"});

    auto cList = df_cols.Take<RVecD>("charge_v");
    auto aList = df_cols.Take<RVecD>("amp_v");
    auto pList = df_cols.Take<RVecD>("pTime_v");
    auto rList = df_cols.Take<RVecD>("rTime_v");
    auto fList = df_cols.Take<RVecD>("fTime_v");
    auto wList = df_cols.Take<RVecD>("fwhm_v");
    auto bList = df_cols.Take<RVecD>("bl_v");
    auto sList = df_cols.Take<RVecD>("blRMS_v");
    auto iList = df_cols.Take<ROOT::VecOps::RVec<int>>("sat_v");
    auto eList = df_cols.Take<UInt_t>("ev_v");

    std::cout << "Calculating and filling tree..." << std::endl;
    TString outName = "production_output_" + TString(gSystem->BaseName(fileName.c_str()));
    TFile *fOut = new TFile(outName, "RECREATE");
    TTree *TOut = new TTree("T_Charge", "Analysis");
    
    std::vector<double> charge_pC(numChannels), amplitude_mV(numChannels), peakTime_ns(numChannels), riseTime_ns(numChannels), fallTime_ns(numChannels), fwhm_ns(numChannels), baseline_ADC(numChannels), baselineRMS_ADC(numChannels);
    std::vector<int> isSaturated(numChannels); UInt_t eventNum;
    TOut->Branch("eventNum", &eventNum, "eventNum/i");
    TOut->Branch("charge_pC", charge_pC.data(), Form("charge_pC[%d]/D", numChannels));
    TOut->Branch("amplitude_mV", amplitude_mV.data(), Form("amplitude_mV[%d]/D", numChannels));
    TOut->Branch("peakTime_ns", peakTime_ns.data(), Form("peakTime_ns[%d]/D", numChannels));
    TOut->Branch("riseTime_ns", riseTime_ns.data(), Form("riseTime_ns[%d]/D", numChannels));
    TOut->Branch("fallTime_ns", fallTime_ns.data(), Form("fallTime_ns[%d]/D", numChannels));
    TOut->Branch("fwhm_ns", fwhm_ns.data(), Form("fwhm_ns[%d]/D", numChannels));
    TOut->Branch("baseline_ADC", baseline_ADC.data(), Form("baseline_ADC[%d]/D", numChannels));
    TOut->Branch("baselineRMS_ADC", baselineRMS_ADC.data(), Form("baselineRMS_ADC[%d]/D", numChannels));
    TOut->Branch("isSaturated", isSaturated.data(), Form("isSaturated[%d]/I", numChannels));

    auto& cv = *cList; auto& av = *aList; auto& pv = *pList; auto& rv = *rList; auto& fv = *fList; auto& wv = *wList; auto& bv = *bList; auto& sv = *sList; auto& iv = *iList; auto& ev = *eList;
    for (size_t i = 0; i < cv.size(); ++i) {
        eventNum = ev[i];
        for (int ch = 0; ch < numChannels; ++ch) {
            charge_pC[ch] = cv[i][ch]; amplitude_mV[ch] = av[i][ch]; peakTime_ns[ch] = pv[i][ch];
            riseTime_ns[ch] = rv[i][ch]; fallTime_ns[ch] = fv[i][ch]; fwhm_ns[ch] = wv[i][ch];
            baseline_ADC[ch] = bv[i][ch]; baselineRMS_ADC[ch] = sv[i][ch]; isSaturated[ch] = iv[i][ch];
        }
        TOut->Fill();
    }
    TOut->Write();

    auto allAvg = *avgWPtr;
    for (int ch = 0; ch < numChannels; ++ch) {
        hPtrs[ch]->Write();
        TH1D hAvg(Form("hAvgWave_Ch%d", ch), "Avg Wave;ns;ADC", bins, 0, bins*ts_ns);
        for(int s=0; s<bins; ++s) hAvg.SetBinContent(s+1, allAvg[ch][s]/nEntries);
        hAvg.SetEntries(nEntries * bins); hAvg.Write();
    }

    // Save Plots (Matching original naming)
    TCanvas *c1 = new TCanvas("c1", "Charge Spectra", 1200, 800);
    int cols = (numChannels <= 4) ? 2 : (numChannels <= 6) ? 3 : 4;
    int rows = (numChannels <= 1) ? 1 : 2;
    c1->Divide(cols, rows);
    for (int ch = 0; ch < numChannels; ++ch) {
        c1->cd(ch+1); gPad->SetLogy();
        hPtrs[ch]->SetLineColor(kBlue + (ch%4));
        hPtrs[ch]->Draw();
    }
    TString base = gSystem->BaseName(fileName.c_str());
    c1->SaveAs(base.ReplaceAll(".root", Form("_ch%d_charge.png", numChannels)));

    TCanvas *c2 = new TCanvas("c2", "Average Waveforms", 1200, 800);
    c2->Divide(cols, rows);
    for (int ch = 0; ch < numChannels; ++ch) {
        c2->cd(ch+1);
        TH1D *hAvg = (TH1D*)fOut->Get(Form("hAvgWave_Ch%d", ch));
        if(hAvg) {
            hAvg->SetLineColor(kRed + (ch%4));
            hAvg->Draw("HIST");
        }
    }
    base = gSystem->BaseName(fileName.c_str());
    c2->SaveAs(base.ReplaceAll(".root", Form("_ch%d_avg_wave.png", numChannels)));

    fOut->Close();
    std::cout << "Done." << std::endl;
}
