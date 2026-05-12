#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

double get_tree_quantile(TTree *tree, const TString &expr, double xmin, double xmax, double prob, const char *name) {
    if (!tree || xmax <= xmin) return xmax;

    TString histName = Form("hQuant_%s", name);
    TString drawExpr = Form("%s>>%s(5000,%g,%g)", expr.Data(), histName.Data(), xmin, xmax);
    tree->Draw(drawExpr, "", "goff");
    TH1D *hist = (TH1D*)gDirectory->Get(histName);
    if (!hist || hist->GetEntries() == 0) return xmax;

    double quantile = xmax;
    double probability = prob;
    hist->GetQuantiles(1, &quantile, &probability);
    delete hist;
    return quantile;
}

void style_source_bg(TH1D *hSrc, TH1D *hBG) {
    hSrc->SetLineColor(kBlue + 1);
    hSrc->SetLineWidth(2);
    hSrc->SetFillStyle(0);

    hBG->SetLineColor(kRed + 1);
    hBG->SetLineWidth(2);
    hBG->SetLineStyle(2);
    hBG->SetFillStyle(0);
}

void configure_tree_cache(TTree *tree, int nChannels) {
    if (!tree) return;
    tree->SetCacheSize(64 * 1024 * 1024);
    if (tree->GetBranch("SyncTime_TTT")) tree->AddBranchToCache("SyncTime_TTT", kTRUE);
    for (int i = 0; i < nChannels; ++i) {
        TString branch = Form("Charge_CH%d", i);
        if (tree->GetBranch(branch)) tree->AddBranchToCache(branch, kTRUE);
    }
}

double get_live_seconds(TTree *tree, double tttClockHz) {
    if (!tree || tttClockHz <= 0) return 0.0;
    const double span = tree->GetMaximum("SyncTime_TTT") - tree->GetMinimum("SyncTime_TTT");
    return span > 0 ? span / tttClockHz : 0.0;
}

TH1D *make_rate_hist(const TH1D *hist, double liveSeconds, const char *name) {
    TH1D *rateHist = (TH1D*)hist->Clone(name);
    if (liveSeconds > 0) rateHist->Scale(1.0 / liveSeconds, "width");
    rateHist->GetYaxis()->SetTitle("Rate [Hz / NPE]");
    return rateHist;
}

double adc_integral_to_pC(double dynamicRangeV, double samplingTimeNs, double resistanceOhm, int adcBits) {
    const double maxADC = std::pow(2.0, adcBits) - 1.0;
    return (dynamicRangeV / maxADC) * (samplingTimeNs * 1.0e-9) / resistanceOhm * 1.0e12;
}

TString source_label_from_file(const char *sourceFile) {
    std::string label(sourceFile ? sourceFile : "");
    const size_t slash = label.find_last_of("/\\");
    if (slash != std::string::npos) label = label.substr(slash + 1);
    const size_t dot = label.rfind(".root");
    if (dot != std::string::npos) label = label.substr(0, dot);
    if (label.empty()) label = "Source";
    return TString(label.c_str());
}

/**
 * @brief 소스 데이터에서 백그라운드를 차감한 NPE 분포를 그리는 매크로
 * 
 * @param sourceFile  Source data ROOT file
 * @param bgFile      백그라운드 데이터 파일 (예: background_1hr_prod.root)
 * @param gain        PMT Gain 값 (기본값: 10^7)
 * @param bgScale     BG 스케일. 0보다 작으면 SyncTime_TTT span 비율로 자동 계산
 * @param xQuantile   x축 상한에 사용할 quantile. 기본값 1은 전체 범위 사용
 * @param nBins       히스토그램 bin 수
 * @param xMinUser    x축 하한. xMaxUser보다 작을 때만 수동 범위 적용
 * @param xMaxUser    x축 상한. xMinUser보다 클 때만 수동 범위 적용
 * @param tttClockHz  SyncTime_TTT clock. DT5730 TTT는 보통 125 MHz
 * @param outPrefix   출력 PNG 파일 prefix
 * @param dynamicRangeV DT5730 ADC full scale voltage range
 * @param samplingTimeNs Sampling period in ns
 * @param resistanceOhm Digitizer input impedance
 * @param adcBits      ADC bit depth
 * @param sourceLabel  Legend label for source data. Empty string derives it from sourceFile.
 */
void plot_npe_subtracted(const char* sourceFile = "source.root",
                         const char* bgFile = "background_1hr_prod.root", 
                         double gain = 1.0e7,
                         double bgScale = -1.0,
                         double xQuantile = 1.0,
                         int nBins = 400,
                         double xMinUser = 0.0,
                         double xMaxUser = -1.0,
                         double tttClockHz = 125.0e6,
                         const char* outPrefix = "npe",
                         double dynamicRangeV = 2.0,
                         double samplingTimeNs = 2.0,
                         double resistanceOhm = 50.0,
                         int adcBits = 14,
                         const char* sourceLabel = "") {
    const auto startTime = std::chrono::steady_clock::now();
    gStyle->SetOptStat(0); // 차감 후에는 통계 박스가 부정확할 수 있어 끔
    TString prefix(outPrefix);
    TString srcLabel(sourceLabel);
    if (srcLabel.Length() == 0) srcLabel = source_label_from_file(sourceFile);

    // 1. 파일 열기
    TFile *fSrc = TFile::Open(sourceFile);
    TFile *fBG = TFile::Open(bgFile);

    if (!fSrc || fSrc->IsZombie() || !fBG || fBG->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    TTree *tSrc = (TTree*)fSrc->Get("phys_tree");
    TTree *tBG = (TTree*)fBG->Get("phys_tree");

    if (!tSrc || !tBG) {
        std::cerr << "Error: 'phys_tree' not found in one of the files!" << std::endl;
        return;
    }

    const double e_charge_pC = 1.60217663e-7;
    const double adcIntegralToPC = adc_integral_to_pC(dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
    const double adcIntegralToNPE = adcIntegralToPC / (gain * e_charge_pC);
    TString npe_conv = Form("* %e", adcIntegralToNPE);
    const int nChannels = 8;
    configure_tree_cache(tSrc, nChannels);
    configure_tree_cache(tBG, nChannels);
    const double srcLiveSeconds = get_live_seconds(tSrc, tttClockHz);
    const double bgLiveSeconds = get_live_seconds(tBG, tttClockHz);
    const double srcRate = srcLiveSeconds > 0 ? tSrc->GetEntries() / srcLiveSeconds : 0.0;
    const double bgRate = bgLiveSeconds > 0 ? tBG->GetEntries() / bgLiveSeconds : 0.0;

    if (bgScale < 0) {
        if (srcLiveSeconds > 0 && bgLiveSeconds > 0) {
            bgScale = srcLiveSeconds / bgLiveSeconds;
        } else {
            bgScale = static_cast<double>(tSrc->GetEntries()) / static_cast<double>(tBG->GetEntries());
        }
    }
    if (nBins <= 0) nBins = 400;

    std::vector<int> activeChannels;
    for (int i = 0; i < nChannels; ++i) {
        TString branch = Form("Charge_CH%d", i);
        if (!tSrc->GetBranch(branch) || !tBG->GetBranch(branch)) continue;
        if (tSrc->GetMaximum(branch) != 0 || tBG->GetMaximum(branch) != 0) {
            activeChannels.push_back(i);
        }
    }

    if (activeChannels.empty()) {
        std::cerr << "Error: no active Charge_CH branches found!" << std::endl;
        return;
    }

    // --- 채널별 차감 분포 ---
    TCanvas *c_sub = new TCanvas("c_sub", Form("Background Subtracted NPE (Gain=%.1e)", gain), 1200, 600);
    c_sub->Divide(static_cast<int>(activeChannels.size()), 1);

    TCanvas *c_overlay_ch = new TCanvas("c_overlay_ch", Form("Source vs Background per Channel (Gain=%.1e)", gain), 1200, 600);
    c_overlay_ch->Divide(static_cast<int>(activeChannels.size()), 1);

    TString total_expr_raw = "(";
    std::vector<double> rawSrcMins(activeChannels.size(), 0.0);
    std::vector<double> rawSrcMaxs(activeChannels.size(), 0.0);
    std::vector<double> rawBgMins(activeChannels.size(), 0.0);
    std::vector<double> rawBgMaxs(activeChannels.size(), 0.0);

    for (size_t idx = 0; idx < activeChannels.size(); ++idx) {
        int i = activeChannels[idx];
        TString branch = Form("Charge_CH%d", i);
        TString npeExpr = Form("(%s) %s", branch.Data(), npe_conv.Data());

        rawSrcMins[idx] = tSrc->GetMinimum(branch);
        rawSrcMaxs[idx] = tSrc->GetMaximum(branch);
        rawBgMins[idx] = tBG->GetMinimum(branch);
        rawBgMaxs[idx] = tBG->GetMaximum(branch);
        const double srcMin = rawSrcMins[idx] * adcIntegralToNPE;
        const double srcMax = rawSrcMaxs[idx] * adcIntegralToNPE;
        const double bgMin = rawBgMins[idx] * adcIntegralToNPE;
        const double bgMax = rawBgMaxs[idx] * adcIntegralToNPE;
        double xmin = std::min(srcMin, bgMin);
        double xmax = std::max(srcMax, bgMax);
        if (xmin == xmax) xmax = xmin + 1.0;
        if (xQuantile > 0 && xQuantile < 1) {
            const double srcQ = get_tree_quantile(tSrc, npeExpr, xmin, xmax, xQuantile, Form("src_ch%d", i));
            const double bgQ = get_tree_quantile(tBG, npeExpr, xmin, xmax, xQuantile, Form("bg_ch%d", i));
            xmax = std::max(srcQ, bgQ);
            if (xmax <= xmin) xmax = std::max(srcMax, bgMax);
        }
        if (xMaxUser > xMinUser) {
            xmin = xMinUser;
            xmax = xMaxUser;
        }

        // 히스토그램 생성 (Source, BG, Subtracted)
        TH1D *hSrc = new TH1D(Form("hSrc_ch%d", i), "", nBins, xmin, xmax);
        TH1D *hBG = new TH1D(Form("hBG_ch%d", i), "", nBins, xmin, xmax);
        
        tSrc->Project(hSrc->GetName(), npeExpr);
        tBG->Project(hBG->GetName(), npeExpr);
        hBG->Scale(bgScale);
        style_source_bg(hSrc, hBG);

        // 차감용 히스토그램 (Source 복사 후 BG 차감)
        TH1D *hSub = (TH1D*)hSrc->Clone(Form("hSub_ch%d", i));
        hSub->SetTitle(Form("CH%d: Source - %.3g #times BG;NPE;Counts", i, bgScale));
        hSub->Add(hBG, -1.0);

        c_overlay_ch->cd(static_cast<int>(idx) + 1);
        gPad->SetLogy();
        hSrc->SetTitle(Form("CH%d: Source vs scaled BG;NPE;Counts", i));
        hSrc->SetMinimum(0.5);
        hSrc->SetMaximum(std::max(hSrc->GetMaximum(), hBG->GetMaximum()) * 2.0);
        hSrc->Draw("HIST");
        hBG->Draw("SAME HIST");
        TLegend *legCh = new TLegend(0.55, 0.72, 0.88, 0.88);
        legCh->AddEntry(hSrc, srcLabel.Data(), "l");
        legCh->AddEntry(hBG, Form("BG #times %.3g", bgScale), "l");
        legCh->Draw();

        const double yMin = std::min(0.0, hSub->GetMinimum());
        const double yMax = std::max(0.0, hSub->GetMaximum());
        const double yPad = std::max(1.0, (yMax - yMin) * 0.15);

        c_sub->cd(static_cast<int>(idx) + 1);
        hSub->SetLineColor(kBlack);
        hSub->SetLineWidth(2);
        hSub->SetFillColor(kGray);
        hSub->SetMinimum(yMin - yPad);
        hSub->SetMaximum(yMax + yPad);
        hSub->Draw("HIST");
        TLine *zero = new TLine(xmin, 0, xmax, 0);
        zero->SetLineStyle(2);
        zero->Draw("SAME");

        if (idx > 0) total_expr_raw += " + ";
        total_expr_raw += branch;
    }
    total_expr_raw += ")";
    TString total_npe_expr = Form("%s %s", total_expr_raw.Data(), npe_conv.Data());

    c_overlay_ch->SaveAs(Form("%s_overlay_channels_log.png", prefix.Data()));
    c_sub->SaveAs(Form("%s_subtracted_channels.png", prefix.Data()));

    // --- 전체 총합 차감 분포 ---
    TCanvas *c_total_sub = new TCanvas("c_total_sub", "Total NPE Subtraction", 800, 600);
    
    double totalRawSrcMin = 0.0;
    double totalRawSrcMax = 0.0;
    double totalRawBgMin = 0.0;
    double totalRawBgMax = 0.0;
    for (size_t idx = 0; idx < activeChannels.size(); ++idx) {
        totalRawSrcMin += rawSrcMins[idx];
        totalRawSrcMax += rawSrcMaxs[idx];
        totalRawBgMin += rawBgMins[idx];
        totalRawBgMax += rawBgMaxs[idx];
    }
    double t_xmin = std::min(totalRawSrcMin, totalRawBgMin) * adcIntegralToNPE;
    double t_xmax = std::max(totalRawSrcMax, totalRawBgMax) * adcIntegralToNPE;
    if (t_xmin == t_xmax) t_xmax = t_xmin + 1.0;
    if (xQuantile > 0 && xQuantile < 1) {
        const double srcQ = get_tree_quantile(tSrc, total_npe_expr, t_xmin, t_xmax, xQuantile, "total_src");
        const double bgQ = get_tree_quantile(tBG, total_npe_expr, t_xmin, t_xmax, xQuantile, "total_bg");
        t_xmax = std::max(srcQ, bgQ);
        if (t_xmax <= t_xmin) t_xmax = std::max(totalRawSrcMax, totalRawBgMax) * adcIntegralToNPE;
    }
    if (xMaxUser > xMinUser) {
        t_xmin = xMinUser;
        t_xmax = xMaxUser;
    }

    TH1D *hTotalSrc = new TH1D("hTotalSrc", "Total NPE;NPE;Counts", nBins, t_xmin, t_xmax);
    TH1D *hTotalBG = new TH1D("hTotalBG", "", nBins, t_xmin, t_xmax);
    
    tSrc->Project("hTotalSrc", total_npe_expr);
    tBG->Project("hTotalBG", total_npe_expr);
    hTotalBG->Scale(bgScale);
    style_source_bg(hTotalSrc, hTotalBG);

    TH1D *hTotalSub = (TH1D*)hTotalSrc->Clone("hTotalSub");
    hTotalSub->SetTitle(Form("Total NPE (Source - %.3g #times BG, Gain=%.1e);NPE;Counts", bgScale, gain));
    hTotalSub->Add(hTotalBG, -1.0);

    TCanvas *c_total_overlay = new TCanvas("c_total_overlay", "Total NPE Source vs Background", 900, 650);
    c_total_overlay->cd();
    hTotalSrc->SetTitle(Form("Total NPE: Source vs scaled BG (Gain=%.1e);NPE;Counts", gain));
    hTotalSrc->SetMinimum(0.0);
    hTotalSrc->SetMaximum(std::max(hTotalSrc->GetMaximum(), hTotalBG->GetMaximum()) * 1.15);
    hTotalSrc->Draw("HIST");
    hTotalBG->Draw("SAME HIST");
    TLegend *legOverlay = new TLegend(0.58, 0.72, 0.88, 0.88);
    legOverlay->AddEntry(hTotalSrc, srcLabel.Data(), "l");
    legOverlay->AddEntry(hTotalBG, Form("BG #times %.3g", bgScale), "l");
    legOverlay->Draw();
    c_total_overlay->SaveAs(Form("%s_overlay_total_linear.png", prefix.Data()));

    TCanvas *c_total_overlay_log = new TCanvas("c_total_overlay_log", "Total NPE Source vs Background Log", 900, 650);
    c_total_overlay_log->cd();
    gPad->SetLogy();
    hTotalSrc->SetMinimum(0.5);
    hTotalSrc->SetMaximum(std::max(hTotalSrc->GetMaximum(), hTotalBG->GetMaximum()) * 2.0);
    hTotalSrc->Draw("HIST");
    hTotalBG->Draw("SAME HIST");
    TLegend *legOverlayLog = new TLegend(0.58, 0.72, 0.88, 0.88);
    legOverlayLog->AddEntry(hTotalSrc, srcLabel.Data(), "l");
    legOverlayLog->AddEntry(hTotalBG, Form("BG #times %.3g", bgScale), "l");
    legOverlayLog->Draw();
    c_total_overlay_log->SaveAs(Form("%s_overlay_total_log.png", prefix.Data()));

    TH1D *hTotalSrcRate = make_rate_hist(hTotalSrc, srcLiveSeconds, "hTotalSrcRate");
    TH1D *hTotalBGRate = make_rate_hist(hTotalBG, srcLiveSeconds, "hTotalBGRate");
    style_source_bg(hTotalSrcRate, hTotalBGRate);

    TCanvas *c_total_rate_log = new TCanvas("c_total_rate_log", "Total NPE Rate-Normalized Source vs Background", 900, 650);
    c_total_rate_log->cd();
    gPad->SetLogy();
    hTotalSrcRate->SetTitle(Form("Total NPE rate density: Source vs scaled BG (Gain=%.1e);NPE;Rate [Hz / NPE]", gain));
    hTotalSrcRate->SetMinimum(1e-5);
    hTotalSrcRate->SetMaximum(std::max(hTotalSrcRate->GetMaximum(), hTotalBGRate->GetMaximum()) * 2.0);
    hTotalSrcRate->Draw("HIST");
    hTotalBGRate->Draw("SAME HIST");
    TLegend *legRate = new TLegend(0.54, 0.72, 0.88, 0.88);
    legRate->AddEntry(hTotalSrcRate, Form("%s %.1f Hz", srcLabel.Data(), srcRate), "l");
    legRate->AddEntry(hTotalBGRate, Form("BG scaled to source live time %.1f Hz", bgRate), "l");
    legRate->Draw();
    c_total_rate_log->SaveAs(Form("%s_overlay_total_rate_log.png", prefix.Data()));

    c_total_sub->cd();
    const double totalYMin = std::min(0.0, hTotalSub->GetMinimum());
    const double totalYMax = std::max(0.0, hTotalSub->GetMaximum());
    const double totalYPad = std::max(1.0, (totalYMax - totalYMin) * 0.15);

    hTotalSub->SetLineColor(kBlack);
    hTotalSub->SetLineWidth(2);
    hTotalSub->SetFillColor(kGray);
    hTotalSub->SetMinimum(totalYMin - totalYPad);
    hTotalSub->SetMaximum(totalYMax + totalYPad);
    hTotalSub->Draw("HIST");
    TLine *zeroTotal = new TLine(t_xmin, 0, t_xmax, 0);
    zeroTotal->SetLineStyle(2);
    zeroTotal->Draw("SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(hTotalSub, "Source - scaled BG", "l");
    leg->Draw();

    c_total_sub->SaveAs(Form("%s_subtracted_total.png", prefix.Data()));

    std::cout << "\n--- Background Subtraction Complete ---" << std::endl;
    std::cout << "Source: " << sourceFile << std::endl;
    std::cout << "Source label: " << srcLabel << std::endl;
    std::cout << "BG: " << bgFile << std::endl;
    std::cout << "Charge_CH unit: ADC-count sample integral from DAQ production_dt5730.cpp" << std::endl;
    std::cout << "ADC integral to pC: " << adcIntegralToPC << " pC / (ADC count * sample)" << std::endl;
    std::cout << "Conversion: NPE = Charge_CH " << npe_conv << std::endl;
    std::cout << "Conversion constants: dynamicRange=" << dynamicRangeV << " V, sampling="
              << samplingTimeNs << " ns, R=" << resistanceOhm << " ohm, ADC bits=" << adcBits
              << ", gain=" << gain << std::endl;
    std::cout << "Entries: source=" << tSrc->GetEntries() << ", BG=" << tBG->GetEntries() << std::endl;
    std::cout << "Live time [s]: source=" << srcLiveSeconds << ", BG=" << bgLiveSeconds << std::endl;
    std::cout << "Trigger rate [Hz]: source=" << srcRate << ", BG=" << bgRate
              << ", excess=" << (srcRate - bgRate) << std::endl;
    std::cout << "Active channels:";
    for (int ch : activeChannels) std::cout << " CH" << ch;
    std::cout << std::endl;
    std::cout << "BG scale: " << bgScale << std::endl;
    std::cout << "X-axis quantile: " << xQuantile << std::endl;
    std::cout << "Bins: " << nBins << ", X range: [" << t_xmin << ", " << t_xmax << "]" << std::endl;
    std::cout << "Output prefix: " << prefix << std::endl;
    const auto endTime = std::chrono::steady_clock::now();
    const double elapsedSeconds = std::chrono::duration<double>(endTime - startTime).count();
    std::cout << "Elapsed wall time [s]: " << elapsedSeconds << std::endl;
    std::cout << "Results saved to: " << prefix << "_subtracted_channels.png, "
              << prefix << "_subtracted_total.png" << std::endl;
    std::cout << "Overlay plots saved to: " << prefix << "_overlay_channels_log.png, "
              << prefix << "_overlay_total_linear.png, "
              << prefix << "_overlay_total_log.png, "
              << prefix << "_overlay_total_rate_log.png" << std::endl;
}
