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
#include <cstdlib>
#include <string>
#include <chrono>
#include <sstream>

struct ChannelThreshold {
    int channel;
    double minNpe;
};

double get_tree_quantile(TTree *tree, const TString &expr, double xmin, double xmax, double prob, const char *name, const TString &selection = "") {
    if (!tree || xmax <= xmin) return xmax;

    TString histName = Form("hQuant_%s", name);
    TString drawExpr = Form("%s>>%s(5000,%g,%g)", expr.Data(), histName.Data(), xmin, xmax);
    tree->Draw(drawExpr, selection, "goff");
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

std::vector<int> parse_channel_list(const char *channels, int nChannels, bool *ok = nullptr) {
    if (ok) *ok = true;
    std::vector<int> selected;
    TString text(channels ? channels : "");
    text.ReplaceAll(" ", "");
    if (text.Length() == 0) return selected;

    std::stringstream stream(text.Data());
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) continue;
        if (token.rfind("CH", 0) == 0 || token.rfind("ch", 0) == 0) token = token.substr(2);

        char *end = nullptr;
        const long channel = std::strtol(token.c_str(), &end, 10);
        if (!end || *end != '\0' || channel < 0 || channel >= nChannels) {
            std::cerr << "Error: invalid channel selection '" << token << "'" << std::endl;
            if (ok) *ok = false;
            continue;
        }
        if (std::find(selected.begin(), selected.end(), static_cast<int>(channel)) == selected.end()) {
            selected.push_back(static_cast<int>(channel));
        }
    }
    return selected;
}

std::vector<ChannelThreshold> parse_threshold_list(const char *thresholds, int nChannels, bool *ok = nullptr) {
    if (ok) *ok = true;
    std::vector<ChannelThreshold> parsed;
    TString text(thresholds ? thresholds : "");
    text.ReplaceAll(" ", "");
    if (text.Length() == 0) return parsed;

    std::stringstream stream(text.Data());
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) continue;
        const size_t colon = token.find(':');
        if (colon == std::string::npos) {
            std::cerr << "Error: invalid threshold token '" << token << "'" << std::endl;
            if (ok) *ok = false;
            continue;
        }

        std::string channelToken = token.substr(0, colon);
        std::string thresholdToken = token.substr(colon + 1);
        if (channelToken.rfind("CH", 0) == 0 || channelToken.rfind("ch", 0) == 0) channelToken = channelToken.substr(2);

        char *channelEnd = nullptr;
        const long channel = std::strtol(channelToken.c_str(), &channelEnd, 10);
        char *thresholdEnd = nullptr;
        const double minNpe = std::strtod(thresholdToken.c_str(), &thresholdEnd);
        if (!channelEnd || *channelEnd != '\0' || channel < 0 || channel >= nChannels ||
            !thresholdEnd || *thresholdEnd != '\0' || !std::isfinite(minNpe)) {
            std::cerr << "Error: invalid threshold token '" << token << "'" << std::endl;
            if (ok) *ok = false;
            continue;
        }

        auto existing = std::find_if(parsed.begin(), parsed.end(), [channel](const ChannelThreshold &threshold) {
            return threshold.channel == channel;
        });
        if (existing == parsed.end()) {
            parsed.push_back({static_cast<int>(channel), minNpe});
        } else {
            existing->minNpe = minNpe;
        }
    }
    return parsed;
}

TString build_threshold_selection(const std::vector<ChannelThreshold> &thresholds, double adcIntegralToNPE) {
    TString selection;
    for (size_t idx = 0; idx < thresholds.size(); ++idx) {
        if (idx > 0) selection += " && ";
        const ChannelThreshold &threshold = thresholds[idx];
        selection += Form("(Charge_CH%d * %.12g >= %.12g)", threshold.channel, adcIntegralToNPE, threshold.minNpe);
    }
    return selection;
}

/**
 * @brief 소스 데이터에서 백그라운드를 차감한 NPE 분포를 그리는 매크로
 * 
 * @param sourceFile  Source data ROOT file
 * @param bgFile      백그라운드 데이터 파일. Empty string runs source-only mode.
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
 * @param selectedChannels Comma-separated channel list such as "0,1,3". Empty string auto-detects active channels.
 * @param channelThresholds Comma-separated NPE cuts such as "0:5,1:5". All cuts must pass for an event.
 */
void plot_compton_edge(const char* sourceFile = "source.root",
                         const char* bgFile = "background_1hr_prod.root", 
                         double gain = 1.0e7,
                         double bgScale = -1.0,
                         double xQuantile = 1.0,
                         int nBins = 400,
                         double xMinUser = 0.0,
                         double xMaxUser = -1.0,
                         double tttClockHz = 125.0e6,
                         const char* outPrefix = "compton_edge",
                         double dynamicRangeV = 2.0,
                         double samplingTimeNs = 2.0,
                         double resistanceOhm = 50.0,
                         int adcBits = 14,
                         const char* sourceLabel = "",
                         const char* selectedChannels = "",
                         const char* channelThresholds = "") {
    const auto startTime = std::chrono::steady_clock::now();
    gStyle->SetOptStat(0); // 차감 후에는 통계 박스가 부정확할 수 있어 끔
    TString prefix(outPrefix);
    TString srcLabel(sourceLabel);
    if (srcLabel.Length() == 0) srcLabel = source_label_from_file(sourceFile);
    TString bgPath(bgFile ? bgFile : "");
    const bool hasBackground = bgPath.Length() > 0;

    // 1. 파일 열기
    TFile *fSrc = TFile::Open(sourceFile);
    TFile *fBG = hasBackground ? TFile::Open(bgFile) : nullptr;

    if (!fSrc || fSrc->IsZombie() || (hasBackground && (!fBG || fBG->IsZombie()))) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    TTree *tSrc = (TTree*)fSrc->Get("phys_tree");
    TTree *tBG = hasBackground ? (TTree*)fBG->Get("phys_tree") : nullptr;

    if (!tSrc || (hasBackground && !tBG)) {
        std::cerr << "Error: 'phys_tree' not found in one of the files!" << std::endl;
        return;
    }

    const double e_charge_pC = 1.60217663e-7;
    const double adcIntegralToPC = adc_integral_to_pC(dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
    const double adcIntegralToNPE = adcIntegralToPC / (gain * e_charge_pC);
    TString npe_conv = Form("* %e", adcIntegralToNPE);
    const int nChannels = 8;
    configure_tree_cache(tSrc, nChannels);
    if (hasBackground) configure_tree_cache(tBG, nChannels);
    const double srcLiveSeconds = get_live_seconds(tSrc, tttClockHz);
    const double bgLiveSeconds = hasBackground ? get_live_seconds(tBG, tttClockHz) : 0.0;
    const double srcRate = srcLiveSeconds > 0 ? tSrc->GetEntries() / srcLiveSeconds : 0.0;
    const double bgRate = bgLiveSeconds > 0 ? tBG->GetEntries() / bgLiveSeconds : 0.0;

    if (!hasBackground) {
        bgScale = 0.0;
    } else if (bgScale < 0) {
        if (srcLiveSeconds > 0 && bgLiveSeconds > 0) {
            bgScale = srcLiveSeconds / bgLiveSeconds;
        } else {
            std::cerr << "Error: automatic background scaling requires valid SyncTime_TTT live time in both files." << std::endl;
            std::cerr << "       Pass an explicit bgScale to plot_compton_edge if live-time scaling is not applicable." << std::endl;
            return;
        }
    }
    if (nBins <= 0) nBins = 400;

    TString thresholdText(channelThresholds ? channelThresholds : "");
    thresholdText.ReplaceAll(" ", "");
    bool thresholdsOk = true;
    const std::vector<ChannelThreshold> thresholds = parse_threshold_list(channelThresholds, nChannels, &thresholdsOk);
    if (!thresholdsOk || (thresholdText.Length() > 0 && thresholds.empty())) {
        std::cerr << "Error: invalid threshold selection. Expected comma-separated CH:MIN pairs, e.g. 0:5,1:5." << std::endl;
        return;
    }
    for (const ChannelThreshold &threshold : thresholds) {
        TString branch = Form("Charge_CH%d", threshold.channel);
        if (!tSrc->GetBranch(branch) || (hasBackground && !tBG->GetBranch(branch))) {
            std::cerr << "Error: threshold branch " << branch << " not found in required input files." << std::endl;
            return;
        }
    }
    const TString eventSelection = build_threshold_selection(thresholds, adcIntegralToNPE);

    TString selectedText(selectedChannels ? selectedChannels : "");
    selectedText.ReplaceAll(" ", "");
    bool channelsOk = true;
    std::vector<int> requestedChannels = parse_channel_list(selectedChannels, nChannels, &channelsOk);
    if (!channelsOk || (selectedText.Length() > 0 && requestedChannels.empty())) {
        std::cerr << "Error: invalid channel selection. Expected comma-separated channels in [0, "
                  << (nChannels - 1) << "], e.g. 0,1,3." << std::endl;
        return;
    }
    std::vector<int> activeChannels;
    const bool useRequestedChannels = selectedText.Length() > 0;
    const std::vector<int> channelsToCheck = useRequestedChannels ? requestedChannels : [&]() {
        std::vector<int> channels;
        for (int i = 0; i < nChannels; ++i) channels.push_back(i);
        return channels;
    }();

    for (int i : channelsToCheck) {
        TString branch = Form("Charge_CH%d", i);
        if (!tSrc->GetBranch(branch) || (hasBackground && !tBG->GetBranch(branch))) {
            if (useRequestedChannels) {
                std::cerr << "Warning: selected branch " << branch << " not found in required input files; skipping." << std::endl;
            }
            continue;
        }
        if (useRequestedChannels || tSrc->GetMaximum(branch) != 0 || (hasBackground && tBG->GetMaximum(branch) != 0)) {
            activeChannels.push_back(i);
        }
    }

    if (activeChannels.empty()) {
        if (useRequestedChannels) {
            std::cerr << "Error: none of the selected Charge_CH branches can be analyzed!" << std::endl;
        } else {
            std::cerr << "Error: no active Charge_CH branches found!" << std::endl;
        }
        return;
    }

    std::vector<TH1D*> histogramsToWrite;

    // --- 채널별 분포 ---
    TCanvas *c_sub = new TCanvas("c_sub", Form("%s NPE (Gain=%.1e)", hasBackground ? "Background Subtracted" : "Source Only", gain), 1200, 600);
    c_sub->Divide(static_cast<int>(activeChannels.size()), 1);

    TCanvas *c_overlay_ch = new TCanvas("c_overlay_ch", Form("%s per Channel (Gain=%.1e)", hasBackground ? "Source vs Background" : "Source", gain), 1200, 600);
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
        rawBgMins[idx] = hasBackground ? tBG->GetMinimum(branch) : rawSrcMins[idx];
        rawBgMaxs[idx] = hasBackground ? tBG->GetMaximum(branch) : rawSrcMaxs[idx];
        const double srcMin = rawSrcMins[idx] * adcIntegralToNPE;
        const double srcMax = rawSrcMaxs[idx] * adcIntegralToNPE;
        const double bgMin = rawBgMins[idx] * adcIntegralToNPE;
        const double bgMax = rawBgMaxs[idx] * adcIntegralToNPE;
        double xmin = std::min(srcMin, bgMin);
        double xmax = std::max(srcMax, bgMax);
        if (xmin == xmax) xmax = xmin + 1.0;
        if (xQuantile > 0 && xQuantile < 1) {
            const double srcQ = get_tree_quantile(tSrc, npeExpr, xmin, xmax, xQuantile, Form("src_ch%d", i), eventSelection);
            const double bgQ = hasBackground ? get_tree_quantile(tBG, npeExpr, xmin, xmax, xQuantile, Form("bg_ch%d", i), eventSelection) : srcQ;
            xmax = std::max(srcQ, bgQ);
            if (xmax <= xmin) xmax = std::max(srcMax, bgMax);
        }
        if (xMaxUser > xMinUser) {
            xmin = xMinUser;
            xmax = xMaxUser;
        }

        // 히스토그램 생성 (Source, optional BG, display)
        TH1D *hSrc = new TH1D(Form("hSrc_ch%d", i), "", nBins, xmin, xmax);
        TH1D *hBG = hasBackground ? new TH1D(Form("hBG_ch%d", i), "", nBins, xmin, xmax) : nullptr;
        
        tSrc->Project(hSrc->GetName(), npeExpr, eventSelection);
        if (hasBackground) {
            tBG->Project(hBG->GetName(), npeExpr, eventSelection);
            hBG->Scale(bgScale);
            style_source_bg(hSrc, hBG);
        } else {
            hSrc->SetLineColor(kBlue + 1);
            hSrc->SetLineWidth(2);
            hSrc->SetFillStyle(0);
        }

        TH1D *hSub = (TH1D*)hSrc->Clone(hasBackground ? Form("hSub_ch%d", i) : Form("hSourceOnly_ch%d", i));
        hSub->SetTitle(hasBackground ? Form("CH%d: Source - %.3g #times BG;NPE;Counts", i, bgScale)
                                     : Form("CH%d: Source only;NPE;Counts", i));
        if (hasBackground) hSub->Add(hBG, -1.0);
        histogramsToWrite.push_back(hSrc);
        if (hBG) histogramsToWrite.push_back(hBG);
        histogramsToWrite.push_back(hSub);

        c_overlay_ch->cd(static_cast<int>(idx) + 1);
        gPad->SetLogy();
        hSrc->SetTitle(hasBackground ? Form("CH%d: Source vs scaled BG;NPE;Counts", i)
                                     : Form("CH%d: Source;NPE;Counts", i));
        hSrc->SetMinimum(0.5);
        hSrc->SetMaximum((hasBackground ? std::max(hSrc->GetMaximum(), hBG->GetMaximum()) : hSrc->GetMaximum()) * 2.0);
        hSrc->Draw("HIST");
        if (hasBackground) hBG->Draw("SAME HIST");
        TLegend *legCh = new TLegend(0.55, 0.72, 0.88, 0.88);
        legCh->AddEntry(hSrc, srcLabel.Data(), "l");
        if (hasBackground) legCh->AddEntry(hBG, Form("BG #times %.3g", bgScale), "l");
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

    // --- 전체 총합 분포 ---
    TCanvas *c_total_sub = new TCanvas("c_total_sub", hasBackground ? "Total NPE Subtraction" : "Total NPE Source Only", 800, 600);
    
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
        const double srcQ = get_tree_quantile(tSrc, total_npe_expr, t_xmin, t_xmax, xQuantile, "total_src", eventSelection);
        const double bgQ = hasBackground ? get_tree_quantile(tBG, total_npe_expr, t_xmin, t_xmax, xQuantile, "total_bg", eventSelection) : srcQ;
        t_xmax = std::max(srcQ, bgQ);
        if (t_xmax <= t_xmin) t_xmax = std::max(totalRawSrcMax, totalRawBgMax) * adcIntegralToNPE;
    }
    if (xMaxUser > xMinUser) {
        t_xmin = xMinUser;
        t_xmax = xMaxUser;
    }

    TH1D *hTotalSrc = new TH1D("hTotalSrc", "Total NPE;NPE;Counts", nBins, t_xmin, t_xmax);
    TH1D *hTotalBG = hasBackground ? new TH1D("hTotalBG", "", nBins, t_xmin, t_xmax) : nullptr;
    
    tSrc->Project("hTotalSrc", total_npe_expr, eventSelection);
    if (hasBackground) {
        tBG->Project("hTotalBG", total_npe_expr, eventSelection);
        hTotalBG->Scale(bgScale);
        style_source_bg(hTotalSrc, hTotalBG);
    } else {
        hTotalSrc->SetLineColor(kBlue + 1);
        hTotalSrc->SetLineWidth(2);
        hTotalSrc->SetFillStyle(0);
    }

    TH1D *hTotalSub = (TH1D*)hTotalSrc->Clone(hasBackground ? "hTotalSub" : "hTotalSourceOnly");
    hTotalSub->SetTitle(hasBackground ? Form("Total NPE (Source - %.3g #times BG, Gain=%.1e);NPE;Counts", bgScale, gain)
                                      : Form("Total NPE (Source only, Gain=%.1e);NPE;Counts", gain));
    if (hasBackground) hTotalSub->Add(hTotalBG, -1.0);
    histogramsToWrite.push_back(hTotalSrc);
    if (hTotalBG) histogramsToWrite.push_back(hTotalBG);
    histogramsToWrite.push_back(hTotalSub);

    TCanvas *c_total_overlay = new TCanvas("c_total_overlay", hasBackground ? "Total NPE Source vs Background" : "Total NPE Source", 900, 650);
    c_total_overlay->cd();
    hTotalSrc->SetTitle(hasBackground ? Form("Total NPE: Source vs scaled BG (Gain=%.1e);NPE;Counts", gain)
                                      : Form("Total NPE: Source (Gain=%.1e);NPE;Counts", gain));
    hTotalSrc->SetMinimum(0.0);
    hTotalSrc->SetMaximum((hasBackground ? std::max(hTotalSrc->GetMaximum(), hTotalBG->GetMaximum()) : hTotalSrc->GetMaximum()) * 1.15);
    hTotalSrc->Draw("HIST");
    if (hasBackground) hTotalBG->Draw("SAME HIST");
    TLegend *legOverlay = new TLegend(0.58, 0.72, 0.88, 0.88);
    legOverlay->AddEntry(hTotalSrc, srcLabel.Data(), "l");
    if (hasBackground) legOverlay->AddEntry(hTotalBG, Form("BG #times %.3g", bgScale), "l");
    legOverlay->Draw();
    c_total_overlay->SaveAs(Form("%s_overlay_total_linear.png", prefix.Data()));

    TCanvas *c_total_overlay_log = new TCanvas("c_total_overlay_log", hasBackground ? "Total NPE Source vs Background Log" : "Total NPE Source Log", 900, 650);
    c_total_overlay_log->cd();
    gPad->SetLogy();
    hTotalSrc->SetMinimum(0.5);
    hTotalSrc->SetMaximum((hasBackground ? std::max(hTotalSrc->GetMaximum(), hTotalBG->GetMaximum()) : hTotalSrc->GetMaximum()) * 2.0);
    hTotalSrc->Draw("HIST");
    if (hasBackground) hTotalBG->Draw("SAME HIST");
    TLegend *legOverlayLog = new TLegend(0.58, 0.72, 0.88, 0.88);
    legOverlayLog->AddEntry(hTotalSrc, srcLabel.Data(), "l");
    if (hasBackground) legOverlayLog->AddEntry(hTotalBG, Form("BG #times %.3g", bgScale), "l");
    legOverlayLog->Draw();
    c_total_overlay_log->SaveAs(Form("%s_overlay_total_log.png", prefix.Data()));

    const bool canMakeRatePlot = srcLiveSeconds > 0 && (!hasBackground || bgLiveSeconds > 0);
    if (canMakeRatePlot) {
        TH1D *hTotalSrcRate = make_rate_hist(hTotalSrc, srcLiveSeconds, "hTotalSrcRate");
        TH1D *hTotalBGRate = hasBackground ? make_rate_hist(hTotalBG, srcLiveSeconds, "hTotalBGRate") : nullptr;
        if (hasBackground) {
            style_source_bg(hTotalSrcRate, hTotalBGRate);
        } else {
            hTotalSrcRate->SetLineColor(kBlue + 1);
            hTotalSrcRate->SetLineWidth(2);
        }
        histogramsToWrite.push_back(hTotalSrcRate);
        if (hTotalBGRate) histogramsToWrite.push_back(hTotalBGRate);

        TCanvas *c_total_rate_log = new TCanvas("c_total_rate_log", hasBackground ? "Total NPE Rate-Normalized Source vs Background" : "Total NPE Rate-Normalized Source", 900, 650);
        c_total_rate_log->cd();
        gPad->SetLogy();
        hTotalSrcRate->SetTitle(hasBackground ? Form("Total NPE rate density: Source vs scaled BG (Gain=%.1e);NPE;Rate [Hz / NPE]", gain)
                                              : Form("Total NPE rate density: Source (Gain=%.1e);NPE;Rate [Hz / NPE]", gain));
        hTotalSrcRate->SetMinimum(1e-5);
        hTotalSrcRate->SetMaximum((hasBackground ? std::max(hTotalSrcRate->GetMaximum(), hTotalBGRate->GetMaximum()) : hTotalSrcRate->GetMaximum()) * 2.0);
        hTotalSrcRate->Draw("HIST");
        if (hasBackground) hTotalBGRate->Draw("SAME HIST");
        TLegend *legRate = new TLegend(0.54, 0.72, 0.88, 0.88);
        legRate->AddEntry(hTotalSrcRate, Form("%s %.1f Hz", srcLabel.Data(), srcRate), "l");
        if (hasBackground) legRate->AddEntry(hTotalBGRate, Form("BG scaled to source live time %.1f Hz", bgRate), "l");
        legRate->Draw();
        c_total_rate_log->SaveAs(Form("%s_overlay_total_rate_log.png", prefix.Data()));
    } else {
        std::cerr << "Warning: skipping rate-normalized overlay because SyncTime_TTT live time is unavailable." << std::endl;
    }

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
    leg->AddEntry(hTotalSub, hasBackground ? "Source - scaled BG" : srcLabel.Data(), "l");
    leg->Draw();

    c_total_sub->SaveAs(Form("%s_subtracted_total.png", prefix.Data()));

    TString histOutName = Form("%s_histograms.root", prefix.Data());
    TFile *fHistOut = TFile::Open(histOutName, "RECREATE");
    if (!fHistOut || fHistOut->IsZombie()) {
        std::cerr << "Error creating histogram output file: " << histOutName << std::endl;
    } else {
        fHistOut->cd();
        for (TH1D *hist : histogramsToWrite) {
            if (hist) hist->Write();
        }
        fHistOut->Close();
    }

    std::cout << "\n--- " << (hasBackground ? "Background Subtraction" : "Source-Only Histogram") << " Complete ---" << std::endl;
    std::cout << "Source: " << sourceFile << std::endl;
    std::cout << "Source label: " << srcLabel << std::endl;
    std::cout << "BG: " << (hasBackground ? bgFile : "none") << std::endl;
    std::cout << "Charge_CH unit: ADC-count sample integral from DAQ production_dt5730.cpp" << std::endl;
    std::cout << "ADC integral to pC: " << adcIntegralToPC << " pC / (ADC count * sample)" << std::endl;
    std::cout << "Conversion: NPE = Charge_CH " << npe_conv << std::endl;
    std::cout << "Conversion constants: dynamicRange=" << dynamicRangeV << " V, sampling="
              << samplingTimeNs << " ns, R=" << resistanceOhm << " ohm, ADC bits=" << adcBits
              << ", gain=" << gain << std::endl;
    std::cout << "Entries: source=" << tSrc->GetEntries();
    if (hasBackground) std::cout << ", BG=" << tBG->GetEntries();
    std::cout << std::endl;
    std::cout << "Live time [s]: source=" << srcLiveSeconds;
    if (hasBackground) std::cout << ", BG=" << bgLiveSeconds;
    std::cout << std::endl;
    std::cout << "Trigger rate [Hz]: source=" << srcRate;
    if (hasBackground) std::cout << ", BG=" << bgRate << ", excess=" << (srcRate - bgRate);
    std::cout << std::endl;
    std::cout << "Active channels:";
    for (int ch : activeChannels) std::cout << " CH" << ch;
    std::cout << std::endl;
    if (useRequestedChannels) std::cout << "Channel selection: " << selectedChannels << std::endl;
    if (!thresholds.empty()) {
        std::cout << "Event NPE thresholds:";
        for (const ChannelThreshold &threshold : thresholds) {
            std::cout << " CH" << threshold.channel << ">=" << threshold.minNpe;
        }
        std::cout << std::endl;
        std::cout << "ROOT event selection: " << eventSelection << std::endl;
    }
    if (hasBackground) std::cout << "BG scale: " << bgScale << std::endl;
    std::cout << "X-axis quantile: " << xQuantile << std::endl;
    std::cout << "Bins: " << nBins << ", X range: [" << t_xmin << ", " << t_xmax << "]" << std::endl;
    std::cout << "Output prefix: " << prefix << std::endl;
    std::cout << "Histogram ROOT file: " << histOutName << std::endl;
    const auto endTime = std::chrono::steady_clock::now();
    const double elapsedSeconds = std::chrono::duration<double>(endTime - startTime).count();
    std::cout << "Elapsed wall time [s]: " << elapsedSeconds << std::endl;
    std::cout << "Results saved to: " << prefix << "_subtracted_channels.png, "
              << prefix << "_subtracted_total.png" << std::endl;
    std::cout << "Overlay plots saved to: " << prefix << "_overlay_channels_log.png, "
              << prefix << "_overlay_total_linear.png, "
              << prefix << "_overlay_total_log.png";
    if (canMakeRatePlot) std::cout << ", " << prefix << "_overlay_total_rate_log.png";
    std::cout << std::endl;
    std::cout << "Histograms saved to: " << histOutName << std::endl;
}
