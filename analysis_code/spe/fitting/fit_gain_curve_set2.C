#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>

#include "spe_fit_common.h"

namespace {

struct SpeFitPoint {
    int channel = -1;
    int voltage = 0;
    double q0 = 0.0;
    double q0Err = 0.0;
    double qSpe = 0.0;
    double qSpeErr = 0.0;
    double gain = 0.0;
    double gainErr = 0.0;
    double mu = 0.0;
    double chi2 = 0.0;
    int ndf = 0;
    double ledPedestalBias = 0.0;
    double ledPedestalBiasErr = 0.0;
    int status = -1;
};

Double_t spe_model_set2(Double_t *x, Double_t *par) {
    const Double_t xx = x[0];
    const Double_t nTotal = par[0];
    const Double_t mu = par[1];
    const Double_t q0 = par[2];
    const Double_t s0 = par[3];
    const Double_t q1 = par[4];
    const Double_t s1 = par[5];
    const Double_t bg = par[6];
    const Double_t expAmp = par[7];
    const Double_t expTau = par[8];

    Double_t sum = 0.0;
    for (int n = 0; n <= 8; ++n) {
        const Double_t pn = TMath::Poisson(n, mu);
        const Double_t qn = q0 + n * q1;
        const Double_t sn = TMath::Sqrt(s0 * s0 + n * s1 * s1);
        const Double_t gn = TMath::Gaus(xx, qn, sn, true);
        sum += pn * gn;
    }

    Double_t lowChargeTail = 0.0;
    if (xx > q0 && expTau > 0.0) {
        lowChargeTail = expAmp * TMath::Exp(-(xx - q0) / expTau);
    }

    return nTotal * sum + bg + lowChargeTail;
}

std::vector<int> parse_channels(const char *channels) {
    std::vector<int> parsed;
    std::stringstream ss(channels ? channels : "0,1");
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(std::remove_if(item.begin(), item.end(), [](unsigned char c) {
            return std::isspace(c);
        }), item.end());
        if (!item.empty()) parsed.push_back(std::stoi(item));
    }
    if (parsed.empty()) parsed = {0, 1};
    return parsed;
}

bool valid_root_file(const TString &path) {
    if (gSystem->AccessPathName(path)) return false;

    Long_t id = 0;
    Long64_t size = 0;
    Long_t flags = 0;
    Long_t modtime = 0;
    gSystem->GetPathInfo(path, &id, &size, &flags, &modtime);
    if (size <= 0) return false;

    TFile file(path, "READ");
    return !file.IsZombie();
}

TString join_path(const char *dir, const char *name) {
    TString out(dir);
    if (!out.EndsWith("/")) out += "/";
    out += name;
    return out;
}

bool estimate_pedestal(TH1D *hist, double &q0, double &s0, double &q0Err) {
    if (!hist || hist->GetEntries() <= 0) return false;

    hist->GetXaxis()->SetRangeUser(-0.3, 0.5);
    const double peak = hist->GetBinCenter(hist->GetMaximumBin());
    hist->GetXaxis()->SetRange(0, 0);

    TF1 ped("ped_set2", "gaus", peak - 0.12, peak + 0.12);
    const int status = hist->Fit(&ped, "RQN");
    if (status != 0 || ped.GetParameter(2) <= 0.0) {
        q0 = peak;
        s0 = 0.05;
        q0Err = 0.0;
        return true;
    }

    q0 = ped.GetParameter(1);
    s0 = std::max(std::abs(ped.GetParameter(2)), 0.005);
    q0Err = ped.GetParError(1);
    return true;
}

bool estimate_dark_bias(TH1D *hist, double &bias, double &biasErr) {
    if (!hist || hist->GetEntries() <= 0) return false;

    hist->GetXaxis()->SetRangeUser(-1.0, 10.0);
    const double peak = hist->GetBinCenter(hist->GetMaximumBin());
    hist->GetXaxis()->SetRange(0, 0);

    TF1 fit("dark_bias_set2", "gaus", peak - 0.35, peak + 0.35);
    const int status = hist->Fit(&fit, "RQN");
    if (status != 0 || fit.GetParameter(2) <= 0.0) {
        bias = peak;
        biasErr = 0.0;
        return true;
    }

    bias = fit.GetParameter(1);
    biasErr = fit.GetParError(1);
    return true;
}

bool estimate_led_pedestal_bias(TH1D *hist, double &bias, double &biasErr) {
    if (!hist || hist->GetEntries() <= 0) return false;

    hist->GetXaxis()->SetRangeUser(0.3, std::min(10.0, hist->GetXaxis()->GetXmax()));
    const double peak = hist->GetBinCenter(hist->GetMaximumBin());
    hist->GetXaxis()->SetRange(0, 0);

    TF1 fit("led_pedestal_bias_set2", "gaus", peak - 0.45, peak + 0.45);
    const int status = hist->Fit(&fit, "RQN");
    if (status != 0 || fit.GetParameter(2) <= 0.0) {
        bias = peak;
        biasErr = 0.0;
        return true;
    }

    bias = fit.GetParameter(1);
    biasErr = fit.GetParError(1);
    return true;
}

double auto_xmax(int voltage, double requested) {
    if (requested > 0.0) return requested;
    return (voltage < 2000) ? 8.0 : 45.0;
}

void style_debug_hist(TH1D *hist, int index) {
    const Color_t colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2,
                              kMagenta + 1, kOrange + 7, kCyan + 1, kViolet};
    if (!hist) return;
    hist->SetDirectory(nullptr);
    hist->SetLineColor(colors[index % 8]);
    hist->SetLineWidth(2);
    hist->SetStats(kFALSE);
}

TH1D *make_pull_hist(const TH1D *hist, TF1 *fit, const char *name) {
    TH1D *pull = static_cast<TH1D*>(hist->Clone(name));
    pull->Reset("ICES");
    pull->SetTitle(";Charge - LED pedestal [pC];Pull");
    pull->SetStats(kFALSE);
    pull->SetLineColor(kBlack);
    pull->SetMarkerStyle(20);
    pull->SetMarkerSize(0.55);
    pull->SetMarkerColor(kBlack);

    double fitMin = 0.0;
    double fitMax = 0.0;
    fit->GetRange(fitMin, fitMax);

    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        const double x = hist->GetBinCenter(bin);
        if (x < fitMin || x > fitMax) continue;
        const double observed = hist->GetBinContent(bin);
        const double expected = fit->Eval(x);
        const double sigma = std::sqrt(std::max(1.0, observed));
        pull->SetBinContent(bin, (observed - expected) / sigma);
        pull->SetBinError(bin, 1.0);
    }
    return pull;
}

void save_fit_diagnostic_canvas(TH1D *hist,
                                TF1 *fit,
                                const TString &fileName,
                                int ch,
                                int voltage,
                                double ledPedestalBias,
                                double qSpe,
                                double gain,
                                double chi2,
                                int ndf,
                                int fitStatus,
                                const char *modelName) {
    TCanvas canvas(TString::Format("c_fit_diag_ch%d_%d", ch, voltage), "", 950, 820);

    TPad top("top", "top", 0.0, 0.30, 1.0, 1.0);
    TPad bottom("bottom", "bottom", 0.0, 0.0, 1.0, 0.30);
    top.SetBottomMargin(0.02);
    bottom.SetTopMargin(0.03);
    bottom.SetBottomMargin(0.30);
    top.Draw();
    bottom.Draw();

    top.cd();
    top.SetLogy();
    hist->SetTitle(TString::Format("Ch%d %dV SPE Fit;Charge - LED pedestal [pC];Counts", ch, voltage));
    hist->GetXaxis()->SetLabelSize(0.0);
    hist->Draw();
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    fit->Draw("same");

    TLatex text;
    text.SetNDC(kTRUE);
    text.SetTextSize(0.035);
    const double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
    text.DrawLatex(0.58, 0.84, TString::Format("#chi^{2}/NDF = %.1f/%d = %.2f", chi2, ndf, chi2ndf));
    text.DrawLatex(0.58, 0.79, TString::Format("status = %d", fitStatus));
    text.DrawLatex(0.58, 0.74, TString::Format("model = %s", modelName));
    text.DrawLatex(0.58, 0.69, TString::Format("LED ped. bias = %.3f pC", ledPedestalBias));
    text.DrawLatex(0.58, 0.64, TString::Format("Q_{SPE} = %.3f pC", qSpe));
    text.DrawLatex(0.58, 0.59, TString::Format("Gain = %.3g", gain));

    bottom.cd();
    TH1D *pull = make_pull_hist(hist, fit, TString::Format("pull_ch%d_%d", ch, voltage));
    pull->SetMinimum(-6.0);
    pull->SetMaximum(6.0);
    pull->GetYaxis()->SetTitleSize(0.09);
    pull->GetYaxis()->SetTitleOffset(0.45);
    pull->GetYaxis()->SetLabelSize(0.08);
    pull->GetXaxis()->SetTitleSize(0.10);
    pull->GetXaxis()->SetLabelSize(0.08);
    pull->Draw("E1");

    double xmin = 0.0;
    double xmax = 0.0;
    fit->GetRange(xmin, xmax);
    TLine zero(xmin, 0.0, xmax, 0.0);
    zero.SetLineColor(kRed + 1);
    zero.SetLineWidth(2);
    zero.Draw("same");
    TLine plus(xmin, 3.0, xmax, 3.0);
    plus.SetLineColor(kGray + 2);
    plus.SetLineStyle(2);
    plus.Draw("same");
    TLine minus(xmin, -3.0, xmax, -3.0);
    minus.SetLineColor(kGray + 2);
    minus.SetLineStyle(2);
    minus.Draw("same");

    canvas.SaveAs(fileName);
    delete pull;
}

} // namespace

void make_set2_debug_charge_plots(const char *dataDir = "/home/sgwon/MD/data/SPE/5inch_set2",
                                  const char *outDir = "set2_spe_results",
                                  const char *channelsText = "0,1",
                                  int vStart = 1600,
                                  int vStop = 2200,
                                  int vStep = 100,
                                  int bins = 700,
                                  double xMin = -1.0,
                                  double xMax = 45.0,
                                  double dynamicRangeV = 2.0,
                                  double samplingTimeNs = 2.0,
                                  double resistanceOhm = 50.0,
                                  int adcBits = 14) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(outDir, kTRUE);

    const std::vector<int> channels = parse_channels(channelsText);

    for (int ch : channels) {
        std::vector<TH1D*> ledHists;
        std::vector<TH1D*> darkHists;
        std::vector<TH1D*> subHists;
        std::vector<int> voltages;

        for (int v = vStart; v <= vStop; v += vStep) {
            const TString ledFile = join_path(dataDir, TString::Format("LED_RUN_%dV_prod.root", v));
            const TString darkFile = join_path(dataDir, TString::Format("DARK_RUN_%dV_prod.root", v));
            if (!valid_root_file(ledFile) || !valid_root_file(darkFile)) {
                std::cout << "skip debug ch" << ch << " " << v << "V, missing LED or DARK" << std::endl;
                continue;
            }

            TH1D *hLed = spe_make_charge_hist(
                ledFile, ch, TString::Format("debug_led_ch%d_%d", ch, v),
                TString::Format("Raw LED Ch%d;Charge [pC];Counts", ch),
                bins, xMin, xMax, dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
            TH1D *hDark = spe_make_charge_hist(
                darkFile, ch, TString::Format("debug_dark_ch%d_%d", ch, v),
                TString::Format("Raw DARK Ch%d;Charge [pC];Counts", ch),
                bins, xMin, xMax, dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
            if (!hLed || !hDark) {
                delete hLed;
                delete hDark;
                continue;
            }

            TH1D *hSub = static_cast<TH1D*>(hLed->Clone(TString::Format("debug_led_minus_dark_ch%d_%d", ch, v)));
            hSub->SetTitle(TString::Format("Raw LED - DARK Ch%d;Charge [pC];LED counts - scaled DARK counts", ch));
            if (hDark->Integral() > 0.0) {
                hSub->Add(hDark, -hLed->Integral() / hDark->Integral());
            }

            const int index = static_cast<int>(voltages.size());
            style_debug_hist(hLed, index);
            style_debug_hist(hDark, index);
            style_debug_hist(hSub, index);

            ledHists.push_back(hLed);
            darkHists.push_back(hDark);
            subHists.push_back(hSub);
            voltages.push_back(v);
        }

        auto draw_overlay = [&](const std::vector<TH1D*> &hists, const char *canvasName,
                                const char *title, const char *fileName, bool logY) {
            if (hists.empty()) return;

            TCanvas canvas(canvasName, title, 980, 720);
            if (logY) canvas.SetLogy();
            canvas.SetGrid();

            double ymax = 0.0;
            double ymin = logY ? 0.5 : hists.front()->GetMinimum();
            for (TH1D *hist : hists) {
                ymax = std::max(ymax, hist->GetMaximum());
                if (!logY) ymin = std::min(ymin, hist->GetMinimum());
            }
            if (logY) {
                hists.front()->SetMinimum(ymin);
                hists.front()->SetMaximum(std::max(10.0, ymax * 4.0));
            } else {
                hists.front()->SetMinimum(std::min(0.0, ymin * 1.2));
                hists.front()->SetMaximum(std::max(1.0, ymax * 1.25));
            }

            hists.front()->SetTitle(title);
            hists.front()->Draw("hist");
            for (size_t i = 1; i < hists.size(); ++i) hists[i]->Draw("hist same");

            TLegend legend(0.70, 0.66, 0.90, 0.90);
            legend.SetBorderSize(1);
            for (size_t i = 0; i < hists.size(); ++i) {
                legend.AddEntry(hists[i], TString::Format("%d V", voltages[i]), "L");
            }
            legend.Draw();
            canvas.SaveAs(join_path(outDir, fileName));
        };

        draw_overlay(ledHists,
                     TString::Format("c_debug_raw_led_ch%d", ch),
                     TString::Format("Raw LED Charge Ch%d;Charge [pC];Counts", ch),
                     TString::Format("debug_raw_led_ch%d.png", ch),
                     true);
        draw_overlay(darkHists,
                     TString::Format("c_debug_raw_dark_ch%d", ch),
                     TString::Format("Raw DARK Charge Ch%d;Charge [pC];Counts", ch),
                     TString::Format("debug_raw_dark_ch%d.png", ch),
                     true);
        draw_overlay(subHists,
                     TString::Format("c_debug_led_minus_dark_ch%d", ch),
                     TString::Format("Raw LED - DARK Charge Ch%d;Charge [pC];LED counts - scaled DARK counts", ch),
                     TString::Format("debug_led_minus_dark_ch%d.png", ch),
                     false);

        for (TH1D *hist : ledHists) delete hist;
        for (TH1D *hist : darkHists) delete hist;
        for (TH1D *hist : subHists) delete hist;
    }
}

void fit_gain_curve_set2(const char *dataDir = "/home/sgwon/MD/data/SPE/5inch_set2",
                         const char *outDir = "set2_spe_results",
                         const char *channelsText = "0,1",
                         int vStart = 1600,
                         int vStop = 2200,
                         int vStep = 100,
                         int bins = 600,
                         double xMin = -1.0,
                         double xMax = -1.0,
                         double dynamicRangeV = 2.0,
                         double samplingTimeNs = 2.0,
                         double resistanceOhm = 50.0,
                         int adcBits = 14,
                         bool alignLedPeakToZero = false) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptFit(0);
    gSystem->mkdir(outDir, kTRUE);

    const std::vector<int> channels = parse_channels(channelsText);
    const double electronCharge = 1.602176634e-19;
    std::vector<SpeFitPoint> results;

    std::ofstream csv(join_path(outDir, "gain_results.csv").Data());
    csv << "channel,voltage,led_pedestal_bias_pC,led_pedestal_bias_err_pC,"
        << "q0_pC,q0_err_pC,"
        << "q_spe_pC,q_spe_err_pC,gain,gain_err,mu,chi2,ndf,chi2_ndf,model,fit_status\n";

    for (int ch : channels) {
        std::cout << "\n>>> Channel " << ch << std::endl;

        for (int v = vStart; v <= vStop; v += vStep) {
            SpeFitPoint point;
            point.channel = ch;
            point.voltage = v;

            const TString ledFile = join_path(dataDir, TString::Format("LED_RUN_%dV_prod.root", v));
            const TString darkFile = join_path(dataDir, TString::Format("DARK_RUN_%dV_prod.root", v));
            std::cout << "  " << v << " V: ";

            if (!valid_root_file(ledFile)) {
                std::cout << "skip, invalid LED file: " << ledFile << std::endl;
                point.status = -10;
                results.push_back(point);
                csv << ch << "," << v << ",,,,,,,,,,,,," << point.status << "\n";
                continue;
            }

            const double xmaxThis = auto_xmax(v, xMax);
            TH1D *hLedRaw = spe_make_charge_hist(
                ledFile, ch, TString::Format("h_led_raw_for_bias_ch%d_%d", ch, v),
                "Raw LED;Charge [pC];Counts",
                bins, xMin, xmaxThis, dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
            double ledPedestalBias = 0.0;
            double ledPedestalBiasErr = 0.0;
            if (!estimate_led_pedestal_bias(hLedRaw, ledPedestalBias, ledPedestalBiasErr)) {
                std::cout << "skip, LED pedestal estimate failed" << std::endl;
                point.status = -13;
                delete hLedRaw;
                results.push_back(point);
                csv << ch << "," << v << ",,,,,,,,,,,,," << point.status << "\n";
                continue;
            }
            delete hLedRaw;

            point.ledPedestalBias = ledPedestalBias;
            point.ledPedestalBiasErr = ledPedestalBiasErr;

            TH1D *hLed = spe_make_charge_hist(
                ledFile, ch, TString::Format("h_led_ch%d_%d", ch, v),
                TString::Format("Set2 Ch%d %dV (LED pedestal corrected);Charge - LED pedestal [pC];Counts", ch, v),
                bins, xMin, xmaxThis, dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits, ledPedestalBias);

            if (!hLed || hLed->GetEntries() <= 0) {
                std::cout << "skip, empty LED histogram" << std::endl;
                point.status = -11;
                delete hLed;
                results.push_back(point);
                csv << ch << "," << v << "," << ledPedestalBias << "," << ledPedestalBiasErr
                    << ",,,,,,,,,," << point.status << "\n";
                continue;
            }

            double q0 = 0.0;
            double s0 = 0.05;
            double q0Err = 0.0;
            bool havePed = estimate_pedestal(hLed, q0, s0, q0Err);
            if (!havePed) {
                std::cout << "skip, pedestal estimate failed" << std::endl;
                point.status = -12;
                delete hLed;
                results.push_back(point);
                csv << ch << "," << v << "," << ledPedestalBias << "," << ledPedestalBiasErr
                    << ",,,,,,,,,," << point.status << "\n";
                continue;
            }

            double searchMin = q0 + std::max(5.0 * s0, 0.8);
            hLed->GetXaxis()->SetRangeUser(searchMin, xmaxThis * 0.8);
            double q1Guess = hLed->GetBinCenter(hLed->GetMaximumBin()) - q0;
            hLed->GetXaxis()->SetRange(0, 0);
            if (q1Guess < 0.5) q1Guess = std::max(0.8, (v - 1500.0) * 0.001 + 0.6);

            const double fitMax = std::min(xmaxThis, q0 + q1Guess * 5.0);
            auto configureFit = [&](TF1 &candidate, bool enableTail) {
                candidate.SetParNames("Norm", "mu", "Q_ped", "s_ped", "Q_spe", "s_spe", "BG", "TailAmp", "TailTau");
                candidate.SetParameters(hLed->Integral() * hLed->GetBinWidth(1), 0.5, q0, s0, q1Guess, q1Guess * 0.4,
                                        1.0, enableTail ? hLed->GetMaximum() * 0.05 : 0.0, std::max(0.5, q1Guess * 0.8));
                candidate.SetParLimits(1, 0.01, 4.0);
                candidate.SetParLimits(2, q0 - 0.15, q0 + 0.15);
                candidate.SetParLimits(3, std::max(0.003, s0 * 0.5), std::max(0.3, s0 * 1.8));
                candidate.SetParLimits(4, std::max(0.5, q1Guess * 0.6), q1Guess * 1.4);
                candidate.SetParLimits(5, q1Guess * 0.08, q1Guess * 1.2);
                candidate.SetParLimits(6, 0.0, std::max(1.0, hLed->GetMaximum()));
                if (enableTail) {
                    candidate.SetParLimits(7, 0.0, std::max(1.0, hLed->GetMaximum()));
                    candidate.SetParLimits(8, 0.1, std::max(2.0, fitMax - q0));
                } else {
                    candidate.FixParameter(7, 0.0);
                    candidate.FixParameter(8, 1.0);
                }
            };

            TF1 fitBase(TString::Format("f_spe_set2_base_ch%d_%d", ch, v), spe_model_set2, xMin, fitMax, 9);
            TF1 fitTail(TString::Format("f_spe_set2_tail_ch%d_%d", ch, v), spe_model_set2, xMin, fitMax, 9);
            configureFit(fitBase, false);
            configureFit(fitTail, true);

            const int baseStatus = hLed->Fit(&fitBase, "RQN");
            const int tailStatus = hLed->Fit(&fitTail, "RQN");
            const double baseChi2Ndf = (fitBase.GetNDF() > 0) ? fitBase.GetChisquare() / fitBase.GetNDF() : 1.0e99;
            const double tailChi2Ndf = (fitTail.GetNDF() > 0) ? fitTail.GetChisquare() / fitTail.GetNDF() : 1.0e99;
            const bool useTail = (tailStatus == 0 && tailChi2Ndf < baseChi2Ndf * 0.85);
            TF1 *fit = useTail ? &fitTail : &fitBase;
            const int fitStatus = useTail ? tailStatus : baseStatus;
            const char *modelName = useTail ? "ped+spe+tail" : "ped+spe";
            point.status = fitStatus;
            point.q0 = fit->GetParameter(2);
            point.q0Err = fit->GetParError(2);
            point.qSpe = fit->GetParameter(4);
            point.qSpeErr = fit->GetParError(4);
            point.mu = fit->GetParameter(1);
            point.gain = point.qSpe * 1.0e-12 / electronCharge;
            point.gainErr = point.qSpeErr * 1.0e-12 / electronCharge;
            point.chi2 = fit->GetChisquare();
            point.ndf = fit->GetNDF();
            const double chi2ndf = (point.ndf > 0) ? point.chi2 / point.ndf : 0.0;

            csv << point.channel << "," << point.voltage << ","
                << point.ledPedestalBias << "," << point.ledPedestalBiasErr << ","
                << point.q0 << "," << point.q0Err << ","
                << point.qSpe << "," << point.qSpeErr << ","
                << point.gain << "," << point.gainErr << ","
                << point.mu << "," << point.chi2 << "," << point.ndf << ","
                << chi2ndf << "," << modelName << "," << point.status << "\n";
            results.push_back(point);

            save_fit_diagnostic_canvas(
                hLed, fit, join_path(outDir, TString::Format("fit_ch%d_%dV.png", ch, v)),
                ch, v, ledPedestalBias, point.qSpe, point.gain, point.chi2, point.ndf, fitStatus, modelName);
            std::cout << "led_pedestal_bias=" << ledPedestalBias
                      << " pC, gain=" << point.gain << ", q_spe=" << point.qSpe
                      << " pC, chi2/ndf=" << chi2ndf
                      << ", model=" << modelName
                      << ", status=" << fitStatus << std::endl;

            delete hLed;
        }
    }

    TCanvas cg("c_gain_set2", "Set2 Gain Curve", 950, 700);
    cg.SetLogy();
    cg.SetGrid();
    const double maxGainFitChi2Ndf = 5.0;
    TH1F *frame = cg.DrawFrame(vStart - 50, 5.0e5, vStop + 50, 5.0e8);
    frame->SetTitle(TString::Format("Set2 Gain vs Voltage (log-linear fit: #chi^{2}/NDF < %.1f);Voltage [V];Gain",
                                    maxGainFitChi2Ndf));

    TLegend leg(0.14, 0.72, 0.44, 0.88);
    leg.SetBorderSize(1);

    const Color_t colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                              kOrange + 7, kCyan + 1, kViolet, kGray + 2};
    bool hasPoint = false;

    for (size_t i = 0; i < channels.size(); ++i) {
        const int ch = channels[i];
        TGraphErrors *grAll = new TGraphErrors();
        grAll->SetName(TString::Format("gr_gain_set2_all_ch%d", ch));
        grAll->SetMarkerStyle(24 + (i % 10));
        grAll->SetMarkerColor(colors[i % 8]);
        grAll->SetLineColor(colors[i % 8]);

        TGraphErrors *grFit = new TGraphErrors();
        grFit->SetName(TString::Format("gr_gain_set2_fit_ch%d", ch));
        grFit->SetMarkerStyle(20 + (i % 10));
        grFit->SetMarkerColor(colors[i % 8]);
        grFit->SetLineColor(colors[i % 8]);

        for (const auto &p : results) {
            if (p.channel != ch || p.status != 0 || p.gain <= 0.0) continue;
            if (!std::isfinite(p.gain) || !std::isfinite(p.gainErr)) continue;

            const int nAll = grAll->GetN();
            grAll->SetPoint(nAll, p.voltage, p.gain);
            grAll->SetPointError(nAll, 0.0, p.gainErr);

            const double chi2ndf = (p.ndf > 0) ? p.chi2 / p.ndf : 1.0e99;
            if (chi2ndf >= maxGainFitChi2Ndf) continue;

            const int nFit = grFit->GetN();
            grFit->SetPoint(nFit, p.voltage, p.gain);
            grFit->SetPointError(nFit, 0.0, p.gainErr);
        }

        if (grAll->GetN() == 0) {
            delete grAll;
            delete grFit;
            continue;
        }
        hasPoint = true;
        grAll->Draw("PE same");
        if (grFit->GetN() > 0) grFit->Draw("PE same");

        if (grFit->GetN() >= 3) {
            TF1 *gainFit = new TF1(TString::Format("f_gain_set2_ch%d", ch),
                                   "TMath::Power(10.0, [0] + [1]*x)", vStart - 50, vStop + 50);
            gainFit->SetLineColor(colors[i % 8]);
            gainFit->SetParameters(-2.0, 0.005);
            grFit->Fit(gainFit, "RQ");
            leg.AddEntry(grFit, TString::Format("Ch%d used", ch), "LP");
            leg.AddEntry(grAll, TString::Format("Ch%d excluded", ch), "P");
        } else {
            leg.AddEntry(grAll, TString::Format("Ch%d all points", ch), "P");
        }
    }

    if (hasPoint) {
        leg.Draw();
        cg.SaveAs(join_path(outDir, "gain_curve_set2.png"));
        std::cout << "\nWrote " << join_path(outDir, "gain_curve_set2.png") << std::endl;
    } else {
        std::cout << "\nNo valid gain points. Check input ROOT files and "
                  << join_path(outDir, "gain_results.csv") << std::endl;
    }

    std::cout << "Wrote " << join_path(outDir, "gain_results.csv") << std::endl;
}
