#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace optical_gel_compare {

struct Metrics {
    double signedRate = 0.0;
    double signedRateError = 0.0;
    double positiveRate = 0.0;
    double meanNPE = 0.0;
    double p90NPE = 0.0;
    double p95NPE = 0.0;
};

TH1D *load_hist(TFile *file, const char *name, const char *suffix) {
    if (!file || file->IsZombie()) return nullptr;
    TH1D *hist = dynamic_cast<TH1D*>(file->Get(name));
    if (!hist) {
        std::cerr << "Error: missing histogram '" << name << "' in " << file->GetName() << std::endl;
        return nullptr;
    }
    TH1D *clone = dynamic_cast<TH1D*>(hist->Clone(Form("%s_%s", name, suffix)));
    if (clone) clone->SetDirectory(nullptr);
    return clone;
}

bool same_binning(const TH1D *a, const TH1D *b) {
    if (!a || !b || a->GetNbinsX() != b->GetNbinsX()) return false;
    const double scale = std::max({1.0,
                                   std::abs(a->GetXaxis()->GetXmin()),
                                   std::abs(a->GetXaxis()->GetXmax())});
    const double tolerance = 1e-12 * scale;
    return std::abs(a->GetXaxis()->GetXmin() - b->GetXaxis()->GetXmin()) <= tolerance &&
           std::abs(a->GetXaxis()->GetXmax() - b->GetXaxis()->GetXmax()) <= tolerance;
}

TH1D *positive_roi_hist(const TH1D *input,
                        const char *name,
                        double roiMin,
                        double roiMax) {
    if (!input || roiMax <= roiMin) return nullptr;
    TH1D *hist = dynamic_cast<TH1D*>(input->Clone(name));
    if (!hist) return nullptr;
    hist->SetDirectory(nullptr);
    for (int bin = 0; bin <= hist->GetNbinsX() + 1; ++bin) {
        const bool regularBin = bin >= 1 && bin <= hist->GetNbinsX();
        const double x = regularBin ? hist->GetBinCenter(bin) : 0.0;
        if (!regularBin || x <= roiMin || x >= roiMax || hist->GetBinContent(bin) < 0.0) {
            hist->SetBinContent(bin, 0.0);
            hist->SetBinError(bin, 0.0);
        }
    }
    return hist;
}

double quantile(const TH1D *positive, double probability) {
    if (!positive || positive->Integral(1, positive->GetNbinsX()) <= 0.0) return 0.0;
    double value = 0.0;
    double p = probability;
    const_cast<TH1D*>(positive)->GetQuantiles(1, &value, &p);
    return value;
}

Metrics compute_metrics(const TH1D *input,
                        double elapsedSeconds,
                        double roiMin,
                        double roiMax,
                        const char *suffix) {
    Metrics result;
    if (!input || elapsedSeconds <= 0.0 || roiMax <= roiMin) return result;

    int firstBin = input->FindFixBin(roiMin);
    int lastBin = input->FindFixBin(roiMax);
    firstBin = std::max(1, firstBin);
    lastBin = std::min(input->GetNbinsX(), lastBin);
    while (firstBin <= lastBin && input->GetBinCenter(firstBin) <= roiMin) ++firstBin;
    while (lastBin >= firstBin && input->GetBinCenter(lastBin) >= roiMax) --lastBin;

    double signedError = 0.0;
    if (firstBin <= lastBin) {
        result.signedRate = input->IntegralAndError(firstBin, lastBin, signedError) / elapsedSeconds;
        result.signedRateError = signedError / elapsedSeconds;
    }

    TH1D *positive = positive_roi_hist(input, Form("hPositive_%s", suffix), roiMin, roiMax);
    if (!positive) return result;
    const double positiveCounts = positive->Integral(1, positive->GetNbinsX());
    if (positiveCounts > 0.0) {
        result.positiveRate = positiveCounts / elapsedSeconds;
        result.meanNPE = positive->GetMean();
        result.p90NPE = quantile(positive, 0.90);
        result.p95NPE = quantile(positive, 0.95);
    }
    delete positive;
    return result;
}

double safe_ratio(double withValue, double withoutValue) {
    return withoutValue != 0.0 ? withValue / withoutValue : 0.0;
}

double ratio_error(double withValue,
                   double withError,
                   double withoutValue,
                   double withoutError) {
    if (withValue == 0.0 || withoutValue == 0.0) return 0.0;
    const double value = withValue / withoutValue;
    return std::abs(value) * std::sqrt(std::pow(withError / withValue, 2) +
                                       std::pow(withoutError / withoutValue, 2));
}

void style_hist(TH1D *hist, Color_t color, Style_t lineStyle) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    hist->SetLineStyle(lineStyle);
    hist->SetFillStyle(0);
}

void scale_to_rate_density(TH1D *hist, double elapsedSeconds) {
    if (hist && elapsedSeconds > 0.0) hist->Scale(1.0 / elapsedSeconds, "width");
}

void set_common_y_range(TH1D *a, TH1D *b) {
    if (!a || !b) return;
    const double maximum = std::max(a->GetMaximum(), b->GetMaximum());
    const double minimum = std::min(a->GetMinimum(), b->GetMinimum());
    const double span = std::max(1e-12, maximum - minimum);
    a->SetMaximum(maximum + 0.15 * span);
    a->SetMinimum(std::min(0.0, minimum - 0.08 * span));
}

void save_total_overlay(TH1D *without,
                        TH1D *with,
                        double withoutElapsed,
                        double withElapsed,
                        const TString &sourceLabel,
                        const TString &labelWithout,
                        const TString &labelWith,
                        const TString &outPrefix) {
    TH1D *withoutRate = dynamic_cast<TH1D*>(without->Clone("hTotalRateDensityWithoutGel"));
    TH1D *withRate = dynamic_cast<TH1D*>(with->Clone("hTotalRateDensityWithGel"));
    withoutRate->SetDirectory(nullptr);
    withRate->SetDirectory(nullptr);
    scale_to_rate_density(withoutRate, withoutElapsed);
    scale_to_rate_density(withRate, withElapsed);
    style_hist(withoutRate, kBlue + 1, 1);
    style_hist(withRate, kRed + 1, 2);
    set_common_y_range(withoutRate, withRate);
    withoutRate->SetTitle(Form("%s optical-gel comparison;NPE;Elapsed-time-normalized rate [Hz / NPE]",
                               sourceLabel.Data()));

    TCanvas canvas("cOpticalGelTotal", "Optical-gel total comparison", 950, 680);
    withoutRate->Draw("HIST");
    withRate->Draw("SAME HIST");
    TLegend legend(0.61, 0.74, 0.88, 0.88);
    legend.AddEntry(withoutRate, labelWithout, "l");
    legend.AddEntry(withRate, labelWith, "l");
    legend.Draw();
    canvas.SaveAs(Form("%s_total_overlay.png", outPrefix.Data()));
    delete withoutRate;
    delete withRate;
}

void save_channel_overlay(TFile *withoutFile,
                          TFile *withFile,
                          double withoutElapsed,
                          double withElapsed,
                          const TString &sourceLabel,
                          const TString &labelWithout,
                          const TString &labelWith,
                          const TString &outPrefix,
                          const std::vector<int> &channels) {
    if (channels.empty()) return;
    TCanvas canvas("cOpticalGelChannels", "Optical-gel channel comparison", 650 * channels.size(), 560);
    canvas.Divide(static_cast<int>(channels.size()), 1);
    std::vector<TH1D*> owned;
    for (size_t index = 0; index < channels.size(); ++index) {
        const int channel = channels[index];
        TH1D *without = load_hist(withoutFile, Form("hSub_ch%d", channel), Form("without_ch%d", channel));
        TH1D *with = load_hist(withFile, Form("hSub_ch%d", channel), Form("with_ch%d", channel));
        if (!without || !with || !same_binning(without, with)) {
            delete without;
            delete with;
            continue;
        }
        owned.push_back(without);
        owned.push_back(with);
        scale_to_rate_density(without, withoutElapsed);
        scale_to_rate_density(with, withElapsed);
        style_hist(without, kBlue + 1, 1);
        style_hist(with, kRed + 1, 2);
        set_common_y_range(without, with);
        without->SetTitle(Form("%s CH%d background-subtracted spectrum;NPE;Rate [Hz / NPE]",
                               sourceLabel.Data(), channel));
        canvas.cd(static_cast<int>(index) + 1);
        without->Draw("HIST");
        with->Draw("SAME HIST");
        TLegend *legend = new TLegend(0.57, 0.76, 0.88, 0.88);
        legend->AddEntry(without, labelWithout, "l");
        legend->AddEntry(with, labelWith, "l");
        legend->Draw();
    }
    canvas.SaveAs(Form("%s_channels_overlay.png", outPrefix.Data()));
    for (TH1D *hist : owned) delete hist;
}

void write_metrics(std::ostream &out,
                   const char *scope,
                   const Metrics &without,
                   const Metrics &with) {
    out << scope << " without"
        << " signed_rate=" << without.signedRate
        << " signed_rate_error=" << without.signedRateError
        << " positive_rate=" << without.positiveRate
        << " mean_npe=" << without.meanNPE
        << " p90_npe=" << without.p90NPE
        << " p95_npe=" << without.p95NPE << '\n';
    out << scope << " with"
        << " signed_rate=" << with.signedRate
        << " signed_rate_error=" << with.signedRateError
        << " positive_rate=" << with.positiveRate
        << " mean_npe=" << with.meanNPE
        << " p90_npe=" << with.p90NPE
        << " p95_npe=" << with.p95NPE << '\n';
    out << scope << " ratio_with_over_without"
        << " signed_rate=" << safe_ratio(with.signedRate, without.signedRate)
        << " signed_rate_error=" << ratio_error(with.signedRate, with.signedRateError,
                                                  without.signedRate, without.signedRateError)
        << " positive_rate=" << safe_ratio(with.positiveRate, without.positiveRate)
        << " mean_npe=" << safe_ratio(with.meanNPE, without.meanNPE)
        << " p90_npe=" << safe_ratio(with.p90NPE, without.p90NPE)
        << " p95_npe=" << safe_ratio(with.p95NPE, without.p95NPE) << '\n';
}

}  // namespace optical_gel_compare

void compare_optical_gel_histograms(
        const char *withoutFileName = "results/optical_gel_comparison/Co60_without_gel_wide_histograms.root",
        const char *withFileName = "results/optical_gel_comparison/Co60_with_gel_wide_histograms.root",
        const char *labelWithout = "without gel",
        const char *labelWith = "with gel",
        const char *outPrefix = "results/optical_gel_comparison/Co60_gel_comparison",
        double withoutElapsedSeconds = 3599.960126,
        double withElapsedSeconds = 3599.947047,
        double totalRoiMin = 50.0,
        double totalRoiMax = 1000.0,
        double channelRoiMin = 50.0,
        double channelRoiMax = 600.0,
        const char *sourceLabel = "Co-60") {
    using namespace optical_gel_compare;
    gStyle->SetOptStat(0);

    if (withoutElapsedSeconds <= 0.0 || withElapsedSeconds <= 0.0) {
        std::cerr << "Error: elapsed spans must be positive." << std::endl;
        return;
    }

    TFile *withoutFile = TFile::Open(withoutFileName);
    TFile *withFile = TFile::Open(withFileName);
    if (!withoutFile || withoutFile->IsZombie() || !withFile || withFile->IsZombie()) {
        std::cerr << "Error: could not open one or both histogram files." << std::endl;
        delete withoutFile;
        delete withFile;
        return;
    }

    TH1D *withoutTotal = load_hist(withoutFile, "hTotalSub", "without_total");
    TH1D *withTotal = load_hist(withFile, "hTotalSub", "with_total");
    if (!withoutTotal || !withTotal || !same_binning(withoutTotal, withTotal)) {
        std::cerr << "Error: total histograms are missing or have different binning." << std::endl;
        delete withoutTotal;
        delete withTotal;
        withoutFile->Close();
        withFile->Close();
        delete withoutFile;
        delete withFile;
        return;
    }

    std::vector<int> channels;
    for (int channel = 0; channel < 8; ++channel) {
        if (withoutFile->Get(Form("hSub_ch%d", channel)) && withFile->Get(Form("hSub_ch%d", channel))) {
            channels.push_back(channel);
        }
    }

    const double fullMin = withoutTotal->GetXaxis()->GetXmin() - 0.5 * withoutTotal->GetBinWidth(1);
    const double fullMax = withoutTotal->GetXaxis()->GetXmax() + 0.5 * withoutTotal->GetBinWidth(withoutTotal->GetNbinsX());
    const Metrics fullWithout = compute_metrics(withoutTotal, withoutElapsedSeconds,
                                                fullMin, fullMax, "full_without");
    const Metrics fullWith = compute_metrics(withTotal, withElapsedSeconds,
                                             fullMin, fullMax, "full_with");
    const Metrics totalWithout = compute_metrics(withoutTotal, withoutElapsedSeconds,
                                                 totalRoiMin, totalRoiMax, "total_without");
    const Metrics totalWith = compute_metrics(withTotal, withElapsedSeconds,
                                              totalRoiMin, totalRoiMax, "total_with");

    save_total_overlay(withoutTotal, withTotal,
                       withoutElapsedSeconds, withElapsedSeconds,
                       sourceLabel,
                       labelWithout, labelWith, outPrefix);
    save_channel_overlay(withoutFile, withFile,
                         withoutElapsedSeconds, withElapsedSeconds,
                         sourceLabel,
                         labelWithout, labelWith, outPrefix, channels);

    std::ofstream report(Form("%s_summary.txt", outPrefix));
    if (!report.is_open()) {
        std::cerr << "Error: could not create summary file for prefix " << outPrefix << std::endl;
    } else {
        report << std::setprecision(9);
        report << "Optical-gel diagnostic spectrum comparison\n";
        report << "source_label " << sourceLabel << '\n';
        report << "ratio_direction with_over_without\n";
        report << "normalization elapsed_span_seconds (not dead-time corrected)\n";
        report << "without_file " << withoutFileName << '\n';
        report << "with_file " << withFileName << '\n';
        report << "without_elapsed_seconds " << withoutElapsedSeconds << '\n';
        report << "with_elapsed_seconds " << withElapsedSeconds << '\n';
        report << "diagnostic_note mean_and_quantiles_clip_negative_background_subtracted_bins_to_zero\n";
        report << "full_roi " << fullMin << ' ' << fullMax << '\n';
        write_metrics(report, "total_full", fullWithout, fullWith);
        report << "total_roi " << totalRoiMin << ' ' << totalRoiMax << '\n';
        write_metrics(report, "total_roi", totalWithout, totalWith);

        for (int channel : channels) {
            TH1D *withoutChannel = load_hist(withoutFile, Form("hSub_ch%d", channel),
                                             Form("without_summary_ch%d", channel));
            TH1D *withChannel = load_hist(withFile, Form("hSub_ch%d", channel),
                                          Form("with_summary_ch%d", channel));
            if (!withoutChannel || !withChannel || !same_binning(withoutChannel, withChannel)) {
                delete withoutChannel;
                delete withChannel;
                continue;
            }
            const Metrics channelWithout = compute_metrics(withoutChannel, withoutElapsedSeconds,
                                                           channelRoiMin, channelRoiMax,
                                                           Form("ch%d_without", channel));
            const Metrics channelWith = compute_metrics(withChannel, withElapsedSeconds,
                                                        channelRoiMin, channelRoiMax,
                                                        Form("ch%d_with", channel));
            report << "channel_roi " << channel << ' ' << channelRoiMin << ' ' << channelRoiMax << '\n';
            write_metrics(report, Form("channel_%d_roi", channel), channelWithout, channelWith);
            delete withoutChannel;
            delete withChannel;
        }
        report.close();
    }

    std::cout << "OPTICAL_GEL_COMPARE_RESULT"
              << " edge_independent_mean_ratio=" << safe_ratio(totalWith.meanNPE, totalWithout.meanNPE)
              << " p90_ratio=" << safe_ratio(totalWith.p90NPE, totalWithout.p90NPE)
              << " p95_ratio=" << safe_ratio(totalWith.p95NPE, totalWithout.p95NPE)
              << " signed_rate_ratio=" << safe_ratio(totalWith.signedRate, totalWithout.signedRate)
              << std::endl;
    std::cout << "Saved: " << outPrefix << "_summary.txt, "
              << outPrefix << "_total_overlay.png, "
              << outPrefix << "_channels_overlay.png" << std::endl;

    delete withoutTotal;
    delete withTotal;
    withoutFile->Close();
    withFile->Close();
    delete withoutFile;
    delete withFile;
}
