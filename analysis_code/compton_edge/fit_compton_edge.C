#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

double compton_edge_erfc(double *x, double *p) {
    const double sigma = std::max(std::abs(p[3]), 1.0e-9);
    return p[0] + p[1] * x[0] + 0.5 * p[2] * TMath::Erfc((x[0] - p[4]) / (TMath::Sqrt2() * sigma));
}

double compton_edge_erfc_gaus(double *x, double *p) {
    const double base = compton_edge_erfc(x, p);
    const double gausSigma = std::max(std::abs(p[7]), 1.0e-9);
    const double arg = (x[0] - p[6]) / gausSigma;
    return base + p[5] * std::exp(-0.5 * arg * arg);
}

double gaussian_component(double *x, double *p) {
    const double gausSigma = std::max(std::abs(p[2]), 1.0e-9);
    const double arg = (x[0] - p[1]) / gausSigma;
    return p[0] * std::exp(-0.5 * arg * arg);
}

TH1D *find_total_histogram(TFile *file, const char *histName) {
    if (!file) return nullptr;
    TH1D *hist = nullptr;
    if (histName && TString(histName).Length() > 0) {
        hist = dynamic_cast<TH1D*>(file->Get(histName));
        if (hist) return hist;
    }
    hist = dynamic_cast<TH1D*>(file->Get("hTotalSub"));
    if (hist) return hist;
    return dynamic_cast<TH1D*>(file->Get("hTotalSourceOnly"));
}

double estimate_edge_x(TH1D *hist, double xmin, double xmax) {
    if (!hist) return 0.5 * (xmin + xmax);

    int firstBin = std::max(2, hist->FindBin(xmin));
    int lastBin = std::min(hist->GetNbinsX() - 1, hist->FindBin(xmax));
    double bestSlope = 0.0;
    double bestX = 0.5 * (xmin + xmax);

    for (int bin = firstBin; bin <= lastBin; ++bin) {
        const double dx = hist->GetBinCenter(bin + 1) - hist->GetBinCenter(bin - 1);
        if (dx <= 0) continue;
        const double slope = (hist->GetBinContent(bin + 1) - hist->GetBinContent(bin - 1)) / dx;
        if (slope < bestSlope) {
            bestSlope = slope;
            bestX = hist->GetBinCenter(bin);
        }
    }
    return bestX;
}

bool use_gaussian_model(const char *modelName) {
    TString model(modelName ? modelName : "erfc_linear");
    model.ToLower();
    return model == "erfc_gaus" || model == "erfc_gaussian" || model == "erfc+gaus" || model == "erfc+gaussian";
}

struct EdgeFitMetrics {
    double xmin = 0.0;
    double xmax = 0.0;
    int status = -1;
    double edge = 0.0;
    double edgeErr = 0.0;
    double sigma = 0.0;
    double chi2ndf = 0.0;
};

EdgeFitMetrics fit_edge_metrics(TH1D *hist, double fitXMin, double fitXMax, const char *modelName, int index) {
    EdgeFitMetrics metrics;
    metrics.xmin = fitXMin;
    metrics.xmax = fitXMax;
    if (!hist || fitXMax <= fitXMin) return metrics;

    const int leftBin = hist->FindBin(fitXMin);
    const int rightBin = hist->FindBin(fitXMax);
    double yLeft = hist->GetBinContent(leftBin);
    const double yRight = hist->GetBinContent(rightBin);
    const double yMax = hist->GetMaximum();
    if (yLeft == 0 && yMax > 0) yLeft = yMax;
    const double amplitude = std::max(yLeft - yRight, 1.0);
    const double slopeEdge = estimate_edge_x(hist, fitXMin, fitXMax);
    const double sigmaGuess = std::max((fitXMax - fitXMin) / 20.0, hist->GetBinWidth(1));
    const bool useGaussian = use_gaussian_model(modelName);

    TF1 *fit = useGaussian
        ? new TF1(Form("fScanComptonEdgeErfcGaussian_%d", index), compton_edge_erfc_gaus, fitXMin, fitXMax, 8)
        : new TF1(Form("fScanComptonEdgeErfc_%d", index), compton_edge_erfc, fitXMin, fitXMax, 5);
    fit->SetParameters(yRight, 0.0, amplitude, sigmaGuess, slopeEdge);
    fit->SetParLimits(2, 0.0, std::max(yMax * 20.0, amplitude * 20.0));
    fit->SetParLimits(3, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    fit->SetParLimits(4, fitXMin, fitXMax);
    if (useGaussian) {
        fit->SetParameter(5, 0.2 * amplitude);
        fit->SetParameter(6, slopeEdge);
        fit->SetParameter(7, std::max((fitXMax - fitXMin) / 10.0, hist->GetBinWidth(1)));
        fit->SetParLimits(5, -std::max(yMax * 20.0, amplitude * 20.0), std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(6, fitXMin, fitXMax);
        fit->SetParLimits(7, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    }

    auto result = hist->Fit(fit, "RQ0SN");
    metrics.status = static_cast<int>(result);
    metrics.edge = fit->GetParameter(4);
    metrics.edgeErr = fit->GetParError(4);
    metrics.sigma = fit->GetParameter(3);
    const int ndf = fit->GetNDF();
    metrics.chi2ndf = ndf > 0 ? fit->GetChisquare() / ndf : 0.0;
    delete fit;
    return metrics;
}

void fit_compton_edge(const char *histFile = "compton_edge_histograms.root",
                      double fitXMin = 0.0,
                      double fitXMax = -1.0,
                      const char *outPrefix = "compton_edge_fit",
                      const char *histName = "",
                      const char *modelName = "erfc_linear") {
    TFile *input = TFile::Open(histFile);
    if (!input || input->IsZombie()) {
        std::cerr << "Error opening histogram ROOT file: " << histFile << std::endl;
        return;
    }

    TH1D *sourceHist = find_total_histogram(input, histName);
    if (!sourceHist) {
        std::cerr << "Error: cannot find hTotalSub or hTotalSourceOnly in " << histFile << std::endl;
        return;
    }

    TH1D *hist = dynamic_cast<TH1D*>(sourceHist->Clone("hComptonEdgeFitInput"));
    hist->SetDirectory(nullptr);
    input->Close();

    if (fitXMax <= fitXMin) {
        fitXMin = hist->GetXaxis()->GetXmin();
        fitXMax = hist->GetXaxis()->GetXmax();
    }
    if (fitXMax <= fitXMin) {
        std::cerr << "Error: invalid fit range." << std::endl;
        return;
    }

    const int leftBin = hist->FindBin(fitXMin);
    const int rightBin = hist->FindBin(fitXMax);
    double yLeft = hist->GetBinContent(leftBin);
    double yRight = hist->GetBinContent(rightBin);
    const double yMax = hist->GetMaximum();
    if (yLeft == 0 && yMax > 0) yLeft = yMax;
    const double amplitude = std::max(yLeft - yRight, 1.0);
    const double slopeEdge = estimate_edge_x(hist, fitXMin, fitXMax);
    const double sigmaGuess = std::max((fitXMax - fitXMin) / 20.0, hist->GetBinWidth(1));

    TString model(modelName ? modelName : "erfc_linear");
    model.ToLower();
    const bool useGaussian = model == "erfc_gaus" || model == "erfc_gaussian" || model == "erfc+gaus" || model == "erfc+gaussian";
    if (!useGaussian && model != "erfc_linear") {
        std::cerr << "Warning: unknown model '" << modelName << "', using erfc_linear." << std::endl;
        model = "erfc_linear";
    }

    TF1 *fit = useGaussian
        ? new TF1("fComptonEdgeErfcGaussian", compton_edge_erfc_gaus, fitXMin, fitXMax, 8)
        : new TF1("fComptonEdgeErfc", compton_edge_erfc, fitXMin, fitXMax, 5);
    if (useGaussian) {
        fit->SetParNames("offset", "slope", "amplitude", "sigma", "edge", "gaus_amp", "gaus_mean", "gaus_sigma");
    } else {
        fit->SetParNames("offset", "slope", "amplitude", "sigma", "edge");
    }
    fit->SetParameters(yRight, 0.0, amplitude, sigmaGuess, slopeEdge);
    fit->SetParLimits(2, 0.0, std::max(yMax * 20.0, amplitude * 20.0));
    fit->SetParLimits(3, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    fit->SetParLimits(4, fitXMin, fitXMax);
    if (useGaussian) {
        fit->SetParameter(5, 0.2 * amplitude);
        fit->SetParameter(6, slopeEdge);
        fit->SetParameter(7, std::max((fitXMax - fitXMin) / 10.0, hist->GetBinWidth(1)));
        fit->SetParLimits(5, -std::max(yMax * 20.0, amplitude * 20.0), std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(6, fitXMin, fitXMax);
        fit->SetParLimits(7, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    }

    auto result = hist->Fit(fit, "RQ0S");
    const int status = static_cast<int>(result);

    const double edge = fit->GetParameter(4);
    const double edgeErr = fit->GetParError(4);
    const double sigma = fit->GetParameter(3);
    const double sigmaErr = fit->GetParError(3);
    const double amplitudeFit = fit->GetParameter(2);
    const double amplitudeErr = fit->GetParError(2);
    const double gausAmp = useGaussian ? fit->GetParameter(5) : 0.0;
    const double gausAmpErr = useGaussian ? fit->GetParError(5) : 0.0;
    const double gausMean = useGaussian ? fit->GetParameter(6) : 0.0;
    const double gausMeanErr = useGaussian ? fit->GetParError(6) : 0.0;
    const double gausSigma = useGaussian ? fit->GetParameter(7) : 0.0;
    const double gausSigmaErr = useGaussian ? fit->GetParError(7) : 0.0;
    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    const double chi2ndf = ndf > 0 ? chi2 / ndf : 0.0;
    const double fitWidth = fitXMax - fitXMin;
    std::vector<TString> warnings;
    if (status != 0) warnings.push_back(Form("fit_status_%d", status));
    if (edge - fitXMin < 0.05 * fitWidth || fitXMax - edge < 0.05 * fitWidth) {
        warnings.push_back("edge_near_fit_boundary");
    }
    if (edgeErr <= 0.0 || edgeErr > 0.5 * fitWidth) {
        warnings.push_back("edge_error_large_or_invalid");
    }
    if (sigma <= hist->GetBinWidth(1) * 0.25 || sigma > 0.5 * fitWidth) {
        warnings.push_back("sigma_at_limit_or_too_large");
    }
    if (ndf <= 0) {
        warnings.push_back("nonpositive_ndf");
    } else if (chi2ndf > 5.0) {
        warnings.push_back("large_chi2_ndf");
    }
    if (useGaussian && (gausSigma <= hist->GetBinWidth(1) * 0.25 || gausSigma > fitWidth)) {
        warnings.push_back("gaussian_sigma_at_limit_or_too_large");
    }

    TString prefix(outPrefix);
    TString pngName = Form("%s_compton_edge_fit.png", prefix.Data());
    TString pdfName = Form("%s_compton_edge_fit.pdf", prefix.Data());
    TString rootName = Form("%s_compton_edge_fit.root", prefix.Data());
    TString txtName = Form("%s_compton_edge_fit.txt", prefix.Data());

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    TCanvas *canvas = new TCanvas("c_compton_edge_fit", "Compton Edge Fit", 1200, 820);
    TPad *topPad = new TPad("fit_top_pad", "fit_top_pad", 0.0, 0.30, 0.76, 1.0);
    TPad *infoPad = new TPad("fit_info_pad", "fit_info_pad", 0.76, 0.30, 1.0, 1.0);
    TPad *pullPad = new TPad("fit_pull_pad", "fit_pull_pad", 0.0, 0.0, 0.76, 0.30);
    topPad->SetBottomMargin(0.03);
    topPad->SetLeftMargin(0.13);
    topPad->SetRightMargin(0.03);
    topPad->SetTicks(1, 1);
    infoPad->SetLeftMargin(0.02);
    infoPad->SetRightMargin(0.06);
    infoPad->SetTopMargin(0.08);
    infoPad->SetBottomMargin(0.06);
    infoPad->SetFrameBorderMode(0);
    pullPad->SetTopMargin(0.04);
    pullPad->SetBottomMargin(0.30);
    pullPad->SetLeftMargin(0.13);
    pullPad->SetRightMargin(0.03);
    pullPad->SetTicks(1, 1);
    topPad->Draw();
    infoPad->Draw();
    pullPad->Draw();

    topPad->cd();
    hist->SetTitle(";NPE;Counts");
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.65);
    hist->SetLineWidth(2);
    hist->SetMaximum(hist->GetMaximum() * 1.18);
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.048);
    hist->GetYaxis()->SetTitleOffset(1.05);
    hist->Draw("E");

    TF1 *erfcComponent = nullptr;
    TF1 *gausComponent = nullptr;
    TLine *gausMeanLine = nullptr;
    if (useGaussian) {
        erfcComponent = new TF1("fComptonEdgeErfcComponent", compton_edge_erfc, fitXMin, fitXMax, 5);
        for (int i = 0; i < 5; ++i) erfcComponent->SetParameter(i, fit->GetParameter(i));
        erfcComponent->SetLineColor(kOrange + 7);
        erfcComponent->SetLineStyle(7);
        erfcComponent->SetLineWidth(2);

        gausComponent = new TF1("fComptonEdgeGaussianComponent", gaussian_component, fitXMin, fitXMax, 3);
        gausComponent->SetParameters(gausAmp, gausMean, gausSigma);
        gausComponent->SetLineColor(kMagenta + 1);
        gausComponent->SetLineStyle(5);
        gausComponent->SetLineWidth(2);
    }

    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(3);
    fit->Draw("SAME");
    if (useGaussian) {
        erfcComponent->Draw("SAME");
        gausComponent->Draw("SAME");
    }

    TLine *edgeLine = new TLine(edge, 0.0, edge, hist->GetMaximum() * 1.05);
    edgeLine->SetLineColor(kBlue + 1);
    edgeLine->SetLineStyle(2);
    edgeLine->SetLineWidth(2);
    edgeLine->Draw("SAME");

    TLine *slopeLine = new TLine(slopeEdge, 0.0, slopeEdge, hist->GetMaximum() * 1.05);
    slopeLine->SetLineColor(kGreen + 2);
    slopeLine->SetLineStyle(3);
    slopeLine->SetLineWidth(2);
    slopeLine->Draw("SAME");

    if (useGaussian) {
        gausMeanLine = new TLine(gausMean, hist->GetMinimum(), gausMean, hist->GetMaximum() * 1.02);
        gausMeanLine->SetLineColor(kMagenta + 1);
        gausMeanLine->SetLineStyle(9);
        gausMeanLine->SetLineWidth(2);
        gausMeanLine->Draw("SAME");
    }

    infoPad->cd();
    TLegend *legend = new TLegend(0.05, 0.55, 0.95, 0.95);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(useGaussian ? 0.043 : 0.060);
    legend->AddEntry(hist, "Histogram", "lep");
    legend->AddEntry(fit, useGaussian ? "linear + erfc edge + Gaussian" : "linear + erfc edge", "l");
    if (useGaussian) {
        legend->AddEntry(erfcComponent, "linear + erfc component", "l");
        legend->AddEntry(gausComponent, "Gaussian contribution", "l");
        legend->AddEntry(gausMeanLine, "Gaussian mean", "l");
    }
    legend->AddEntry(edgeLine, "edge", "l");
    legend->AddEntry(slopeLine, "max falling slope", "l");
    legend->Draw();

    TH1D *pullHist = dynamic_cast<TH1D*>(hist->Clone("hComptonEdgeFitPull"));
    pullHist->Reset("ICES");
    pullHist->SetDirectory(nullptr);
    pullHist->SetTitle(";NPE;(data-fit)/error");
    pullHist->SetLineColor(kBlack);
    pullHist->SetMarkerStyle(20);
    pullHist->SetMarkerSize(0.6);

    int pullCount = 0;
    double pullSum2 = 0.0;
    const int firstFitBin = hist->FindBin(fitXMin);
    const int lastFitBin = hist->FindBin(fitXMax);
    for (int bin = firstFitBin; bin <= lastFitBin; ++bin) {
        const double x = hist->GetBinCenter(bin);
        if (x < fitXMin || x > fitXMax) continue;
        const double y = hist->GetBinContent(bin);
        double err = hist->GetBinError(bin);
        if (err <= 0.0) err = std::sqrt(std::max(std::abs(y), 1.0));
        const double pull = (y - fit->Eval(x)) / err;
        pullHist->SetBinContent(bin, pull);
        pullHist->SetBinError(bin, 1.0);
        pullSum2 += pull * pull;
        ++pullCount;
    }
    const double pullRms = pullCount > 0 ? std::sqrt(pullSum2 / pullCount) : 0.0;
    if (pullRms > 3.0) warnings.push_back("large_pull_rms");

    infoPad->cd();
    TPaveText *text = new TPaveText(0.05, 0.08, 0.95, 0.50, "NDC");
    text->SetFillColor(0);
    text->SetFillStyle(0);
    text->SetBorderSize(0);
    text->SetTextFont(42);
    text->SetTextSize(0.055);
    text->SetTextAlign(12);
    text->AddText(Form("Edge = %.6g #pm %.3g NPE", edge, edgeErr));
    text->AddText(Form("Slope seed = %.6g NPE", slopeEdge));
    text->AddText(Form("#sigma = %.6g #pm %.3g NPE", sigma, sigmaErr));
    if (useGaussian) text->AddText(Form("Gaus #mu = %.6g #pm %.3g NPE", gausMean, gausMeanErr));
    text->AddText(Form("#chi^{2}/ndf = %.3g", chi2ndf));
    text->AddText(Form("Pull RMS = %.3g", pullRms));
    if (!warnings.empty()) text->AddText(Form("Warnings: %zu", warnings.size()));
    text->AddText(Form("Range: %.6g to %.6g NPE", fitXMin, fitXMax));
    text->Draw();

    pullPad->cd();
    pullHist->GetXaxis()->SetTitleSize(0.11);
    pullHist->GetXaxis()->SetLabelSize(0.09);
    pullHist->GetYaxis()->SetTitleSize(0.10);
    pullHist->GetYaxis()->SetLabelSize(0.08);
    pullHist->GetYaxis()->SetTitleOffset(0.45);
    pullHist->SetMinimum(-5.0);
    pullHist->SetMaximum(5.0);
    pullHist->Draw("E");
    TLine *zeroPull = new TLine(fitXMin, 0.0, fitXMax, 0.0);
    zeroPull->SetLineStyle(2);
    zeroPull->Draw("SAME");
    TLine *plusThree = new TLine(fitXMin, 3.0, fitXMax, 3.0);
    plusThree->SetLineStyle(3);
    plusThree->SetLineColor(kGray + 2);
    plusThree->Draw("SAME");
    TLine *minusThree = new TLine(fitXMin, -3.0, fitXMax, -3.0);
    minusThree->SetLineStyle(3);
    minusThree->SetLineColor(kGray + 2);
    minusThree->Draw("SAME");

    canvas->cd();

    canvas->SaveAs(pngName);
    canvas->SaveAs(pdfName);

    TFile *output = TFile::Open(rootName, "RECREATE");
    if (output && !output->IsZombie()) {
        hist->Write("hComptonEdgeFitInput");
        pullHist->Write("hComptonEdgeFitPull");
        fit->Write("fComptonEdgeErfc");
        if (erfcComponent) erfcComponent->Write("fComptonEdgeErfcComponent");
        if (gausComponent) gausComponent->Write("fComptonEdgeGaussianComponent");
        canvas->Write("c_compton_edge_fit");
        output->Close();
    }

    std::ofstream txt(txtName.Data());
    txt << "model " << (useGaussian ? "erfc_gaussian" : "erfc_linear") << "\n";
    txt << "status " << status << "\n";
    txt << "hist_file " << histFile << "\n";
    txt << "hist_name " << (histName && TString(histName).Length() > 0 ? histName : hist->GetName()) << "\n";
    txt << "fit_xmin " << fitXMin << "\n";
    txt << "fit_xmax " << fitXMax << "\n";
    txt << "edge " << edge << "\n";
    txt << "edge_error " << edgeErr << "\n";
    txt << "slope_edge " << slopeEdge << "\n";
    txt << "sigma " << sigma << "\n";
    txt << "sigma_error " << sigmaErr << "\n";
    txt << "amplitude " << amplitudeFit << "\n";
    txt << "amplitude_error " << amplitudeErr << "\n";
    if (useGaussian) {
        txt << "gaus_amp " << gausAmp << "\n";
        txt << "gaus_amp_error " << gausAmpErr << "\n";
        txt << "gaus_mean " << gausMean << "\n";
        txt << "gaus_mean_error " << gausMeanErr << "\n";
        txt << "gaus_sigma " << gausSigma << "\n";
        txt << "gaus_sigma_error " << gausSigmaErr << "\n";
    }
    txt << "chi2 " << chi2 << "\n";
    txt << "ndf " << ndf << "\n";
    txt << "chi2_ndf " << chi2ndf << "\n";
    txt << "pull_rms " << pullRms << "\n";
    txt << "warnings";
    if (warnings.empty()) {
        txt << " none";
    } else {
        for (const TString &warning : warnings) txt << " " << warning;
    }
    txt << "\n";
    txt << "png " << pngName << "\n";
    txt << "pdf " << pdfName << "\n";
    txt << "root " << rootName << "\n";
    txt.close();

    std::cout << "\n--- Compton Edge Fit Complete ---" << std::endl;
    std::cout << "Model: " << (useGaussian ? "erfc_gaussian" : "erfc_linear") << std::endl;
    std::cout << "Histogram ROOT file: " << histFile << std::endl;
    std::cout << "Fit range [NPE]: [" << fitXMin << ", " << fitXMax << "]" << std::endl;
    std::cout << "Fit status: " << status << std::endl;
    std::cout << "Edge [NPE]: " << edge << " +/- " << edgeErr << std::endl;
    std::cout << "Max falling slope edge seed [NPE]: " << slopeEdge << std::endl;
    std::cout << "Sigma [NPE]: " << sigma << " +/- " << sigmaErr << std::endl;
    std::cout << "Amplitude: " << amplitudeFit << " +/- " << amplitudeErr << std::endl;
    if (useGaussian) {
        std::cout << "Gaussian amp: " << gausAmp << " +/- " << gausAmpErr << std::endl;
        std::cout << "Gaussian mean [NPE]: " << gausMean << " +/- " << gausMeanErr << std::endl;
        std::cout << "Gaussian sigma [NPE]: " << gausSigma << " +/- " << gausSigmaErr << std::endl;
    }
    std::cout << "Chi2/NDF: " << chi2 << "/" << ndf << " = " << chi2ndf << std::endl;
    std::cout << "Pull RMS: " << pullRms << std::endl;
    if (warnings.empty()) {
        std::cout << "Fit warnings: none" << std::endl;
    } else {
        std::cout << "Fit warnings:";
        for (const TString &warning : warnings) std::cout << " " << warning;
        std::cout << std::endl;
    }
    std::cout << "Fit plot: " << pngName << std::endl;
    std::cout << "Fit PDF: " << pdfName << std::endl;
    std::cout << "Fit ROOT file: " << rootName << std::endl;
    std::cout << "Fit summary: " << txtName << std::endl;
    std::cout << "FIT_RESULT model=" << (useGaussian ? "erfc_gaussian" : "erfc_linear")
              << " edge=" << edge
              << " edge_error=" << edgeErr
              << " slope_edge=" << slopeEdge
              << " sigma=" << sigma
              << " sigma_error=" << sigmaErr
              << " chi2_ndf=" << chi2ndf
              << " pull_rms=" << pullRms;
    if (useGaussian) {
        std::cout << " gaus_mean=" << gausMean
                  << " gaus_mean_error=" << gausMeanErr
                  << " gaus_sigma=" << gausSigma
                  << " gaus_sigma_error=" << gausSigmaErr;
    }
    std::cout << " status=" << status;
    if (!warnings.empty()) {
        std::cout << " warnings=";
        for (size_t i = 0; i < warnings.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << warnings[i];
        }
    } else {
        std::cout << " warnings=none";
    }
    std::cout << std::endl;
}

void scan_compton_edge_fit(const char *histFile = "compton_edge_histograms.root",
                           double fitXMin = 0.0,
                           double fitXMax = -1.0,
                           const char *outPrefix = "compton_edge_fit",
                           const char *histName = "",
                           const char *modelName = "erfc_linear",
                           double scanFraction = 0.10) {
    TFile *input = TFile::Open(histFile);
    if (!input || input->IsZombie()) {
        std::cerr << "Error opening histogram ROOT file: " << histFile << std::endl;
        return;
    }

    TH1D *sourceHist = find_total_histogram(input, histName);
    if (!sourceHist) {
        std::cerr << "Error: cannot find hTotalSub or hTotalSourceOnly in " << histFile << std::endl;
        return;
    }

    TH1D *hist = dynamic_cast<TH1D*>(sourceHist->Clone("hComptonEdgeScanInput"));
    hist->SetDirectory(nullptr);
    input->Close();

    if (fitXMax <= fitXMin) {
        fitXMin = hist->GetXaxis()->GetXmin();
        fitXMax = hist->GetXaxis()->GetXmax();
    }
    if (fitXMax <= fitXMin) {
        std::cerr << "Error: invalid fit range." << std::endl;
        return;
    }

    const double axisMin = hist->GetXaxis()->GetXmin();
    const double axisMax = hist->GetXaxis()->GetXmax();
    const double baseWidth = fitXMax - fitXMin;
    const double delta = std::max(0.0, scanFraction) * baseWidth;
    const std::vector<double> shifts = {-delta, 0.0, delta};
    std::vector<EdgeFitMetrics> scans;

    int scanIndex = 0;
    for (double lowShift : shifts) {
        for (double highShift : shifts) {
            double xmin = std::max(axisMin, fitXMin + lowShift);
            double xmax = std::min(axisMax, fitXMax + highShift);
            if (xmax <= xmin || xmax - xmin < 0.25 * baseWidth) continue;
            scans.push_back(fit_edge_metrics(hist, xmin, xmax, modelName, scanIndex++));
        }
    }

    std::vector<EdgeFitMetrics> goodScans;
    for (const EdgeFitMetrics &scan : scans) {
        if (scan.status == 0 && std::isfinite(scan.edge)) goodScans.push_back(scan);
    }

    double edgeMean = 0.0;
    double edgeRms = 0.0;
    double edgeMin = 0.0;
    double edgeMax = 0.0;
    if (!goodScans.empty()) {
        edgeMin = goodScans.front().edge;
        edgeMax = goodScans.front().edge;
        for (const EdgeFitMetrics &scan : goodScans) {
            edgeMean += scan.edge;
            edgeMin = std::min(edgeMin, scan.edge);
            edgeMax = std::max(edgeMax, scan.edge);
        }
        edgeMean /= goodScans.size();
        for (const EdgeFitMetrics &scan : goodScans) {
            const double diff = scan.edge - edgeMean;
            edgeRms += diff * diff;
        }
        edgeRms = std::sqrt(edgeRms / goodScans.size());
    }

    TString prefix(outPrefix);
    TString pngName = Form("%s_compton_edge_scan.png", prefix.Data());
    TString pdfName = Form("%s_compton_edge_scan.pdf", prefix.Data());
    TString rootName = Form("%s_compton_edge_scan.root", prefix.Data());
    TString txtName = Form("%s_compton_edge_scan.txt", prefix.Data());

    TGraphErrors *graph = new TGraphErrors(scans.size());
    graph->SetName("gComptonEdgeRangeScan");
    graph->SetTitle("Compton edge range scan;Scan index;Fitted edge [NPE]");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->SetLineColor(kBlue + 1);
    graph->SetMarkerColor(kBlue + 1);
    for (size_t i = 0; i < scans.size(); ++i) {
        graph->SetPoint(i, static_cast<double>(i + 1), scans[i].edge);
        graph->SetPointError(i, 0.0, std::max(scans[i].edgeErr, 0.0));
    }

    TCanvas *canvas = new TCanvas("c_compton_edge_scan", "Compton Edge Range Scan", 900, 650);
    graph->Draw("AP");
    if (!goodScans.empty()) {
        TLine *meanLine = new TLine(0.5, edgeMean, scans.size() + 0.5, edgeMean);
        meanLine->SetLineColor(kRed + 1);
        meanLine->SetLineStyle(2);
        meanLine->SetLineWidth(2);
        meanLine->Draw("SAME");

        TPaveText *text = new TPaveText(0.14, 0.72, 0.52, 0.88, "NDC");
        text->SetFillColor(0);
        text->SetBorderSize(1);
        text->AddText(Form("Mean edge = %.6g NPE", edgeMean));
        text->AddText(Form("Range RMS = %.3g NPE", edgeRms));
        text->AddText(Form("Min/Max = %.6g / %.6g", edgeMin, edgeMax));
        text->AddText(Form("Good fits = %zu / %zu", goodScans.size(), scans.size()));
        text->Draw();
    }
    canvas->SaveAs(pngName);
    canvas->SaveAs(pdfName);

    TFile *output = TFile::Open(rootName, "RECREATE");
    if (output && !output->IsZombie()) {
        hist->Write("hComptonEdgeScanInput");
        graph->Write("gComptonEdgeRangeScan");
        canvas->Write("c_compton_edge_scan");
        output->Close();
    }

    std::ofstream txt(txtName.Data());
    txt << "model " << (use_gaussian_model(modelName) ? "erfc_gaussian" : "erfc_linear") << "\n";
    txt << "hist_file " << histFile << "\n";
    txt << "hist_name " << (histName && TString(histName).Length() > 0 ? histName : hist->GetName()) << "\n";
    txt << "base_fit_xmin " << fitXMin << "\n";
    txt << "base_fit_xmax " << fitXMax << "\n";
    txt << "scan_fraction " << scanFraction << "\n";
    txt << "n_good " << goodScans.size() << "\n";
    txt << "n_total " << scans.size() << "\n";
    txt << "edge_mean " << edgeMean << "\n";
    txt << "edge_rms " << edgeRms << "\n";
    txt << "edge_min " << edgeMin << "\n";
    txt << "edge_max " << edgeMax << "\n";
    txt << "index xmin xmax status edge edge_error sigma chi2_ndf\n";
    for (size_t i = 0; i < scans.size(); ++i) {
        txt << (i + 1) << " "
            << scans[i].xmin << " "
            << scans[i].xmax << " "
            << scans[i].status << " "
            << scans[i].edge << " "
            << scans[i].edgeErr << " "
            << scans[i].sigma << " "
            << scans[i].chi2ndf << "\n";
    }
    txt << "png " << pngName << "\n";
    txt << "pdf " << pdfName << "\n";
    txt << "root " << rootName << "\n";
    txt.close();

    std::cout << "\n--- Compton Edge Range Scan Complete ---" << std::endl;
    std::cout << "Histogram ROOT file: " << histFile << std::endl;
    std::cout << "Base fit range [NPE]: [" << fitXMin << ", " << fitXMax << "]" << std::endl;
    std::cout << "Scan fraction: " << scanFraction << std::endl;
    std::cout << "Good fits: " << goodScans.size() << " / " << scans.size() << std::endl;
    std::cout << "Edge mean [NPE]: " << edgeMean << std::endl;
    std::cout << "Edge RMS from range scan [NPE]: " << edgeRms << std::endl;
    std::cout << "Edge min/max [NPE]: " << edgeMin << " / " << edgeMax << std::endl;
    std::cout << "Range scan plot: " << pngName << std::endl;
    std::cout << "Range scan PDF: " << pdfName << std::endl;
    std::cout << "Range scan ROOT file: " << rootName << std::endl;
    std::cout << "Range scan summary: " << txtName << std::endl;
    std::cout << "SCAN_RESULT model=" << (use_gaussian_model(modelName) ? "erfc_gaussian" : "erfc_linear")
              << " edge_mean=" << edgeMean
              << " edge_rms=" << edgeRms
              << " edge_min=" << edgeMin
              << " edge_max=" << edgeMax
              << " n_good=" << goodScans.size()
              << " n_total=" << scans.size();
    if (goodScans.size() < 3) {
        std::cout << " warnings=too_few_good_fits";
    } else {
        std::cout << " warnings=none";
    }
    std::cout << std::endl;
}
