#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace source_edge_compare {

struct EdgeResult {
    double edge = 0.0;
    double error = 0.0;
    int status = -1;
};

bool load_edge_result(const char *path, EdgeResult &result) {
    std::ifstream input(path);
    if (!input.is_open()) {
        std::cerr << "Error: cannot open fit summary: " << path << std::endl;
        return false;
    }
    std::string key;
    std::string value;
    while (input >> key >> value) {
        try {
            if (key == "edge") result.edge = std::stod(value);
            else if (key == "edge_error") result.error = std::stod(value);
            else if (key == "status") result.status = std::stoi(value);
        } catch (...) {
        }
        std::string remainder;
        std::getline(input, remainder);
    }
    if (result.status != 0 || result.edge <= 0.0 || result.error <= 0.0) {
        std::cerr << "Error: invalid fit result in " << path
                  << " (status=" << result.status
                  << ", edge=" << result.edge
                  << ", error=" << result.error << ")" << std::endl;
        return false;
    }
    return true;
}

double ratio_error(double numerator, double numeratorError,
                   double denominator, double denominatorError) {
    if (numerator == 0.0 || denominator == 0.0) return 0.0;
    const double value = numerator / denominator;
    return std::abs(value) * std::sqrt(std::pow(numeratorError / numerator, 2) +
                                       std::pow(denominatorError / denominator, 2));
}

}  // namespace source_edge_compare

void compare_source_edge_linearity(
        const char *csWithoutSummary = "results/optical_gel_comparison/Cs137_without_gel_total_40_180_compton_edge_fit.txt",
        const char *csWithSummary = "results/optical_gel_comparison/Cs137_with_gel_total_85_385_compton_edge_fit.txt",
        const char *coWithoutSummary = "results/optical_gel_comparison/Co60_without_gel_total_100_400_compton_edge_fit.txt",
        const char *coWithSummary = "results/optical_gel_comparison/Co60_with_gel_total_250_1000_compton_edge_fit.txt",
        const char *outPrefix = "results/optical_gel_comparison/Cs137_Co60_edge_linearity") {
    using namespace source_edge_compare;
    gStyle->SetOptStat(0);

    EdgeResult csWithout;
    EdgeResult csWith;
    EdgeResult coWithout;
    EdgeResult coWith;
    if (!load_edge_result(csWithoutSummary, csWithout) ||
        !load_edge_result(csWithSummary, csWith) ||
        !load_edge_result(coWithoutSummary, coWithout) ||
        !load_edge_result(coWithSummary, coWith)) {
        return;
    }

    constexpr double csEdgeKeV = 477.334018;
    constexpr double coUpperEdgeKeV = 1118.101047;
    const double energyRatio = csEdgeKeV / coUpperEdgeKeV;

    const double csGain = csWith.edge / csWithout.edge;
    const double csGainError = ratio_error(csWith.edge, csWith.error,
                                          csWithout.edge, csWithout.error);
    const double coGain = coWith.edge / coWithout.edge;
    const double coGainError = ratio_error(coWith.edge, coWith.error,
                                          coWithout.edge, coWithout.error);
    const double doubleRatio = csGain / coGain;
    const double doubleRatioError = ratio_error(csGain, csGainError, coGain, coGainError);

    const double csWithoutYield = csWithout.edge / csEdgeKeV;
    const double csWithYield = csWith.edge / csEdgeKeV;
    const double coWithoutYield = coWithout.edge / coUpperEdgeKeV;
    const double coWithYield = coWith.edge / coUpperEdgeKeV;
    const double withoutLinearity = csWithoutYield / coWithoutYield;
    const double withLinearity = csWithYield / coWithYield;
    const double predictedCsWith = coWith.edge * energyRatio;
    const double predictedDifferencePercent = 100.0 * (csWith.edge / predictedCsWith - 1.0);

    double energy[2] = {csEdgeKeV, coUpperEdgeKeV};
    double energyError[2] = {0.0, 0.0};
    double withoutEdge[2] = {csWithout.edge, coWithout.edge};
    double withoutError[2] = {csWithout.error, coWithout.error};
    double withEdge[2] = {csWith.edge, coWith.edge};
    double withError[2] = {csWith.error, coWith.error};

    TGraphErrors graphWithout(2, energy, withoutEdge, energyError, withoutError);
    TGraphErrors graphWith(2, energy, withEdge, energyError, withError);
    graphWithout.SetName("gEdgeWithoutGel");
    graphWith.SetName("gEdgeWithGel");
    graphWithout.SetMarkerStyle(20);
    graphWithout.SetMarkerSize(1.3);
    graphWithout.SetMarkerColor(kBlue + 1);
    graphWithout.SetLineColor(kBlue + 1);
    graphWith.SetMarkerStyle(21);
    graphWith.SetMarkerSize(1.3);
    graphWith.SetMarkerColor(kRed + 1);
    graphWith.SetLineColor(kRed + 1);

    const double withoutGuideSlope = 0.5 * (csWithoutYield + coWithoutYield);
    const double withGuideSlope = 0.5 * (csWithYield + coWithYield);
    TLine withoutGuide(0.0, 0.0, 1200.0, 1200.0 * withoutGuideSlope);
    TLine withGuide(0.0, 0.0, 1200.0, 1200.0 * withGuideSlope);
    withoutGuide.SetLineColor(kBlue + 1);
    withoutGuide.SetLineStyle(2);
    withoutGuide.SetLineWidth(2);
    withGuide.SetLineColor(kRed + 1);
    withGuide.SetLineStyle(2);
    withGuide.SetLineWidth(2);

    TCanvas canvas("cSourceEdgeLinearity", "Cs-137 and Co-60 edge comparison", 950, 700);
    graphWithout.SetTitle("Cs-137 vs Co-60 effective Compton edge;Physical Compton-edge energy [keV];Fitted effective edge [NPE]");
    graphWithout.GetXaxis()->SetLimits(0.0, 1200.0);
    graphWithout.SetMinimum(0.0);
    graphWithout.SetMaximum(560.0);
    graphWithout.Draw("AP");
    graphWith.Draw("P SAME");
    withoutGuide.Draw("SAME");
    withGuide.Draw("SAME");

    TLegend legend(0.16, 0.70, 0.45, 0.88);
    legend.AddEntry(&graphWithout, "without gel", "p");
    legend.AddEntry(&graphWith, "with gel", "p");
    legend.AddEntry(&withoutGuide, "through-origin guide", "l");
    legend.Draw();

    TPaveText note(0.50, 0.16, 0.88, 0.39, "NDC");
    note.SetFillColor(0);
    note.SetBorderSize(1);
    note.AddText(Form("gel gain: Cs-137 = %.3f, Co-60 = %.3f", csGain, coGain));
    note.AddText(Form("double ratio = %.4f #pm %.4f (formal)", doubleRatio, doubleRatioError));
    note.AddText(Form("NPE/keV agreement: %.2f%% (no gel), %.2f%% (gel)",
                      100.0 * (withoutLinearity - 1.0),
                      100.0 * (withLinearity - 1.0)));
    note.AddText("Co-60 point uses the 1332-keV upper edge");
    note.Draw();

    canvas.SaveAs(Form("%s.png", outPrefix));
    canvas.SaveAs(Form("%s.pdf", outPrefix));

    TFile output(Form("%s.root", outPrefix), "RECREATE");
    if (!output.IsZombie()) {
        graphWithout.Write();
        graphWith.Write();
        canvas.Write("cSourceEdgeLinearity");
        output.Close();
    }

    std::ofstream report(Form("%s.txt", outPrefix));
    if (report.is_open()) {
        report << std::setprecision(10);
        report << "cs_edge_energy_kev " << csEdgeKeV << '\n';
        report << "co_upper_edge_energy_kev " << coUpperEdgeKeV << '\n';
        report << "physical_energy_ratio_cs_over_co " << energyRatio << '\n';
        report << "cs_without_edge " << csWithout.edge << " error " << csWithout.error << '\n';
        report << "cs_with_edge " << csWith.edge << " error " << csWith.error << '\n';
        report << "co_without_edge " << coWithout.edge << " error " << coWithout.error << '\n';
        report << "co_with_edge " << coWith.edge << " error " << coWith.error << '\n';
        report << "cs_gel_gain " << csGain << " formal_error " << csGainError << '\n';
        report << "co_gel_gain " << coGain << " formal_error " << coGainError << '\n';
        report << "gain_double_ratio_cs_over_co " << doubleRatio
               << " formal_error " << doubleRatioError << '\n';
        report << "cs_without_npe_per_kev " << csWithoutYield << '\n';
        report << "co_without_npe_per_kev " << coWithoutYield << '\n';
        report << "without_linearity_ratio_cs_over_co " << withoutLinearity << '\n';
        report << "cs_with_npe_per_kev " << csWithYield << '\n';
        report << "co_with_npe_per_kev " << coWithYield << '\n';
        report << "with_linearity_ratio_cs_over_co " << withLinearity << '\n';
        report << "cs_with_predicted_from_co " << predictedCsWith << '\n';
        report << "cs_with_minus_prediction_percent " << predictedDifferencePercent << '\n';
        report << "note Co60 effective edge is mapped to the 1332.492-keV upper physical Compton edge\n";
        report << "note errors are formal fit errors and exclude fit-range, model, and run systematics\n";
        report.close();
    }

    std::cout << "SOURCE_EDGE_COMPARE_RESULT"
              << " cs_gain=" << csGain
              << " co_gain=" << coGain
              << " double_ratio=" << doubleRatio
              << " double_ratio_error=" << doubleRatioError
              << " no_gel_linearity=" << withoutLinearity
              << " with_gel_linearity=" << withLinearity
              << " cs_prediction_delta_percent=" << predictedDifferencePercent
              << std::endl;
}
