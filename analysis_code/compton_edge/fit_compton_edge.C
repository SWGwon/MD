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

struct GammaLine {
    TString label;
    double energyKeV;
};

std::vector<GammaLine> gResponseGammaLines;
constexpr int kDerivativeSmoothHalfWindow = 3;

double compton_edge_energy_kev(double gammaKeV);

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

double smeared_edge_response_model(double *x, double *p) {
    const double npePerKeV = std::max(p[3], 1.0e-9);
    const double sigmaNPE = std::max(std::abs(p[4]), 1.0e-9);
    double value = p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
    for (size_t i = 0; i < gResponseGammaLines.size(); ++i) {
        const double edgeNPE = npePerKeV * compton_edge_energy_kev(gResponseGammaLines[i].energyKeV);
        value += 0.5 * p[5 + i] * TMath::Erfc((x[0] - edgeNPE) / (TMath::Sqrt2() * sigmaNPE));
    }
    return value;
}

double gaussian_component(double *x, double *p) {
    const double gausSigma = std::max(std::abs(p[2]), 1.0e-9);
    const double arg = (x[0] - p[1]) / gausSigma;
    return p[0] * std::exp(-0.5 * arg * arg);
}

double linear_gaussian_model(double *x, double *p) {
    const double sigma = std::max(std::abs(p[4]), 1.0e-9);
    const double arg = (x[0] - p[3]) / sigma;
    return p[0] + p[1] * x[0] + p[2] * std::exp(-0.5 * arg * arg);
}

double compton_edge_energy_kev(double gammaKeV) {
    constexpr double electronRestKeV = 510.99895;
    return gammaKeV * (1.0 - 1.0 / (1.0 + 2.0 * gammaKeV / electronRestKeV));
}

double klein_nishina_recoil_shape(double electronKeV, double gammaKeV) {
    constexpr double electronRestKeV = 510.99895;
    const double edgeKeV = compton_edge_energy_kev(gammaKeV);
    if (electronKeV < 0.0 || electronKeV > edgeKeV || electronKeV >= gammaKeV) return 0.0;

    const double scatteredRatio = (gammaKeV - electronKeV) / gammaKeV;
    if (scatteredRatio <= 0.0) return 0.0;

    const double cosTheta = 1.0 - electronRestKeV * electronKeV / (gammaKeV * (gammaKeV - electronKeV));
    const double sin2Theta = std::max(0.0, 1.0 - cosTheta * cosTheta);
    return std::max(0.0, scatteredRatio + 1.0 / scatteredRatio - sin2Theta);
}

double gaussian_pdf(double x, double mean, double sigma) {
    constexpr double sqrtTwoPi = 2.5066282746310002;
    const double width = std::max(std::abs(sigma), 1.0e-9);
    const double arg = (x - mean) / width;
    return std::exp(-0.5 * arg * arg) / (sqrtTwoPi * width);
}

double convolved_compton_line(double npe, double gammaKeV, double npePerKeV, double sigmaNPE) {
    const double edgeKeV = compton_edge_energy_kev(gammaKeV);
    if (edgeKeV <= 0.0 || npePerKeV <= 0.0) return 0.0;

    constexpr int nSteps = 160;
    const double step = edgeKeV / nSteps;
    double integral = 0.0;
    double norm = 0.0;
    for (int i = 0; i < nSteps; ++i) {
        const double electronKeV = (i + 0.5) * step;
        const double shape = klein_nishina_recoil_shape(electronKeV, gammaKeV);
        if (shape <= 0.0) continue;
        integral += shape * gaussian_pdf(npe, npePerKeV * electronKeV, sigmaNPE) * step;
        norm += shape * step;
    }
    return norm > 0.0 ? integral / norm : 0.0;
}

double compton_response_model(double *x, double *p) {
    const double npePerKeV = std::max(p[2], 1.0e-9);
    const double sigmaNPE = std::max(std::abs(p[3]), 1.0e-9);
    double value = p[0] + p[1] * x[0];
    for (size_t i = 0; i < gResponseGammaLines.size(); ++i) {
        value += p[4 + i] * convolved_compton_line(x[0], gResponseGammaLines[i].energyKeV, npePerKeV, sigmaNPE);
    }
    return value;
}

double derivative_response_model(double *x, double *p) {
    const double npePerKeV = std::max(p[2], 1.0e-9);
    const double sigmaNPE = std::max(std::abs(p[3]), 1.0e-9);
    double value = p[0] + p[1] * x[0];
    for (size_t i = 0; i < gResponseGammaLines.size(); ++i) {
        const double edgeNPE = npePerKeV * compton_edge_energy_kev(gResponseGammaLines[i].energyKeV);
        const double arg = (x[0] - edgeNPE) / sigmaNPE;
        value += p[4 + i] * std::exp(-0.5 * arg * arg);
    }
    return value;
}

bool use_response_model(const char *modelName) {
    TString model(modelName ? modelName : "");
    model.ToLower();
    return model == "compton_response" || model == "kn_convolution" ||
           model == "kn_response" || model == "response";
}

bool use_edge_response_model(const char *modelName) {
    TString model(modelName ? modelName : "");
    model.ToLower();
    return model == "edge_response" || model == "smeared_edge" ||
           model == "erfc_response" || model == "multi_erfc";
}

bool use_derivative_response_model(const char *modelName) {
    TString model(modelName ? modelName : "");
    model.ToLower();
    return model == "derivative_response" || model == "derivative_edge" ||
           model == "diff_response" || model == "diff_edge";
}

TString normalized_source_name(const char *sourceName,
                               const char *histFile = "",
                               const char *outPrefix = "") {
    TString source(sourceName ? sourceName : "");
    if (source.Length() == 0) source = Form("%s %s", histFile ? histFile : "", outPrefix ? outPrefix : "");
    source.ToLower();
    source.ReplaceAll("-", "");
    source.ReplaceAll("_", "");
    source.ReplaceAll(" ", "");
    if (source.Contains("cs137") || source.Contains("137cs")) return "Cs137";
    if (source.Contains("co60") || source.Contains("60co")) return "Co60";
    if (source.Contains("na22") || source.Contains("22na")) return "Na22";
    if (source.Contains("mn54") || source.Contains("54mn")) return "Mn54";
    return "";
}

std::vector<GammaLine> gamma_lines_for_source(const char *sourceName,
                                              const char *histFile = "",
                                              const char *outPrefix = "") {
    const TString source = normalized_source_name(sourceName, histFile, outPrefix);
    std::vector<GammaLine> lines;
    if (source == "Cs137") {
        lines.push_back({"Cs137_661.657", 661.657});
    } else if (source == "Co60") {
        lines.push_back({"Co60_1173.228", 1173.228});
        lines.push_back({"Co60_1332.492", 1332.492});
    } else if (source == "Na22") {
        lines.push_back({"Na22_511.000", 511.000});
        lines.push_back({"Na22_1274.537", 1274.537});
    } else if (source == "Mn54") {
        lines.push_back({"Mn54_834.848", 834.848});
    }
    return lines;
}

double highest_response_edge_kev(const std::vector<GammaLine> &lines) {
    double edge = 0.0;
    for (const GammaLine &line : lines) edge = std::max(edge, compton_edge_energy_kev(line.energyKeV));
    return edge;
}

double lowest_response_edge_kev(const std::vector<GammaLine> &lines) {
    double edge = 0.0;
    for (const GammaLine &line : lines) {
        const double lineEdge = compton_edge_energy_kev(line.energyKeV);
        if (edge <= 0.0 || lineEdge < edge) edge = lineEdge;
    }
    return edge;
}

int find_inrange_bin(const TH1D *hist, double x, bool upperEdge = false) {
    if (!hist) return 0;
    const int nBins = hist->GetNbinsX();
    if (upperEdge && x >= hist->GetXaxis()->GetXmax()) return nBins;
    const int bin = hist->GetXaxis()->FindFixBin(x);
    return std::min(std::max(bin, 1), nBins);
}

double local_maximum_in_window(const TH1D *hist, double x, double halfWidth) {
    if (!hist) return 0.0;
    const int firstBin = find_inrange_bin(hist, x - halfWidth);
    const int lastBin = find_inrange_bin(hist, x + halfWidth, true);
    double best = -1.0e300;
    for (int bin = firstBin; bin <= lastBin; ++bin) {
        best = std::max(best, hist->GetBinContent(bin));
    }
    return best > -1.0e299 ? best : 0.0;
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

    int firstBin = std::max(2, find_inrange_bin(hist, xmin));
    int lastBin = std::min(hist->GetNbinsX() - 1, find_inrange_bin(hist, xmax, true));
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
    return model == "erfc_gaus" || model == "erfc_gaussian" || model == "erfc+gaus" ||
           model == "erfc+gaussian" || model == "recommended" || model == "best_current";
}

TString canonical_model_name(const char *modelName) {
    if (use_derivative_response_model(modelName)) return "derivative_response";
    if (use_edge_response_model(modelName)) return "edge_response";
    if (use_response_model(modelName)) return "compton_response";
    if (use_gaussian_model(modelName)) return "erfc_gaussian";
    return "erfc_linear";
}

TH1D *make_derivative_histogram(const TH1D *hist, int smoothHalfWindow, const char *name) {
    if (!hist) return nullptr;
    const int halfWindow = std::max(1, smoothHalfWindow);
    TH1D *derivative = static_cast<TH1D*>(hist->Clone(name));
    derivative->Reset("ICES");
    derivative->SetDirectory(nullptr);
    derivative->SetTitle(";NPE;-dCounts/dNPE");
    for (int bin = halfWindow + 2; bin <= hist->GetNbinsX() - halfWindow - 1; ++bin) {
        double left = 0.0;
        double right = 0.0;
        double leftErr2 = 0.0;
        double rightErr2 = 0.0;
        for (int offset = 1; offset <= halfWindow; ++offset) {
            left += hist->GetBinContent(bin - offset);
            right += hist->GetBinContent(bin + offset);
            leftErr2 += std::pow(hist->GetBinError(bin - offset), 2);
            rightErr2 += std::pow(hist->GetBinError(bin + offset), 2);
        }
        left /= halfWindow;
        right /= halfWindow;
        const double dx = hist->GetBinCenter(bin + halfWindow) - hist->GetBinCenter(bin - halfWindow);
        if (dx <= 0.0) continue;
        const double value = -(right - left) / dx;
        const double error = std::sqrt(leftErr2 + rightErr2) / (halfWindow * dx);
        derivative->SetBinContent(bin, value);
        derivative->SetBinError(bin, error > 0.0 ? error : 1.0);
    }
    return derivative;
}

double maximum_bin_center_in_range(const TH1D *hist, double xmin, double xmax) {
    if (!hist) return 0.5 * (xmin + xmax);
    double bestY = -1.0e300;
    double bestX = 0.5 * (xmin + xmax);
    for (int bin = find_inrange_bin(hist, xmin); bin <= find_inrange_bin(hist, xmax, true); ++bin) {
        const double x = hist->GetBinCenter(bin);
        if (x < xmin || x > xmax) continue;
        const double y = hist->GetBinContent(bin);
        if (y > bestY) {
            bestY = y;
            bestX = x;
        }
    }
    return bestX;
}

double response_scale_seed_from_peak(const TH1D *fitHist,
                                     const std::vector<GammaLine> &lines,
                                     double fitXMin,
                                     double fitXMax,
                                     double fallbackEdgeX) {
    const double highestEdgeKeV = highest_response_edge_kev(lines);
    if (lines.empty() || highestEdgeKeV <= 0.0) return 1.0e-4;

    const double peakX = maximum_bin_center_in_range(fitHist, fitXMin, fitXMax);
    const double fallbackScale = std::max(fallbackEdgeX / highestEdgeKeV, 1.0e-4);
    if (peakX <= 0.0) return fallbackScale;

    const double binWidth = fitHist ? fitHist->GetBinWidth(1) : (fitXMax - fitXMin) / 100.0;
    const double supportWindow = std::max(3.0 * binWidth, (fitXMax - fitXMin) / 80.0);
    const double clampedFallback = std::min(std::max(fallbackEdgeX, fitXMin), fitXMax);
    double bestScale = fallbackScale;
    double bestScore = -1.0e300;

    for (const GammaLine &seedLine : lines) {
        const double seedEdgeKeV = compton_edge_energy_kev(seedLine.energyKeV);
        if (seedEdgeKeV <= 0.0) continue;
        const double scale = peakX / seedEdgeKeV;
        if (scale <= 0.0 || !std::isfinite(scale)) continue;

        const double highEdgeX = scale * highestEdgeKeV;
        double score = 0.0;
        if (highEdgeX < fitXMin || highEdgeX > fitXMax) {
            const double nearest = std::min(std::max(highEdgeX, fitXMin), fitXMax);
            score -= 1.0e6 + std::abs(highEdgeX - nearest);
        } else {
            score -= 0.01 * std::abs(highEdgeX - clampedFallback);
        }

        for (const GammaLine &line : lines) {
            const double edgeX = scale * compton_edge_energy_kev(line.energyKeV);
            if (edgeX < fitXMin || edgeX > fitXMax) {
                score -= 0.05 * std::max(1.0, std::abs(edgeX - std::min(std::max(edgeX, fitXMin), fitXMax)));
                continue;
            }
            score += std::max(0.0, local_maximum_in_window(fitHist, edgeX, supportWindow));
        }

        if (score > bestScore) {
            bestScore = score;
            bestScale = scale;
        }
    }

    return std::max(bestScale, 1.0e-4);
}

size_t response_line_index_nearest_peak(const std::vector<GammaLine> &lines,
                                        double npePerKeV,
                                        double peakX) {
    size_t bestIndex = 0;
    double bestDistance = 1.0e300;
    for (size_t i = 0; i < lines.size(); ++i) {
        const double edgeX = npePerKeV * compton_edge_energy_kev(lines[i].energyKeV);
        const double distance = std::abs(edgeX - peakX);
        if (distance < bestDistance) {
            bestDistance = distance;
            bestIndex = i;
        }
    }
    return bestIndex;
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

EdgeFitMetrics fit_edge_metrics(TH1D *hist,
                                double fitXMin,
                                double fitXMax,
                                const char *modelName,
                                const char *sourceName,
                                int index) {
    EdgeFitMetrics metrics;
    metrics.xmin = fitXMin;
    metrics.xmax = fitXMax;
    if (!hist || fitXMax <= fitXMin) return metrics;

    const bool useGaussian = use_gaussian_model(modelName);
    const bool useResponse = use_response_model(modelName);
    const bool useEdgeResponse = use_edge_response_model(modelName);
    const bool useDerivative = use_derivative_response_model(modelName);
    const bool useResponseLike = useResponse || useEdgeResponse || useDerivative;

    TH1D *derivativeHist = nullptr;
    TH1D *fitHist = hist;
    if (useDerivative) {
        derivativeHist = make_derivative_histogram(hist, kDerivativeSmoothHalfWindow, Form("hScanDerivative_%d", index));
        if (!derivativeHist) return metrics;
        fitHist = derivativeHist;
    }

    const int leftBin = find_inrange_bin(fitHist, fitXMin);
    const int rightBin = find_inrange_bin(fitHist, fitXMax, true);
    double yLeft = fitHist->GetBinContent(leftBin);
    const double yRight = fitHist->GetBinContent(rightBin);
    const double yMax = fitHist->GetMaximum();
    if (yLeft == 0 && yMax > 0) yLeft = yMax;
    const double amplitude = std::max(yLeft - yRight, 1.0);
    const double slopeEdge = estimate_edge_x(hist, fitXMin, fitXMax);
    const double derivativePeak = useDerivative ? maximum_bin_center_in_range(fitHist, fitXMin, fitXMax) : slopeEdge;
    const double sigmaGuess = std::max((fitXMax - fitXMin) / 20.0, fitHist->GetBinWidth(1));
    std::vector<GammaLine> responseLines;
    double responseHighestEdgeKeV = 0.0;
    if (useResponseLike) {
        responseLines = gamma_lines_for_source(sourceName);
        responseHighestEdgeKeV = highest_response_edge_kev(responseLines);
        if (responseLines.empty() || responseHighestEdgeKeV <= 0.0) {
            delete derivativeHist;
            metrics.status = -2;
            return metrics;
        }
        gResponseGammaLines = responseLines;
    }

    TF1 *fit = nullptr;
    size_t derivativeSeedLineIndex = 0;
    double derivativeSeedLineEdgeKeV = 0.0;
    double derivativeFitXMin = fitXMin;
    double derivativeFitXMax = fitXMax;
    if (useDerivative) {
        const double npePerKeVGuess = response_scale_seed_from_peak(fitHist, responseLines, fitXMin, fitXMax, slopeEdge);
        derivativeSeedLineIndex = response_line_index_nearest_peak(responseLines, npePerKeVGuess, derivativePeak);
        derivativeSeedLineEdgeKeV = compton_edge_energy_kev(responseLines[derivativeSeedLineIndex].energyKeV);
        const double localHalfWidth = std::max(8.0 * fitHist->GetBinWidth(1), 0.25 * (fitXMax - fitXMin));
        derivativeFitXMin = std::max(fitXMin, derivativePeak - localHalfWidth);
        derivativeFitXMax = std::min(fitXMax, derivativePeak + localHalfWidth);
        if (derivativeFitXMax <= derivativeFitXMin) {
            derivativeFitXMin = fitXMin;
            derivativeFitXMax = fitXMax;
        }

        fit = new TF1(Form("fScanDerivativePeak_%d", index), linear_gaussian_model, derivativeFitXMin, derivativeFitXMax, 5);
        const double baselineGuess = std::min(yLeft, yRight);
        const double peakAmpGuess = std::max(yMax - baselineGuess, 1.0);
        fit->SetParameters(baselineGuess, 0.0, peakAmpGuess, derivativePeak, sigmaGuess);
        fit->SetParLimits(2, 0.0, std::max(peakAmpGuess * 100.0, 1.0));
        fit->SetParLimits(3, derivativeFitXMin, derivativeFitXMax);
        fit->SetParLimits(4, std::max(2.0 * fitHist->GetBinWidth(1), 0.03 * (derivativeFitXMax - derivativeFitXMin)),
                          derivativeFitXMax - derivativeFitXMin);
    } else if (useEdgeResponse) {
        fit = new TF1(Form("fScanEdgeResponse_%d", index), smeared_edge_response_model, fitXMin, fitXMax, 5 + responseLines.size());
        const double observedEdge = (slopeEdge > fitXMin && slopeEdge < fitXMax) ? slopeEdge : derivativePeak;
        const double seedEdgeKeV = responseLines.size() > 1 ? lowest_response_edge_kev(responseLines) : responseHighestEdgeKeV;
        double npePerKeVGuess = seedEdgeKeV > 0.0 ? observedEdge / seedEdgeKeV : 1.0e-4;
        if (npePerKeVGuess * responseHighestEdgeKeV < fitXMin || npePerKeVGuess * responseHighestEdgeKeV > fitXMax) {
            npePerKeVGuess = std::max(observedEdge / responseHighestEdgeKeV, 1.0e-4);
        }
        fit->SetParameters(yRight, 0.0, 0.0, npePerKeVGuess, sigmaGuess);
        fit->SetParLimits(3, std::max(1.0e-6, 0.5 * npePerKeVGuess),
                          std::max(1.0e-5, 1.5 * npePerKeVGuess));
        fit->SetParLimits(4, fitHist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        const double responseAmpGuess = std::max(amplitude, yMax);
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParameter(5 + i, responseAmpGuess / responseLines.size());
            fit->SetParLimits(5 + i, 0.0, std::max(responseAmpGuess * 50.0, 1.0));
        }
    } else if (useResponse) {
        fit = new TF1(Form("fScanComptonResponse_%d", index), compton_response_model, fitXMin, fitXMax, 4 + responseLines.size());
        const double edgeSeed = (slopeEdge > fitXMin && slopeEdge < fitXMax) ? slopeEdge : 0.5 * (fitXMin + fitXMax);
        const double npePerKeVGuess = std::max(edgeSeed / responseHighestEdgeKeV, 1.0e-4);
        fit->SetParameters(yRight, 0.0, npePerKeVGuess, sigmaGuess);
        fit->SetParLimits(2, std::max(1.0e-6, 0.2 * fitXMin / responseHighestEdgeKeV),
                          std::max(1.0e-5, 2.0 * fitXMax / responseHighestEdgeKeV));
        fit->SetParLimits(3, fitHist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        const double responseAmpGuess = std::max(yMax * std::max(fitXMax - fitXMin, 1.0), amplitude);
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParameter(4 + i, responseAmpGuess / responseLines.size());
            fit->SetParLimits(4 + i, 0.0, responseAmpGuess * 100.0);
        }
    } else {
        fit = useGaussian
            ? new TF1(Form("fScanComptonEdgeErfcGaussian_%d", index), compton_edge_erfc_gaus, fitXMin, fitXMax, 8)
            : new TF1(Form("fScanComptonEdgeErfc_%d", index), compton_edge_erfc, fitXMin, fitXMax, 5);
        fit->SetParameters(yRight, 0.0, amplitude, sigmaGuess, slopeEdge);
        fit->SetParLimits(2, 0.0, std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(3, fitHist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        fit->SetParLimits(4, fitXMin, fitXMax);
    }
    if (useGaussian && !useResponseLike) {
        fit->SetParameter(5, 0.2 * amplitude);
        fit->SetParameter(6, slopeEdge);
        fit->SetParameter(7, std::max((fitXMax - fitXMin) / 10.0, hist->GetBinWidth(1)));
        fit->SetParLimits(5, -std::max(yMax * 20.0, amplitude * 20.0), std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(6, fitXMin, fitXMax);
        fit->SetParLimits(7, fitHist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    }

    auto result = fitHist->Fit(fit, "RQ0SN");
    metrics.status = static_cast<int>(result);
    if (useDerivative) {
        const double peakCenter = fit->GetParameter(3);
        const double peakCenterErr = fit->GetParError(3);
        const double npePerKeV = derivativeSeedLineEdgeKeV > 0.0 ? peakCenter / derivativeSeedLineEdgeKeV : 0.0;
        const double npePerKeVErr = derivativeSeedLineEdgeKeV > 0.0 ? peakCenterErr / derivativeSeedLineEdgeKeV : 0.0;
        metrics.edge = npePerKeV * responseHighestEdgeKeV;
        metrics.edgeErr = npePerKeVErr * responseHighestEdgeKeV;
        metrics.sigma = fit->GetParameter(4);
    } else if (useEdgeResponse) {
        metrics.edge = fit->GetParameter(3) * responseHighestEdgeKeV;
        metrics.edgeErr = fit->GetParError(3) * responseHighestEdgeKeV;
        metrics.sigma = fit->GetParameter(4);
    } else if (useResponse) {
        metrics.edge = fit->GetParameter(2) * responseHighestEdgeKeV;
        metrics.edgeErr = fit->GetParError(2) * responseHighestEdgeKeV;
        metrics.sigma = fit->GetParameter(3);
    } else {
        metrics.edge = fit->GetParameter(4);
        metrics.edgeErr = fit->GetParError(4);
        metrics.sigma = fit->GetParameter(3);
    }
    const int ndf = fit->GetNDF();
    metrics.chi2ndf = ndf > 0 ? fit->GetChisquare() / ndf : 0.0;
    delete fit;
    delete derivativeHist;
    return metrics;
}

void fit_compton_edge(const char *histFile = "compton_edge_histograms.root",
                      double fitXMin = 0.0,
                      double fitXMax = -1.0,
                      const char *outPrefix = "compton_edge_fit",
                      const char *histName = "",
                      const char *modelName = "erfc_linear",
                      const char *sourceName = "") {
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
    TH1D *originalHist = hist;
    input->Close();

    if (fitXMax <= fitXMin) {
        fitXMin = hist->GetXaxis()->GetXmin();
        fitXMax = hist->GetXaxis()->GetXmax();
    }
    if (fitXMax <= fitXMin) {
        std::cerr << "Error: invalid fit range." << std::endl;
        return;
    }

    TString model(modelName ? modelName : "erfc_linear");
    model.ToLower();
    const bool useResponse = use_response_model(modelName);
    const bool useEdgeResponse = use_edge_response_model(modelName);
    const bool useDerivative = use_derivative_response_model(modelName);
    const bool useResponseLike = useResponse || useEdgeResponse || useDerivative;
    const bool useGaussian = use_gaussian_model(modelName);

    TH1D *derivativeHist = nullptr;
    if (useDerivative) {
        derivativeHist = make_derivative_histogram(originalHist, kDerivativeSmoothHalfWindow, "hComptonEdgeDerivativeInput");
        if (!derivativeHist) {
            std::cerr << "Error: failed to build derivative histogram for derivative_response model." << std::endl;
            return;
        }
        hist = derivativeHist;
        model = "derivative_response";
    }

    const int leftBin = find_inrange_bin(hist, fitXMin);
    const int rightBin = find_inrange_bin(hist, fitXMax, true);
    double yLeft = hist->GetBinContent(leftBin);
    double yRight = hist->GetBinContent(rightBin);
    const double yMax = hist->GetMaximum();
    if (yLeft == 0 && yMax > 0) yLeft = yMax;
    const double amplitude = std::max(yLeft - yRight, 1.0);
    const double slopeEdge = estimate_edge_x(originalHist, fitXMin, fitXMax);
    const double derivativePeak = useDerivative ? maximum_bin_center_in_range(hist, fitXMin, fitXMax) : slopeEdge;
    const double sigmaGuess = std::max((fitXMax - fitXMin) / 20.0, hist->GetBinWidth(1));

    TString responseSource;
    std::vector<GammaLine> responseLines;
    double responseHighestEdgeKeV = 0.0;
    if (useResponseLike) {
        responseSource = normalized_source_name(sourceName, histFile, outPrefix);
        responseLines = gamma_lines_for_source(responseSource.Data(), histFile, outPrefix);
        responseHighestEdgeKeV = highest_response_edge_kev(responseLines);
        if (responseLines.empty() || responseHighestEdgeKeV <= 0.0) {
            std::cerr << "Error: cannot infer gamma line list for " << model << " model." << std::endl;
            std::cerr << "       Pass -S Cs137, Co60, Na22, or Mn54 to the fit wrapper." << std::endl;
            return;
        }
        gResponseGammaLines = responseLines;
        if (useResponse) model = "compton_response";
        if (useEdgeResponse) model = "edge_response";
    } else if (useGaussian) {
        model = "erfc_gaussian";
    } else if (!useGaussian && model != "erfc_linear") {
        std::cerr << "Warning: unknown model '" << modelName << "', using erfc_linear." << std::endl;
        model = "erfc_linear";
    }

    TF1 *fit = nullptr;
    size_t derivativeSeedLineIndex = 0;
    double derivativeSeedLineEdgeKeV = 0.0;
    double derivativeFitXMin = fitXMin;
    double derivativeFitXMax = fitXMax;
    if (useDerivative) {
        const double npePerKeVGuess = response_scale_seed_from_peak(hist, responseLines, fitXMin, fitXMax, slopeEdge);
        derivativeSeedLineIndex = response_line_index_nearest_peak(responseLines, npePerKeVGuess, derivativePeak);
        derivativeSeedLineEdgeKeV = compton_edge_energy_kev(responseLines[derivativeSeedLineIndex].energyKeV);
        const double localHalfWidth = std::max(8.0 * hist->GetBinWidth(1), 0.25 * (fitXMax - fitXMin));
        derivativeFitXMin = std::max(fitXMin, derivativePeak - localHalfWidth);
        derivativeFitXMax = std::min(fitXMax, derivativePeak + localHalfWidth);
        if (derivativeFitXMax <= derivativeFitXMin) {
            derivativeFitXMin = fitXMin;
            derivativeFitXMax = fitXMax;
        }

        fit = new TF1("fDerivativePeak", linear_gaussian_model, derivativeFitXMin, derivativeFitXMax, 5);
        fit->SetParNames("offset", "slope", "amplitude", "peak_npe", "sigma_npe");
    } else if (useEdgeResponse) {
        fit = new TF1("fEdgeResponse", smeared_edge_response_model, fitXMin, fitXMax, 5 + responseLines.size());
        fit->SetParName(0, "offset");
        fit->SetParName(1, "slope");
        fit->SetParName(2, "curvature");
        fit->SetParName(3, "npe_per_kev");
        fit->SetParName(4, "sigma_npe");
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParName(5 + i, Form("amp_%s", responseLines[i].label.Data()));
        }
    } else if (useResponse) {
        fit = new TF1("fComptonResponse", compton_response_model, fitXMin, fitXMax, 4 + responseLines.size());
        fit->SetParName(0, "offset");
        fit->SetParName(1, "slope");
        fit->SetParName(2, "npe_per_kev");
        fit->SetParName(3, "sigma_npe");
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParName(4 + i, Form("amp_%s", responseLines[i].label.Data()));
        }
    } else if (useGaussian) {
        fit = new TF1("fComptonEdgeErfcGaussian", compton_edge_erfc_gaus, fitXMin, fitXMax, 8);
        fit->SetParNames("offset", "slope", "amplitude", "sigma", "edge", "gaus_amp", "gaus_mean", "gaus_sigma");
    } else {
        fit = new TF1("fComptonEdgeErfc", compton_edge_erfc, fitXMin, fitXMax, 5);
        fit->SetParNames("offset", "slope", "amplitude", "sigma", "edge");
    }
    if (useDerivative) {
        const double baselineGuess = std::min(yLeft, yRight);
        const double peakAmpGuess = std::max(yMax - baselineGuess, 1.0);
        fit->SetParameters(baselineGuess, 0.0, peakAmpGuess, derivativePeak, sigmaGuess);
        fit->SetParLimits(2, 0.0, std::max(peakAmpGuess * 100.0, 1.0));
        fit->SetParLimits(3, derivativeFitXMin, derivativeFitXMax);
        fit->SetParLimits(4, std::max(2.0 * hist->GetBinWidth(1), 0.03 * (derivativeFitXMax - derivativeFitXMin)),
                          derivativeFitXMax - derivativeFitXMin);
    } else if (useEdgeResponse) {
        const double observedEdge = (slopeEdge > fitXMin && slopeEdge < fitXMax) ? slopeEdge : 0.5 * (fitXMin + fitXMax);
        const double seedEdgeKeV = responseLines.size() > 1 ? lowest_response_edge_kev(responseLines) : responseHighestEdgeKeV;
        double npePerKeVGuess = seedEdgeKeV > 0.0 ? observedEdge / seedEdgeKeV : 1.0e-4;
        if (npePerKeVGuess * responseHighestEdgeKeV < fitXMin || npePerKeVGuess * responseHighestEdgeKeV > fitXMax) {
            npePerKeVGuess = std::max(observedEdge / responseHighestEdgeKeV, 1.0e-4);
        }
        fit->SetParameters(yRight, 0.0, 0.0, npePerKeVGuess, sigmaGuess);
        fit->SetParLimits(3, std::max(1.0e-6, 0.5 * npePerKeVGuess),
                          std::max(1.0e-5, 1.5 * npePerKeVGuess));
        fit->SetParLimits(4, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        const double responseAmpGuess = std::max(amplitude, yMax);
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParameter(5 + i, responseAmpGuess / responseLines.size());
            fit->SetParLimits(5 + i, 0.0, std::max(responseAmpGuess * 50.0, 1.0));
        }
    } else if (useResponse) {
        const double edgeSeed = (slopeEdge > fitXMin && slopeEdge < fitXMax) ? slopeEdge : 0.5 * (fitXMin + fitXMax);
        const double npePerKeVGuess = std::max(edgeSeed / responseHighestEdgeKeV, 1.0e-4);
        fit->SetParameters(yRight, 0.0, npePerKeVGuess, sigmaGuess);
        fit->SetParLimits(2, std::max(1.0e-6, 0.2 * fitXMin / responseHighestEdgeKeV),
                          std::max(1.0e-5, 2.0 * fitXMax / responseHighestEdgeKeV));
        fit->SetParLimits(3, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        const double responseAmpGuess = std::max(yMax * std::max(fitXMax - fitXMin, 1.0), amplitude);
        for (size_t i = 0; i < responseLines.size(); ++i) {
            fit->SetParameter(4 + i, responseAmpGuess / responseLines.size());
            fit->SetParLimits(4 + i, 0.0, responseAmpGuess * 100.0);
        }
    } else {
        fit->SetParameters(yRight, 0.0, amplitude, sigmaGuess, slopeEdge);
        fit->SetParLimits(2, 0.0, std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(3, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
        fit->SetParLimits(4, fitXMin, fitXMax);
    }
    if (useGaussian && !useResponseLike) {
        fit->SetParameter(5, 0.2 * amplitude);
        fit->SetParameter(6, slopeEdge);
        fit->SetParameter(7, std::max((fitXMax - fitXMin) / 10.0, hist->GetBinWidth(1)));
        fit->SetParLimits(5, -std::max(yMax * 20.0, amplitude * 20.0), std::max(yMax * 20.0, amplitude * 20.0));
        fit->SetParLimits(6, fitXMin, fitXMax);
        fit->SetParLimits(7, hist->GetBinWidth(1) * 0.25, fitXMax - fitXMin);
    }

    auto result = hist->Fit(fit, "RQ0S");
    const int status = static_cast<int>(result);

    double npePerKeV = 0.0;
    double npePerKeVErr = 0.0;
    double edge = 0.0;
    double edgeErr = 0.0;
    double sigma = 0.0;
    double sigmaErr = 0.0;
    if (useDerivative) {
        const double peakCenter = fit->GetParameter(3);
        const double peakCenterErr = fit->GetParError(3);
        npePerKeV = derivativeSeedLineEdgeKeV > 0.0 ? peakCenter / derivativeSeedLineEdgeKeV : 0.0;
        npePerKeVErr = derivativeSeedLineEdgeKeV > 0.0 ? peakCenterErr / derivativeSeedLineEdgeKeV : 0.0;
        edge = npePerKeV * responseHighestEdgeKeV;
        edgeErr = npePerKeVErr * responseHighestEdgeKeV;
        sigma = fit->GetParameter(4);
        sigmaErr = fit->GetParError(4);
    } else if (useEdgeResponse) {
        npePerKeV = fit->GetParameter(3);
        npePerKeVErr = fit->GetParError(3);
        edge = npePerKeV * responseHighestEdgeKeV;
        edgeErr = npePerKeVErr * responseHighestEdgeKeV;
        sigma = fit->GetParameter(4);
        sigmaErr = fit->GetParError(4);
    } else if (useResponse) {
        npePerKeV = fit->GetParameter(2);
        npePerKeVErr = fit->GetParError(2);
        edge = npePerKeV * responseHighestEdgeKeV;
        edgeErr = npePerKeVErr * responseHighestEdgeKeV;
        sigma = fit->GetParameter(3);
        sigmaErr = fit->GetParError(3);
    } else {
        edge = fit->GetParameter(4);
        edgeErr = fit->GetParError(4);
        sigma = fit->GetParameter(3);
        sigmaErr = fit->GetParError(3);
    }
    double amplitudeFit = 0.0;
    double amplitudeErr = 0.0;
    if (useDerivative) {
        amplitudeFit = fit->GetParameter(2);
        amplitudeErr = fit->GetParError(2);
    } else if (useEdgeResponse) {
        for (size_t i = 0; i < responseLines.size(); ++i) amplitudeFit += fit->GetParameter(5 + i);
    } else if (useResponse) {
        for (size_t i = 0; i < responseLines.size(); ++i) amplitudeFit += fit->GetParameter(4 + i);
    } else {
        amplitudeFit = fit->GetParameter(2);
        amplitudeErr = fit->GetParError(2);
    }
    const double gausAmp = (useGaussian && !useResponseLike) ? fit->GetParameter(5) : 0.0;
    const double gausAmpErr = (useGaussian && !useResponseLike) ? fit->GetParError(5) : 0.0;
    const double gausMean = (useGaussian && !useResponseLike) ? fit->GetParameter(6) : 0.0;
    const double gausMeanErr = (useGaussian && !useResponseLike) ? fit->GetParError(6) : 0.0;
    const double gausSigma = (useGaussian && !useResponseLike) ? fit->GetParameter(7) : 0.0;
    const double gausSigmaErr = (useGaussian && !useResponseLike) ? fit->GetParError(7) : 0.0;
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
    hist->SetTitle(useDerivative ? ";NPE;-dCounts/dNPE" : ";NPE;Counts");
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.65);
    hist->SetLineWidth(2);
    hist->SetMaximum(hist->GetMaximum() * 1.18);
    if (useDerivative) hist->SetMinimum(std::min(0.0, hist->GetMinimum() * 1.15));
    hist->GetXaxis()->SetLabelSize(0);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.048);
    hist->GetYaxis()->SetTitleOffset(1.05);
    hist->Draw("E");
    const double plotYMin = hist->GetMinimum();
    const double plotYMax = hist->GetMaximum();

    TF1 *erfcComponent = nullptr;
    TF1 *gausComponent = nullptr;
    TLine *gausMeanLine = nullptr;
    std::vector<TLine*> responseEdgeLines;
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

    TLine *edgeLine = new TLine(edge, plotYMin, edge, plotYMax * 1.05);
    edgeLine->SetLineColor(kBlue + 1);
    edgeLine->SetLineStyle(2);
    edgeLine->SetLineWidth(2);
    edgeLine->Draw("SAME");

    if (useResponseLike) {
        for (size_t i = 0; i < responseLines.size(); ++i) {
            const double lineEdge = npePerKeV * compton_edge_energy_kev(responseLines[i].energyKeV);
            if (std::abs(lineEdge - edge) < 1.0e-9) continue;
            TLine *line = new TLine(lineEdge, plotYMin, lineEdge, plotYMax * 1.02);
            line->SetLineColor(kBlue + 2);
            line->SetLineStyle(7);
            line->SetLineWidth(2);
            line->Draw("SAME");
            responseEdgeLines.push_back(line);
        }
    }

    TLine *slopeLine = new TLine(slopeEdge, plotYMin, slopeEdge, plotYMax * 1.05);
    slopeLine->SetLineColor(kGreen + 2);
    slopeLine->SetLineStyle(3);
    slopeLine->SetLineWidth(2);
    slopeLine->Draw("SAME");

    if (useGaussian) {
        gausMeanLine = new TLine(gausMean, plotYMin, gausMean, plotYMax * 1.02);
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
    legend->SetTextSize((useGaussian || useResponseLike) ? 0.043 : 0.060);
    legend->AddEntry(hist, useDerivative ? "Derivative histogram" : "Histogram", "lep");
    legend->AddEntry(fit, useDerivative ? "local derivative peak response"
                         : (useEdgeResponse ? "line-tied smeared edge response"
                                            : (useResponse ? "Klein-Nishina convolution response"
                                                           : (useGaussian ? "linear + erfc edge + Gaussian" : "linear + erfc edge"))), "l");
    if (useGaussian) {
        legend->AddEntry(erfcComponent, "linear + erfc component", "l");
        legend->AddEntry(gausComponent, "Gaussian contribution", "l");
        legend->AddEntry(gausMeanLine, "Gaussian mean", "l");
    }
    legend->AddEntry(edgeLine, useResponseLike ? "highest-#gamma Compton edge" : "edge", "l");
    if (useResponseLike && !responseEdgeLines.empty()) legend->AddEntry(responseEdgeLines.front(), "other line edge", "l");
    legend->AddEntry(slopeLine, useDerivative ? "raw-spectrum max falling slope" : "max falling slope", "l");
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
    const double pullXMin = useDerivative ? derivativeFitXMin : fitXMin;
    const double pullXMax = useDerivative ? derivativeFitXMax : fitXMax;
    const int firstFitBin = find_inrange_bin(hist, pullXMin);
    const int lastFitBin = find_inrange_bin(hist, pullXMax, true);
    for (int bin = firstFitBin; bin <= lastFitBin; ++bin) {
        const double x = hist->GetBinCenter(bin);
        if (x < pullXMin || x > pullXMax) continue;
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
    if (useDerivative) text->AddText(Form("Derivative peak seed = %.6g NPE", derivativePeak));
    text->AddText(Form("#sigma = %.6g #pm %.3g NPE", sigma, sigmaErr));
    if (useGaussian) text->AddText(Form("Gaus #mu = %.6g #pm %.3g NPE", gausMean, gausMeanErr));
    if (useResponseLike) text->AddText(Form("NPE/keV = %.6g #pm %.3g", npePerKeV, npePerKeVErr));
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
    TLine *zeroPull = new TLine(pullXMin, 0.0, pullXMax, 0.0);
    zeroPull->SetLineStyle(2);
    zeroPull->Draw("SAME");
    TLine *plusThree = new TLine(pullXMin, 3.0, pullXMax, 3.0);
    plusThree->SetLineStyle(3);
    plusThree->SetLineColor(kGray + 2);
    plusThree->Draw("SAME");
    TLine *minusThree = new TLine(pullXMin, -3.0, pullXMax, -3.0);
    minusThree->SetLineStyle(3);
    minusThree->SetLineColor(kGray + 2);
    minusThree->Draw("SAME");

    canvas->cd();

    canvas->SaveAs(pngName);
    canvas->SaveAs(pdfName);

    TFile *output = TFile::Open(rootName, "RECREATE");
    if (output && !output->IsZombie()) {
        if (useDerivative) originalHist->Write("hComptonEdgeOriginalInput");
        hist->Write("hComptonEdgeFitInput");
        pullHist->Write("hComptonEdgeFitPull");
        fit->Write(useDerivative ? "fDerivativePeak" : (useEdgeResponse ? "fEdgeResponse" : (useResponse ? "fComptonResponse" : "fComptonEdgeErfc")));
        if (erfcComponent) erfcComponent->Write("fComptonEdgeErfcComponent");
        if (gausComponent) gausComponent->Write("fComptonEdgeGaussianComponent");
        canvas->Write("c_compton_edge_fit");
        output->Close();
    }

    std::ofstream txt(txtName.Data());
    txt << "model " << model << "\n";
    txt << "status " << status << "\n";
    txt << "hist_file " << histFile << "\n";
    txt << "hist_name " << (histName && TString(histName).Length() > 0 ? histName : hist->GetName()) << "\n";
    txt << "fit_input " << (useDerivative ? "smoothed_negative_derivative" : "histogram") << "\n";
    if (useResponseLike) txt << "response_source " << responseSource << "\n";
    txt << "fit_xmin " << fitXMin << "\n";
    txt << "fit_xmax " << fitXMax << "\n";
    if (useDerivative) {
        txt << "derivative_fit_xmin " << derivativeFitXMin << "\n";
        txt << "derivative_fit_xmax " << derivativeFitXMax << "\n";
        txt << "derivative_seed_line_index " << derivativeSeedLineIndex << "\n";
        txt << "derivative_seed_line_label " << responseLines[derivativeSeedLineIndex].label << "\n";
    }
    txt << "edge " << edge << "\n";
    txt << "edge_error " << edgeErr << "\n";
    txt << "slope_edge " << slopeEdge << "\n";
    if (useDerivative) txt << "derivative_peak_seed " << derivativePeak << "\n";
    txt << "sigma " << sigma << "\n";
    txt << "sigma_error " << sigmaErr << "\n";
    txt << "amplitude " << amplitudeFit << "\n";
    txt << "amplitude_error " << amplitudeErr << "\n";
    if (useResponseLike) {
        txt << "npe_per_kev " << npePerKeV << "\n";
        txt << "npe_per_kev_error " << npePerKeVErr << "\n";
        txt << "response_lines " << responseLines.size() << "\n";
        for (size_t i = 0; i < responseLines.size(); ++i) {
            const double edgeKeV = compton_edge_energy_kev(responseLines[i].energyKeV);
            txt << "line_" << i
                << " label=" << responseLines[i].label
                << " gamma_kev=" << responseLines[i].energyKeV
                << " edge_kev=" << edgeKeV
                << " edge_npe=" << npePerKeV * edgeKeV
                << " edge_npe_error=" << npePerKeVErr * edgeKeV
                << " amplitude=" << (useDerivative ? (i == derivativeSeedLineIndex ? amplitudeFit : 0.0)
                                                     : (useEdgeResponse ? fit->GetParameter(5 + i) : fit->GetParameter(4 + i)))
                << " amplitude_error=" << (useDerivative ? (i == derivativeSeedLineIndex ? amplitudeErr : 0.0)
                                                           : (useEdgeResponse ? fit->GetParError(5 + i) : fit->GetParError(4 + i)))
                << "\n";
        }
    }
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
    std::cout << "Model: " << model << std::endl;
    if (useResponseLike) std::cout << "Response source: " << responseSource << std::endl;
    std::cout << "Histogram ROOT file: " << histFile << std::endl;
    std::cout << "Fit range [NPE]: [" << fitXMin << ", " << fitXMax << "]" << std::endl;
    std::cout << "Fit status: " << status << std::endl;
    std::cout << "Edge [NPE]: " << edge << " +/- " << edgeErr << std::endl;
    std::cout << "Max falling slope edge seed [NPE]: " << slopeEdge << std::endl;
    if (useDerivative) std::cout << "Derivative peak seed [NPE]: " << derivativePeak << std::endl;
    if (useDerivative) {
        std::cout << "Derivative local fit range [NPE]: [" << derivativeFitXMin << ", " << derivativeFitXMax << "]" << std::endl;
        std::cout << "Derivative seed line: " << responseLines[derivativeSeedLineIndex].label << std::endl;
    }
    std::cout << "Sigma [NPE]: " << sigma << " +/- " << sigmaErr << std::endl;
    std::cout << "Amplitude: " << amplitudeFit << " +/- " << amplitudeErr << std::endl;
    if (useResponseLike) {
        std::cout << "NPE/keV: " << npePerKeV << " +/- " << npePerKeVErr << std::endl;
        for (size_t i = 0; i < responseLines.size(); ++i) {
            const double edgeKeV = compton_edge_energy_kev(responseLines[i].energyKeV);
            std::cout << "Line " << responseLines[i].label
                      << ": gamma=" << responseLines[i].energyKeV
                      << " keV, edge=" << edgeKeV
                      << " keV -> " << npePerKeV * edgeKeV
                      << " +/- " << npePerKeVErr * edgeKeV
                      << " NPE, amplitude="
                      << (useDerivative ? (i == derivativeSeedLineIndex ? amplitudeFit : 0.0)
                                        : (useEdgeResponse ? fit->GetParameter(5 + i) : fit->GetParameter(4 + i)))
                      << " +/- "
                      << (useDerivative ? (i == derivativeSeedLineIndex ? amplitudeErr : 0.0)
                                        : (useEdgeResponse ? fit->GetParError(5 + i) : fit->GetParError(4 + i)))
                      << std::endl;
        }
    }
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
    std::cout << "FIT_RESULT model=" << model
              << " edge=" << edge
              << " edge_error=" << edgeErr
              << " slope_edge=" << slopeEdge
              << " sigma=" << sigma
              << " sigma_error=" << sigmaErr
              << " chi2_ndf=" << chi2ndf
              << " pull_rms=" << pullRms;
    if (useDerivative) {
        std::cout << " derivative_peak_seed=" << derivativePeak
                  << " derivative_seed_line=" << responseLines[derivativeSeedLineIndex].label
                  << " derivative_fit_xmin=" << derivativeFitXMin
                  << " derivative_fit_xmax=" << derivativeFitXMax;
    }
    if (useResponseLike) {
        std::cout << " source=" << responseSource
                  << " npe_per_kev=" << npePerKeV
                  << " npe_per_kev_error=" << npePerKeVErr
                  << " response_lines=" << responseLines.size();
    }
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
                           double scanFraction = 0.10,
                           const char *sourceName = "") {
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
    const bool useResponseLike = use_response_model(modelName) ||
                                 use_edge_response_model(modelName) ||
                                 use_derivative_response_model(modelName);
    if (useResponseLike && gamma_lines_for_source(sourceName, histFile, outPrefix).empty()) {
        std::cerr << "Error: cannot infer gamma line list for " << canonical_model_name(modelName) << " model." << std::endl;
        std::cerr << "       Pass -S Cs137, Co60, Na22, or Mn54 to the scan wrapper." << std::endl;
        return;
    }
    const TString modelLabel = canonical_model_name(modelName);
    const TString responseSource = useResponseLike ? normalized_source_name(sourceName, histFile, outPrefix) : "";

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
            scans.push_back(fit_edge_metrics(hist, xmin, xmax, modelName, responseSource.Data(), scanIndex++));
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
    txt << "model " << modelLabel << "\n";
    if (useResponseLike) txt << "response_source " << responseSource << "\n";
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
    std::cout << "Model: " << modelLabel << std::endl;
    if (useResponseLike) std::cout << "Response source: " << responseSource << std::endl;
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
    std::cout << "SCAN_RESULT model=" << modelLabel
              << " edge_mean=" << edgeMean
              << " edge_rms=" << edgeRms
              << " edge_min=" << edgeMin
              << " edge_max=" << edgeMax
              << " n_good=" << goodScans.size()
              << " n_total=" << scans.size();
    if (useResponseLike) std::cout << " source=" << responseSource;
    if (goodScans.size() < 3) {
        std::cout << " warnings=too_few_good_fits";
    } else {
        std::cout << " warnings=none";
    }
    std::cout << std::endl;
}
