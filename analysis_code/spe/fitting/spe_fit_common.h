#ifndef MD_SPE_FIT_COMMON_H
#define MD_SPE_FIT_COMMON_H

#include <cmath>
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TTree.h>

inline double spe_adc_integral_to_pC(double dynamicRangeV = 2.0,
                                     double samplingTimeNs = 2.0,
                                     double resistanceOhm = 50.0,
                                     int adcBits = 14) {
    const double maxADC = std::pow(2.0, adcBits) - 1.0;
    return (dynamicRangeV / maxADC) * (samplingTimeNs * 1.0e-9) / resistanceOhm * 1.0e12;
}

inline TTree *spe_find_charge_tree(TFile *file) {
    if (!file || file->IsZombie()) return nullptr;
    if (TTree *tree = static_cast<TTree*>(file->Get("phys_tree"))) return tree;
    return static_cast<TTree*>(file->Get("T_Charge"));
}

inline TString spe_charge_expr(TTree *tree,
                               int ch,
                               double dynamicRangeV = 2.0,
                               double samplingTimeNs = 2.0,
                               double resistanceOhm = 50.0,
                               int adcBits = 14,
                               double chargeOffsetPC = 0.0) {
    if (!tree) return "";

    TString newBranch = Form("Charge_CH%d", ch);
    if (tree->GetBranch(newBranch)) {
        const double toPC = spe_adc_integral_to_pC(dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits);
        return Form("(%s * %.17g - %.17g)", newBranch.Data(), toPC, chargeOffsetPC);
    }

    if (tree->GetBranch("charge_pC")) {
        return Form("(charge_pC[%d] - %.17g)", ch, chargeOffsetPC);
    }

    return "";
}

inline TH1D *spe_make_charge_hist(const char *fileName,
                                  int ch,
                                  const char *histName,
                                  const char *title,
                                  int bins,
                                  double xmin,
                                  double xmax,
                                  double dynamicRangeV = 2.0,
                                  double samplingTimeNs = 2.0,
                                  double resistanceOhm = 50.0,
                                  int adcBits = 14,
                                  double chargeOffsetPC = 0.0) {
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << fileName << std::endl;
        delete file;
        return nullptr;
    }

    TTree *tree = spe_find_charge_tree(file);
    if (!tree) {
        std::cerr << "Cannot find phys_tree or T_Charge in " << fileName << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    TString expr = spe_charge_expr(tree, ch, dynamicRangeV, samplingTimeNs, resistanceOhm, adcBits, chargeOffsetPC);
    if (expr.IsNull()) {
        std::cerr << "Cannot find Charge_CH" << ch << " or charge_pC[" << ch << "] in " << fileName << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    tree->SetCacheSize(64 * 1024 * 1024);
    TString branch = Form("Charge_CH%d", ch);
    if (tree->GetBranch(branch)) tree->AddBranchToCache(branch, kTRUE);
    if (tree->GetBranch("charge_pC")) tree->AddBranchToCache("charge_pC", kTRUE);

    TH1D *hist = new TH1D(histName, title, bins, xmin, xmax);
    tree->Project(hist->GetName(), expr, "", "goff");
    hist->SetDirectory(nullptr);

    file->Close();
    delete file;
    return hist;
}

inline TString spe_format_file(const char *pattern, int voltage) {
    return TString::Format(pattern, voltage);
}

#endif
