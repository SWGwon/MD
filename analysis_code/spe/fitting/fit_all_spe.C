#include <iostream>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSpectrum.h>
#include "spe_fit_common.h"

Double_t spe_model_final(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t N = par[0];
    Double_t mu = par[1];
    Double_t Q0 = par[2];
    Double_t s0 = par[3];
    Double_t Q1 = par[4];
    Double_t s1 = par[5];
    Double_t bg = par[6];

    Double_t sum = 0.0;
    for (int n = 0; n <= 8; n++) {
        Double_t P_n = TMath::Poisson(n, mu);
        Double_t Q_n = Q0 + n * Q1;
        Double_t sigma_n = TMath::Sqrt(s0 * s0 + n * s1 * s1);
        Double_t G_n = (1.0 / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma_n)) * TMath::Exp(-0.5 * TMath::Power((xx - Q_n) / sigma_n, 2));
        sum += P_n * G_n;
    }
    return N * sum + bg;
}

void fit_all_spe(int ch = 0, const char *ledPattern = "output_charge_run2_%dV.root") {
    gStyle->SetOptFit(1111);
    std::vector<int> voltages = {1600, 1700, 1800, 1900, 2000, 2100, 2200};
    TGraphErrors *gr_gain = new TGraphErrors();
    gr_gain->SetMarkerStyle(21);
    gr_gain->SetMarkerColor(kRed);

    for (int v : voltages) {
        TString led_fn = spe_format_file(ledPattern, v);
        
        double x_max = (v < 2000) ? 8.0 : 45.0;
        TH1D *h = spe_make_charge_hist(
            led_fn, ch, "h",
            TString::Format("%dV Ch%d;Charge [pC];Counts", v, ch),
            600, -1.0, x_max
        );
        if (!h) continue;

        // 1. Pedestal Fit (Focus on small area near 0)
        TF1 *f_ped_pre = new TF1("f_ped_pre", "gaus", -0.2, 0.3);
        h->Fit(f_ped_pre, "RQN");
        double q0 = f_ped_pre->GetParameter(1);
        double s0 = f_ped_pre->GetParameter(2);

        // 2. Search for 1st PE peak excluding pedestal
        double search_min = q0 + 0.3;
        if(v >= 2000) search_min = 2.0;
        
        h->GetXaxis()->SetRangeUser(search_min, x_max * 0.8);
        double q1_peak = h->GetBinCenter(h->GetMaximumBin());
        h->GetXaxis()->SetRange(0, 0); // reset
        
        double q1_guess = q1_peak - q0;
        if (q1_guess < 0.1) q1_guess = 0.3;

        // 3. Main Fit with constraints
        TF1 *f_spe = new TF1("f_spe", spe_model_final, -0.5, q0 + q1_guess * 5.0, 7);
        f_spe->SetParNames("Norm", "mu", "Q_ped", "s_ped", "Q_spe", "s_spe", "BG");
        
        double norm_init = h->Integral() * h->GetBinWidth(1);
        f_spe->SetParameters(norm_init, 0.5, q0, s0, q1_guess, q1_guess*0.4, 1.0);
        
        f_spe->SetParLimits(1, 0.01, 4.0); // mu
        f_spe->SetParLimits(2, q0-0.15, q0+0.15); // Q_ped
        f_spe->SetParLimits(3, 0.01, 0.2); // s_ped
        f_spe->SetParLimits(4, q1_guess*0.7, q1_guess*1.3); // Q_spe
        f_spe->SetParLimits(5, q1_guess*0.1, q1_guess*1.0); // s_spe
        f_spe->SetParLimits(6, 0.0, h->GetMaximum());

        h->Fit(f_spe, "RQ");

        double q_spe = f_spe->GetParameter(4);
        double q_spe_err = f_spe->GetParError(4);
        double gain = q_spe * 1e-12 / 1.602e-19;
        double gain_err = q_spe_err * 1e-12 / 1.602e-19;
        
        int n = gr_gain->GetN();
        gr_gain->SetPoint(n, (double)v, gain);
        gr_gain->SetPointError(n, 0, gain_err);

        TCanvas *c = new TCanvas(); c->SetLogy(); h->Draw(); f_spe->Draw("same");
        c->SaveAs(TString::Format("fit_res_%dV.png", v));
        delete c; delete h;
    }
    if (gr_gain->GetN() == 0) {
        std::cerr << "No valid gain points. Skipping gain curve plot." << std::endl;
        return;
    }
    TCanvas *cg = new TCanvas(); cg->SetLogy(); gr_gain->Draw("AP");
    cg->SaveAs("gain_curve_final.png");
}
