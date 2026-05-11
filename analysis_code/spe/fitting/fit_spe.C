#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

// Model for SPE distribution
// x: charge [pC]
// par[0]: N (Normalization)
// par[1]: mu (Mean of Poisson distribution)
// par[2]: Q0 (Pedestal position)
// par[3]: sigma0 (Pedestal width)
// par[4]: Q1 (Single Photoelectron charge)
// par[5]: sigma1 (Single Photoelectron width)
Double_t spe_model(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t N = par[0];
    Double_t mu = par[1];
    Double_t Q0 = par[2];
    Double_t sigma0 = par[3];
    Double_t Q1 = par[4];
    Double_t sigma1 = par[5];

    Double_t sum = 0.0;
    // Sum over multiple photoelectron peaks (n=0: pedestal, n=1: 1pe, ...)
    for (int n = 0; n <= 10; n++) {
        Double_t P_n = TMath::Poisson(n, mu);
        Double_t Q_n = Q0 + n * Q1;
        Double_t sigma_n = TMath::Sqrt(sigma0 * sigma0 + n * sigma1 * sigma1);
        Double_t G_n = (1.0 / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma_n)) * TMath::Exp(-0.5 * TMath::Power((xx - Q_n) / sigma_n, 2));
        sum += P_n * G_n;
    }
    return N * sum;
}

void fit_spe(int voltage = 1600, int ch = 0) {
    gStyle->SetOptFit(1111);
    
    TString led_filename = TString::Format("output_charge_run2_%dV.root", voltage);
    TString dark_filename = TString::Format("output_charge_run2_dark_%dV.root", voltage);
    
    TFile *f_led = TFile::Open(led_filename);
    if (!f_led || f_led->IsZombie()) {
        std::cerr << "Cannot open " << led_filename << std::endl;
        return;
    }
    
    TH1D *h_led = (TH1D*)f_led->Get(TString::Format("hCharge_Ch%d", ch));
    if (!h_led) {
        std::cerr << "Cannot find hCharge_Ch" << ch << " in LED file" << std::endl;
        return;
    }
    
    // Initial pedestal parameters from dark run if available
    double q0_init = 0.0;
    double s0_init = 0.01;
    
    TFile *f_dark = TFile::Open(dark_filename);
    if (f_dark && !f_dark->IsZombie()) {
        TH1D *h_dark = (TH1D*)f_dark->Get(TString::Format("hCharge_Ch%d", ch));
        if (h_dark) {
            // Fit dark data with a simple Gaussian for pedestal
            TF1 *f_ped = new TF1("f_ped", "gaus", -0.2, 0.2);
            h_dark->Fit(f_ped, "QN");
            q0_init = f_ped->GetParameter(1);
            s0_init = f_ped->GetParameter(2);
        }
        f_dark->Close();
    } else {
        std::cout << "Warning: Dark file for " << voltage << "V not found. Using defaults for pedestal." << std::endl;
        // Simple heuristic for pedestal if dark run is missing
        q0_init = h_led->GetBinCenter(h_led->GetMaximumBin());
        s0_init = 0.015;
    }

    TCanvas *c1 = new TCanvas("c1", "SPE Fit", 800, 600);
    h_led->Draw();
    
    // Estimate mu from the ratio of pedestal counts (Poisson P(0) = e^-mu)
    // N_ped = Total * e^-mu => mu = -ln(N_ped/Total)
    // Here we'll just start with mu = 0.5 as an initial guess.
    
    TF1 *f_spe = new TF1("f_spe", spe_model, h_led->GetXaxis()->GetXmin(), h_led->GetXaxis()->GetXmax(), 6);
    f_spe->SetParNames("Norm", "mu", "Q_ped", "sigma_ped", "Q_spe", "sigma_spe");
    
    double norm_init = h_led->Integral("width");
    double q1_init = (voltage - 1500.0) * 0.0005 + 0.05; // Very rough guess for 1pe charge (in pC)
    
    f_spe->SetParameters(norm_init, 0.5, q0_init, s0_init, q1_init, 0.02);
    f_spe->SetParLimits(1, 0.0, 5.0); // mu
    f_spe->SetParLimits(2, q0_init - 0.05, q0_init + 0.05); // Q0
    f_spe->SetParLimits(3, 0.001, 0.1); // sigma0
    f_spe->SetParLimits(4, 0.01, 1.0); // Q1
    f_spe->SetParLimits(5, 0.005, 0.1); // sigma1
    
    h_led->Fit(f_spe, "R");
    
    c1->SaveAs(TString::Format("fit_run2_%dV_ch%d.png", voltage, ch));
    
    double q_spe = f_spe->GetParameter(4);
    double gain = q_spe * 1e-12 / 1.602e-19;
    std::cout << "Voltage: " << voltage << "V" << std::endl;
    std::cout << "SPE Charge: " << q_spe << " pC" << std::endl;
    std::cout << "Estimated Gain: " << gain << std::endl;
}
