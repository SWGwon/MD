#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

// SPE 분포 모델 (Poisson + Gaussians + Background)
Double_t spe_model_full(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t N = par[0], mu = par[1], Q0 = par[2], s0 = par[3], Q1 = par[4], s1 = par[5], bg = par[6];
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

void fit_channels_dark_led() {
    gStyle->SetOptFit(0); 
    std::vector<int> voltages = {1600, 1700, 1800, 1900, 2000, 2100, 2200};
    int channels[] = {0, 1};
    Color_t colors[] = {kRed, kBlue};
    
    TGraphErrors *gr_gains[2];
    for(int i=0; i<2; i++) {
        gr_gains[i] = new TGraphErrors();
        gr_gains[i]->SetName(TString::Format("gr_gain_full_ch%d", i));
        gr_gains[i]->SetMarkerStyle(20 + i);
        gr_gains[i]->SetMarkerColor(colors[i]);
        gr_gains[i]->SetLineColor(colors[i]);
    }

    for (int ch : channels) {
        std::cout << "\n>>> Processing Channel " << ch << " (Dark + LED) <<<" << std::endl;
        
        for (int v : voltages) {
            TString led_fn = TString::Format("output_charge_run2_%dV.root", v);
            TString dark_fn = TString::Format("output_charge_run2_dark_%dV.root", v);
            
            TFile *f_led = TFile::Open(led_fn);
            if (!f_led || f_led->IsZombie()) continue;
            TTree *t_led = (TTree*)f_led->Get("T_Charge");
            
            double x_max = (v < 2000) ? 8.0 : 45.0;
            TH1D *h_led = new TH1D("h_led", TString::Format("Ch%d %dV Fit (Dark+LED);Charge [pC];Counts", ch, v), 600, -1.0, x_max);
            t_led->Draw(TString::Format("charge_pC[%d]>>h_led", ch), "", "goff");

            // 1. Pedestal 위치 추정 (Dark 데이터 우선 사용)
            double q0 = 0.0, s0 = 0.05;
            TFile *f_dark = TFile::Open(dark_fn);
            if (f_dark && !f_dark->IsZombie()) {
                TTree *t_dark = (TTree*)f_dark->Get("T_Charge");
                if (t_dark) {
                    TH1D *h_dark = new TH1D("h_dark", "Dark", 300, -0.5, 0.5);
                    t_dark->Draw(TString::Format("charge_pC[%d]>>h_dark", ch), "", "goff");
                    TF1 *f_ped_dark = new TF1("f_ped_dark", "gaus", -0.2, 0.2);
                    h_dark->Fit(f_ped_dark, "RQN");
                    q0 = f_ped_dark->GetParameter(1);
                    s0 = f_ped_dark->GetParameter(2);
                    delete h_dark; delete f_ped_dark;
                }
                f_dark->Close();
            } else {
                // Dark 파일이 없는 경우(예: 1900V) LED에서 직접 추정
                h_led->GetXaxis()->SetRangeUser(-0.3, 0.4);
                q0 = h_led->GetBinCenter(h_led->GetMaximumBin());
                h_led->GetXaxis()->SetRange(0, 0);
            }

            // 2. 1-PE 피크 위치 자동 탐색
            double search_min = q0 + 0.3;
            if(v >= 2000) search_min = 2.0;
            h_led->GetXaxis()->SetRangeUser(search_min, x_max * 0.8);
            double q1_guess = h_led->GetBinCenter(h_led->GetMaximumBin()) - q0;
            if (q1_guess < 0.1) q1_guess = 0.3;
            h_led->GetXaxis()->SetRange(0, 0);

            // 3. 메인 SPE 피팅
            TF1 *f_spe = new TF1("f_spe", spe_model_full, -0.6, q0 + q1_guess * 5.0, 7);
            f_spe->SetParameters(h_led->Integral()*h_led->GetBinWidth(1), 0.5, q0, s0, q1_guess, q1_guess*0.4, 1.0);
            f_spe->SetParLimits(1, 0.01, 4.0);
            f_spe->SetParLimits(2, q0-0.1, q0+0.1); // Pedestal 고정력 강화
            f_spe->SetParLimits(4, q1_guess*0.6, q1_guess*1.4);
            f_spe->SetLineColor(colors[ch]);

            h_led->Fit(f_spe, "RQ");

            // Gain 계산 및 저장
            double q_spe = f_spe->GetParameter(4);
            double q_spe_err = f_spe->GetParError(4);
            double gain = q_spe * 1e-12 / 1.602e-19;
            double gain_err = q_spe_err * 1e-12 / 1.602e-19;
            
            int n = gr_gains[ch]->GetN();
            gr_gains[ch]->SetPoint(n, (double)v, gain);
            gr_gains[ch]->SetPointError(n, 0, gain_err);

            // 개별 피팅 결과 저장
            TCanvas *c = new TCanvas("c", "", 800, 600);
            c->SetLogy();
            h_led->Draw();
            f_spe->Draw("same");
            c->SaveAs(TString::Format("fit_full_ch%d_%dV.png", ch, v));
            
            delete c; delete h_led; f_led->Close();
        }
    }

    // --- 두 채널 통합 Gain Curve 그리기 ---
    TCanvas *cg = new TCanvas("cg", "Gain Comparison (Dark+LED)", 900, 700);
    cg->SetLogy();
    cg->SetGrid();

    double y_min = 5e5, y_max = 5e8;
    TH1F *frame = cg->DrawFrame(1550, y_min, 2250, y_max);
    frame->SetTitle("Gain vs Voltage Comparison (Dark+LED);Voltage [V];Gain");

    TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg->SetBorderSize(1);

    for (int i = 0; i < 2; i++) {
        gr_gains[i]->Draw("P same");
        TF1 *f_gain = new TF1(TString::Format("f_gain_full_ch%d", i), "[0]*TMath::Power(x, [1])", 1550, 2250);
        f_gain->SetLineColor(colors[i]);
        f_gain->SetParameters(1e-25, 10.0);
        gr_gains[i]->Fit(f_gain, "RQ");
        
        // 범례에서 k값 제외, 채널 이름만 표시
        leg->AddEntry(gr_gains[i], TString::Format("Channel %d", i), "LP");
    }
    
    leg->Draw();
    cg->SaveAs("gain_comp_full_ch0_ch1.png");
    
    std::cout << "\n*** Analysis Complete. Check 'gain_comp_full_ch0_ch1.png' ***" << std::endl;
}
