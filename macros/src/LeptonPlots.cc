#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(){
  TString save_path = get_save_path();

  TString uhh = "2016/muon/uhh2.AnalysisModuleRunner.MC.";
  vector<TString> name = {"nominal", "1695", "1715", "1735", "1755"};
  vector<int> colors = {kRed, kBlue, kGreen, kYellow, kBlack};

  TFile *f_tt = new TFile(dir+uhh+"TTbar.root");
  TFile *f_1695 = new TFile(dir+uhh+"TTbar_mtop1695_2016v3.root");
  TFile *f_1715 = new TFile(dir+uhh+"TTbar_mtop1715_2016v3.root");
  TFile *f_1735 = new TFile(dir+uhh+"TTbar_mtop1735_2016v3.root");
  TFile *f_1755 = new TFile(dir+uhh+"TTbar_mtop1755_2016v3.root");
  vector<TFile*> files = {f_tt, f_1695, f_1715, f_1735, f_1755};
  TString hist_perp = "LeptonicTop/muon_momentum_perp_full";
  TString hist_boost = "LeptonicTop/boosted_muon_pt_full";
  TString hist_boost_perp = "LeptonicTop/boosted_muon_momentum_perp_full";
  vector<TH1F*> pt_perp, pt_boost, pt_perp_norm, pt_boost_norm;
  vector<TH1F*> pt_perp_rebin_10, pt_boost_rebin_10, pt_perp_rebin_10_norm, pt_boost_rebin_10_norm;
  vector<TH1F*> momentum_perp_boosted, momentum_perp_boosted_norm;
  vector<TH1F*> momentum_perp_boosted_rebin_10, momentum_perp_boosted_norm_rebin_10;

  TH1F *hist, *hist_norm;
  double integral;


  for(unsigned int i=0; i<files.size(); i++){
    // --- Boost ----------------------------
    hist = (TH1F*)files[i]->Get(hist_boost);
    pt_boost.push_back(hist);
    // --- Boost Norm -----------------------
    hist_norm = normalize(hist);
    pt_boost_norm.push_back(hist_norm);
    // --- Rebin 10 -------------------------
    hist->Rebin(10);
    pt_boost_rebin_10.push_back(hist);
    // --- R10 Norm -------------------------
    hist_norm = normalize(hist);
    pt_boost_rebin_10_norm.push_back(hist_norm);

    // --------------------------------------
    // --- Perp boost -----------------------
    hist = (TH1F*)files[i]->Get(hist_boost_perp);
    momentum_perp_boosted.push_back(hist);
    // --- Perp Norm boost ------------------
    hist_norm = normalize(hist);
    momentum_perp_boosted_norm.push_back(hist_norm);
    // --- Rebin 10 -------------------------
    hist->Rebin(10);
    momentum_perp_boosted_rebin_10.push_back(hist);
    // --- R10 Norm -------------------------
    hist_norm = normalize(hist);
    momentum_perp_boosted_norm_rebin_10.push_back(hist_norm);

    // --------------------------------------
    // --- Perp -----------------------------
    hist = (TH1F*)files[i]->Get(hist_perp);
    pt_perp.push_back(hist);
    // --- Perp Norm ------------------------
    hist_norm = normalize(hist);
    pt_perp_norm.push_back(hist_norm);
    // --- Rebin 10 -------------------------
    hist->Rebin(10);
    pt_perp_rebin_10.push_back(hist);
    // --- R10 Norm -------------------------
    hist_norm = normalize(hist);
    pt_perp_rebin_10_norm.push_back(hist_norm);
  };


  vector<vector<TH1F*>> all_vector, all_vector_boost_perp, all_vector_rebin10;
  all_vector = {pt_perp, pt_boost, pt_perp_norm, pt_boost_norm,};
  vector<TString> all_vector_name = {"perpendicular_momentum_to_top", "pt_boost", "perp_momentum_to_top_norm", "pt_boost_norm"};

  all_vector_boost_perp = {momentum_perp_boosted, momentum_perp_boosted_norm, momentum_perp_boosted_rebin_10, momentum_perp_boosted_norm_rebin_10};
  vector<TString> all_vector_boost_perp_name = {"boosted_perpendicular_momentum", "boosted_perpendicular_momentum_norm", "boosted_perpendicular_momentum_rebin_10", "boosted_perpendicular_momentum_norm_rebin_10"};

  all_vector_rebin10 = {pt_perp_rebin_10, pt_boost_rebin_10, pt_perp_rebin_10_norm, pt_boost_rebin_10_norm};
  vector<TString> all_vector_rebin10_name = {"perpendicular_momentum_to_top_rebin_10", "pt_boost_rebin_10", "perpendicular_momentum_to_top_rebin_10_norm", "pt_boost_rebin_10_norm"};

  /*
  ██████  ██       ██████  ████████
  ██   ██ ██      ██    ██    ██
  ██████  ██      ██    ██    ██
  ██      ██      ██    ██    ██
  ██      ███████  ██████     ██
  */

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  TLegend *leg;

  for(unsigned int i=0; i<all_vector.size(); i++){
    all_vector[i].at(0)->SetTitle("");
    if(i == 1 || i == 3) all_vector[i].at(0)->GetXaxis()->SetRangeUser(0, 1000);
    if(i == 0 || i == 2) all_vector[i].at(0)->GetXaxis()->SetRangeUser(0, 100);
    all_vector[i].at(0)->GetYaxis()->SetRangeUser(0, get_highest_peak(all_vector[i])*1.2);
    all_vector[i].at(0)->GetXaxis()->SetNdivisions(505);
    all_vector[i].at(0)->GetYaxis()->SetNdivisions(505);
    all_vector[i].at(0)->GetXaxis()->SetTitleSize(0.05);
    all_vector[i].at(0)->GetYaxis()->SetTitleSize(0.04);
    all_vector[i].at(0)->GetXaxis()->SetTitleOffset(0.9);
    all_vector[i].at(0)->GetYaxis()->SetTitleOffset(1.5);
    if(i == 1 || i == 3) all_vector[i].at(0)->GetXaxis()->SetTitle("p_{T}^{muon*}");
    if(i == 0 || i == 2) all_vector[i].at(0)->GetXaxis()->SetTitle("#||{#vec{p}_{#perp}^{muon}}");
    all_vector[i].at(0)->GetYaxis()->SetTitle("");

    for(unsigned int j=0; j<files.size(); j++){
      all_vector[i].at(j)->SetLineWidth(2);
      all_vector[i].at(j)->SetLineColor(colors[j]);
    }

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    all_vector[i].at(0)->Draw("HIST");
    all_vector[i].at(1)->Draw("HIST SAME");
    all_vector[i].at(4)->Draw("HIST SAME");
    leg = new TLegend(0.65,0.65,0.85,0.85);
    leg->AddEntry(all_vector[i].at(0),"nominal 172.5","l");
    leg->AddEntry(all_vector[i].at(1),"mTop 169.5","l");
    leg->AddEntry(all_vector[i].at(4),"mTop 175.5","l");
    leg->SetTextSize(0.03);
    leg->Draw();
    gPad->RedrawAxis();
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_name[i]+"_main_2016.pdf");
    all_vector[i].at(2)->Draw("SAME HIST");
    all_vector[i].at(3)->Draw("SAME HIST");
    leg->AddEntry(all_vector[i].at(2),"mTop 171.5","l");
    leg->AddEntry(all_vector[i].at(3),"mTop 173.5","l");
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_name[i]+"_full_2016.pdf");
    delete A;
    leg->Clear();
  }

  /*
  ███    ███  ██████  ███    ███ ███████ ███    ██ ████████ ██    ██ ███    ███
  ████  ████ ██    ██ ████  ████ ██      ████   ██    ██    ██    ██ ████  ████
  ██ ████ ██ ██    ██ ██ ████ ██ █████   ██ ██  ██    ██    ██    ██ ██ ████ ██
  ██  ██  ██ ██    ██ ██  ██  ██ ██      ██  ██ ██    ██    ██    ██ ██  ██  ██
  ██      ██  ██████  ██      ██ ███████ ██   ████    ██     ██████  ██      ██
  */

  for(unsigned int i=0; i<all_vector_boost_perp.size(); i++){
    all_vector_boost_perp[i].at(0)->SetTitle("");
    all_vector_boost_perp[i].at(0)->GetXaxis()->SetRangeUser(0, 100);
    all_vector_boost_perp[i].at(0)->GetYaxis()->SetRangeUser(0, get_highest_peak(all_vector_boost_perp[i])*1.2);
    all_vector_boost_perp[i].at(0)->GetXaxis()->SetNdivisions(505);
    all_vector_boost_perp[i].at(0)->GetYaxis()->SetNdivisions(505);
    all_vector_boost_perp[i].at(0)->GetXaxis()->SetTitleSize(0.05);
    all_vector_boost_perp[i].at(0)->GetYaxis()->SetTitleSize(0.04);
    all_vector_boost_perp[i].at(0)->GetXaxis()->SetTitleOffset(0.9);
    all_vector_boost_perp[i].at(0)->GetYaxis()->SetTitleOffset(1.5);
    if(i == 1 || i == 3) all_vector_boost_perp[i].at(0)->GetXaxis()->SetTitle("#||{#vec{p}_{#perp}^{muon*}}");
    if(i == 0 || i == 2) all_vector_boost_perp[i].at(0)->GetXaxis()->SetTitle("#||{#vec{p}_{#perp}^{muon}}");
    all_vector_boost_perp[i].at(0)->GetYaxis()->SetTitle("");

    for(unsigned int j=0; j<files.size(); j++){
      all_vector_boost_perp[i].at(j)->SetLineWidth(2);
      all_vector_boost_perp[i].at(j)->SetLineColor(colors[j]);
    }

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    all_vector_boost_perp[i].at(0)->Draw("HIST");
    all_vector_boost_perp[i].at(1)->Draw("HIST SAME");
    all_vector_boost_perp[i].at(4)->Draw("HIST SAME");
    leg = new TLegend(0.65,0.65,0.85,0.85);
    leg->AddEntry(all_vector_boost_perp[i].at(0),"nominal 172.5","l");
    leg->AddEntry(all_vector_boost_perp[i].at(1),"mTop 169.5","l");
    leg->AddEntry(all_vector_boost_perp[i].at(4),"mTop 175.5","l");
    leg->SetTextSize(0.03);
    leg->Draw();
    gPad->RedrawAxis();
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_boost_perp_name[i]+"_main_2016.pdf");
    all_vector_boost_perp[i].at(2)->Draw("SAME HIST");
    all_vector_boost_perp[i].at(3)->Draw("SAME HIST");
    leg->AddEntry(all_vector_boost_perp[i].at(2),"mTop 171.5","l");
    leg->AddEntry(all_vector_boost_perp[i].at(3),"mTop 173.5","l");
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_boost_perp_name[i]+"_full_2016.pdf");
    delete A;
    leg->Clear();
  }

  /*
  ██████  ███████ ██████  ██ ███    ██
  ██   ██ ██      ██   ██ ██ ████   ██
  ██████  █████   ██████  ██ ██ ██  ██
  ██   ██ ██      ██   ██ ██ ██  ██ ██
  ██   ██ ███████ ██████  ██ ██   ████
  */

  for(unsigned int i=0; i<all_vector_rebin10.size(); i++){
    all_vector_rebin10[i].at(0)->SetTitle("");
    if(i == 1 || i == 3) all_vector_rebin10[i].at(0)->GetXaxis()->SetRangeUser(0, 1000);
    if(i == 0 || i == 2) all_vector_rebin10[i].at(0)->GetXaxis()->SetRangeUser(0, 100);
    all_vector_rebin10[i].at(0)->GetYaxis()->SetRangeUser(0, get_highest_peak(all_vector_rebin10[i])*1.2);
    all_vector_rebin10[i].at(0)->GetXaxis()->SetNdivisions(505);
    all_vector_rebin10[i].at(0)->GetYaxis()->SetNdivisions(505);
    all_vector_rebin10[i].at(0)->GetXaxis()->SetTitleSize(0.05);
    all_vector_rebin10[i].at(0)->GetYaxis()->SetTitleSize(0.04);
    all_vector_rebin10[i].at(0)->GetXaxis()->SetTitleOffset(0.9);
    all_vector_rebin10[i].at(0)->GetYaxis()->SetTitleOffset(1.5);
    if(i == 1 || i == 3) all_vector_rebin10[i].at(0)->GetXaxis()->SetTitle("p_{T}^{muon*}");
    if(i == 0 || i == 2) all_vector_rebin10[i].at(0)->GetXaxis()->SetTitle("#||{#vec{p}_{#perp}^{muon}}");
    all_vector_rebin10[i].at(0)->GetYaxis()->SetTitle("");

    for(unsigned int j=0; j<files.size(); j++){
      all_vector_rebin10[i].at(j)->SetLineWidth(2);
      all_vector_rebin10[i].at(j)->SetLineColor(colors[j]);
    }

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    all_vector_rebin10[i].at(0)->Draw("HIST");
    all_vector_rebin10[i].at(1)->Draw("HIST SAME");
    all_vector_rebin10[i].at(4)->Draw("HIST SAME");
    leg = new TLegend(0.65,0.65,0.85,0.85);
    leg->AddEntry(all_vector_rebin10[i].at(0),"nominal 172.5","l");
    leg->AddEntry(all_vector_rebin10[i].at(1),"mTop 169.5","l");
    leg->AddEntry(all_vector_rebin10[i].at(4),"mTop 175.5","l");
    leg->SetTextSize(0.03);
    leg->Draw();
    gPad->RedrawAxis();
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_rebin10_name[i]+"_main_2016.pdf");
    all_vector_rebin10[i].at(2)->Draw("SAME HIST");
    all_vector_rebin10[i].at(3)->Draw("SAME HIST");
    leg->AddEntry(all_vector_rebin10[i].at(2),"mTop 171.5","l");
    leg->AddEntry(all_vector_rebin10[i].at(3),"mTop 173.5","l");
    A->SaveAs(save_path+"/Plots/LeptonPlots/"+all_vector_rebin10_name[i]+"_full_2016.pdf");
    delete A;
    leg->Clear();
  }

}
