#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  vector<TString> mass_sample = {"1695", "1715", "1735", "1735"};

  TFile *file_mass_1695 = new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695_2016v3.root");
  TFile *file_mass_1715 = new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1715_2016v3.root");
  TFile *file_mass_1735 = new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1735_2016v3.root");
  TFile *file_mass_1755 = new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755_2016v3.root");

  TH1F *hist_1695 = (TH1F*)file_mass_1695->Get("LeptonicTop/MuonPtPerpFull");
  TH1F *hist_1715 = (TH1F*)file_mass_1715->Get("LeptonicTop/MuonPtPerpFull");
  TH1F *hist_1735 = (TH1F*)file_mass_1735->Get("LeptonicTop/MuonPtPerpFull");
  TH1F *hist_1755 = (TH1F*)file_mass_1755->Get("LeptonicTop/MuonPtPerpFull");


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

  hist_1755->GetXaxis()->SetRangeUser(0, 120);
  hist_1755->GetYaxis()->SetRangeUser(0, (hist_1755->GetMaximum())*1.2);
  hist_1755->GetXaxis()->SetNdivisions(505);
  hist_1755->GetYaxis()->SetNdivisions(505);
  hist_1755->GetXaxis()->SetTitleSize(0.05);
  hist_1755->GetYaxis()->SetTitleSize(0.04);
  hist_1755->GetXaxis()->SetTitleOffset(0.9);
  hist_1755->GetYaxis()->SetTitleOffset(1.5);
  hist_1755->GetXaxis()->SetTitle("P_{T,#perp}");
  hist_1755->GetYaxis()->SetTitle("Events");
  hist_1755->SetLineWidth(3);
  hist_1755->SetLineColor(kRed);
  hist_1735->SetLineWidth(3);
  hist_1735->SetLineColor(kBlue);
  hist_1715->SetLineWidth(3);
  hist_1715->SetLineColor(kGreen);
  hist_1695->SetLineWidth(3);
  hist_1695->SetLineColor(kBlack);

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *c = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hist_1755->Draw("HIST");
  hist_1695->Draw("SAME HIST");
  hist_1735->Draw("SAME HIST");
  hist_1715->Draw("SAME HIST");
  TLegend *leg = new TLegend(0.7,0.60,0.9,0.80);
  leg->AddEntry(hist_1755,"1755","l");
  leg->AddEntry(hist_1735,"1735","l");
  leg->AddEntry(hist_1715,"1715","l");
  leg->AddEntry(hist_1695,"1695","l");
  leg->SetTextSize(0.05);
  leg->Draw();
  gPad->RedrawAxis();
  c->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_2016.pdf");


}
