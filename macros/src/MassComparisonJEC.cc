#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  vector<TString> year = {"2016", "2017", "2018"};
  TFile *file_nominal =  new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *file_JECup =  new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
  TFile *file_JECdown =  new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");

  TH1F *mass_nominal = (TH1F*)file_nominal->Get("XCone_cor/M_jet1");
  TH1F *mass_JECup = (TH1F*)file_JECup->Get("XCone_cor/M_jet1");
  TH1F *mass_JECdown = (TH1F*)file_JECdown->Get("XCone_cor/M_jet1");

  std::vector<int> color = {kRed, kBlue, kGreen};
  std::vector<TString> mass_name = {"nominal", "JECup", "JECdown"};

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

  mass_nominal->GetXaxis()->SetRangeUser(50, 300);
  mass_nominal->GetYaxis()->SetRangeUser(0, (mass_nominal->GetMaximum())*1.2);
  mass_nominal->GetXaxis()->SetNdivisions(505);
  mass_nominal->GetYaxis()->SetNdivisions(505);
  mass_nominal->GetXaxis()->SetTitleSize(0.05);
  mass_nominal->GetYaxis()->SetTitleSize(0.04);
  mass_nominal->GetXaxis()->SetTitleOffset(0.9);
  mass_nominal->GetYaxis()->SetTitleOffset(1.5);
  mass_nominal->GetXaxis()->SetTitle("M_{jet}");
  mass_nominal->GetYaxis()->SetTitle("Compare MassDistribution");
  mass_nominal->SetLineWidth(4);
  mass_nominal->SetLineColor(kRed);
  mass_JECup->SetLineWidth(4);
  mass_JECup->SetLineColor(kBlue);
  mass_JECdown->SetLineWidth(4);
  mass_JECdown->SetLineColor(kGreen);

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *c = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  mass_nominal->Draw("HIST");
  // zero_line->Draw("SAME");
  mass_nominal->Draw("SAME HIST");
  mass_JECup->Draw("SAME HIST");
  mass_JECdown->Draw("SAME HIST");
  TLegend *leg = new TLegend(0.05,0.05,0.25,0.25);
  leg->AddEntry(mass_nominal,"nominal","l");
  leg->AddEntry(mass_JECup,"JEC up","l");
  leg->AddEntry(mass_JECdown,"JEC down","l");
  leg->SetTextSize(0.05);
  leg->Draw();
  gPad->RedrawAxis();
  c->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/MassComparison_2016.pdf");


}
