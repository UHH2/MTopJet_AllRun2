#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  vector<TString> year = {"2016", "2017", "2018"};
  vector<TFile*> file_nominal, file_JECup, file_JECdown;

  TFile *file_filing;

  for(unsigned int i=0; i<year.size(); i++){
    file_filing = new TFile(dir+year[i]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    file_nominal.push_back(file_filing);
    file_filing = new TFile(dir+year[i]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
    file_JECup.push_back(file_filing);
    file_filing = new TFile(dir+year[i]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");
    file_JECdown.push_back(file_filing);
  }

  std::vector<int> color = {kRed, kBlue, kGreen};
  std::vector<TString> mass_name = {"nominal", "JECup", "JECdown"};

  vector<TH1F*> mass_nominal, mass_JECup, mass_JECdown;
  TH1F *mass_filing;

  for(unsigned int i=0; i<year.size(); i++){
    mass_filing = (TH1F*)file_nominal[i]->Get("XCone_cor/M_jet1");
    mass_nominal.push_back(mass_filing);
    mass_filing = (TH1F*)file_JECup[i]->Get("XCone_cor/M_jet1");
    mass_JECup.push_back(mass_filing);
    mass_filing = (TH1F*)file_JECdown[i]->Get("XCone_cor/M_jet1");
    mass_JECdown.push_back(mass_filing);
  }

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

  for(unsigned int i=0; i<year.size(); i++){
    mass_nominal[i]->GetXaxis()->SetRangeUser(50, 300);
    mass_nominal[i]->GetYaxis()->SetRangeUser(0, (mass_nominal[i]->GetMaximum())*1.2);
    mass_nominal[i]->GetXaxis()->SetNdivisions(505);
    mass_nominal[i]->GetYaxis()->SetNdivisions(505);
    mass_nominal[i]->GetXaxis()->SetTitleSize(0.05);
    mass_nominal[i]->GetYaxis()->SetTitleSize(0.04);
    mass_nominal[i]->GetXaxis()->SetTitleOffset(0.9);
    mass_nominal[i]->GetYaxis()->SetTitleOffset(1.5);
    mass_nominal[i]->GetXaxis()->SetTitle("M_{jet}");
    mass_nominal[i]->GetYaxis()->SetTitle("Events");
    mass_nominal[i]->SetLineWidth(4);
    mass_nominal[i]->SetLineColor(kRed);
    mass_JECup[i]->SetLineWidth(4);
    mass_JECup[i]->SetLineColor(kBlue);
    mass_JECdown[i]->SetLineWidth(4);
    mass_JECdown[i]->SetLineColor(kGreen);

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    TCanvas *c = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    mass_nominal[i]->Draw("HIST");
    mass_nominal[i]->Draw("SAME HIST");
    mass_JECup[i]->Draw("SAME HIST");
    mass_JECdown[i]->Draw("SAME HIST");
    TLegend *leg = new TLegend(0.2,0.60,0.4,0.80);
    leg->AddEntry(mass_nominal[i],"nominal","l");
    leg->AddEntry(mass_JECup[i],"JEC up","l");
    leg->AddEntry(mass_JECdown[i],"JEC down","l");
    leg->SetTextSize(0.05);
    leg->Draw();
    gPad->RedrawAxis();
    c->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/MassComparison_"+year[i]+".pdf");
  }

}
