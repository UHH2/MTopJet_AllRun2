#include "../include/CentralInclude.h"


using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);


int main(int argc, char* argv[]){
  TFile* fileELEC = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* fileMUON = new TFile(dir_elec+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  vector<TString> dirnames = {"WeightHists01_noSF", "WeightHists02_PU", "WeightHists03_Lepton", "WeightHists04_BTag"};
  vector<TH1F*> h_muon, h_elec;
  for(auto name: dirnames){
    h_muon.push_back((TH1F*)fileMUON->Get(name+"/weight_1"));
    h_elec.push_back((TH1F*)fileELEC->Get(name+"/weight_1"));
  }

  for(auto h: h_muon){
    double integral = h->Integral();
    h->Scale(1/integral);
    h->SetLineWidth(4);
    h->SetLineColor(kAzure+7);
  }
  for(auto h: h_elec){
    double integral = h->Integral();
    h->Scale(1/integral);
    h->SetLineWidth(4);
    h->SetLineColor(kRed+1);
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<h_muon.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_muon[i]->GetXaxis()->SetRangeUser(0, 0.5);
    h_muon[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_muon[i]->GetMaximum());
    h_muon[i]->GetXaxis()->SetNdivisions(505);
    h_muon[i]->GetYaxis()->SetNdivisions(505);
    h_muon[i]->GetYaxis()->SetTitleOffset(1.6);
    h_muon[i]->GetXaxis()->SetTitleOffset(1.3);
    h_muon[i]->SetTitle(" ");
    h_muon[i]->GetYaxis()->SetTitle("events");
    h_muon[i]->Draw("HIST");
    h_elec[i]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.53,0.63,0.88,0.88);
    leg->AddEntry(h_muon[i], "muon channel", "l");
    leg->AddEntry(h_elec[i], "elec channel", "l");
    leg->Draw();
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/WeightHists/"+dirnames[i]+".pdf");
    delete a;
  }

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
