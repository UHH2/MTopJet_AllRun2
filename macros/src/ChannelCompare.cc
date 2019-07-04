#include "../include/CentralInclude.h"

using namespace std;

void Plot(TH1F* elec, TH1F* elec_stat, TH1F* muon, TH1F* muon_stat);
/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/


int main(int argc, char* argv[]){
  TFile* f_elec = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Results_data_elec.root");
  TFile* f_muon = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Results_data_muon.root");

  TH1F * elec = (TH1F*)f_elec->Get("Unfold_XS_totuncert");
  TH1F * elec_stat = (TH1F*)f_elec->Get("Unfold_XS_statuncert");
  TH1F * muon = (TH1F*)f_muon->Get("Unfold_XS_totuncert");
  TH1F * muon_stat = (TH1F*)f_muon->Get("Unfold_XS_statuncert");
  Plot(elec, elec_stat, muon, muon_stat);
}

void Plot(TH1F* elec, TH1F* elec_stat, TH1F* muon, TH1F* muon_stat){

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  double max = elec->GetMaximum();
  double ymax = 1.5 * max;

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  elec->SetTitle(" ");
  elec->GetYaxis()->SetRangeUser(0., ymax);
  elec->GetXaxis()->SetTitle("Leading-jet mass [GeV]");
  elec->GetYaxis()->SetTitle("#frac{d#sigma}{dm_{jet}} [#frac{fb}{GeV}]");
  elec->GetYaxis()->SetTitleOffset(1.1);
  elec->GetXaxis()->SetTitleOffset(0.9);
  elec->GetYaxis()->SetTitleSize(0.05);
  elec->GetXaxis()->SetTitleSize(0.05);
  elec->GetYaxis()->SetNdivisions(505);
  elec->SetLineColor(17);
  elec->SetFillColor(17);
  elec->SetMarkerColor(kBlack);
  elec->SetMarkerStyle(8);
  elec->SetMarkerSize(1);
  elec->Draw("E2");
  elec_stat->SetFillColor(15);
  elec_stat->SetMarkerColor(kBlack);
  elec_stat->SetMarkerStyle(8);
  elec_stat->SetMarkerSize(1);
  elec_stat->Draw("E2 SAME");

  muon->SetLineColor(kRed+1);
  muon->SetMarkerColor(kRed+1);
  muon->SetMarkerStyle(8);
  muon->SetMarkerSize(1);
  muon->Draw("E1 SAME");
  muon_stat->SetLineColor(kRed+1);
  muon_stat->SetMarkerColor(kRed+1);
  muon_stat->SetMarkerStyle(8);
  muon_stat->SetMarkerSize(1);
  muon_stat->Draw("E1 SAME");
  muon->Draw("E1 SAME");

  gStyle->SetEndErrorSize(5);

  CMSLabel(true, 0.2, 0.85);

  TLegend *l=new TLegend(0.55,0.68,0.85,0.87);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(elec,"elec channel","pf");
  l->AddEntry(muon,"muon channel","ple");
  l->SetTextSize(0.03);
  l->Draw();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/ChannelComparison.pdf");
  delete c;
}
