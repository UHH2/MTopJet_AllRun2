#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(){
  TString uhh = "2016/muon/uhh2.AnalysisModuleRunner.MC.";
  vector<TString> name = {"nominal", "1695", "1715", "1735", "1755"};
  vector<int> colors = {kRed, kBlue, kGreen, kYellow, kBlack};

  TFile *f_tt = new TFile(dir+uhh+"TTbar.root");
  TFile *f_1695 = new TFile(dir+uhh+"TTbar_mtop1695_2016v3.root");
  TFile *f_1715 = new TFile(dir+uhh+"TTbar_mtop1715_2016v3.root");
  TFile *f_1735 = new TFile(dir+uhh+"TTbar_mtop1735_2016v3.root");
  TFile *f_1755 = new TFile(dir+uhh+"TTbar_mtop1755_2016v3.root");
  vector<TFile*> files = {f_tt, f_1695, f_1715, f_1735, f_1755};
  vector<TString> hists_name = {"nominal 172.5", "mTop 1695", "mTop 171.5", "mTop 173.5", "mTop 175.5"};

  TString hist_jetmass = "XCone_cor/M_jet1";

  vector<TH1F*> hists;
  TH1F *hist;
  for(unsigned int i=0; i<files.size(); i++) hists.push_back((TH1F*)files[i]->Get(hist_jetmass));


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

  hists[0]->SetTitle("");
  hists[0]->GetXaxis()->SetRangeUser(0, 400);
  hists[0]->GetYaxis()->SetRangeUser(0, get_highest_peak(hists)*1.2);
  hists[0]->GetXaxis()->SetNdivisions(505);
  hists[0]->GetYaxis()->SetNdivisions(505);
  hists[0]->GetXaxis()->SetTitleSize(0.05);
  hists[0]->GetYaxis()->SetTitleSize(0.04);
  hists[0]->GetXaxis()->SetTitleOffset(0.9);
  hists[0]->GetYaxis()->SetTitleOffset(1.5);
  hists[0]->GetXaxis()->SetTitle("M_{jet}^{had}");
  hists[0]->GetYaxis()->SetTitle("");

  for(unsigned int i=0; i<hists.size(); i++){
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colors[i]);
  }

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hists[0]->Draw("HIST");
  hists[1]->Draw("HIST SAME");
  hists[4]->Draw("HIST SAME");
  leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->AddEntry(hists[0],"nominal 172.5","l");
  leg->AddEntry(hists[1],"mTop 169.5","l");
  leg->AddEntry(hists[4],"mTop 175.5","l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/MassSensitivity/HadjetMass_main_2016.pdf");
  hists[2]->Draw("SAME HIST");
  hists[3]->Draw("SAME HIST");
  leg->AddEntry(hists[2],"mTop 171.5","l");
  leg->AddEntry(hists[3],"mTop 173.5","l");
  A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/MassSensitivity/HadjetMass_full_2016.pdf");
  delete A;
  leg->Clear();

}
