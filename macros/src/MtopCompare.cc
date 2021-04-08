#include "../include/CentralInclude.h"

using namespace std;

TH1F* GetRun2Plot(TString histname, TString dir, TString prefix, TString process);
void Plot(TH1F* h_data, vector<TH1F*> h_ttbar);
TH1F* GetRatio(TH1F* h1, TH1F* h2);
TH1F* GetRatioUncert(TH1F* mc);
/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  TString directory = "/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/";
  TString prefix_mc = "uhh2.AnalysisModuleRunner.MC.";
  TString prefix_data = "uhh2.AnalysisModuleRunner.DATA.";

  vector<TString> backgrounds = {"SingleTop", "WJets", "other"};
  vector<TString> ttbars = {"TTbar_mtop1695", "TTbar_mtop1715", "TTbar", "TTbar_mtop1735", "TTbar_mtop1755"};
  TString histname = "JetMassScaleHists/hadjet_jms_mass";

  TH1F* h_data = GetRun2Plot(histname, directory, prefix_data, "DATA");
  vector<TH1F*> h_bkg, h_ttbar;
  for(auto bkg: backgrounds) h_bkg.push_back(GetRun2Plot(histname, directory, prefix_mc, bkg));
  for(auto tt: ttbars) h_ttbar.push_back(GetRun2Plot(histname, directory, prefix_mc, tt));

  // Subtract Background from data
  for(auto h: h_bkg) h_data->Add(h, -1);

  // Rebin
  int rebin = 2;
  h_data->Rebin(rebin);
  for(auto h: h_ttbar) h->Rebin(rebin);

  // Normalize
  h_data->Scale(1/h_data->Integral());
  for(auto h: h_ttbar) h->Scale(1/h->Integral());

  Plot(h_data, h_ttbar);

}

TH1F* GetRun2Plot(TString histname, TString dir, TString prefix, TString process){
  bool isfirst = true;
  vector<TString> channels = {"muon", "elec"};
  vector<TString> years = {"2016v3", "2017v2", "2018"};
  TH1F* h_all;
  for(auto channel: channels){
    for(auto year: years){
      TFile* file = new TFile(dir+"/"+channel+"/"+prefix+process+"_"+year+".root");
      TH1F* h = (TH1F*) file->Get(histname);
      if(isfirst){
        h_all = (TH1F*) h->Clone();
      }
      else{
        h_all->Add(h);
        isfirst = false;
      }
    }
  }
  return h_all;
}

void Plot(TH1F* h_data, vector<TH1F*> h_ttbar){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  TCanvas* c = new TCanvas("c","c",600,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  double xmin = 100;
  double ymax = 0.4;
  int titlefont = 43;

  h_data->SetTitle("");
  h_data->GetXaxis()->SetTitle("#it{m}_{jet}");
  h_data->GetYaxis()->SetTitle("a.u.");
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetNdivisions(505);
  h_data->GetXaxis()->SetRangeUser(xmin, 340);
  h_data->GetYaxis()->SetRangeUser(0, ymax);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerColor(kBlack);
  h_data->SetLineColor(kBlack);
  h_data->Draw("E1");
  Int_t color[] = {kRed-4, kRed-4, 14, kAzure-7, kAzure-7};
  Int_t style[] = {3, 1, 1, 1, 2};
  for(unsigned int i=0; i<h_ttbar.size();i++){
    h_ttbar[i]->SetLineWidth(3);
    h_ttbar[i]->SetLineColor(color[i]);
    h_ttbar[i]->SetLineStyle(style[i]);
    h_ttbar[i]->Draw("HIST SAME");
  }
  h_data->Draw("E1 SAME");
  vector<TString> legnames = {"#it{m}_{t} = 169.5 GeV", "#it{m}_{t} = 171.5 GeV", "#it{m}_{t} = 172.5 GeV", "#it{m}_{t} = 173.5 GeV","#it{m}_{t} = 175.5 GeV"};
  TLegend * leg = new TLegend(0.5, 0.5, 0.85, 0.85);
  leg->AddEntry(h_data, "Data - backgrounds", "pe");
  for(unsigned int i=0; i<h_ttbar.size();i++) leg->AddEntry(h_ttbar[i], legnames[i], "l");
  leg->Draw();

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  h_data->GetXaxis()->SetLabelSize(0.);
  h_data->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( xmin, 0, xmin, ymax, 0, ymax, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(titlefont);
  axis->SetLabelSize(21);
  axis->Draw();

  // Ratio Plot kommt in unteres pad
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.38);
  pad2->Draw();
  pad2->cd();

  TH1F * h_ratio_data = GetRatioUncert(h_data);
  vector<TH1F*> h_ratio_ttbar;
  for(auto h: h_ttbar) h_ratio_ttbar.push_back(GetRatio(h, h_data));
  h_ratio_data->GetXaxis()->SetTickLength(0.07);
  h_ratio_data->GetXaxis()->SetTitleSize(25);
  h_ratio_data->GetXaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetXaxis()->SetTitleOffset(4.0);
  h_ratio_data->GetXaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetXaxis()->SetLabelSize(21);
  h_ratio_data->GetXaxis()->SetLabelOffset(0.035);
  h_ratio_data->GetYaxis()->SetTitle("#frac{Theory}{Data}");
  h_ratio_data->GetYaxis()->CenterTitle();
  h_ratio_data->GetYaxis()->SetTitleSize(22);
  h_ratio_data->GetYaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetYaxis()->SetTitleOffset(2.2);
  h_ratio_data->GetYaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetYaxis()->SetLabelSize(19);
  h_ratio_data->GetYaxis()->SetLabelOffset(0.009);
  h_ratio_data->GetXaxis()->SetNdivisions(505);
  h_ratio_data->GetYaxis()->SetNdivisions(505);

  h_ratio_data->SetTitle(" ");
  h_ratio_data->GetYaxis()->SetRangeUser(0.4, 1.6);
  h_ratio_data->GetXaxis()->SetTitle("#it{m}_{jet} #left[GeV#right]");
  h_ratio_data->SetLineColor(kBlack);
  h_ratio_data->SetMarkerColor(kBlack);
  h_ratio_data->SetMarkerStyle(8);
  h_ratio_data->SetMarkerSize(1);
  h_ratio_data->Draw("E1");
  for(unsigned int i=0; i<h_ratio_ttbar.size();i++){
    h_ratio_ttbar[i]->SetLineWidth(3);
    h_ratio_ttbar[i]->SetLineColor(color[i]);
    h_ratio_ttbar[i]->SetLineStyle(style[i]);
    h_ratio_ttbar[i]->Draw("HIST SAME");
  }
  h_ratio_data->Draw("E1 SAME");

  gStyle->SetEndErrorSize(5);

  gPad->RedrawAxis();

  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MtopComparison/MtopComparison.pdf");
  delete c;
  return;
}

TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      // double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      double error = E1/N2; // only consider uncert from MC
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}

TH1F* GetRatioUncert(TH1F* mc){
  TH1F* mc_uncert = (TH1F*) mc->Clone();
  int Nbins = mc->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    mc_uncert->SetBinContent(i, 1.0);
    double central = mc->GetBinContent(i);
    double error = mc->GetBinError(i);
    double error_ratio = error/central;
    mc_uncert->SetBinError(i, error_ratio);
  }
  return mc_uncert;
}
