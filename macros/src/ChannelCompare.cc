#include "../include/CentralInclude.h"

using namespace std;

void Plot(TGraphAsymmErrors* elec, TGraphAsymmErrors* elec_stat, TGraphAsymmErrors* muon, TGraphAsymmErrors* muon_stat);
TGraphAsymmErrors* ConvertToGraph(TH1F*, TString);
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

  TH1F * h_elec = (TH1F*)f_elec->Get("Unfold_XS_totuncert");
  TH1F * h_elec_stat = (TH1F*)f_elec->Get("Unfold_XS_statuncert");
  TH1F * h_muon = (TH1F*)f_muon->Get("Unfold_XS_totuncert");
  TH1F * h_muon_stat = (TH1F*)f_muon->Get("Unfold_XS_statuncert");

  TGraphAsymmErrors* elec = ConvertToGraph(h_elec, "left");
  TGraphAsymmErrors* elec_stat = ConvertToGraph(h_elec_stat, "left");
  TGraphAsymmErrors* muon = ConvertToGraph(h_muon, "right");
  TGraphAsymmErrors* muon_stat = ConvertToGraph(h_muon_stat, "right");
  Plot(elec, elec_stat, muon, muon_stat);
}

void Plot(TGraphAsymmErrors* elec, TGraphAsymmErrors* elec_stat, TGraphAsymmErrors* muon, TGraphAsymmErrors* muon_stat){

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  elec->Draw("AP");
  elec->SetTitle(" ");
  elec->GetXaxis()->SetRangeUser(112, 232);
  elec->GetYaxis()->SetRangeUser(0, 11);
  elec->GetXaxis()->SetTitle("m_{jet} [GeV]");
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
  elec_stat->SetFillColor(15);
  elec_stat->SetMarkerColor(kBlack);
  elec_stat->SetMarkerStyle(8);
  elec_stat->SetMarkerSize(1);
  elec->Draw("E2 SAME");
  elec_stat->Draw("E2 SAME");
  elec_stat->Draw("PX SAME");

  muon->SetLineColor(kRed+1);
  muon->SetMarkerColor(kRed+1);
  muon->SetMarkerStyle(8);
  muon->SetMarkerSize(1);
  muon->Draw("P SAME");
  muon_stat->SetLineColor(kRed+1);
  muon_stat->SetMarkerColor(kRed+1);
  muon_stat->SetMarkerStyle(8);
  muon_stat->SetMarkerSize(1);
  muon_stat->Draw("P SAME");
  muon->Draw("P SAME");

  gStyle->SetEndErrorSize(5);
  gPad->RedrawAxis();

  CMSLabel(true, 0.2, 0.85);

  TLegend *l=new TLegend(0.55,0.68,0.85,0.87);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(elec,"elec channel","pf");
  l->AddEntry(muon,"muon channel","pe");
  l->SetTextSize(0.03);
  l->Draw();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/ChannelComparison.pdf");
  delete c;
}


TGraphAsymmErrors* ConvertToGraph(TH1F* hist, TString leftright){
  vector<double> x,y,xel,xer,ye;
  int nbins = hist->GetXaxis()->GetNbins();
  for(int bin = 1; bin<=nbins; bin++){
    double x_ = hist->GetBinCenter(bin);
    if(leftright == "left") x_ -= 2.0;
    if(leftright == "right") x_ += 2.0;
    // cout << x_ << ", " << hist->GetBinContent(bin) << endl;
    x.push_back(x_);
    double xel_ = hist->GetBinWidth(bin)/2;
    double xer_ = hist->GetBinWidth(bin)/2;
    if(leftright == "left"){
      xel_ -= 2.0;
      xer_ += 2.0;
    }
    if(leftright == "right"){
      // xel_ += 2.0;
      // xer_ -= 2.0;
      xel_ = 0.0;
      xer_ = 0.0;
    }
    xel.push_back(xel_);
    xer.push_back(xer_);
    y.push_back(hist->GetBinContent(bin));
    ye.push_back(hist->GetBinError(bin));
  }
  TGraphAsymmErrors *g = new TGraphAsymmErrors(nbins, &x[0], &y[0], &xel[0], &xer[0], &ye[0], &ye[0]);
  return g;
}
