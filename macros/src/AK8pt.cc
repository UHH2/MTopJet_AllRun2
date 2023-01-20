#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"
#include "../include/normalise.h"
using namespace std;

void DrawHistogram(TString path, VecH hists, TH1F* rData, TH1F* rMC, VecI colors, VecTS names, TString xaxis, double xmin = 0, TString binwidth = "");
TH1F* GetHist(TFile* file, TString name, int rebin);

Double_t lumi_plot = 138;

TString save = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/";

bool debug = false;

int main(int argc, char* argv[]){

  gErrorIgnoreLevel = kWarning;

  TFile* infile = new TFile("FSR_hists_pt0toInf.root");

  cout << "Get Histograms" << endl;
  int rebin = 20;
  TH1F *data = GetHist(infile, "data", rebin);
  TH1F *ttbar = GetHist(infile, "nominal", rebin);
  TH1F *wjet = GetHist(infile, "bgr_wj", rebin);
  TH1F *single = GetHist(infile, "bgr_st", rebin);
  TH1F *other = GetHist(infile, "bgr_ot", rebin);

  TH1F *MC = (TH1F*) ttbar->Clone();
  MC->Add(wjet, 1);
  MC->Add(single, 1);
  MC->Add(other, 1);

  TH1F* ratio_data = GetRatio(data, MC, true);
  TH1F* ratio_mc = GetRatio(MC, MC, true);

  TString binwidth = "/"+to_string(rebin)+" GeV";
  DrawHistogram("AK8pt_fullRun2", {data, other, wjet, single, ttbar}, ratio_data, ratio_mc, {kBlack, kOrange+1, kGreen+2, kBlue+1, kRed+1}, {"Data", "Other", "WJets", "Single Top", "t#bar{t}"}, "#it{p}_{T}", 300, binwidth);

  return 0;
}

TH1F* GetHist(TFile* file, TString name, int rebin)
{
  TH1F *h16 = (TH1F*) file->Get(name+"_pt_2016");
  TH1F *h17 = (TH1F*) file->Get(name+"_pt_2017");
  TH1F *h18 = (TH1F*) file->Get(name+"_pt_2018");
  TH1F *h = (TH1F*) h16->Clone();
  h->Add(h17, 1);
  h->Add(h18, 1);
  h->Rebin(rebin);
  return h;
}

void DrawHistogram(TString path, VecH hists, TH1F* rData, TH1F* rMC, VecI colors, VecTS names, TString xaxis, double xmin, TString binwidth)
{
  TString unique_path = save+"AK8pt_fullRun2.pdf";
  double nHists = hists.size();
  double bins = hists[0]->GetNbinsX();

  if(debug) cout << "Draw Control Hists ..." << endl;
  hists[0]->SetTitle("");
  hists[0]->GetXaxis()->SetTitle("");
  hists[0]->GetYaxis()->SetTitle("Events"+binwidth);
  hists[0]->SetMarkerStyle(8);
  hists[0]->SetMarkerColor(colors[0]);
  hists[0]->SetLineColor(colors[0]);

  THStack* stack = new THStack(unique_path, "");
  for(unsigned int c=1; c<colors.size(); c++){
    hists[c]->SetFillColor(colors[c]);
    hists[c]->SetLineColor(colors[c]);
    stack->Add(hists[c]);
  }
  double maximum = stack->GetMaximum();

  TLegend *leg = new TLegend(0.7,0.6,0.85,0.85);
  TCanvas *A = new TCanvas(unique_path, unique_path, 600, 600);

  gStyle->SetOptTitle(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  TH1F* errorUp = (TH1F*) rMC->Clone(); TH1F* errorDown = (TH1F*) rMC->Clone();
  for(unsigned int i=1; i<=rMC->GetNbinsX(); i++){
    errorUp->SetBinContent(i, rMC->GetBinContent(i)+rMC->GetBinError(i));
    errorDown->SetBinContent(i, rMC->GetBinContent(i)-rMC->GetBinError(i));
  }
  errorUp->SetFillColorAlpha(kGray+2, 0.5);
  errorUp->SetLineWidth(0);
  errorDown->SetFillColorAlpha(kWhite, 1);
  errorDown->SetLineWidth(0);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  hists[0]->GetXaxis()->SetLabelSize(0);
  hists[0]->SetStats(false);
  hists[0]->GetXaxis()->SetRangeUser(xmin, 1000);
  hists[0]->GetYaxis()->SetRangeUser(0, maximum*1.3);
  hists[0]->GetYaxis()->SetTitle("Events");
  hists[0]->GetYaxis()->SetLabelSize(0.05);
  hists[0]->GetYaxis()->SetTitleSize(0.06);
  hists[0]->GetYaxis()->SetTitleOffset(1.15);
  hists[0]->Draw("P");
  stack->Draw("same hist");
  hists[0]->Draw("same PE"); // redraw
  gPad->RedrawAxis();

  leg->AddEntry(hists[0], "Data", "pl");
  for(int i=nHists-1; 1<=i; i--) leg->AddEntry(hists[i], names[i], "f");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();

  CMSLabel(0.25, 0.83, "Work in Progress");
  LumiInfo(lumi_plot, false, 0.9, 0.96);

  // Ratio Plot kommt in unteres pad
  A->cd();
  if(debug) cout << "\t ... Pad2" << endl;
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  A->SetTickx(0);
  A->SetTicky(0);

  TLine* lMC = new TLine(-2.5, 1, 2.5, 1);
  lMC->SetLineWidth(1);
  lMC->SetLineColor(kBlack);

  rData->SetMarkerStyle(8);
  rData->SetLineColor(colors[0]);
  rMC->SetFillStyle(0);
  rMC->SetLineStyle(kSolid);
  rMC->SetLineWidth(0);
  rMC->SetStats(false);
  rMC->GetXaxis()->SetRangeUser(xmin, 1000);
  rMC->GetXaxis()->SetLabelSize(0.15);
  rMC->GetYaxis()->SetLabelSize(0.13);
  rMC->GetXaxis()->SetTitleSize(0.1);
  rMC->GetYaxis()->SetRangeUser(0.4, 1.6);
  rMC->GetXaxis()->SetTitle(xaxis);
  rMC->GetYaxis()->SetTitle("#frac{DATA}{MC}");
  rMC->GetXaxis()->SetTitleOffset(0.8);
  rMC->GetXaxis()->SetTitleSize(0.2);
  rMC->GetYaxis()->SetTitleSize(0.15);
  rMC->GetYaxis()->SetNdivisions(505);
  rMC->GetYaxis()->SetTitleOffset(0.4);
  rMC->GetYaxis()->CenterTitle();
  rMC->Draw("Axis");
  errorUp->Draw("Hist same");
  errorDown->Draw("Hist same");
  lMC->Draw("same");
  rData->Draw("SAME PE");

  gPad->RedrawAxis();
  A->SaveAs(unique_path);

  delete A;
  leg->Clear();
}
