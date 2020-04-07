#pragma once
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <iostream>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <TColor.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>

using namespace std;

TString dir = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/";

// -------------------------------------------------------------------------------------------
void SetupGlobalStyle()
{
  // general appearance and style

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");

  gStyle->UseCurrentStyle();

}

// -------------------------------------------------------------------------------------------
void CMSLabel(bool prelim, double x=0.25, double y=0.83){
  TString cmstext = "CMS";
  TLatex *text = new TLatex(3.5, 24, cmstext);
  text->SetNDC();
  text->SetTextAlign(13);
  text->SetX(x);
  text->SetTextFont(62);
  text->SetTextSize(0.07);
  text->SetY(y);
  text->Draw();

  if(prelim){
    TString simtext = "Preliminary";
    TLatex *text3 = new TLatex(3.5, 24, simtext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(x);
    text3->SetTextFont(52);
    text3->SetTextSize(0.05);
    text3->SetY(y-0.06);
    text3->Draw();
  }
}

// -------------------------------------------------------------------------------------------
void CMSSimLabel(bool prelim, double x=0.24, double y=0.9){
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(x);
  text2->SetTextFont(62);
  text2->SetTextSize(0.06);
  text2->SetY(y);
  text2->Draw();

  TString simtext = "Simulation";
  TLatex *text3 = new TLatex(3.5, 24, simtext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(x);
  text3->SetTextFont(52);
  text3->SetTextSize(0.05);
  text3->SetY(y-0.06);
  text3->Draw();

  if(prelim){
    TString simtext = "Preliminary";
    TLatex *text4 = new TLatex(3.5, 24, simtext);
    text4->SetNDC();
    text4->SetTextAlign(13);
    text4->SetX(x);
    text4->SetTextFont(52);
    text4->SetTextSize(0.05);
    text4->SetY(y-0.12);
    text4->Draw();
  }
}

// -------------------------------------------------------------------------------------------
void LumiInfo(double lumi = 35.9, bool bratio = false, double x=0.9, double y=0.945){
  TString infotext = TString::Format("%3.1f fb^{-1} (13 TeV)", lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(x);
  text1->SetY(y);
  text1->Draw();
}

// -------------------------------------------------------------------------------------------
TH1F *normalize(TH1F *hist){
  TH1F *hist_norm = (TH1F*) hist->Clone();
  hist_norm->Scale(1/hist->Integral());
  return hist_norm;
}

// -------------------------------------------------------------------------------------------
double get_highest_peak(vector<TH1F*> hists){
  vector<double> max;
  double max_value;
  for(unsigned int i = 0; i<hists.size(); i++) max.push_back(hists[i]->GetMaximum());
  sort(max.begin(), max.end());
  double top_index = max.size()-1;
  return max[top_index];
}

// -------------------------------------------------------------------------------------------
