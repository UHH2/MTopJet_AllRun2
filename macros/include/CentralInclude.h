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



TString dir = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/muon/";
TString dir_elec = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/elec/";
TString dir_combine = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/combine/";

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

void CMSLabel(bool bratio = false, double xl=0., bool ontop=false)
{
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.24+xl);
  text2->SetTextFont(62);
  if (bratio){
    text2->SetTextSize(0.08);
    if (ontop){
      text2->SetY(0.97);
    } else {
      text2->SetY(0.90);
    }
  } else {
    text2->SetTextSize(0.05);
    if (ontop){
      text2->SetY(0.995);
    } else {
      text2->SetY(0.90);
    }
  }
  text2->Draw();
}

void CMSSimLabel(bool bratio = false)
{
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.24);
  text2->SetTextFont(62);
  if (bratio){
    text2->SetTextSize(0.09);
    text2->SetY(0.90);
  } else {
    text2->SetTextSize(0.06);
    text2->SetY(0.90);
  }
  text2->Draw();

  TString simtext = "Simulation";
  TLatex *text3 = new TLatex(3.5, 24, simtext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.24);
  text3->SetTextFont(52);
  if (bratio){
    text3->SetTextSize(0.08);
    text3->SetY(0.84);
  } else {
    text3->SetTextSize(0.05);
    text3->SetY(0.84);
  }
  text3->Draw();

}

void SetHist(TH1F *hist, bool bPlotRatio = false){

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(1.15);
    hist->GetYaxis()->SetTickLength(0.02);

    // only this histogram
  } else {

    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.045);
    // hist->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetXaxis()->SetTickLength(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetNdivisions(505);


  // hist->GetYaxis()->SetNdivisions(505);
  
  hist->SetLineWidth(2);
}


void LumiInfo(double lumi = 19.7, bool bratio = false)
{
  TString infotext = TString::Format("%3.1f fb^{-1} (8 TeV)", lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetX(0.95);
  text1->SetTextFont(42);
  if (bratio){
    text1->SetTextSize(0.06);
    text1->SetY(1.);
  } else {
    text1->SetTextSize(0.045);
    text1->SetY(1.);
  }
  text1->Draw();

}
