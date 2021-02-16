
#include <TGraph.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TROOT.h>

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

void SetGraph(TGraph *hist){
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetXaxis()->SetTickLength(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);
  // hist->GetXaxis()->SetNdivisions(505);

  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.045);
  // hist->SetLabelSize(0.045);
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

  // hist->GetYaxis()->SetNdivisions(505);

   hist->SetLineWidth(2);
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

void CMSPrelimSimLabel(bool bratio = false)
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

    TString pretext = "Preliminary";
    TLatex *text4 = new TLatex(3.5, 24, pretext);
    text4->SetNDC();
    text4->SetTextAlign(13);
    text4->SetX(0.24);
    text4->SetTextFont(52);
    text4->SetTextSize(0.05);
    text4->SetY(0.78);
    text4->Draw();
}


void LumiInfo(double lumi = 35.9, bool bratio = false)
{
  TString infotext = TString::Format("%3.1f fb^{-1} (13 TeV)", lumi);
  if(lumi > 100) infotext = TString::Format("%3.0f fb^{-1} (13 TeV)", lumi);

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



void ScaletoPercent(TH1F *hist)
{
  for (int i=1; i<hist->GetNbinsX()+1; ++i){
    hist->SetBinContent(i, 100*hist->GetBinContent(i));
    hist->SetBinError(i, 0.1);
  }
}

void ModStatUnc(TH1F *hist)
{
  for (int i=1; i<hist->GetNbinsX()+1; ++i){
    hist->SetBinContent(i, hist->GetBinError(i));
    hist->SetBinError(i, 0.001);
  }
}

void CalcTotalUnc(vector<TH1F*>& hists)
{
  TH1F* hsum = NULL;
  if (hists.size()>0){
    hsum = (TH1F*) hists.at(0);
  }

  for (int i=1; i<hsum->GetNbinsX()+1; ++i){
    double sum2 = 0;
    for (unsigned int j=1; j<hists.size(); ++j){
      //cout << "hist " << j << " name = " << hists.at(j)->GetName()
      //     << " error in bin " << i << " = " << hists.at(j)->GetBinContent(i) << endl;
      sum2 += TMath::Power(hists.at(j)->GetBinContent(i),2);
    }

    hsum->SetBinContent(i, TMath::Sqrt(sum2));
    //cout << "bin " << i << " total = " << TMath::Sqrt(sum2) << endl;

  }

  //hists.push_back(hsum);

}


void SingleEPSRatioCosmetics(TH1* hist)
{

  hist->GetYaxis()->SetRangeUser(0.3, 1.7);
  //hist->GetYaxis()->SetRangeUser(0.05, 1.95);
  hist->SetMarkerSize(0.7);

  // cosmetics for portrait mode

    hist->SetTitle("");

    // x-axis
    hist->GetXaxis()->SetLabelSize(0.12);
    hist->GetXaxis()->SetTickLength(0.08);
    hist->GetXaxis()->SetTitleSize(0.12);
    hist->GetXaxis()->SetTitleOffset(1.25);

    // y-axis
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.12);
    hist->GetYaxis()->SetTitleOffset(0.66);
    hist->GetYaxis()->SetLabelSize(0.11);
    //hist->GetYaxis()->SetNdivisions(210);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    // cosmetics for landscape mode
    /* } else {

    hist->SetTitle("");
    hist->SetTitleOffset(1.1, "X");
    hist->SetTitleOffset(0.5, "Y");
    hist->SetLabelOffset(0.02, "X");
    hist->SetLabelOffset(0.01, "Y");

    hist->GetXaxis()->SetLabelSize(0.14);
    hist->GetXaxis()->SetTickLength(0.07);
    hist->GetXaxis()->SetTitleSize(0.15);

    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(0.11);
    hist->GetYaxis()->SetLabelSize(0.12);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetTickLength(0.03);

  }
    */
}

void XConeLabel(double x, double y){

  TString xconetext = "XCone, R_{jet}_{ }=_{ }1.2";
  TLatex *t = new TLatex(3.5, 24, xconetext);
  t->SetNDC();
  t->SetTextAlign(13);
  t->SetTextFont(42);
  t->SetTextSize(0.04);
  t->SetX(x);
  t->SetY(y);
  t->Draw();

  TString subtext = "N_{sub}_{ }=_{ }3, R_{sub}_{ }=_{ }0.4";
  TLatex *t2 = new TLatex(3.5, 24, subtext);
  t2->SetNDC();
  t2->SetTextAlign(13);
  t2->SetX(x);
  t2->SetTextFont(42);
  t2->SetTextSize(0.04);
  t2->SetY(y-0.055);
  t2->Draw();

  TString pttext = "p_{T}_{ }>_{ } 400 GeV";
  TLatex *t3 = new TLatex(3.5, 24, pttext);
  t3->SetNDC();
  t3->SetTextAlign(13);
  t3->SetX(x);
  t3->SetTextFont(42);
  t3->SetTextSize(0.04);
  t3->SetY(y-0.11);
  t3->Draw();

}
