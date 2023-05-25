#pragma once
#include "CentralInclude.h"

void SetupCanvas(bool bPlotRatio=false);
void HistCosmetics(TH1* hist);
void PortraitCosmetics(TH1* hist, bool bPlotRatio);
void YieldCosmetics(TH1* hist);
void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width);
void Cosmetics(TH1* hist, TString xtitle, int color, int style, int marker, int width, bool bPlotRatio);
void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre, TString name="");
void DrawLegend(vector<TH1F*> hists, TString hname);


TCanvas* m_can;
TPad* m_pad1;
TPad* m_pad2;

TPad* m_rp1_top;
TPad* m_rp1;
TPad* m_rp2_top;
TPad* m_rp2;

void SetupCanvas(bool bPlotRatio) // Copied from SPlotter.cxx
{
  // set up a canvas for single EPS files
  // optimised plots for including in theses or publications and documents
  // different possibilities
  // ratio/no ratio plots

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 400;
  CanHeight = 400;

  // set up the canvas
  m_can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //
  Float_t y1, y2, y3;                           //  y3 +-------------+
  y3 = 0.99;                                    //     |             |
  y2 = y3-yplot;                                //     |     pad1    |
  y1 = y2-yratio;                               //  y2 |-------------|
  Float_t x1, x2;                               //     |     rp1     |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
  //
  // No Pad 2!


  m_rp1_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp1 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad1 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);

  m_rp2_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp2 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad2 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);


  m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
  m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);

  m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.04);  m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
  m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.04);  m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);

  m_rp1->SetTopMargin(0.0);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19);  m_rp1->SetRightMargin(0.05);
  m_rp2->SetTopMargin(0.0);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19);  m_rp2->SetRightMargin(0.05);

  m_pad1->Draw();
  m_pad2->Draw();

  if (bPlotRatio){
    m_rp1_top->Draw();
    m_rp1->Draw();
  }
  else{
    m_pad1->Draw();
  }

  return;

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void YieldCosmetics(TH1* hist)
{
  // cosmetics for the lumi yield histogram
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetXaxis()->SetTickLength(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);

  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

  hist->GetXaxis()->SetTitle("integrated luminosity [fb^{-1}]");
  double dlum = hist->GetXaxis()->GetBinWidth(1);
  TString xtit = TString::Format("events per %3.1f fb^{-1}", dlum);
  hist->GetYaxis()->SetTitle(xtit);
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void PortraitCosmetics(TH1* hist, bool bPlotRatio)
{

  // top histogram of the ratio plot
  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);
    hist->GetXaxis()->SetLabelSize(0.00); //remove labels from Xaxis of upper plot
    hist->GetXaxis()->SetTitleSize(0.00); // remove title from Xaxis

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTickLength(0.02);

    // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.045);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
{
  PortraitCosmetics(hist, bPlotRatio);

  hist->SetLineColor(color);
  hist->SetLineStyle(style);
  hist->SetLineWidth(width);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.7);

  // set Y-axis title
  hist->GetYaxis()->SetTitle("a.u.");
  if(xtitle.Contains("JEC")) hist->GetYaxis()->SetTitle("#it{f}^{XCone}");

  // set X-axis title
  hist->GetXaxis()->SetTitle(xtitle);
  hist->SetTitle("");

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  // top histogram of the ratio plot
  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetNdivisions(505);

    TString hname = hist->GetName();
    hist->GetXaxis()->SetNdivisions(505);

    // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.045);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
{
  hist->SetLineColor(color);
  hist->SetLineStyle(style);
  hist->SetLineWidth(width);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.7);

  hist->GetYaxis()->SetRangeUser(0.3, 1.7);
  if(xtitle.Contains("#it{m}_{W}")) hist->GetYaxis()->SetRangeUser(0.7, 1.3);
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerSize(0.7);

  hist->SetTitle("");
  // x-axis
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetLabelSize(0.12);
  hist->GetXaxis()->SetTickLength(0.08);
  hist->GetXaxis()->SetTitleSize(0.15);
  hist->GetXaxis()->SetTitleOffset(1.15);
  if(xtitle.Contains("#it{m}_{W}")) hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetXaxis()->SetLabelOffset(0.02);
  TString hname = hist->GetName();
  hist->GetXaxis()->SetNdivisions(505);

  // y-axis
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle("Ratio");
  hist->GetYaxis()->SetTitleSize(0.12);
  hist->GetYaxis()->SetTitleOffset(0.78); //0.66
  hist->GetYaxis()->SetLabelSize(0.11);
  //hist->GetYaxis()->SetNdivisions(210);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre, TString name)
{
  TString infotext = TString::Format("%3.0f fb^{-1} (13 TeV)", m_lumi);
  if(m_lumi==36.3) infotext = TString::Format("%2.1f fb^{-1} (13 TeV)", m_lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetX(0.95);
  text1->SetTextFont(42);
  if (bPlotRatio){
    text1->SetTextSize(0.06);
    text1->SetY(1.);
  } else {
    text1->SetTextSize(0.045);
    text1->SetX(0.90);
    text1->SetY(1.);
  }
  if(name.Contains("JMS")) text1->SetX(0.82);
  text1->Draw();

  if (DrawCMS || forPre){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.22); // standard was 0.24
    text2->SetTextFont(62);
    if (bPlotRatio){
      text2->SetTextSize(0.08);
      text2->SetY(0.87);
    } else {
      text2->SetTextSize(0.05);
      text2->SetY(0.895);
    }
    if(name.Contains("JMS")) text2->SetX(0.15);
    text2->Draw();
  }

  if (forPre){
    // TString preltext = "Preliminary";
    TString preltext = "Work in Progress";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.22); // standard was 0.24
    text3->SetTextFont(52);
    if (bPlotRatio){
      //  text3->SetTextSize(0.06);
      text3->SetTextSize(0.065);
      text3->SetY(0.78);
    } else {
      // text3->SetTextSize(0.035);
      text3->SetTextSize(0.04);
      text3->SetX(0.33); // standard was 0.24
      if(preltext.Contains("Work")) text3->SetX(0.32);
      text3->SetY(0.89);
    }
    if(name.Contains("JMS")) text3->SetX(0.26);
    text3->Draw();
  }

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void DrawLegend(vector<TH1F*> hists, TString hname)
{
  // draw a legend
  cout << "\t\t ... In Legend " << endl;
  int narr = hists.size();
  float yfrac = 0.07;

  float top = 0.85;

  float ysize = yfrac*narr;
  float xleft = 0.64;
  float xright = 0.92;

  if (hname.Contains("tau32")){
    xleft = 0.21;
    xright = 0.39;
    top = 0.76;
    ysize = 0.112*narr;
  }
  if (hname.Contains("JMS")){
    xleft = 0.22;
    xright = 0.4;
    top = 0.68;
    ysize = 0.09*narr;
  }
  if (hname.Contains("wmass")){
    if(hname.Contains("ll")){
      xleft = 0.5;
      xright = 0.77;
      top = 0.54;
      ysize = 0.08*narr;
    }
    else if(hname.Contains("lh")){
      xleft = 0.5;
      xright = 0.77;
      top = 0.54;
      ysize = 0.08*narr;
    }
    else if(hname.Contains("jer")){
      xleft = 0.65;
      xright = 0.90;
      top = 0.85;
      ysize = 0.08*narr;
    }
    else{
      xleft = 0.45;
      xright = 0.7;
      top = 0.54;
      ysize = 0.08*narr;
    }
  }
  if (hname.Contains("ECF")){
    xleft = 0.6;
    xright = 0.78;
    top = 0.76;
    ysize = 0.112*narr;
  }
  if (hname.Contains("r21")){
    xleft = 0.67;
    xright = 0.85;
    top = 0.76;
    ysize = 0.112*narr;
  }
  if (hname.Contains("r31")){
    xleft = 0.7;
    xright = 0.88;
    top = 0.76;
    ysize = 0.112*narr;
  }
  if (hname.Contains("r") && hname.Contains("gen")) {
    xleft = 0.5;
    xright = 0.68;
    top = 0.55;
    ysize = 0.112*narr;
  }

  xright = xleft + 0.22;
  top += 0.02;

  TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.067);
  if(hname.Contains("tau32")) leg->SetTextSize(0.0645);

  for (Int_t i=0; i<narr; ++i){
    TH1F* sh = hists[i];

    TString legname = TString::Format("leg_entry_%i",i);
    TString legtitle;
    if(hname.Contains("tau32_combine")){
      if(i==0) legtitle = "Data";
      if(i==1) legtitle = "t#bar{t}";
      if(i==2) legtitle = "t#bar{t} #it{f}_{FSR} = #frac{1}{4}";
      if(i==3) legtitle = "t#bar{t} #it{f}_{FSR} = 4";
    }
    else if(hname.Contains("tau32_2016") || hname.Contains("ECF") || hname.Contains("r")){
      if(i==0) legtitle = "Data";
      if(i==1) legtitle = "t#bar{t}";
      if(i==2) legtitle = "t#bar{t} #it{f}_{FSR} = #frac{1}{2}";
      if(i==3) legtitle = "t#bar{t} #it{f}_{FSR} = 2";
    }
    else if(hname.Contains("wmass")){
      if(hname.Contains("jer")){
        if(i==0) legtitle = "Data #minus Bkg";
        if(i==1) legtitle = "t#bar{t}";
        if(i==2) legtitle = "t#bar{t} #it{f}^{JER} = #plus 1";
        if(i==3) legtitle = "t#bar{t} #it{f}^{JER} = #minus 1";
      }
      else{
        if(i==0) legtitle = "Data";
        if(i==1) legtitle = "t#bar{t}";
        if(i==2) legtitle = "t#bar{t} #it{f}^{JEC} = #plus 1";
        if(i==3) legtitle = "t#bar{t} #it{f}^{JEC} = #minus 1";
        if(i==4) legtitle = "t#bar{t} #it{f}^{XCone} = #plus 1";
        if(i==5) legtitle = "t#bar{t} #it{f}^{XCone} = #minus 1";
      }
    }
    else if(hname.Contains("JMS")){
      if(i==0) legtitle = "Minimum";
      if(i==1) legtitle = "68#% CL";
      if(i==2) legtitle = "95#% CL";
    }

    TLegendEntry* entry = NULL;
    int marker = sh->GetMarkerStyle();
    int lstyle = sh->GetLineStyle();
    if (marker>0){
      entry = leg->AddEntry(sh, legtitle, "pe");
      entry->SetLineWidth(1);
      entry->SetLineColor(sh->GetLineColor());
      entry->SetMarkerColor(sh->GetLineColor());
      entry->SetMarkerStyle(marker);
      entry->SetMarkerSize(1.0);
      entry->SetMarkerSize(0.8);
    } else {

      entry = leg->AddEntry(sh, legtitle, "l");
      entry->SetLineColor(sh->GetLineColor());
      entry->SetMarkerStyle(0);
      entry->SetMarkerSize(0);
      entry->SetMarkerColor(sh->GetLineColor());
      entry->SetLineWidth(2.7);
      entry->SetLineStyle(lstyle);

      entry->SetTextAlign(12);
      //entry->SetTextColor(fSampleColors.At(i));
    }

  }

  leg->Draw();

}
