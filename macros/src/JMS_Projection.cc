#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
// #include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/tdrstyle_all.h"

// #include "../CovMatrices/JMS/CollectCovHeaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>
#include "TSystem.h"

using namespace std;

bool debug = true;
bool prelim = false;

void SetupCanvas(bool bPlotRatio=false);
void HistCosmetics(TH1* hist);
void PortraitCosmetics(TH1* hist, bool bPlotRatio);
void YieldCosmetics(TH1* hist);
void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width);
void Cosmetics(TH1* hist, TString xtitle, int color, int style, int marker, int width, bool bPlotRatio);
void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre, TString name="");
void DrawLegend(vector<TH1F*> hists, TString hname);
void DrawLegend(vector<TPolyMarker3D*> hists, TString hname);
void CompareHistStructureLocal(TH1F* h1, TH1F* h2);
void DrawJMS(TF2* chi2, VecDD &points, VecDD &extrema, VecD &minimum, TString hist);
TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool equal, bool isData = false);

TCanvas* m_can;
TPad* m_pad1;
TPad* m_pad2;

TPad* m_rp1_top;
TPad* m_rp1;
TPad* m_rp2_top;
TPad* m_rp2;

TString save_afs = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/PaperPlots/";
TString save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/Projections/";

int main(int argc, char* argv[]){

  TFile *file = new TFile("files/PaperPlots_Peak.root", "read");

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();
  gStyle->SetErrorX(0);

  // prelim = stob(argv[1]);
  // if(prelim) save_nfs.ReplaceAll("PaperPlots/", "PaperPlots/preliminary/");

  // ===========================================================================
  // === JMS 2D Chi2
  cout << "Start with JMS plots ..." << endl;

  TF2* chi2 = (TF2*) file->Get("Functions/JMS_Chi2");
  TString formula = chi2->GetExpFormula();

  TString hold = "x";
  vector<double> range = {0.0,1.2};
  double step=0.1;
  TH1F* hist = new TH1F("hist", "hist", abs(range[0]-range[1])/step, range[0], range[1]);
  for(unsigned int c=1;c<=hist->GetNbinsX();c++){
    cout << c<<setw(15)<<hist->GetXaxis()->GetBinCenter(c)<<setw(15)<<hist->GetXaxis()->GetBinLowEdge(c)<<setw(15)<<hist->GetXaxis()->GetBinUpEdge(c)<<endl;
  }
  double start=range[0];
  for(double i=range[0]; i<range[1]; i+=step){
    double val = i+step/2;
    string value = to_string(val);
    TString f_temp = formula;
    f_temp.ReplaceAll(hold, value);
    if(hold.EqualTo("x")) f_temp.ReplaceAll('y', 'x');
    TF1 *func = new TF1("func", f_temp, -1.5, 1.5);
    double minimum = func->GetMinimum();
    double xmin = func->GetX(minimum, -1.5, 1.5);
    double eu = func->GetX(minimum+1, xmin, 1.5)-xmin;
    double ed = xmin-func->GetX(minimum+1, -1.5, xmin);
    int bin = hist->FindBin(val);
    hist->SetBinContent(bin, minimum);
    hist->SetBinError(bin, (eu+ed)/2);
    cout << value << "\t" << bin << "\t" <<  minimum << "\t" << xmin << "\t" << eu << "\t" << ed << endl;
  }
  TF1 *fit = new TF1("parabel", "[0]*x*x+[1]*x+[2]", range[0], range[1]);
  hist->Fit(fit, "MQ", "", range[0], range[1]);
  cout<< fit->Eval(0.6) << "   " << fit->GetParameter(0) << endl;

  double minimum = fit->GetMinimum();
  double ymin = fit->GetX(minimum, range[0], range[1]);
  double eu = fit->GetX(minimum+1, ymin, range[1])-ymin;
  double ed = ymin-fit->GetX(minimum+1, range[0], ymin);
  cout << minimum << " " << ymin  << " +" << eu  << " -" << ed << endl;

  TCanvas *canvas = new TCanvas("JEC");
  hist->Draw("L");
  TLine *lmin=new TLine(ymin, 128, ymin, minimum);
  TLine *lerr=new TLine(range[0], minimum+1, range[0], minimum+1);
  TLine *lerrU=new TLine(ymin+eu, 128, ymin+eu, minimum+1);
  TLine *lerrD=new TLine(ymin-ed, 128, ymin-ed, minimum+1);
  lmin->SetLineColor(kGray+2);
  lerr->SetLineColor(kGray+2);
  lerrU->SetLineColor(kGray+2);
  lerrD->SetLineColor(kGray+2);
  lmin->SetLineStyle(2);
  lerr->SetLineStyle(1);
  lerrU->SetLineStyle(1);
  lerrD->SetLineStyle(1);
  lmin->Draw("same");
  lerr->Draw("same");
  lerrU->Draw("same");
  lerrD->Draw("same");
  // TLegend *leg = new TLegend(0.3, 0.7, 0.7, 0.8);
  // leg->AddEntry(lmin, "Minimum ('+str("%.2f"%min_x)+','+str("%.2f"%min_y)+')', 'l');
  // leg->AddEntry(lerr, "1#sigma ('+str("%.2f"%(min_x-ed))+','+str("%.2f"%(min_x+eu))+')', 'l');
  // leg->Draw();
  canvas->Print(save_nfs+"1D_JEC_with_c.pdf");
  canvas->Close();

  hold = "y";
  range = {-0.8,0.8};
  step=0.1;
  hist = new TH1F("hist", "hist", abs(range[0]-range[1])/step, range[0], range[1]);
  for(unsigned int c=1;c<=hist->GetNbinsX();c++){
    cout << c<<setw(15)<<hist->GetXaxis()->GetBinCenter(c)<<setw(15)<<hist->GetXaxis()->GetBinLowEdge(c)<<setw(15)<<hist->GetXaxis()->GetBinUpEdge(c)<<endl;
  }
  start=range[0];
  for(double i=range[0]; i<range[1]; i+=step){
    double val = i+step/2;
    string value = to_string(val);
    TString f_temp = formula;
    f_temp.ReplaceAll(hold, value);
    if(hold.EqualTo("x")) f_temp.ReplaceAll('y', 'x');
    TF1 *func = new TF1("func", f_temp, -1.5, 1.5);
    double minimum = func->GetMinimum();
    double xmin = func->GetX(minimum, -1.5, 1.5);
    double eu = func->GetX(minimum+1, xmin, 1.5)-xmin;
    double ed = xmin-func->GetX(minimum+1, -1.5, xmin);
    int bin = hist->FindBin(val);
    hist->SetBinContent(bin, minimum);
    hist->SetBinError(bin, (eu+ed)/2);
    cout << value << "\t" << bin << "\t" <<  minimum << "\t" << xmin << "\t" << eu << "\t" << ed << endl;
  }
  fit = new TF1("parabel", "[0]*x*x+[1]*x+[2]", range[0], range[1]);
  hist->Fit(fit, "MQ", "", range[0], range[1]);
  cout<< fit->Eval(0.6) << "   " << fit->GetParameter(0) << endl;

  minimum = fit->GetMinimum();
  ymin = fit->GetX(minimum, range[0], range[1]);
  eu = fit->GetX(minimum+1, ymin, range[1])-ymin;
  ed = ymin-fit->GetX(minimum+1, range[0], ymin);
  cout << minimum << " " << ymin  << " +" << eu  << " -" << ed << endl;

  canvas = new TCanvas("JEC");
  hist->Draw("L");
  lmin=new TLine(ymin, 128, ymin, minimum);
  lerr=new TLine(range[0], minimum+1, range[0], minimum+1);
  lerrU=new TLine(ymin+eu, 128, ymin+eu, minimum+1);
  lerrD=new TLine(ymin-ed, 128, ymin-ed, minimum+1);
  lmin->SetLineColor(kGray+2);
  lerr->SetLineColor(kGray+2);
  lerrU->SetLineColor(kGray+2);
  lerrD->SetLineColor(kGray+2);
  lmin->SetLineStyle(2);
  lerr->SetLineStyle(1);
  lerrU->SetLineStyle(1);
  lerrD->SetLineStyle(1);
  lmin->Draw("same");
  lerr->Draw("same");
  lerrU->Draw("same");
  lerrD->Draw("same");
  // TLegend *leg = new TLegend(0.3, 0.7, 0.7, 0.8);
  // leg->AddEntry(lmin, "Minimum ('+str("%.2f"%min_x)+','+str("%.2f"%min_y)+')', 'l');
  // leg->AddEntry(lerr, "1#sigma ('+str("%.2f"%(min_x-ed))+','+str("%.2f"%(min_x+eu))+')', 'l');
  // leg->Draw();
  canvas->Print(save_nfs+"1D_XCone_with_c.pdf");
  canvas->Close();

}

// ===========================================================================
// ===========================================================================
// ===========================================================================

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

  m_rp1_top->Draw();
  // m_rp2_top->Draw();

  if (bPlotRatio){
    m_rp1->Draw();
    // m_rp2->Draw();
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
  if(xtitle.Contains("JEC")) hist->GetYaxis()->SetTitle("f^{XCone}");

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
    hist->GetYaxis()->SetTitleOffset(1.35);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetNdivisions(505);

    TString hname = hist->GetName();
    hist->GetXaxis()->SetNdivisions(505);

    // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleSize(0.05);
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
  hist->GetXaxis()->SetTitleSize(0.12);
  hist->GetXaxis()->SetTitleOffset(1.25);
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
    text1->SetY(0.95);
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
      text2->SetY(0.945);
    }
    if(name.Contains("JMS")) text2->SetX(0.15);
    text2->Draw();
  }

  if (forPre){
    TString preltext = "Preliminary";
    // TString preltext = "Work in Progress";
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
      text3->SetY(0.94);
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

  int narr = hists.size();
  float yfrac = 0.07;

  float top = 0.85;

  float ysize = yfrac*narr;
  float xleft = 0.64;
  float xright = 0.92;

  if (hname.Contains("tau32")){
    xleft = 0.22;
    xright = 0.4;
    top = 0.68;
    ysize = 0.09*narr;
  }
  if (hname.Contains("JMS")){
    xleft = 0.22;
    xright = 0.4;
    top = 0.68;
    ysize = 0.09*narr;
  }
  if (hname.Contains("wmass")){
    if(hname.Contains("ll")){
      xleft = 0.60;
      xright = 0.75;
      top = 0.5;
      ysize = 0.07*narr;
    }
    else{
      xleft = 0.45;
      xright = 0.7;
      top = 0.5;
      ysize = 0.07*narr;
    }
  }

  xright = xleft + 0.22;
  top += 0.02;

  TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);

  for (Int_t i=0; i<narr; ++i){
    TH1F* sh = hists[i];

    TString legname = TString::Format("leg_entry_%i",i);
    TString legtitle;
    if(hname.Contains("tau32_combine")){
      if(i==0) legtitle = "Data";
      if(i==1) legtitle = "t#bar{t}";
      if(i==2) legtitle = "t#bar{t} f_{FSR} = #frac{1}{4}";
      if(i==3) legtitle = "t#bar{t} f_{FSR} = 4";
    }
    else if(hname.Contains("tau32_2016")){
      if(i==0) legtitle = "Data";
      if(i==1) legtitle = "t#bar{t}";
      if(i==2) legtitle = "t#bar{t} f_{FSR} = #frac{1}{2}";
      if(i==3) legtitle = "t#bar{t} f_{FSR} = 2";
    }
    else if(hname.Contains("wmass")){
      if(i==0) legtitle = "Data";
      if(i==1) legtitle = "t#bar{t}";
      if(i==2) legtitle = "t#bar{t} f^{JEC} = +1";
      if(i==3) legtitle = "t#bar{t} f^{JEC} = -1";
      if(i==4) legtitle = "t#bar{t} f^{XCone} = +1";
      if(i==5) legtitle = "t#bar{t} f^{XCone} = -1";
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
      entry->SetLineWidth(2);
      entry->SetLineStyle(lstyle);

      entry->SetTextAlign(12);
      //entry->SetTextColor(fSampleColors.At(i));
    }

  }

  // static Int_t LightGray     = TColor::GetColor( "#aaaaaa" );
  // TLegendEntry* entry = leg->AddEntry("unc", "Total unc.", "f");
  // entry->SetLineWidth(0);
  // entry->SetLineColor(kWhite);
  // entry->SetFillColor(LightGray);
  // entry->SetFillStyle(3245);

  leg->Draw();


  // if (hname == "Mass_Leading_Reco_Pt400to500"){
  //   TString infotext = TString::Format("400 < p_{T} < 500 GeV");
  //   TLatex *text1 = new TLatex(3.5, 24, infotext);
  //   text1->SetNDC();
  //   text1->SetTextAlign(13);
  //   text1->SetX(0.65);
  //   text1->SetTextFont(42);
  //   text1->SetTextSize(0.05);
  //   text1->SetY(0.46);
  //   text1->Draw();
  // }

}

void DrawLegend(vector<TPolyMarker3D*> hists, TString hname)
{
  // draw a legend

  int narr = hists.size();
  float yfrac = 0.07;

  float top = 0.85;

  float ysize = yfrac*narr;
  float xleft = 0.64;
  float xright = 0.92;

  if (hname.Contains("JMS")){
    xleft = 0.22;
    xright = 0.4;
    top = 0.8;
    ysize = 0.05*narr;
  }

  xright = xleft + 0.22;
  top += 0.02;

  TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);

  for (Int_t i=0; i<narr; ++i){
    TPolyMarker3D* sh = hists[i];
    TLine *line = new TLine(0,0,0,0);
    line->SetLineColor(sh->GetMarkerColor());
    line->SetLineStyle(kSolid);

    TString legname = TString::Format("leg_entry_%i",i);
    TString legtitle;
    if(hname.Contains("JMS")){
      if(i==0) legtitle = "Best fit";
      if(i==1) legtitle = "68% CL";
      if(i==2) legtitle = "95% CL";
    }

    TLegendEntry* entry = NULL;
    int marker = sh->GetMarkerStyle();
    int lstyle = kSolid;

    if (i>0){
      entry = leg->AddEntry(sh, legtitle, "l");
      entry->SetLineWidth(1);
      entry->SetLineColor(sh->GetMarkerColor());
      entry->SetMarkerColor(sh->GetMarkerColor());
    } else {

      entry = leg->AddEntry(sh, legtitle, "p");
      entry->SetLineColor(sh->GetMarkerColor());
      entry->SetMarkerStyle(0);
      entry->SetMarkerSize(0);
      entry->SetMarkerColor(sh->GetMarkerColor());
      entry->SetLineWidth(2);
      entry->SetLineStyle(lstyle);

      entry->SetTextAlign(12);
      //entry->SetTextColor(fSampleColors.At(i));
    }

    leg->Draw();

  }
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void CompareHistStructureLocal(TH1F* h1, TH1F* h2){
  // Check if histograms have the same structure
  bool sameBins = h1->GetNbinsX() == h2->GetNbinsX();
  bool sameCenter = true;
  for(unsigned int i=1; i<h1->GetNbinsX(); i++){
    sameCenter = h1->GetBinCenter(i) == h2->GetBinCenter(i);
    if(!sameCenter) break;
  }
  if(!sameBins||!sameCenter) throw runtime_error("GetRatioLocal: Histograms do not have the same structure!");
}

TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool equal, bool isData){
  CompareHistStructureLocal(h1, h2);
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetNbinsX();
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, isData?10:1);
      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
