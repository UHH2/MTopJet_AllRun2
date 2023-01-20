#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/Plotting.h"
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

bool drawW = true;
bool drawT = true;
bool drawC = true;

void DrawLegend(vector<TPolyMarker3D*> hists, TString hname);
void CompareHistStructureLocal(TH1F* h1, TH1F* h2);
void DrawWMass(TString bin);
void DrawTau(TString year);
void WMassLabel(double x, double y, TString bin);

TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool equal, bool isData = false);

TString save_afs = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/PaperPlots/";
// TString save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/";
TString save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/peak/";

TFile* file;

int main(int argc, char* argv[]){
  cout << "Start" << endl;
  // file = new TFile("files/PaperPlots_Peaks.root", "read");

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();
  gStyle->SetErrorX(0);

  if(argc==1){
    cout << "No info on Preliminary and specific plots are given: Draw all wo. preliminiary." << endl;
  }
  else if(argc==2){
    cout << "Draw all specific plots with preliminary";
    prelim = stob(argv[1]);
    if(prelim) save_nfs.ReplaceAll("PaperPlots/", "PaperPlots/preliminary/");
  }
  else if(argc==4){
    cout << "Draw only specific plots without preliminary" << endl;
    drawW = stob(argv[1]);
    drawT = stob(argv[2]);
    drawC = stob(argv[3]);
  }
  else if(argc==5){
    cout << "Draw only specific plots with preliminary" << endl;
    prelim = stob(argv[1]);
    drawW = stob(argv[2]);
    drawT = stob(argv[3]);
    drawC = stob(argv[4]);
    if(prelim) save_nfs.ReplaceAll("PaperPlots/", "PaperPlots/preliminary/");
  }
  else{
    throw runtime_error("<E> Weird arguments ... ");
  }

  // ===========================================================================
  // === Plot WMass plots
  if(drawW){
    cout << "Start with WMass plots ..." << endl;
    file = new TFile("files/PaperPlots_Peak.root", "read");
    DrawWMass("hh");
    DrawWMass("hl");
    DrawWMass("lh");
    DrawWMass("ll");
    file->Close();
  }

  // ===========================================================================
  // === Plot tau32 plots
  if(drawT){
    cout << "Start with WMass plots ..." << endl;
    file = new TFile("files/PaperPlots.root", "read");
    DrawTau("2016");
    DrawTau("combine");
    file->Close();
  }

  // ===========================================================================
  // === JMS 2D Chi2
  if(drawC){
    cout << "Start with JMS plots ..." << endl;

    file = new TFile("files/PaperPlots_Peak.root", "read");
    TF2* chi2 = (TF2*) file->Get("Functions/JMS_Chi2");
    // cout << chi2->Eval(0,0) << endl;
    double twoD_minX, twoD_minY;
    double twoD_minZ = chi2->GetMinimumXY(twoD_minX,twoD_minY);

    TPolyMarker3D* minimum = (TPolyMarker3D*) file->Get("Graphs/JMS_nominal");
    minimum->SetMarkerColor(kBlack);
    // minimum->SetMarkerStyle(kFullCircle);
    minimum->SetMarkerStyle(34);
    minimum->SetMarkerSize(0.5);
    double x,y,z;
    minimum->GetPoint(0, x,y,z);
    cout << x << "\t" << y << "\t" << z << endl;

    TEllipse* ell = new TEllipse(twoD_minX, twoD_minY, 0.85, 0.14113, 0, 360, 111.966);
    cout << ell->GetX1() << endl;
    ell->SetLineColor(kYellow+2);
    ell->SetLineWidth(2);
    ell->SetFillStyle(0);
    // ell->SetFillColorAlpha(kBlack, 0.0);

    TPolyMarker3D* ellipse = (TPolyMarker3D*) file->Get("Graphs/JMS_ellipse");
    ellipse->SetMarkerColor(kRed-7);
    ellipse->SetMarkerStyle(kFullCircle);
    ellipse->SetMarkerSize(0.15);

    TPolyMarker3D* ellipse2 = (TPolyMarker3D*) file->Get("Graphs/JMS_ellipse_2sigma");
    ellipse2->SetMarkerColor(kRed+2);
    ellipse2->SetMarkerStyle(kFullCircle);
    ellipse2->SetMarkerSize(0.15);

    SetupCanvas(false);
    m_can->cd();
    m_can->SetLeftMargin(0.15);
    m_can->SetBottomMargin(0.15);
    m_can->SetRightMargin(0.19);

    cout << "\t ... Set main function" << endl;
    Int_t nb = 50;
    chi2->SetTitle("");
    chi2->SetFillStyle(1000);
    chi2->SetLineWidth(4);
    chi2->SetRange(-2, -2, 2, 2);
    chi2->SetContour(nb); // Contours

    cout << "\t ... Define contours" << endl;

    const Int_t Number = 4;
    Double_t Red[Number]    = { 0.99, 0.49, 0.10, 0.00};
    Double_t Green[Number]  = { 0.99, 0.80, 0.40, 0.00};
    Double_t Blue[Number]   = { 0.99, 0.88, 0.66, 0.33};
    Double_t Length[Number] = { 0.00, 0.33, 0.66, 1.00};

    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

    // m_can->SetRightMargin(0.12);
    cout << "\t ... Start plotting" << endl;

    // m_can->SetLogz();
    chi2->Draw("cont4z");
    // ell->Draw("SAME");
    cout << "\t ... Cosmetics" << endl;
    // ell->Paint("SAME");
    Cosmetics(chi2->GetHistogram(), "#it{f}^{ JEC}", 0, 0, 0, 0, false);

    cout << "\t ... Last settings" << endl;
    chi2->GetHistogram()->GetYaxis()->SetTitle("#it{f}^{ XCone}");
    chi2->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
    chi2->GetHistogram()->GetZaxis()->SetTitleOffset(1.1);
    chi2->GetHistogram()->GetZaxis()->SetTitleSize(0.06);
    chi2->GetHistogram()->GetZaxis()->SetLabelSize(0.04);
    chi2->GetHistogram()->GetZaxis()->SetTickLength(0.02);
    chi2->GetHistogram()->GetZaxis()->SetLabelOffset(0.011);
    chi2->GetHistogram()->GetZaxis()->SetTitle("#chi^{ 2}");
    chi2->GetHistogram()->GetZaxis()->CenterTitle();

    chi2->SetMinimum(100);
    gPad->RedrawAxis();
    m_can->SetTheta(90);
    m_can->SetPhi(0);
    cout << "\t ... Draw rest" << endl;
    ellipse->Draw("SAME P");
    ellipse2->Draw("SAME P");
    minimum->Draw("SAME P");
    // ell->Draw("SAME");
    // ell->Paint("SAME");
    DrawLumi(138., false, true, prelim, "JMS");

    cout << "\t ... Draw Legend" << endl;
    DrawLegend({minimum, ellipse, ellipse2}, "chi2_JMS");

    m_can->SaveAs(save_afs+"/chi2_JMS.pdf");
    m_can->SaveAs(save_nfs+"/chi2_JMS.pdf");
    file->Close();
  }
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

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
      entry->SetLineWidth(2.5);
      entry->SetLineColor(sh->GetMarkerColor());
      entry->SetMarkerColor(sh->GetMarkerColor());
      entry->SetLineWidth(3.5);
    } else {

      entry = leg->AddEntry(sh, legtitle, "p");
      entry->SetLineColor(sh->GetMarkerColor());
      entry->SetMarkerStyle(0);
      entry->SetMarkerSize(2);
      entry->SetMarkerColor(sh->GetMarkerColor());
      entry->SetLineWidth(1);
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

void DrawWMass(TString bin){
  cout << "\t ... " << bin << ": get plots" << endl;
  TH1F *wmass_data = (TH1F*) file->Get("wmass_"+bin+"__DATA");
  TH1F *wmass_ttbar = (TH1F*) file->Get("wmass_"+bin+"__TTbar");
  TH1F *wmass_JECup = (TH1F*) file->Get("wmass_"+bin+"__TTbarJECup");
  TH1F *wmass_JECdown = (TH1F*) file->Get("wmass_"+bin+"__TTbarJECdown");
  TH1F *wmass_CORup = (TH1F*) file->Get("wmass_"+bin+"__TTbarCORup");
  TH1F *wmass_CORdown = (TH1F*) file->Get("wmass_"+bin+"__TTbarCORdown");

  TH1F* ratio_data    = GetRatioLocal(wmass_data,    wmass_ttbar, false, true);
  TH1F* ratio_ttbar   = GetRatioLocal(wmass_ttbar,   wmass_ttbar, true);
  TH1F* ratio_JECup   = GetRatioLocal(wmass_JECup,   wmass_ttbar, false);
  TH1F* ratio_JECdown = GetRatioLocal(wmass_JECdown, wmass_ttbar, false);
  TH1F* ratio_CORup   = GetRatioLocal(wmass_CORup,   wmass_ttbar, false);
  TH1F* ratio_CORdown = GetRatioLocal(wmass_CORdown, wmass_ttbar, false);

  cout << "\t ... " << bin << ": draw plots" << endl;
  SetupCanvas(true);
  // void CosmeticsLocal(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_ttbar,   "#it{m}_{W} [GeV]", kBlack,  kSolid,  0, 2, true);
  Cosmetics(wmass_JECup,   "#it{m}_{W} [GeV]", kAzure-3, kSolid,  0, 2, true);
  Cosmetics(wmass_JECdown, "#it{m}_{W} [GeV]", kAzure-3, kDashed,  0, 2, true);
  Cosmetics(wmass_CORup,   "#it{m}_{W} [GeV]", kRed+1, kSolid, 0, 2, true);
  Cosmetics(wmass_CORdown, "#it{m}_{W} [GeV]", kRed+1, kDashed, 0, 2, true);
  Cosmetics(wmass_data,    "#it{m}_{W} [GeV]", kBlack,   kSolid,  8, 2, true);

  // void RatioCosmeticsLocal(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ttbar,   "#it{m}_{W} [GeV]", kBlack,  kSolid,  2, 0);
  RatioCosmetics(ratio_JECup,   "#it{m}_{W} [GeV]", kAzure-3, kSolid,  2, 0);
  RatioCosmetics(ratio_JECdown, "#it{m}_{W} [GeV]", kAzure-3, kDashed,  2, 0);
  RatioCosmetics(ratio_CORup,   "#it{m}_{W} [GeV]", kRed+1, kSolid, 2, 0);
  RatioCosmetics(ratio_CORdown, "#it{m}_{W} [GeV]", kRed+1, kDashed, 2, 0);
  RatioCosmetics(ratio_data,    "#it{m}_{W} [GeV]", kBlack,   kSolid,  2, 8);

  wmass_ttbar->GetYaxis()->SetTitleSize(0.08);
  wmass_ttbar->GetYaxis()->SetTitleOffset(1.05);
  wmass_ttbar->GetYaxis()->SetLabelOffset(0.01);
  ratio_ttbar->GetYaxis()->SetTitleSize(0.15);
  ratio_ttbar->GetYaxis()->SetTitleOffset(0.58);
  ratio_ttbar->GetYaxis()->SetLabelOffset(0.011);

  m_rp1_top->cd();
  wmass_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_ttbar->GetMaximum()*1.2);
  wmass_ttbar->Draw("hist ][");
  wmass_JECup->Draw("hist same ][");
  wmass_JECdown->Draw("hist same ][");
  wmass_CORup->Draw("hist same ][");
  wmass_CORdown->Draw("hist same ][");
  wmass_ttbar->Draw("axis same ][");
  wmass_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_data, wmass_ttbar, wmass_JECup, wmass_JECdown, wmass_CORup, wmass_CORdown}, "wmass_"+bin);
  // leg->SetTextSize(0.065);

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, prelim, "wmass_"+bin);
  if(bin.Contains("hl")||bin.Contains("hh")) WMassLabel(0.68, 0.87, bin);
  else                                       WMassLabel(0.22, 0.78, bin);

  m_rp1->cd();
  ratio_ttbar->Draw("hist ][");
  ratio_JECup->Draw("hist same ][");
  ratio_CORup->Draw("hist same ][");
  ratio_JECdown->Draw("hist same ][");
  ratio_CORdown->Draw("hist same ][");
  ratio_data->Draw("pe same ][");
  gPad->RedrawAxis();

  cout << save_nfs+"wmass_"+bin+".pdf" << endl;
  m_can->Print(save_afs+"wmass_"+bin+".pdf");
  m_can->Print(save_nfs+"wmass_"+bin+".pdf");
}

void DrawTau(TString year){
  if(!(year.EqualTo("2016") || year.EqualTo("combine"))) throw runtime_error("Wrong year for Tau32 plot");
  bool is16 = year.EqualTo("2016");
  cout << "\t ... "+year+": get plots" << endl;
  TH1F *tau32_data, *tau32_ttbar, *tau32_upsqrt2, *tau32_up2, *tau32_up4, *tau32_downsqrt2, *tau32_down2, *tau32_down4;
  tau32_data      = (TH1F*) file->Get("tau32__DATA__"+year);
  tau32_ttbar     = (TH1F*) file->Get("tau32__TTbar__"+year);
  tau32_up2       = (TH1F*) file->Get("tau32__TTbarFSRup2__"+year);
  tau32_down2     = (TH1F*) file->Get("tau32__TTbarFSRdown2__"+year);
  if(!is16){
    tau32_upsqrt2   = (TH1F*) file->Get("tau32__TTbarFSRupsqrt2__"+year);
    tau32_downsqrt2 = (TH1F*) file->Get("tau32__TTbarFSRdownsqrt2__"+year);
    tau32_up4       = (TH1F*) file->Get("tau32__TTbarFSRup4__"+year);
    tau32_down4     = (TH1F*) file->Get("tau32__TTbarFSRdown4__"+year);
  }

  TH1F *ratio_data, *ratio_ttbar, *ratio_upsqrt2, *ratio_up2, *ratio_up4, *ratio_downsqrt2, *ratio_down2, *ratio_down4;
  ratio_data      = GetRatioLocal(tau32_data,      tau32_ttbar, false, true);
  ratio_ttbar     = GetRatioLocal(tau32_ttbar,     tau32_ttbar, true);
  ratio_up2       = GetRatioLocal(tau32_up2,       tau32_ttbar, false);
  ratio_down2     = GetRatioLocal(tau32_down2,     tau32_ttbar, false);
  if(!is16){
    ratio_upsqrt2   = GetRatioLocal(tau32_upsqrt2,   tau32_ttbar, false);
    ratio_downsqrt2 = GetRatioLocal(tau32_downsqrt2, tau32_ttbar, false);
    ratio_up4       = GetRatioLocal(tau32_up4,       tau32_ttbar, false);
    ratio_down4     = GetRatioLocal(tau32_down4,     tau32_ttbar, false);
  }

  cout << "\t ... combine: draw plots" << endl;
  SetupCanvas(true);
  int styleNom = kSolid;
  int styleUp = kSolid;
  int styleDown = kDashed;
  Cosmetics(tau32_ttbar,      "#tau_{32}", kBlack,  styleNom,  0, 3, true);
  Cosmetics(is16?tau32_up2:tau32_up4,     "#tau_{32}", is16?kAzure+7:kRed+1, styleUp,  0, 3, true);
  Cosmetics(is16?tau32_down2:tau32_down4, "#tau_{32}", is16?kAzure+7:kRed+1, styleDown, 0, 3, true);
  Cosmetics(tau32_data,       "#tau_{32}", kBlack,   kSolid,  8, 2, true);

  RatioCosmetics(ratio_ttbar, "#tau_{32}", kBlack,  styleNom,  3, 0);
  RatioCosmetics(is16?ratio_up2:ratio_up4,     "#tau_{32}", is16?kAzure+7:kRed+1, styleUp,  3, 0);
  RatioCosmetics(is16?ratio_down2:ratio_down4, "#tau_{32}", is16?kAzure+7:kRed+1, styleDown, 3, 0);
  RatioCosmetics(ratio_data,  "#tau_{32}", kBlack,   kSolid,  2, 8);

  tau32_ttbar->GetYaxis()->SetTitleSize(0.08);
  tau32_ttbar->GetYaxis()->SetTitleOffset(1.05);
  tau32_ttbar->GetYaxis()->SetLabelOffset(0.01);
  ratio_ttbar->GetYaxis()->SetTitleSize(0.15);
  ratio_ttbar->GetYaxis()->SetTitleOffset(0.58);
  ratio_ttbar->GetYaxis()->SetLabelOffset(0.011);

  m_rp1_top->cd();
  tau32_ttbar->GetYaxis()->SetRangeUser(0, tau32_ttbar->GetMaximum()*1.25);
  tau32_ttbar->Draw("hist ][");
  if(is16){
    tau32_up2->Draw("hist same ][");
    tau32_down2->Draw("hist same ][");
  }
  else{
    tau32_up4->Draw("hist same ][");
    tau32_down4->Draw("hist same ][");
  }
  tau32_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend(
    {
      tau32_data,
      tau32_ttbar,
      is16?tau32_up2:tau32_up4,
      is16?tau32_down2:tau32_down4
    },
    "tau32_"+year
  );

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(year.EqualTo("2016")?36.3:101, true, true, prelim);

  m_rp1->cd();
  ratio_ttbar->Draw("hist ][");
  if(is16){
    ratio_up2->Draw("hist same ][");
    ratio_down2->Draw("hist same ][");
  }
  else{
    ratio_up4->Draw("hist same ][");
    ratio_down4->Draw("hist same ][");
  }
  ratio_data->Draw("pe same ][");
  gPad->RedrawAxis();
  m_can->SaveAs(save_afs+"tau32_"+year+".pdf");
  m_can->SaveAs(save_nfs+"tau32_"+year+".pdf");
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

// =============================================================================
// ===                                                                       ===
// =============================================================================

void WMassLabel(double x, double y, TString bin){
  double tsize = 0.06;
  TString symbol = (bin.BeginsWith("h"))?">":"<";
  TString wmasslabel = "#it{p}_{T}^{W} "+symbol+" 300 GeV";
  TLatex *t = new TLatex(3.5, 24, wmasslabel);
  t->SetNDC();
  t->SetTextAlign(13);
  t->SetTextFont(42);
  t->SetTextSize(tsize);
  t->SetX(x);
  t->SetY(y);
  t->Draw();

  symbol = (bin.EndsWith("h"))?">":"<";
  TString subtext = "#it{r}_{p_{T}}  "+symbol+" 0.7";
  TLatex *t2 = new TLatex(3.5, 24, subtext);
  t2->SetNDC();
  t2->SetTextAlign(13);
  t2->SetX(x);
  t2->SetTextFont(42);
  t2->SetTextSize(tsize);
  t2->SetY(y-1.4*tsize);
  t2->Draw();

}
