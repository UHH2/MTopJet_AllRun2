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

void SetupCanvas(bool bPlotRatio=false);
void HistCosmetics(TH1* hist);
void PortraitCosmetics(TH1* hist, bool bPlotRatio);
void YieldCosmetics(TH1* hist);
void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width);
void Cosmetics(TH1* hist, TString xtitle, int color, int style, int marker, int width, bool bPlotRatio);
void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
void DrawLegend(vector<TH1F*> hists, TString hname);
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
TString save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/";

int main(int argc, char* argv[]){

  TFile *file = new TFile("files/PaperPlotsSort.root", "read");

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();
  gStyle->SetErrorX(0);

  // ===========================================================================
  // === Plot WMass plots
  cout << "Start with WMass plots ..." << endl;

  // === hh ===
  cout << "\t ... hh: get plots" << endl;

  TH1F *wmass_hh_data = (TH1F*) file->Get("Hists/wmass_hh__DATA");
  TH1F *wmass_hh_ttbar = (TH1F*) file->Get("Hists/wmass_hh__TTbar");
  TH1F *wmass_hh_JECup = (TH1F*) file->Get("Hists/wmass_hh__TTbarJECup");
  TH1F *wmass_hh_JECdown = (TH1F*) file->Get("Hists/wmass_hh__TTbarJECdown");
  TH1F *wmass_hh_CORup = (TH1F*) file->Get("Hists/wmass_hh__TTbarCORup");
  TH1F *wmass_hh_CORdown = (TH1F*) file->Get("Hists/wmass_hh__TTbarCORdown");

  TH1F* ratio_hh_data    = GetRatioLocal(wmass_hh_data,    wmass_hh_ttbar, false, true);
  TH1F* ratio_hh_ttbar   = GetRatioLocal(wmass_hh_ttbar,   wmass_hh_ttbar, true);
  TH1F* ratio_hh_JECup   = GetRatioLocal(wmass_hh_JECup,   wmass_hh_ttbar, false);
  TH1F* ratio_hh_JECdown = GetRatioLocal(wmass_hh_JECdown, wmass_hh_ttbar, false);
  TH1F* ratio_hh_CORup   = GetRatioLocal(wmass_hh_CORup,   wmass_hh_ttbar, false);
  TH1F* ratio_hh_CORdown = GetRatioLocal(wmass_hh_CORdown, wmass_hh_ttbar, false);

  cout << "\t ... hh: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_hh_ttbar,      "#it{m}_{W}", 14,  kSolid,  0, 2, true);
  Cosmetics(wmass_hh_JECup,        "#it{m}_{W}", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(wmass_hh_JECdown,        "#it{m}_{W}", kAzure+7, kDashed,  0, 2, true);
  Cosmetics(wmass_hh_CORup,      "#it{m}_{W}", 810, kSolid, 0, 2, true);
  Cosmetics(wmass_hh_CORdown,      "#it{m}_{W}", 810, kDashed, 0, 2, true);
  Cosmetics(wmass_hh_data,       "#it{m}_{W}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_hh_ttbar, "#it{m}_{W}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_hh_JECup,   "#it{m}_{W}", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_hh_JECdown,   "#it{m}_{W}", kAzure+7, kDashed,  2, 0);
  RatioCosmetics(ratio_hh_CORup, "#it{m}_{W}", 810, kSolid, 2, 0);
  RatioCosmetics(ratio_hh_CORdown, "#it{m}_{W}", 810, kDashed, 2, 0);
  RatioCosmetics(ratio_hh_data,  "#it{m}_{W}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  wmass_hh_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_hh_ttbar->GetMaximum()*1.2);
  wmass_hh_ttbar->Draw("hist ][");
  wmass_hh_JECup->Draw("hist same ][");
  wmass_hh_JECdown->Draw("hist same ][");
  wmass_hh_CORup->Draw("hist same ][");
  wmass_hh_CORdown->Draw("hist same ][");
  wmass_hh_ttbar->Draw("axis same ][");
  wmass_hh_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_hh_data, wmass_hh_ttbar, wmass_hh_JECup, wmass_hh_JECdown, wmass_hh_CORup, wmass_hh_CORdown}, "wmass_hh");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, true);

  m_rp1->cd();
  ratio_hh_ttbar->Draw("hist ][");
  ratio_hh_JECup->Draw("hist same ][");
  ratio_hh_CORup->Draw("hist same ][");
  ratio_hh_JECdown->Draw("hist same ][");
  ratio_hh_CORdown->Draw("hist same ][");
  ratio_hh_data->Draw("pe same ][");
  gPad->RedrawAxis();

  m_can->Print(save_afs+"wmass_hh.pdf");
  m_can->Print(save_nfs+"wmass_hh.pdf");

  // === hl ===
  cout << "\t ... hl: get plots" << endl;

  TH1F *wmass_hl_data = (TH1F*) file->Get("Hists/wmass_hl__DATA");
  TH1F *wmass_hl_ttbar = (TH1F*) file->Get("Hists/wmass_hl__TTbar");
  TH1F *wmass_hl_JECup = (TH1F*) file->Get("Hists/wmass_hl__TTbarJECup");
  TH1F *wmass_hl_JECdown = (TH1F*) file->Get("Hists/wmass_hl__TTbarJECdown");
  TH1F *wmass_hl_CORup = (TH1F*) file->Get("Hists/wmass_hl__TTbarCORup");
  TH1F *wmass_hl_CORdown = (TH1F*) file->Get("Hists/wmass_hl__TTbarCORdown");

  TH1F* ratio_hl_data    = GetRatioLocal(wmass_hl_data,    wmass_hl_ttbar, false, true);
  TH1F* ratio_hl_ttbar   = GetRatioLocal(wmass_hl_ttbar,   wmass_hl_ttbar, true);
  TH1F* ratio_hl_JECup   = GetRatioLocal(wmass_hl_JECup,   wmass_hl_ttbar, false);
  TH1F* ratio_hl_JECdown = GetRatioLocal(wmass_hl_JECdown, wmass_hl_ttbar, false);
  TH1F* ratio_hl_CORup   = GetRatioLocal(wmass_hl_CORup,   wmass_hl_ttbar, false);
  TH1F* ratio_hl_CORdown = GetRatioLocal(wmass_hl_CORdown, wmass_hl_ttbar, false);

  cout << "\t ... hl: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_hl_ttbar,      "#it{m}_{W}", 14,  kSolid,  0, 2, true);
  Cosmetics(wmass_hl_JECup,        "#it{m}_{W}", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(wmass_hl_JECdown,        "#it{m}_{W}", kAzure+7, kDashed,  0, 2, true);
  Cosmetics(wmass_hl_CORup,      "#it{m}_{W}", 810, kSolid, 0, 2, true);
  Cosmetics(wmass_hl_CORdown,      "#it{m}_{W}", 810, kDashed, 0, 2, true);
  Cosmetics(wmass_hl_data,       "#it{m}_{W}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_hl_ttbar, "#it{m}_{W}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_hl_JECup,   "#it{m}_{W}", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_hl_JECdown,   "#it{m}_{W}", kAzure+7, kDashed,  2, 0);
  RatioCosmetics(ratio_hl_CORup, "#it{m}_{W}", 810, kSolid, 2, 0);
  RatioCosmetics(ratio_hl_CORdown, "#it{m}_{W}", 810, kDashed, 2, 0);
  RatioCosmetics(ratio_hl_data,  "#it{m}_{W}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  wmass_hl_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_hl_ttbar->GetMaximum()*1.2);
  wmass_hl_ttbar->Draw("hist ][");
  wmass_hl_JECup->Draw("hist same ][");
  wmass_hl_JECdown->Draw("hist same ][");
  wmass_hl_CORup->Draw("hist same ][");
  wmass_hl_CORdown->Draw("hist same ][");
  wmass_hl_ttbar->Draw("axis same ][");
  wmass_hl_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_hl_data, wmass_hl_ttbar, wmass_hl_JECup, wmass_hl_JECdown, wmass_hl_CORup, wmass_hl_CORdown}, "wmass_hl");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, true);

  m_rp1->cd();
  ratio_hl_ttbar->Draw("hist ][");
  ratio_hl_JECup->Draw("hist same ][");
  ratio_hl_CORup->Draw("hist same ][");
  ratio_hl_JECdown->Draw("hist same ][");
  ratio_hl_CORdown->Draw("hist same ][");
  ratio_hl_data->Draw("pe same ][");
  gPad->RedrawAxis();

  m_can->SaveAs(save_afs+"wmass_hl.pdf");
  m_can->SaveAs(save_nfs+"wmass_hl.pdf");

  // === lh ===
  cout << "\t ... lh: get plots" << endl;

  TH1F *wmass_lh_data = (TH1F*) file->Get("Hists/wmass_lh__DATA");
  TH1F *wmass_lh_ttbar = (TH1F*) file->Get("Hists/wmass_lh__TTbar");
  TH1F *wmass_lh_JECup = (TH1F*) file->Get("Hists/wmass_lh__TTbarJECup");
  TH1F *wmass_lh_JECdown = (TH1F*) file->Get("Hists/wmass_lh__TTbarJECdown");
  TH1F *wmass_lh_CORup = (TH1F*) file->Get("Hists/wmass_lh__TTbarCORup");
  TH1F *wmass_lh_CORdown = (TH1F*) file->Get("Hists/wmass_lh__TTbarCORdown");

  TH1F* ratio_lh_data    = GetRatioLocal(wmass_lh_data,    wmass_lh_ttbar, false, true);
  TH1F* ratio_lh_ttbar   = GetRatioLocal(wmass_lh_ttbar,   wmass_lh_ttbar, true);
  TH1F* ratio_lh_JECup   = GetRatioLocal(wmass_lh_JECup,   wmass_lh_ttbar, false);
  TH1F* ratio_lh_JECdown = GetRatioLocal(wmass_lh_JECdown, wmass_lh_ttbar, false);
  TH1F* ratio_lh_CORup   = GetRatioLocal(wmass_lh_CORup,   wmass_lh_ttbar, false);
  TH1F* ratio_lh_CORdown = GetRatioLocal(wmass_lh_CORdown, wmass_lh_ttbar, false);

  cout << "\t ... lh: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_lh_ttbar,      "#it{m}_{W}", 14,  kSolid,  0, 2, true);
  Cosmetics(wmass_lh_JECup,        "#it{m}_{W}", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(wmass_lh_JECdown,        "#it{m}_{W}", kAzure+7, kDashed,  0, 2, true);
  Cosmetics(wmass_lh_CORup,      "#it{m}_{W}", 810, kSolid, 0, 2, true);
  Cosmetics(wmass_lh_CORdown,      "#it{m}_{W}", 810, kDashed, 0, 2, true);
  Cosmetics(wmass_lh_data,       "#it{m}_{W}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_lh_ttbar, "#it{m}_{W}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_lh_JECup,   "#it{m}_{W}", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_lh_JECdown,   "#it{m}_{W}", kAzure+7, kDashed,  2, 0);
  RatioCosmetics(ratio_lh_CORup, "#it{m}_{W}", 810, kSolid, 2, 0);
  RatioCosmetics(ratio_lh_CORdown, "#it{m}_{W}", 810, kDashed, 2, 0);
  RatioCosmetics(ratio_lh_data,  "#it{m}_{W}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  wmass_lh_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_lh_ttbar->GetMaximum()*1.2);
  wmass_lh_ttbar->Draw("hist ][");
  wmass_lh_JECup->Draw("hist same ][");
  wmass_lh_JECdown->Draw("hist same ][");
  wmass_lh_CORup->Draw("hist same ][");
  wmass_lh_CORdown->Draw("hist same ][");
  wmass_lh_ttbar->Draw("axis same ][");
  wmass_lh_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_lh_data, wmass_lh_ttbar, wmass_lh_JECup, wmass_lh_JECdown, wmass_lh_CORup, wmass_lh_CORdown}, "wmass_lh");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, true);

  m_rp1->cd();
  ratio_lh_ttbar->Draw("hist ][");
  ratio_lh_JECup->Draw("hist same ][");
  ratio_lh_CORup->Draw("hist same ][");
  ratio_lh_JECdown->Draw("hist same ][");
  ratio_lh_CORdown->Draw("hist same ][");
  ratio_lh_data->Draw("pe same ][");
  gPad->RedrawAxis();

  m_can->SaveAs(save_afs+"wmass_lh.pdf");
  m_can->SaveAs(save_nfs+"wmass_lh.pdf");

  // === ll ===
  cout << "\t ... ll: get plots" << endl;

  TH1F *wmass_ll_data = (TH1F*) file->Get("Hists/wmass_ll__DATA");
  TH1F *wmass_ll_ttbar = (TH1F*) file->Get("Hists/wmass_ll__TTbar");
  TH1F *wmass_ll_JECup = (TH1F*) file->Get("Hists/wmass_ll__TTbarJECup");
  TH1F *wmass_ll_JECdown = (TH1F*) file->Get("Hists/wmass_ll__TTbarJECdown");
  TH1F *wmass_ll_CORup = (TH1F*) file->Get("Hists/wmass_ll__TTbarCORup");
  TH1F *wmass_ll_CORdown = (TH1F*) file->Get("Hists/wmass_ll__TTbarCORdown");

  TH1F* ratio_ll_data    = GetRatioLocal(wmass_ll_data,    wmass_ll_ttbar, false, true);
  TH1F* ratio_ll_ttbar   = GetRatioLocal(wmass_ll_ttbar,   wmass_ll_ttbar, true);
  TH1F* ratio_ll_JECup   = GetRatioLocal(wmass_ll_JECup,   wmass_ll_ttbar, false);
  TH1F* ratio_ll_JECdown = GetRatioLocal(wmass_ll_JECdown, wmass_ll_ttbar, false);
  TH1F* ratio_ll_CORup   = GetRatioLocal(wmass_ll_CORup,   wmass_ll_ttbar, false);
  TH1F* ratio_ll_CORdown = GetRatioLocal(wmass_ll_CORdown, wmass_ll_ttbar, false);

  cout << "\t ... ll: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_ll_ttbar,      "#it{m}_{W}", 14,  kSolid,  0, 2, true);
  Cosmetics(wmass_ll_JECup,        "#it{m}_{W}", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(wmass_ll_JECdown,        "#it{m}_{W}", kAzure+7, kDashed,  0, 2, true);
  Cosmetics(wmass_ll_CORup,      "#it{m}_{W}", 810, kSolid, 0, 2, true);
  Cosmetics(wmass_ll_CORdown,      "#it{m}_{W}", 810, kDashed, 0, 2, true);
  Cosmetics(wmass_ll_data,       "#it{m}_{W}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ll_ttbar, "#it{m}_{W}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_ll_JECup,   "#it{m}_{W}", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_ll_JECdown,   "#it{m}_{W}", kAzure+7, kDashed,  2, 0);
  RatioCosmetics(ratio_ll_CORup, "#it{m}_{W}", 810, kSolid, 2, 0);
  RatioCosmetics(ratio_ll_CORdown, "#it{m}_{W}", 810, kDashed, 2, 0);
  RatioCosmetics(ratio_ll_data,  "#it{m}_{W}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  wmass_ll_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_ll_ttbar->GetMaximum()*1.2);
  wmass_ll_ttbar->Draw("hist ][");
  wmass_ll_JECup->Draw("hist same ][");
  wmass_ll_JECdown->Draw("hist same ][");
  wmass_ll_CORup->Draw("hist same ][");
  wmass_ll_CORdown->Draw("hist same ][");
  wmass_ll_ttbar->Draw("axis same ][");
  wmass_ll_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_ll_data, wmass_ll_ttbar, wmass_ll_JECup, wmass_ll_JECdown, wmass_ll_CORup, wmass_ll_CORdown}, "wmass_ll");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, true);

  m_rp1->cd();
  ratio_ll_ttbar->Draw("hist ][");
  ratio_ll_JECup->Draw("hist same ][");
  ratio_ll_CORup->Draw("hist same ][");
  ratio_ll_JECdown->Draw("hist same ][");
  ratio_ll_CORdown->Draw("hist same ][");
  ratio_ll_data->Draw("SAME PE");
  gPad->RedrawAxis();

  m_can->SaveAs(save_afs+"wmass_ll.pdf");
  m_can->SaveAs(save_nfs+"wmass_ll.pdf");

  // ===========================================================================
  // === Plot tau32 plots
  cout << "Start with WMass plots ..." << endl;

  // === 2016
  cout << "\t ... 2016: get plots" << endl;

  TH1F *tau32_data_2016  = (TH1F*) file->Get("Hists/tau32__DATA__2016");
  TH1F *tau32_ttbar_2016 = (TH1F*) file->Get("Hists/tau32__TTbar__2016");
  TH1F *tau32_up2_2016   = (TH1F*) file->Get("Hists/tau32__TTbarFSRup2__2016");
  TH1F *tau32_down2_2016 = (TH1F*) file->Get("Hists/tau32__TTbarFSRdown2__2016");

  TH1F *ratio_data_2016  = GetRatioLocal(tau32_data_2016,  tau32_ttbar_2016, false, true);
  TH1F *ratio_ttbar_2016 = GetRatioLocal(tau32_ttbar_2016, tau32_ttbar_2016, true);
  TH1F *ratio_up2_2016   = GetRatioLocal(tau32_up2_2016,   tau32_ttbar_2016, false);
  TH1F *ratio_down2_2016 = GetRatioLocal(tau32_down2_2016, tau32_ttbar_2016, false);

  cout << "\t ... 2016: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(tau32_ttbar_2016,      "#tau_{32}", 14,  kSolid,  0, 2, true);
  Cosmetics(tau32_up2_2016,        "#tau_{32}", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(tau32_down2_2016,      "#tau_{32}", kAzure+7, kDashed, 0, 2, true);
  Cosmetics(tau32_data_2016,       "#tau_{32}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ttbar_2016, "#tau_{32}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_up2_2016,   "#tau_{32}", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_down2_2016, "#tau_{32}", kAzure+7, kDashed, 2, 0);
  RatioCosmetics(ratio_data_2016,  "#tau_{32}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  tau32_ttbar_2016->GetYaxis()->SetRangeUser(0., tau32_ttbar_2016->GetMaximum()*1.25);
  tau32_ttbar_2016->Draw("hist ][");
  tau32_up2_2016->Draw("hist same ][");
  tau32_down2_2016->Draw("hist same ][");
  tau32_data_2016->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({tau32_data_2016, tau32_ttbar_2016, tau32_up2_2016, tau32_down2_2016}, "tau32_2016");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(36., true, true, true);

  m_rp1->cd();
  ratio_ttbar_2016->Draw("hist ][");
  ratio_up2_2016->Draw("hist same ][");
  ratio_down2_2016->Draw("hist same ][");
  ratio_data_2016->Draw("pe same ][");
  gPad->RedrawAxis();
  m_can->SaveAs(save_afs+"tau32_2016.pdf");
  m_can->SaveAs(save_nfs+"tau32_2016.pdf");

  // === combine
  cout << "\t ... combine: get plots" << endl;

  TH1F *tau32_data_combine      = (TH1F*) file->Get("Hists/tau32__DATA__combine");
  TH1F *tau32_ttbar_combine     = (TH1F*) file->Get("Hists/tau32__TTbar__combine");
  TH1F *tau32_upsqrt2_combine   = (TH1F*) file->Get("Hists/tau32__TTbarFSRupsqrt2__combine");
  TH1F *tau32_up2_combine       = (TH1F*) file->Get("Hists/tau32__TTbarFSRup2__combine");
  TH1F *tau32_up4_combine       = (TH1F*) file->Get("Hists/tau32__TTbarFSRup4__combine");
  TH1F *tau32_downsqrt2_combine = (TH1F*) file->Get("Hists/tau32__TTbarFSRdownsqrt2__combine");
  TH1F *tau32_down2_combine     = (TH1F*) file->Get("Hists/tau32__TTbarFSRdown2__combine");
  TH1F *tau32_down4_combine     = (TH1F*) file->Get("Hists/tau32__TTbarFSRdown4__combine");

  TH1F *ratio_data_combine      = GetRatioLocal(tau32_data_combine,      tau32_ttbar_combine, false, true);
  TH1F *ratio_ttbar_combine     = GetRatioLocal(tau32_ttbar_combine,     tau32_ttbar_combine, true);
  TH1F *ratio_upsqrt2_combine   = GetRatioLocal(tau32_upsqrt2_combine,   tau32_ttbar_combine, false);
  TH1F *ratio_up2_combine       = GetRatioLocal(tau32_up2_combine,       tau32_ttbar_combine, false);
  TH1F *ratio_up4_combine       = GetRatioLocal(tau32_up4_combine,       tau32_ttbar_combine, false);
  TH1F *ratio_downsqrt2_combine = GetRatioLocal(tau32_downsqrt2_combine, tau32_ttbar_combine, false);
  TH1F *ratio_down2_combine     = GetRatioLocal(tau32_down2_combine,     tau32_ttbar_combine, false);
  TH1F *ratio_down4_combine     = GetRatioLocal(tau32_down4_combine,     tau32_ttbar_combine, false);

  cout << "\t ... combine: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(tau32_ttbar_combine,      "#tau_{32}", 14,  kSolid,  0, 2, true);
  Cosmetics(tau32_up4_combine,        "#tau_{32}", 810, kSolid,  0, 2, true);
  Cosmetics(tau32_down4_combine,      "#tau_{32}", 810, kDashed, 0, 2, true);
  Cosmetics(tau32_data_combine,       "#tau_{32}", 1,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ttbar_combine, "#tau_{32}", 14,  kSolid,  2, 0);
  RatioCosmetics(ratio_up4_combine,   "#tau_{32}", 810, kSolid,  2, 0);
  RatioCosmetics(ratio_down4_combine, "#tau_{32}", 810, kDashed, 2, 0);
  RatioCosmetics(ratio_data_combine,  "#tau_{32}", 1,   kSolid,  2, 8);

  m_rp1_top->cd();
  tau32_ttbar_combine->GetYaxis()->SetRangeUser(0, tau32_ttbar_combine->GetMaximum()*1.25);
  tau32_ttbar_combine->Draw("hist ][");
  tau32_up4_combine->Draw("hist same ][");
  tau32_down4_combine->Draw("hist same ][");
  tau32_data_combine->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({tau32_data_combine, tau32_ttbar_combine, tau32_up4_combine, tau32_down4_combine}, "tau32_combine");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(101., true, true, true);

  m_rp1->cd();
  ratio_ttbar_combine->Draw("hist ][");
  ratio_up4_combine->Draw("hist same ][");
  ratio_down4_combine->Draw("hist same ][");
  ratio_data_combine->Draw("pe same ][");
  gPad->RedrawAxis();
  m_can->SaveAs(save_afs+"tau32_combine.pdf");
  m_can->SaveAs(save_nfs+"tau32_combine.pdf");

  // ===========================================================================
  // === JMS 2D Chi2
  cout << "Start with JMS plots ..." << endl;

  TF2* chi2 = (TF2*) file->Get("Functions/JMS_Chi2");
  cout << chi2->Eval(0,0) << endl;

  TPolyMarker3D* minimum = (TPolyMarker3D*) file->Get("Graphs/JMS_nominal");
  minimum->SetMarkerColor(kBlack);
  minimum->SetMarkerStyle(kFullCircle);
  minimum->SetMarkerSize(0.4);

  TPolyMarker3D* ellipse = (TPolyMarker3D*) file->Get("Graphs/JMS_ellipse");
  ellipse->SetMarkerColor(kRed+2);
  ellipse->SetMarkerStyle(kFullCircle);
  ellipse->SetMarkerSize(0.1);

  SetupCanvas(false);
  m_can->cd();
  m_can->SetLeftMargin(0.20);
  m_can->SetBottomMargin(0.15);
  m_can->SetRightMargin(0.15);

  cout << "\t ... Set Minimum" << endl;
  minimum->SetMarkerColor(kBlack);
  minimum->SetMarkerStyle(kFullCircle);
  minimum->SetMarkerSize(0.4);

  cout << "\t ... Set Ellipse" << endl;
  ellipse->SetMarkerColor(kRed);
  ellipse->SetMarkerStyle(kFullCircle);
  ellipse->SetMarkerSize(0.1);

  cout << "\t ... Set main function" << endl;
  Int_t nb = 50;
  chi2->SetTitle("");
  chi2->SetFillStyle(1000);
  chi2->SetLineWidth(1);
  chi2->SetRange(-2, -2, 2, 2);
  chi2->SetContour(nb); // Contours

  cout << "\t ... Define contours" << endl;

  const Int_t Number = 4;
  Double_t Red[Number]    = { 0.99, 0.49, 0.00, 0.00};
  Double_t Green[Number]  = { 0.99, 0.80, 0.40, 0.00};
  Double_t Blue[Number]   = { 0.99, 0.88, 0.66, 0.33};
  Double_t Length[Number] = { 0.00, 0.30, 0.70, 1.00};

  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  // m_can->SetRightMargin(0.12);
  cout << "\t ... Start plotting" << endl;

  m_can->SetLogz();
  chi2->Draw("cont4z");

  cout << "\t ... Cosmetics" << endl;
  Cosmetics(chi2->GetHistogram(), "f^{JEC}", 0, 0, 0, 0, false);

  cout << "\t ... Last settings" << endl;
  chi2->GetHistogram()->GetXaxis()->SetRangeUser(-2, 2);
  chi2->GetHistogram()->GetYaxis()->SetRangeUser(-2, 2);
  chi2->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
  chi2->GetHistogram()->GetZaxis()->CenterTitle();
  // chi2->SetMinimum(95);
  gPad->RedrawAxis();
  m_can->SetTheta(90);
  m_can->SetPhi(0);
  cout << "\t ... Draw rrest" << endl;
  ellipse->Draw("SAME P");
  minimum->Draw("SAME P");
  DrawLumi(138., false, true, true);
  m_can->SaveAs(save_afs+"/chi2_JMS.pdf");
  m_can->SaveAs(save_nfs+"/chi2_JMS.pdf");
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

void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre)
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
  // if (hname == "wmass"){
  //   top += 0.02;
  //   xleft = 0.73;
  // }

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
      if(i==4) legtitle = "t#bar{t} f^{COR} = +1";
      if(i==5) legtitle = "t#bar{t} f^{COR} = -1";
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
