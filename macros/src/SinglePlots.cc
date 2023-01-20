#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
// #include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"
#include "../include/tdrstyle_all.h"

// #include "../CovMatrices/JMS/CollectCovHeaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>
#include "TSystem.h"

using namespace std;

bool debug = true;

void DrawWMassJER(TString bin);
void CompareHistStructureLocal(TH1F* h1, TH1F* h2);
TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool equal, bool isData = false);

TString save_afs = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/";
TString save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/";
TString innfs = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/";

TFile* file;

int main(int argc, char* argv[]){

  gErrorIgnoreLevel = kError;
  SetupGlobalStyle();

  file = new TFile("files/WMassPlots.root", "read");
  DrawWMassJER("hh");
  DrawWMassJER("hl");
  DrawWMassJER("lh");
  DrawWMassJER("ll");

}

// ===========================================================================
// ===========================================================================
// ===========================================================================

void DrawWMassJER(TString bin){
  cout << "\t ... " << bin << ": get plots" << endl;
  TH1F *wmass_data = (TH1F*) file->Get("wmass_"+bin+"__DATA__combine");
  TH1F *wmass_ttbar = (TH1F*) file->Get("wmass_"+bin+"__TTbar__combine");
  TH1F *wmass_JERup = (TH1F*) file->Get("wmass_"+bin+"__TTbarJERup__combine");
  TH1F *wmass_JERdown = (TH1F*) file->Get("wmass_"+bin+"__TTbarJERdown__combine");

  TH1F* ratio_data    = GetRatioLocal(wmass_data,    wmass_ttbar, false, true);
  TH1F* ratio_ttbar   = GetRatioLocal(wmass_ttbar,   wmass_ttbar, true);
  TH1F* ratio_JERup   = GetRatioLocal(wmass_JERup,   wmass_ttbar, false);
  TH1F* ratio_JERdown = GetRatioLocal(wmass_JERdown, wmass_ttbar, false);

  cout << "\t ... " << bin << ": draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(wmass_ttbar,   "#it{m}_{W} [GeV]", kBlack,  kSolid,  0, 2, true);
  Cosmetics(wmass_JERup,   "#it{m}_{W} [GeV]", kTeal+5, kSolid,  0, 2, true);
  Cosmetics(wmass_JERdown, "#it{m}_{W} [GeV]", kTeal+5, kDashed,  0, 2, true);
  Cosmetics(wmass_data,    "#it{m}_{W} [GeV]", kBlack,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ttbar,   "#it{m}_{W} [GeV]", kBlack,  kSolid,  2, 0);
  RatioCosmetics(ratio_JERup,   "#it{m}_{W} [GeV]", kTeal+5, kSolid,  2, 0);
  RatioCosmetics(ratio_JERdown, "#it{m}_{W} [GeV]", kTeal+5, kDashed,  2, 0);
  RatioCosmetics(ratio_data,    "#it{m}_{W} [GeV]", kBlack,   kSolid,  2, 8);

  m_rp1_top->cd();
  wmass_ttbar->GetYaxis()->SetRangeUser(0.001, wmass_ttbar->GetMaximum()*1.2);
  wmass_ttbar->Draw("hist ][");
  wmass_JERup->Draw("hist same ][");
  wmass_JERdown->Draw("hist same ][");
  wmass_ttbar->Draw("axis same ][");
  wmass_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({wmass_data, wmass_ttbar, wmass_JERup, wmass_JERdown}, "wmass_"+bin);

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, false, "wmass_"+bin);

  m_rp1->cd();
  ratio_ttbar->Draw("hist ][");
  ratio_JERup->Draw("hist same ][");
  ratio_JERdown->Draw("hist same ][");
  ratio_data->Draw("pe same ][");
  gPad->RedrawAxis();

  m_can->Print(save_nfs+"JMS/JER/wmass_"+bin+".pdf");
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

/*
TString name = "mjet_flavor";
TString hist = "JetMassScaleHists/hadjet_jms_mass";

vector<TH1F*> hists_muon_nominal = get_all_hists(ttbar_muon, hist);
vector<TH1F*> hists_elec_nominal = get_all_hists(ttbar_elec, hist);
vector<TH1F*> hists_nominal = AddHists(hists_muon_nominal, hists_elec_nominal, 1);

vector<TH1F*> hists_muon_up = get_all_hists(flavor_up_muon, hist);
vector<TH1F*> hists_elec_up = get_all_hists(flavor_up_elec, hist);
vector<TH1F*> hists_up = AddHists(hists_muon_up, hists_elec_up, 1);

vector<TH1F*> hists_muon_down = get_all_hists(flavor_down_muon, hist);
vector<TH1F*> hists_elec_down = get_all_hists(flavor_down_elec, hist);
vector<TH1F*> hists_down = AddHists(hists_muon_down, hists_elec_down, 1);

gErrorIgnoreLevel = kError;
SetupGlobalStyle();
gStyle->SetErrorX(0);

// ===========================================================================
// === Plot WMass plots
cout << "Start with mjet plots ..." << endl;

vector<TH1F*> ratios_nominal,ratios_up, ratios_down;
for(unsigned int i=0; i<4; i++){
  ratios_nominal.push_back(GetRatioLocal(hists_nominal[i], hists_nominal[i],true));
  ratios_up.push_back(GetRatioLocal(hists_up[i], hists_nominal[i],false));
  ratios_down.push_back(GetRatioLocal(hists_down[i], hists_nominal[i],false));

  SetupCanvas(true);
  Cosmetics(hists_nominal[i], "#it{m}_{had}", kGray+2, kSolid,  0, 2, true);
  Cosmetics(hists_up[i],      "#it{m}_{had}", kRed+2,  kSolid,  0, 2, true);
  Cosmetics(hists_down[i],    "#it{m}_{had}", kRed+2,  kDashed, 0, 2, true);

  RatioCosmetics(ratios_nominal[i], "#it{m}_{jet}", kGray+2, kSolid,  2, 0);
  RatioCosmetics(ratios_up[i],      "#it{m}_{jet}", kRed+2,  kSolid,  2, 0);
  RatioCosmetics(ratios_down[i],    "#it{m}_{jet}", kRed+2,  kDashed, 2, 0);

  m_rp1_top->cd();
  hists_nominal[i]->GetYaxis()->SetRangeUser(0., hists_nominal[i]->GetMaximum()*1.4);
  hists_nominal[i]->Draw("hist ][");
  hists_up[i]->Draw("hist same ][");
  hists_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();
  DrawLegend({hists_nominal[i], hists_up[i], hists_down[i]}, name+"_"+years[i]);
  DrawLumi(lumis[i], true, true, true);

  m_rp1->cd();
  ratios_nominal[i]->Draw("hist ][");
  ratios_up[i]->Draw("hist same ][");
  ratios_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();

  m_can->Print(save_afs+name+"_"+years[i]+".pdf");
  m_can->Print(save_nfs+name+"_"+years[i]+".pdf");
}

vector<TH1F*> ratios_muon_nominal,ratios_muon_up, ratios_muon_down;
for(unsigned int i=0; i<3; i++){
  ratios_muon_nominal.push_back(GetRatioLocal(hists_muon_nominal[i], hists_muon_nominal[i],true));
  ratios_muon_up.push_back(GetRatioLocal(hists_muon_up[i], hists_muon_nominal[i],false));
  ratios_muon_down.push_back(GetRatioLocal(hists_muon_down[i], hists_muon_nominal[i],false));

  SetupCanvas(true);
  Cosmetics(hists_muon_nominal[i], "#it{m}_{had}", kGray+2, kSolid,  0, 2, true);
  Cosmetics(hists_muon_up[i],      "#it{m}_{had}", kRed+2,  kSolid,  0, 2, true);
  Cosmetics(hists_muon_down[i],    "#it{m}_{had}", kRed+2,  kDashed, 0, 2, true);

  RatioCosmetics(ratios_muon_nominal[i], "#it{m}_{jet}", kGray+2, kSolid,  2, 0);
  RatioCosmetics(ratios_muon_up[i],      "#it{m}_{jet}", kRed+2,  kSolid,  2, 0);
  RatioCosmetics(ratios_muon_down[i],    "#it{m}_{jet}", kRed+2,  kDashed, 2, 0);

  m_rp1_top->cd();
  hists_muon_nominal[i]->GetYaxis()->SetRangeUser(0., hists_muon_nominal[i]->GetMaximum()*1.4);
  hists_muon_nominal[i]->Draw("hist ][");
  hists_muon_up[i]->Draw("hist same ][");
  hists_muon_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();
  DrawLegend({hists_muon_nominal[i], hists_muon_up[i], hists_muon_down[i]}, name+"_"+years[i]);
  DrawLumi(lumis[i], true, true, true);

  m_rp1->cd();
  ratios_muon_nominal[i]->Draw("hist ][");
  ratios_muon_up[i]->Draw("hist same ][");
  ratios_muon_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();

  m_can->Print(save_afs+name+"_"+years[i]+"_muon.pdf");
  m_can->Print(save_nfs+name+"_"+years[i]+"_muon.pdf");
}

vector<TH1F*> ratios_elec_nominal,ratios_elec_up, ratios_elec_down;
for(unsigned int i=0; i<3; i++){
  ratios_elec_nominal.push_back(GetRatioLocal(hists_elec_nominal[i], hists_elec_nominal[i],true));
  ratios_elec_up.push_back(GetRatioLocal(hists_elec_up[i], hists_elec_nominal[i],false));
  ratios_elec_down.push_back(GetRatioLocal(hists_elec_down[i], hists_elec_nominal[i],false));

  SetupCanvas(true);
  Cosmetics(hists_elec_nominal[i], "#it{m}_{had}", kGray+2, kSolid,  0, 2, true);
  Cosmetics(hists_elec_up[i],      "#it{m}_{had}", kRed+2,  kSolid,  0, 2, true);
  Cosmetics(hists_elec_down[i],    "#it{m}_{had}", kRed+2,  kDashed, 0, 2, true);

  RatioCosmetics(ratios_elec_nominal[i], "#it{m}_{jet}", kGray+2, kSolid,  2, 0);
  RatioCosmetics(ratios_elec_up[i],      "#it{m}_{jet}", kRed+2,  kSolid,  2, 0);
  RatioCosmetics(ratios_elec_down[i],    "#it{m}_{jet}", kRed+2,  kDashed, 2, 0);

  m_rp1_top->cd();
  hists_elec_nominal[i]->GetYaxis()->SetRangeUser(0., hists_elec_nominal[i]->GetMaximum()*1.4);
  hists_elec_nominal[i]->Draw("hist ][");
  hists_elec_up[i]->Draw("hist same ][");
  hists_elec_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();
  DrawLegend({hists_elec_nominal[i], hists_elec_up[i], hists_elec_down[i]}, name+"_"+years[i]);
  DrawLumi(lumis[i], true, true, true);

  m_rp1->cd();
  ratios_elec_nominal[i]->Draw("hist ][");
  ratios_elec_up[i]->Draw("hist same ][");
  ratios_elec_down[i]->Draw("hist same ][");
  gPad->RedrawAxis();

  m_can->Print(save_afs+name+"_"+years[i]+"_elec.pdf");
  m_can->Print(save_nfs+name+"_"+years[i]+"_elec.pdf");
}

*/
