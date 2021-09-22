#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/tdrstyle_all.h"

using namespace std;

TString indir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/muon/";
TString ttbar = "uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root";
TString outdir = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/Prefire/";
TString dir_up = "Prefire_up/";
TString dir_down = "Prefire_down/";
TString dir_off = "PrefireOFF/";

TString year;

int main(int argc, char* argv[]){

  year = argv[1];

  cout << "Get Files ..." << endl;
  TFile* f_nom = new TFile(indir+ttbar);
  TFile* f_up = new TFile(indir+dir_up+ttbar);
  TFile* f_down = new TFile(indir+dir_down+ttbar);
  TFile* f_off = new TFile(indir+dir_off+ttbar);

  cout << "Get Hists (HT) ..." << endl;
  TString hist = "XCone_cor/M_jet1_";
  TH1F* h_nom = (TH1F*) f_nom->Get(hist);
  TH1F* h_up = (TH1F*) f_up->Get(hist);
  TH1F* h_down = (TH1F*) f_down->Get(hist);
  TH1F* h_off = (TH1F*) f_off->Get(hist);

  cout << "\t ... normalize hists" << endl;
  TH1F* h_nom_norm = Normalize(h_nom);
  TH1F* h_up_norm = Normalize(h_up);
  TH1F* h_down_norm = Normalize(h_down);
  TH1F* h_off_norm = Normalize(h_off);

  // TH1F* GetRatio(TH1F* h1, TH1F* h2, bool equal, bool isEffi){
  TH1F* ratio_nom = GetRatio(h_nom, h_nom, true);
  TH1F* ratio_up = GetRatio(h_up, h_nom, false);
  TH1F* ratio_down = GetRatio(h_down, h_nom, false);
  TH1F* ratio_off = GetRatio(h_off, h_nom, false);

  TH1F* ratio_nom_norm = GetRatio(h_nom_norm, h_nom_norm, true);
  TH1F* ratio_up_norm = GetRatio(h_up_norm, h_nom_norm, false);
  TH1F* ratio_down_norm = GetRatio(h_down_norm, h_nom_norm, false);
  TH1F* ratio_off_norm = GetRatio(h_off_norm, h_nom_norm, false);

  cout << "Draw ..." << endl;
  cout << "\t ... default" << endl;
  TLegend* leg = tdrLeg(0.5, 0.7, 0.9, 0.9, 0.04);
  TCanvas* canv = tdrDiCanvas2("prefire", 0, 500, 0, h_nom->GetMaximum()*1.1, 0.4, 1.6, "#it{m}_{Jet}", "Events", "#frac{variation}{nominal}");
  canv->SetTickx(0);
  canv->SetTicky(0);
  canv->cd(1);
  tdrDraw(h_nom, "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(h_up, "hist same", 1, kBlack, kSolid, kGreen+2, 0);
  tdrDraw(h_down, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  tdrDraw(h_off, "hist same", 1, kBlack, kSolid, kBlue+2, 0);
  leg->AddEntry(h_nom, "nominal", "l");
  leg->AddEntry(h_up, "up", "l");
  leg->AddEntry(h_down, "down", "l");
  leg->AddEntry(h_off, "off", "l");
  leg->Draw("same");
  canv->cd(2);
  tdrDraw(ratio_nom, "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(ratio_up, "hist same", 1, kBlack, kSolid, kGreen+2, 0);
  tdrDraw(ratio_down, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  tdrDraw(ratio_off, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  canv->Update();
  canv->Print(outdir+"Prefire_SYS_"+year+".pdf","pdf");
  leg->Clear();
  delete canv;

  canv = tdrDiCanvas2("prefire_norm", 0, 500, 0, h_nom_norm->GetMaximum()*1.1, 0.4, 1.6, "#it{m}_{Jet}", "a.u.", "#frac{variation}{nominal}");
  canv->SetTickx(0);
  canv->SetTicky(0);
  canv->cd(1);
  tdrDraw(h_nom_norm, "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(h_up_norm, "hist same", 1, kBlack, kSolid, kGreen+2, 0);
  tdrDraw(h_down_norm, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  tdrDraw(h_off_norm, "hist same", 1, kBlack, kSolid, kBlue+2, 0);
  leg->AddEntry(h_nom_norm, "nominal", "l");
  leg->AddEntry(h_up_norm, "up", "l");
  leg->AddEntry(h_down_norm, "down", "l");
  leg->AddEntry(h_off_norm, "off", "l");
  leg->Draw("same");
  canv->cd(2);
  tdrDraw(ratio_nom_norm, "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(ratio_up_norm, "hist same", 1, kBlack, kSolid, kGreen+2, 0);
  tdrDraw(ratio_down_norm, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  tdrDraw(ratio_off_norm, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  canv->Update();
  canv->Print(outdir+"Prefire_SYS_norm_"+year+".pdf","pdf");
  leg->Clear();
  delete canv;

  // cout << "\t ... normalized" << endl;
  // canv = tdrCanvas("prefire_norm", 0, 500, 0, h_nom_norm->GetMaximum()*1.1, "HT", "a.u.");
  // canv->SetTickx(0);
  // canv->SetTicky(0);
  // tdrDraw(h_nom_norm, "hist same", 1, kBlack, kSolid, kRed+2, 0);
  // tdrDraw(h_up_norm, "hist same", 1, kBlack, kSolid, kGreen+2, 0);
  // tdrDraw(h_down_norm, "hist same", 1, kBlack, kDashed, kGreen+2, 0);
  // leg->AddEntry(h_nom, "nominal", "l");
  // leg->AddEntry(h_up, "up", "l");
  // leg->AddEntry(h_down, "down", "l");
  // leg->Draw("same");
  // canv->Update();
  // canv->Print(outdir+"Prefire_SYS_norm_"+year+".pdf","pdf");
  // leg->Clear();

  // TCanvas* tdrDiCanvas2(const char* canvName, double x_min, double x_max, double y_min, double y_max, double y_min2, double y_max2,  const char* nameXaxis, const char* nameYaxis, const char* nameYaxis2, bool square = kRectangular, int iPeriod = 4, int iPos = 11) {


}
