#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"
#include "../include/CreatHists.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

using namespace std;

void PlotTau(TH1F* data, TH1F* ttbar, TH1F* FSRup, TH1F* FSRdown, TString tau, TString year, TString var);
void SetSettings(TH1F* hist, int color, int width, int style, TString title, TString x, TString y);
void DrawHist(TCanvas* c, TH1F* hist, TString option, TLegend* leg, TString ltext, TString loption, int color);

bool debug = true;
bool forTalk = true;
TString save_path;

int main(int argc, char* argv[]){

  if(argc != 2){
    cout << "Usage: ./FSRtau43 <year>" << endl;
    cout << "  year = 2016/2017/2018/combine" << endl;
    return 0;
  }

  SetupGlobalStyle();
  gErrorIgnoreLevel = kWarning;

  TString year = argv[1];
  save_path = get_save_path()+"/Plots/FSRuncertainty";
  save_path = creat_folder_and_path(save_path, "Tau43");

  TString tau43 = "comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau43";
  TString tau32 = "comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau32";

  //////////////////////////////////////////////////////////////////////////////
  ///               Get Histograms                                           ///
  //////////////////////////////////////////////////////////////////////////////

  cout << "Get Hists ..." << endl;
  cout << "\t ... tau43" << endl;
  if(debug) cout << "\t\t ... data" << endl;
  VecH h_43_data_muon          = get_all_hists(data_muon, tau43);
  VecH h_43_data_elec          = get_all_hists(data_elec, tau43);
  VecH h_43_data_combine       = combine_channels(h_43_data_muon, h_43_data_elec);

  if(debug) cout << "\t\t ... ttbar" << endl;
  VecH h_43_ttbar_muon         = get_all_hists(ttbar_muon, tau43);
  VecH h_43_ttbar_elec         = get_all_hists(ttbar_elec, tau43);
  VecH h_43_ttbar_combine      = combine_channels(h_43_ttbar_muon, h_43_ttbar_elec);

  if(debug) cout << "\t\t ... single top" << endl;
  VecH h_43_st_muon            = get_all_hists(st_muon, tau43);
  VecH h_43_st_elec            = get_all_hists(st_elec, tau43);
  VecH h_43_st_combine         = combine_channels(h_43_st_muon, h_43_st_elec);

  if(debug) cout << "\t\t ... wjets" << endl;
  VecH h_43_wjets_muon         = get_all_hists(wjets_muon, tau43);
  VecH h_43_wjets_elec         = get_all_hists(wjets_elec, tau43);
  VecH h_43_wjets_combine      = combine_channels(h_43_wjets_muon, h_43_wjets_elec);

  if(debug) cout << "\t\t ... other" << endl;
  VecH h_43_other_muon         = get_all_hists(other_muon, tau43);
  VecH h_43_other_elec         = get_all_hists(other_elec, tau43);
  VecH h_43_other_combine      = combine_channels(h_43_other_muon, h_43_other_elec);

  if(debug) cout << "\t\t ... add background" << endl;
  VecH h_43_background = AddHists(h_43_st_combine, h_43_wjets_combine, h_43_other_combine, 1);
  VecH h_43_ttbar_= AddHists(h_43_ttbar_combine, h_43_background, 1);

  if(debug) cout << "\t\t ... up 2016" << endl;
  TH1F* h_43_fsr_up_16 = get_hist(fsr_up_muon_16, tau43);
  TH1F* h_43_fsr_down_16 = get_hist(fsr_down_muon_16, tau43);

  if(debug) cout << "\t\t ... up sqrt2" << endl;
  VecH h_43_fsr_up_sqrt2_muon         = get_all_hists(fsr_up_sqrt2_muon, tau43);
  VecH h_43_fsr_up_sqrt2_elec         = get_all_hists(fsr_up_sqrt2_elec, tau43);
  VecH h_43_fsr_up_sqrt2_combine      = combine_channels(h_43_fsr_up_sqrt2_muon, h_43_fsr_up_sqrt2_elec);

  if(debug) cout << "\t\t ... up 2" << endl;
  VecH h_43_fsr_up_2_muon         = get_all_hists(fsr_up_2_muon, tau43);
  VecH h_43_fsr_up_2_elec         = get_all_hists(fsr_up_2_elec, tau43);
  VecH h_43_fsr_up_2_combine      = combine_channels(h_43_fsr_up_2_muon, h_43_fsr_up_2_elec);

  if(debug) cout << "\t\t ... up 4" << endl;
  VecH h_43_fsr_up_4_muon         = get_all_hists(fsr_up_4_muon, tau43);
  VecH h_43_fsr_up_4_elec         = get_all_hists(fsr_up_4_elec, tau43);
  VecH h_43_fsr_up_4_combine      = combine_channels(h_43_fsr_up_4_muon, h_43_fsr_up_4_elec);

  if(debug) cout << "\t\t ... down sqrt2" << endl;
  VecH h_43_fsr_down_sqrt2_muon         = get_all_hists(fsr_down_sqrt2_muon, tau43);
  VecH h_43_fsr_down_sqrt2_elec         = get_all_hists(fsr_down_sqrt2_elec, tau43);
  VecH h_43_fsr_down_sqrt2_combine      = combine_channels(h_43_fsr_down_sqrt2_muon, h_43_fsr_down_sqrt2_elec);

  if(debug) cout << "\t\t ... down 2" << endl;
  VecH h_43_fsr_down_2_muon         = get_all_hists(fsr_down_2_muon, tau43);
  VecH h_43_fsr_down_2_elec         = get_all_hists(fsr_down_2_elec, tau43);
  VecH h_43_fsr_down_2_combine      = combine_channels(h_43_fsr_down_2_muon, h_43_fsr_down_2_elec);

  if(debug) cout << "\t\t ... fsr down 4" << endl;
  VecH h_43_fsr_down_4_muon         = get_all_hists(fsr_down_4_muon, tau43);
  VecH h_43_fsr_down_4_elec         = get_all_hists(fsr_down_4_elec, tau43);
  VecH h_43_fsr_down_4_combine      = combine_channels(h_43_fsr_down_4_muon, h_43_fsr_down_4_elec);

  if(debug) cout << "\t\t ... Rebin ";
  int rebin = 2;
  for(auto hist: h_43_data_combine) hist->Rebin(rebin);
  for(auto hist: h_43_ttbar_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_up_sqrt2_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_up_2_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_up_4_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_down_sqrt2_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_down_2_combine) hist->Rebin(rebin);
  for(auto hist: h_43_fsr_down_4_combine) hist->Rebin(rebin);
  (TH1F*) h_43_fsr_up_16->Rebin(rebin);
  (TH1F*) h_43_fsr_down_16->Rebin(rebin);
  cout << " with factor " << rebin << " - New bin width " <<  h_43_data_combine[0]->GetBinWidth(1) << endl;

  if(debug) cout << "\t\t ... Normalize" << endl;
  TH1F* h_43_data_norm = Normalize(h_43_data_combine[h_43_data_combine.size()-1]); // last entry at v[v.size()-1]
  TH1F* h_43_ttbar_norm = Normalize(h_43_ttbar_combine[h_43_ttbar_combine.size()-1]);
  TH1F* h_43_fsr_up_sqrt2_norm = Normalize(h_43_fsr_up_sqrt2_combine[h_43_fsr_up_sqrt2_combine.size()-1]);
  TH1F* h_43_fsr_up_2_norm = Normalize(h_43_fsr_up_2_combine[h_43_fsr_up_2_combine.size()-1]);
  TH1F* h_43_fsr_up_4_norm = Normalize(h_43_fsr_up_4_combine[h_43_fsr_up_4_combine.size()-1]);
  TH1F* h_43_fsr_down_sqrt2_norm = Normalize(h_43_fsr_down_sqrt2_combine[h_43_fsr_down_sqrt2_combine.size()-1]);
  TH1F* h_43_fsr_down_2_norm = Normalize(h_43_fsr_down_2_combine[h_43_fsr_down_2_combine.size()-1]);
  TH1F* h_43_fsr_down_4_norm = Normalize(h_43_fsr_down_4_combine[h_43_fsr_down_4_combine.size()-1]);
  TH1F* h_43_fsr_up_16_norm = Normalize(h_43_fsr_up_16);
  TH1F* h_43_fsr_down_16_norm = Normalize(h_43_fsr_down_16);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "\t ... tau32" << endl;

  if(debug) cout << "\t\t ... data" << endl;
  VecH h_32_data_muon          = get_all_hists(data_muon, tau32);
  VecH h_32_data_elec          = get_all_hists(data_elec, tau32);
  VecH h_32_data_combine       = combine_channels(h_32_data_muon, h_32_data_elec);

  if(debug) cout << "\t\t ... ttbar" << endl;
  VecH h_32_ttbar_muon         = get_all_hists(ttbar_muon, tau32);
  VecH h_32_ttbar_elec         = get_all_hists(ttbar_elec, tau32);
  VecH h_32_ttbar_combine      = combine_channels(h_32_ttbar_muon, h_32_ttbar_elec);

  if(debug) cout << "\t\t ... single top" << endl;
  VecH h_32_st_muon            = get_all_hists(st_muon, tau32);
  VecH h_32_st_elec            = get_all_hists(st_elec, tau32);
  VecH h_32_st_combine         = combine_channels(h_32_st_muon, h_32_st_elec);

  if(debug) cout << "\t\t ... wjets" << endl;
  VecH h_32_wjets_muon         = get_all_hists(wjets_muon, tau32);
  VecH h_32_wjets_elec         = get_all_hists(wjets_elec, tau32);
  VecH h_32_wjets_combine      = combine_channels(h_32_wjets_muon, h_32_wjets_elec);

  if(debug) cout << "\t\t ... other" << endl;
  VecH h_32_other_muon         = get_all_hists(other_muon, tau32);
  VecH h_32_other_elec         = get_all_hists(other_elec, tau32);
  VecH h_32_other_combine      = combine_channels(h_32_other_muon, h_32_other_elec);

  if(debug) cout << "\t\t ... add background" << endl;
  VecH h_32_background = AddHists(h_32_st_combine, h_32_wjets_combine, h_32_other_combine, 1);
  VecH h_32_ttbar_= AddHists(h_32_ttbar_combine, h_32_background, 1);

  if(debug) cout << "\t\t ... up 2016" << endl;
  TH1F* h_32_fsr_up_16 = get_hist(fsr_up_muon_16, tau32);
  TH1F* h_32_fsr_down_16 = get_hist(fsr_down_muon_16, tau32);

  if(debug) cout << "\t\t ... up sqrt2" << endl;
  VecH h_32_fsr_up_sqrt2_muon         = get_all_hists(fsr_up_sqrt2_muon, tau32);
  VecH h_32_fsr_up_sqrt2_elec         = get_all_hists(fsr_up_sqrt2_elec, tau32);
  VecH h_32_fsr_up_sqrt2_combine      = combine_channels(h_32_fsr_up_sqrt2_muon, h_32_fsr_up_sqrt2_elec);

  if(debug) cout << "\t\t ... up 2" << endl;
  VecH h_32_fsr_up_2_muon         = get_all_hists(fsr_up_2_muon, tau32);
  VecH h_32_fsr_up_2_elec         = get_all_hists(fsr_up_2_elec, tau32);
  VecH h_32_fsr_up_2_combine      = combine_channels(h_32_fsr_up_2_muon, h_32_fsr_up_2_elec);

  if(debug) cout << "\t\t ... up 4" << endl;
  VecH h_32_fsr_up_4_muon         = get_all_hists(fsr_up_4_muon, tau32);
  VecH h_32_fsr_up_4_elec         = get_all_hists(fsr_up_4_elec, tau32);
  VecH h_32_fsr_up_4_combine      = combine_channels(h_32_fsr_up_4_muon, h_32_fsr_up_4_elec);

  if(debug) cout << "\t\t ... down sqrt2" << endl;
  VecH h_32_fsr_down_sqrt2_muon         = get_all_hists(fsr_down_sqrt2_muon, tau32);
  VecH h_32_fsr_down_sqrt2_elec         = get_all_hists(fsr_down_sqrt2_elec, tau32);
  VecH h_32_fsr_down_sqrt2_combine      = combine_channels(h_32_fsr_down_sqrt2_muon, h_32_fsr_down_sqrt2_elec);

  if(debug) cout << "\t\t ... down 2" << endl;
  VecH h_32_fsr_down_2_muon         = get_all_hists(fsr_down_2_muon, tau32);
  VecH h_32_fsr_down_2_elec         = get_all_hists(fsr_down_2_elec, tau32);
  VecH h_32_fsr_down_2_combine      = combine_channels(h_32_fsr_down_2_muon, h_32_fsr_down_2_elec);

  if(debug) cout << "\t\t ... fsr down 4" << endl;
  VecH h_32_fsr_down_4_muon         = get_all_hists(fsr_down_4_muon, tau32);
  VecH h_32_fsr_down_4_elec         = get_all_hists(fsr_down_4_elec, tau32);
  VecH h_32_fsr_down_4_combine      = combine_channels(h_32_fsr_down_4_muon, h_32_fsr_down_4_elec);

  if(debug) cout << "\t\t ... Rebin ";
  for(auto hist: h_32_data_combine) hist->Rebin(rebin);
  for(auto hist: h_32_ttbar_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_up_sqrt2_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_up_2_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_up_4_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_down_sqrt2_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_down_2_combine) hist->Rebin(rebin);
  for(auto hist: h_32_fsr_down_4_combine) hist->Rebin(rebin);
  (TH1F*) h_32_fsr_up_16->Rebin(rebin);
  (TH1F*) h_32_fsr_down_16->Rebin(rebin);
  cout << " with factor " << rebin << " - New bin width " <<  h_32_data_combine[0]->GetBinWidth(1) << endl;

  if(debug) cout << "\t\t ... Normalize" << endl;
  TH1F* h_32_data_norm = Normalize(h_32_data_combine[h_32_data_combine.size()-1]); // last entry at v[v.size()-1]
  TH1F* h_32_ttbar_norm = Normalize(h_32_ttbar_combine[h_32_ttbar_combine.size()-1]);
  TH1F* h_32_fsr_up_sqrt2_norm = Normalize(h_32_fsr_up_sqrt2_combine[h_32_fsr_up_sqrt2_combine.size()-1]);
  TH1F* h_32_fsr_up_2_norm = Normalize(h_32_fsr_up_2_combine[h_32_fsr_up_2_combine.size()-1]);
  TH1F* h_32_fsr_up_4_norm = Normalize(h_32_fsr_up_4_combine[h_32_fsr_up_4_combine.size()-1]);
  TH1F* h_32_fsr_down_sqrt2_norm = Normalize(h_32_fsr_down_sqrt2_combine[h_32_fsr_down_sqrt2_combine.size()-1]);
  TH1F* h_32_fsr_down_2_norm = Normalize(h_32_fsr_down_2_combine[h_32_fsr_down_2_combine.size()-1]);
  TH1F* h_32_fsr_down_4_norm = Normalize(h_32_fsr_down_4_combine[h_32_fsr_down_4_combine.size()-1]);
  TH1F* h_32_fsr_up_16_norm = Normalize(h_32_fsr_up_16);
  TH1F* h_32_fsr_down_16_norm = Normalize(h_32_fsr_down_16);

  //////////////////////////////////////////////////////////////////////////////
  ///               I Plot all Histograms                                    ///
  //////////////////////////////////////////////////////////////////////////////

  // void SetSettings(TH1F* hist, int color, int width, int style, TString title, TString x, TString y){
  cout << "Set Hist Settings ..." << endl;
  if(debug) cout << "\t ... tau43" << endl;
  SetSettings(h_43_data_norm, 1, 1, 8, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_ttbar_norm, 14, 2, 1, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_up_sqrt2_norm, kAzure+7, 2, 1, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_up_2_norm, kGreen-4, 2, 1, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_up_4_norm, 810, 2, 1, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_down_sqrt2_norm, kAzure+7, 2, 2, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_down_2_norm, kGreen-4, 2, 2, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_down_4_norm, 810, 2, 2, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_up_16_norm, kAzure+7, 2, 1, " ", "#tau_{43}", "a.u.");
  SetSettings(h_43_fsr_down_16_norm, kAzure+7, 2, 2, " ", "#tau_{43}", "a.u.");

  if(debug) cout << "\t ... tau32" << endl;
  SetSettings(h_32_data_norm, 1, 1, 8, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_ttbar_norm, 14, 2, 1, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_up_sqrt2_norm, kAzure+7, 2, 1, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_up_2_norm, kGreen-4, 2, 1, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_up_4_norm, 810, 2, 1, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_down_sqrt2_norm, kAzure+7, 2, 2, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_down_2_norm, kGreen-4, 2, 2, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_down_4_norm, 810, 2, 2, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_up_16_norm, kAzure+7, 2, 1, " ", "#tau_{32}", "a.u.");
  SetSettings(h_32_fsr_down_16_norm, kAzure+7, 2, 2, " ", "#tau_{32}", "a.u.");

  // void PlotTau(TH1F* data, TH1F* ttbar, TH1F* FSRup, TH1F* FSRdown, TString tau, TString year, TString var="2"){
  cout << "Start plotting ..." << endl;
  if(debug) cout << "\t ... tau32" << endl;
  PlotTau(h_32_data_norm, h_32_ttbar_norm, h_32_fsr_up_4_norm, h_32_fsr_down_4_norm, "32", "combine", "4");
  PlotTau(h_32_data_norm, h_32_ttbar_norm, h_32_fsr_up_16_norm, h_32_fsr_down_16_norm, "32", "2016", "2");

  if(debug) cout << "\t ... tau43" << endl;
  PlotTau(h_43_data_norm, h_43_ttbar_norm, h_43_fsr_up_4_norm, h_43_fsr_down_4_norm, "43", "combine", "4");
  PlotTau(h_43_data_norm, h_43_ttbar_norm, h_43_fsr_up_16_norm, h_43_fsr_down_16_norm, "43", "2016", "2");

  // void DrawHist(TCanvas* c, TH1F* hist, TString option, TLegend* leg, TString ltext, TString loption, int color)
  if(debug) cout << "\t ... both" << endl;

  if(debug) cout << "\t\t ... 2016" << endl;
  TCanvas *c = new TCanvas("both2016", "both2016", 600, 600);
  gPad->SetLeftMargin(.2);
  gPad->RedrawAxis();

  TLegend * leg = new TLegend(.25, .55, .50, .85);
  leg->SetTextSize(0.03);

  DrawHist(c, h_43_ttbar_norm, "HIST", leg, "t#bar{t} nominal (#tau_{43})", "l", kGray+3);
  DrawHist(c, h_43_fsr_up_16_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = #frac{1}{2} (#tau_{43})", "l", kGreen+2);
  DrawHist(c, h_43_fsr_down_16_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = 2 (#tau_{43})", "l", kGreen+2);

  DrawHist(c, h_32_ttbar_norm, "HIST SAME", leg, "t#bar{t} nominal (#tau_{32})", "l", 14);
  DrawHist(c, h_32_fsr_up_16_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = #frac{1}{2} (#tau_{32})", "l", kAzure+7);
  DrawHist(c, h_32_fsr_down_16_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = 2 (#tau_{32})", "l", kAzure+7);

  leg->Draw();
  if(forTalk) CMSLabel(true, 0.5, 0.84);
  c->SaveAs(save_path+"/TauBoth_2016.pdf");
  leg->Clear();
  delete c;

  if(debug) cout << "\t\t ... combine" << endl;
  c = new TCanvas("both2016", "both2016", 600, 600);
  gPad->SetLeftMargin(.2);
  gPad->RedrawAxis();

  leg = new TLegend(.25, .55, .50, .85);
  leg->SetTextSize(0.03);
  DrawHist(c, h_43_ttbar_norm, "HIST", leg, "t#bar{t} nominal (#tau_{43})", "l", kGray+3);
  DrawHist(c, h_43_fsr_up_4_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = #frac{1}{4} (#tau_{43})", "l", kGreen+2);
  DrawHist(c, h_43_fsr_down_4_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = 4 (#tau_{43})", "l", kGreen+2);

  DrawHist(c, h_32_ttbar_norm, "HIST SAME", leg, "t#bar{t} nominal (#tau_{32})", "l", 14);
  DrawHist(c, h_32_fsr_up_4_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = #frac{1}{4} (#tau_{32})", "l", 810);
  DrawHist(c, h_32_fsr_down_4_norm, "HIST SAME", leg, "t#bar{t} f^{FSR} = 4 (#tau_{32})", "l", 810);

  leg->Draw();
  if(forTalk) CMSLabel(true, 0.5, 0.84);
  c->SaveAs(save_path+"/TauBoth_combine.pdf");
  leg->Clear();
  delete c;


}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetSettings(TH1F* hist, int color, int width, int style, TString title, TString x, TString y){
  hist->SetMarkerStyle(style);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->SetLineWidth(width);
  hist->SetLineStyle(style);
  hist->SetTitle(title);
  hist->GetXaxis()->SetTitle(x);
  hist->GetYaxis()->SetTitle(y);
  hist->GetYaxis()->SetRangeUser(0., 1.4*hist->GetMaximum());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotTau(TH1F* data, TH1F* ttbar, TH1F* FSRup, TH1F* FSRdown, TString tau, TString year, TString var="2"){
  TCanvas *c = new TCanvas(tau+year, tau+year, 600, 600);
  gPad->SetLeftMargin(.2);

  ttbar->Draw("HIST");
  data->Draw("EP SAME");
  FSRup->Draw("HIST SAME");
  FSRdown->Draw("HIST SAME");

  TLegend * leg = new TLegend(.25, .65, .50, .85);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "Data", "pe");
  leg->AddEntry(ttbar, "t#bar{t} nominal", "l");
  leg->AddEntry(FSRdown, "t#bar{t} f^{FSR} = #frac{1}{"+var+"}", "l");
  leg->AddEntry(FSRup, "t#bar{t} f^{FSR} = "+var, "l");
  leg->Draw();

  gPad->RedrawAxis();
  if(forTalk) CMSLabel(true, 0.5, 0.84);

  c->SaveAs(save_path+"/Tau"+tau+"_"+year+".pdf");
  delete c;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DrawHist(TCanvas* c, TH1F* hist, TString option, TLegend* leg, TString ltext, TString loption, int color){
  hist->SetLineColor(color);
  hist->Draw(option);
  leg->AddEntry(hist, ltext, loption);
}
