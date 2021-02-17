#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

int main(int argc, char* argv[]){
  bool debug = true;
  TString path_save = get_save_path()+"/Plots/JEC_SYS";
  TString path_plot = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/";
  vector<TString> years = {"2016", "2017", "2018"};
  vector<TH1F*> hists;
  vector<TH2F*> hists2d;

  TString hist_pt    = "comparison_topjet_xcone_pass_rec/wjet_pt_match_S1divW";
  TString hist_ratio = "comparison_topjet_xcone_pass_rec/wjet_pt_match_S1divW_W";

  if(debug) cout << "Get File\n";
  TFile *file_tt = new TFile(path_plot+"/combined/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  if(debug) cout << "2D Hist\n";
  TH2F* hist2D_tt = (TH2F*) file_tt->Get(hist_ratio);
  hist2D_tt->SetTitle("");
  hist2D_tt->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");
  hist2D_tt->GetYaxis()->SetTitle("p_{T}^{1st subjet}/p_{T}^{W}");
  hist2D_tt->GetXaxis()->SetTitleSize(0.05);
  hist2D_tt->GetYaxis()->SetTitleSize(0.055);
  hist2D_tt->GetXaxis()->SetTitleOffset(0.9);
  hist2D_tt->GetYaxis()->SetTitleOffset(0.88);
  // hist2D_tt->GetXaxis()->SetRangeUser();
  hist2D_tt->GetYaxis()->SetRangeUser(0.5, 1);

  TLine *xaxis = new TLine(0, 0.7, 1000, 0.7);
  TLine *yaxis = new TLine(300, 0.5, 300, 1);
  xaxis->SetLineWidth(1);
  yaxis->SetLineWidth(1);
  xaxis->SetLineColor(kRed);
  yaxis->SetLineColor(kRed);

  TCanvas *A = new TCanvas("A", "A", 700, 500);
  hist2D_tt->Draw("COLZ");
  hist2D_tt->SetStats(kFALSE); // Eliminate Stat-Box
  xaxis->Draw("same");
  yaxis->Draw("same");
  A->SaveAs(path_save+"/S1divW_tt.pdf");
  delete A;

  TFile *file_other = new TFile(path_plot+"/combined/muon/uhh2.AnalysisModuleRunner.MC.other.root");
  TFile *file_st    = new TFile(path_plot+"/combined/muon/uhh2.AnalysisModuleRunner.MC.SingleTop.root");
  TFile *file_wjet  = new TFile(path_plot+"/combined/muon/uhh2.AnalysisModuleRunner.MC.WJets.root");

  TH2F* hist2D_bkg = (TH2F*) file_other->Get(hist_ratio);
  hist2D_bkg->Add((TH2F*) file_st->Get(hist_ratio) ,1);
  hist2D_bkg->Add((TH2F*) file_wjet->Get(hist_ratio) ,1);

  A = new TCanvas("A", "A", 700, 500);
  // A->SetRightMargin(0.15);
  A->SetLeftMargin(0.12);

  // A->SetBottomMargin(0.13);
  // A->SetTopMargin(0.13);

  hist2D_tt->Add(hist2D_bkg, 1);
  hist2D_tt->Draw("COLZ");
  xaxis->Draw("same");
  yaxis->Draw("same");
  gPad->Modified();
  gPad->Update();
  A->SaveAs(path_save+"/S1divW.pdf");
  delete A;

  /*
  ██ ███    ██ ████████ ███████  ██████  ██████   █████  ██
  ██ ████   ██    ██    ██      ██       ██   ██ ██   ██ ██
  ██ ██ ██  ██    ██    █████   ██   ███ ██████  ███████ ██
  ██ ██  ██ ██    ██    ██      ██    ██ ██   ██ ██   ██ ██
  ██ ██   ████    ██    ███████  ██████  ██   ██ ██   ██ ███████
  */

  if(debug) cout << "Get Integral\n";
  int ny=70; // up to 0.7 - Bin 70
  int nx=30;
  cout << hist2D_tt->Integral() << endl; // Only in given Range
  double int_hh = hist2D_tt->Integral(nx+1,   101, ny+1,  101);
  double int_hl = hist2D_tt->Integral(nx+1,   101,    1,   ny);
  double int_lh = hist2D_tt->Integral(   1,    nx, ny+1,  101);
  double int_ll = hist2D_tt->Integral(   1,    nx,    1,   ny);
  cout << "High High: " << int_hh << endl;
  cout << "High Low : " << int_hl << endl;
  cout << "Low  High: " << int_lh << endl;
  cout << "Low  Low : " << int_ll << endl;

  if(debug) cout << "Get Integral - 2\n";
  double integral1 = hist2D_tt->Integral() ;
  double integral2 = hist2D_tt->Integral(0, 101, 0, 101); // All bins outside the range are summarized in one extra bin
  double integral3 = hist2D_tt->Integral(1, 105, 1, 105);
  cout << integral1 << "   " << integral2 << "   " << integral3 << endl;


  // double sum=0; double sum1=0; double sum2=0;
  // for(int bin=0; bin<hist->GetNbinsX(); bin++){
  //   cout << "-------------------------------------- " << bin << endl;
  //   cout << hist->Integral(bin+1, bin+1) << endl;
  //   cout << hist->GetBinContent(bin+1) << endl;
  // }
}
