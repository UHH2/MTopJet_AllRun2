#include "../include/CentralInclude.h"

using namespace std;

// compile with:
// g++ -o resolution_subjets_flavorJEC resolution_subjets_flavorJEC.cc `root-config --cflags --glibs`

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown);
vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error);

int main(int argc, char* argv[]){

  TString channel = "muon";
  bool use_median = false;
  if(argc > 1 && strcmp(argv[1], "median") == 0) use_median = true;


  /*
  .█████  ██      ██          ██    ██ ███████  █████  ██████  ███████
  ██   ██ ██      ██           ██  ██  ██      ██   ██ ██   ██ ██
  ███████ ██      ██            ████   █████   ███████ ██████  ███████
  ██   ██ ██      ██             ██    ██      ██   ██ ██   ██      ██
  ██   ██ ███████ ███████        ██    ███████ ██   ██ ██   ██ ███████
  */

  // Do it more compact ... with vectors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vector<string> dirs = {"RecGenHists_subjets/", "RecGenHists_subjets_noJEC/", "RecGenHists_subjets_corrected/"};


  vector<TFile*> files;
  TFile *file_2016 = new TFile("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *file_2017 = new TFile("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2017/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *file_2018 = new TFile("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2018/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  files.push_back(file_2016);
  files.push_back(file_2017);
  files.push_back(file_2018);

  int n_ptbin = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TString bin_strings[] = {"0", "50", "80", "120", "170", "220", "270", "320", "370", "420", "600"};

  // now get plots from file
  vector<TH1F*> reso_16, reso_17, reso_18;
  vector<TH1F*> reso_cor_16, reso_cor_17, reso_cor_18;
  vector<TH1F*> reso_noJEC_16, reso_noJEC_17, reso_noJEC_18;
  vector<TH1F*> ptrec_16, ptrec_17, ptrec_18;

  std::string dir_jec = "RecGenHists_subjets/";
  std::string dir_raw = "RecGenHists_subjets_noJEC/";
  std::string dir_cor = "RecGenHists_subjets_corrected/";

  std::string name_jec;
  std::string name_raw;
  std::string name_cor;
  std::string name_ptrec;

  vector<vector<TH1F*>> reso_full = {reso_16, reso_17, reso_18};
  vector<vector<TH1F*>> reso_cor_full = {reso_cor_16, reso_cor_17, reso_cor_18};
  vector<vector<TH1F*>> reso_noJEC_full = {reso_noJEC_16, reso_noJEC_17, reso_noJEC_18};
  vector<vector<TH1F*>> ptrec_full = {ptrec_16, ptrec_17, ptrec_18};

  for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){
    name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
    name_raw = dir_raw + "PtResolution_" + std::to_string(ptbin);
    name_cor = dir_cor + "PtResolution_" + std::to_string(ptbin);
    name_ptrec = dir_jec + "PtRec_" + std::to_string(ptbin);


    for(unsigned int i; i < reso_full.size(); i++){
      reso_full[i].push_back( (TH1F*)files[i]->Get(name_jec.c_str()) );
      reso_cor_full[i].push_back( (TH1F*)files[i]->Get(name_cor.c_str()) );
      reso_noJEC_full[i].push_back( (TH1F*)files[i]->Get(name_raw.c_str()) );
      ptrec_full[i].push_back( (TH1F*)files[i]->Get(name_ptrec.c_str()) );
    }
    // //-------------- 16 ------------------------------------------------------
    // reso_16.push_back( (TH1F*)file_2016->Get(name_jec.c_str()) );
    // reso_cor_16.push_back( (TH1F*)file_2016->Get(name_cor.c_str()) );
    // reso_noJEC_16.push_back( (TH1F*)file_2016->Get(name_raw.c_str()) );
    // ptrec_16.push_back( (TH1F*)file_2016->Get(name_ptrec.c_str()) );
    //
    // //-------------- 17 ------------------------------------------------------
    // reso_17.push_back( (TH1F*)file_2017->Get(name_jec.c_str()) );
    // reso_cor_17.push_back( (TH1F*)file_2017->Get(name_cor.c_str()) );
    // reso_noJEC_17.push_back( (TH1F*)file_2018->Get(name_raw.c_str()) );
    // ptrec_17.push_back( (TH1F*)file_2017->Get(name_ptrec.c_str()) );
    //
    // //-------------- 18 ------------------------------------------------------
    // reso_18.push_back( (TH1F*)file_2018->Get(name_jec.c_str()) );
    // reso_cor_18.push_back( (TH1F*)file_2018->Get(name_cor.c_str()) );
    // reso_noJEC_18.push_back( (TH1F*)file_2018->Get(name_raw.c_str()) );
    // ptrec_18.push_back( (TH1F*)file_2018->Get(name_ptrec.c_str()) );
  };


  // ===========================================================================
  TH1F* resolution_16, *resolution_17, *resolution_18;
  TH1F* resolution_cor_16, *resolution_cor_17, *resolution_cor_18;
  TH1F* resolution_noJEC_16, *resolution_noJEC_17, *resolution_noJEC_18;
  vector<double> pt_16, pt_17, pt_18;
  vector<double> pt_err_16, pt_err_17, pt_err_18;
  vector<double> MeanForUncert_16, MeanForUncert_17, MeanForUncert_18;
  vector<double> MeanForUncert_err_16, MeanForUncert_err_17, MeanForUncert_err_18;
  TGraphErrors* non_closure_16, *non_closure_17, *non_closure_18;
  int n_last_16, n_last_17, n_last_18;
  TGraph* area_16, *area_17, *area_18;
  TGraph* upvar_16, *upvar_17, *upvar_18;
  TGraph* downvar_16, *downvar_17, *downvar_18;

  // ---------------- Full Vectors ---------------------------------------------
  vector<TH1F*> resolution_full = {resolution_16, resolution_17, resolution_18};
  vector<TH1F*> resolution_cor_full = {resolution_cor_16, resolution_cor_17, resolution_cor_18};
  vector<TH1F*> resolution_noJEC_full = {resolution_noJEC_16, resolution_noJEC_17, resolution_noJEC_18};

  vector<vector<double>> pt_full = {pt_16, pt_17, pt_18};
  vector<vector<double>> pt_err_full = {pt_err_16, pt_err_17, pt_err_18};

  vector<vector<double>> MeanForUncert_full = {MeanForUncert_16, MeanForUncert_17, MeanForUncert_18};
  vector<vector<double>> MeanForUncert_err_full = {MeanForUncert_err_16, MeanForUncert_err_17, MeanForUncert_err_18};

  vector<TGraphErrors*> non_closure_full = {non_closure_16, non_closure_17, non_closure_18};

  vector<int> n_last_full = {n_last_16, n_last_17, n_last_18};
  vector<TGraph*> area_full = {area_16, area_17, area_18};
  vector<TGraph*> upvar_full = {upvar_16, upvar_17, upvar_18};
  vector<TGraph*> downvar_full= {downvar_16, downvar_17, downvar_18};

  for(unsigned int i=0; i < 3; i++){
    resolution_full[i] = GetResoPlot(reso_full[i], use_median);
    resolution_cor_full[i] = GetResoPlot(reso_cor_full[i], use_median);
    resolution_noJEC_full[i] = GetResoPlot(reso_noJEC_full[i], use_median);

    pt_full[i] = GetMeans(ptrec_16, use_median, true, "mean");
  }

  //-------------- 16 --------------------------------------------------------
  // fit histograms and extract mean/median value
  resolution_16  = GetResoPlot(reso_16, use_median);
  resolution_cor_16 = GetResoPlot(reso_cor_16, use_median);
  resolution_noJEC_16= GetResoPlot(reso_noJEC_16, use_median);

  // calculate non-closure as function of ptrec
  // arguments: (hists, use median?, do ptrec?, mean or error)
  pt_16 = GetMeans(ptrec_16, use_median, true, "mean");
  pt_err_16 = GetMeans(ptrec_16, use_median, true, "error");
  MeanForUncert_16 = GetMeans(reso_cor_16, use_median, false, "mean");
  MeanForUncert_err_16 = GetMeans(reso_cor_16, use_median, false, "error");
  non_closure_16 = new TGraphErrors(n_ptbin, &pt_16[0], &MeanForUncert_16[0], &pt_err_16[0], &MeanForUncert_err_16[0]);

  // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
  // do same to reach ptrec = 600
  non_closure_16->SetPoint(non_closure_16->GetN(), 0, MeanForUncert_16[0]);
  n_last_16 = MeanForUncert_16.size()-1;
  non_closure_16->SetPoint(non_closure_16->GetN(), 600, MeanForUncert_16[n_last_16]);

  area_16 = AreaFromEnvelope(non_closure_16, "area");
  upvar_16 = AreaFromEnvelope(non_closure_16, "up");
  downvar_16 = AreaFromEnvelope(non_closure_16, "down");

  //-------------- 17 --------------------------------------------------------
  // fit histograms and extract mean/median value
  resolution_17  = GetResoPlot(reso_17, use_median);
  resolution_cor_17 = GetResoPlot(reso_cor_17, use_median);
  resolution_noJEC_17 = GetResoPlot(reso_noJEC_17, use_median);

  // calculate non-closure as function of ptrec
  // arguments: (hists, use median?, do ptrec?, mean or error)
  pt_17 = GetMeans(ptrec_17, use_median, true, "mean");
  pt_err_17 = GetMeans(ptrec_17, use_median, true, "error");
  MeanForUncert_17 = GetMeans(reso_cor_17, use_median, false, "mean");
  MeanForUncert_err_17 = GetMeans(reso_cor_17, use_median, false, "error");
  non_closure_17 = new TGraphErrors(n_ptbin, &pt_17[0], &MeanForUncert_17[0], &pt_err_17[0], &MeanForUncert_err_17[0]);

  // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
  // do same to reach ptrec = 600
  non_closure_17->SetPoint(non_closure_17->GetN(), 0, MeanForUncert_17[0]);
  n_last_17 = MeanForUncert_17.size()-1;
  non_closure_17->SetPoint(non_closure_17->GetN(), 600, MeanForUncert_17[n_last_17]);

  area_17 = AreaFromEnvelope(non_closure_17, "area");
  upvar_17 = AreaFromEnvelope(non_closure_17, "up");
  downvar_17 = AreaFromEnvelope(non_closure_17, "down");

  //-------------- 18 --------------------------------------------------------
  // fit histograms and extract mean/median value
  resolution_18  = GetResoPlot(reso_18, use_median);
  resolution_cor_18 = GetResoPlot(reso_cor_18, use_median);
  resolution_noJEC_18 = GetResoPlot(reso_noJEC_18, use_median);

  // calculate non-closure as function of ptrec
  // arguments: (hists, use median?, do ptrec?, mean or error)
  pt_18 = GetMeans(ptrec_18, use_median, true, "mean");
  pt_err_18 = GetMeans(ptrec_18, use_median, true, "error");
  MeanForUncert_18 = GetMeans(reso_cor_18, use_median, false, "mean");
  MeanForUncert_err_18 = GetMeans(reso_cor_18, use_median, false, "error");
  non_closure_18 = new TGraphErrors(n_ptbin, &pt_18[0], &MeanForUncert_18[0], &pt_err_18[0], &MeanForUncert_err_18[0]);

  // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
  // do same to reach ptrec = 600
  non_closure_18->SetPoint(non_closure_18->GetN(), 0, MeanForUncert_18[0]);
  n_last_18 = MeanForUncert_18.size()-1;
  non_closure_18->SetPoint(non_closure_18->GetN(), 600, MeanForUncert_18[n_last_18]);

  area_18 = AreaFromEnvelope(non_closure_18, "area");
  upvar_18 = AreaFromEnvelope(non_closure_18, "up");
  downvar_18 = AreaFromEnvelope(non_closure_18, "down");

  //==========================================================================

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  resolution_16->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_16->GetXaxis()->SetNdivisions(505);
  resolution_16->GetYaxis()->SetNdivisions(505);
  resolution_16->GetXaxis()->SetTitleSize(0.05);
  resolution_16->GetYaxis()->SetTitleSize(0.04);
  resolution_16->GetXaxis()->SetTitleOffset(0.9);
  resolution_16->GetYaxis()->SetTitleOffset(1.5);
  resolution_16->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_16->GetYaxis()->SetTitle("mean #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_16->SetLineWidth(4);
  resolution_16->SetLineColor(kRed);
  resolution_cor_16->SetLineWidth(4);
  resolution_cor_16->SetLineColor(kRed);
  resolution_noJEC_16->SetLineWidth(4);
  resolution_noJEC_16->SetLineColor(kRed);

  non_closure_16->GetXaxis()->SetRangeUser(0, 600);
  non_closure_16->GetYaxis()->SetRangeUser(-5, 5);
  non_closure_16->GetXaxis()->SetNdivisions(505);
  non_closure_16->GetYaxis()->SetNdivisions(505);
  non_closure_16->GetXaxis()->SetTitleSize(0.05);
  non_closure_16->GetYaxis()->SetTitleSize(0.04);
  non_closure_16->GetXaxis()->SetTitleOffset(0.9);
  non_closure_16->GetYaxis()->SetTitleOffset(1.5);
  non_closure_16->GetXaxis()->SetTitle("p_{T}^{rec}");
  non_closure_16->GetYaxis()->SetTitle("non-closure uncertainty [%]");
  non_closure_16->SetTitle(" ");
  non_closure_16->SetMarkerStyle(20);
  non_closure_16->SetMarkerSize(0.8);

  resolution_17->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_17->GetXaxis()->SetNdivisions(505);
  resolution_17->GetYaxis()->SetNdivisions(505);
  resolution_17->GetXaxis()->SetTitleSize(0.05);
  resolution_17->GetYaxis()->SetTitleSize(0.04);
  resolution_17->GetXaxis()->SetTitleOffset(0.9);
  resolution_17->GetYaxis()->SetTitleOffset(1.5);
  resolution_17->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_17->GetYaxis()->SetTitle("mean #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_17->SetLineWidth(4);
  resolution_17->SetLineColor(kGreen);
  resolution_cor_17->SetLineWidth(4);
  resolution_cor_17->SetLineColor(kGreen);
  resolution_noJEC_17->SetLineWidth(4);
  resolution_noJEC_17->SetLineColor(kGreen);

  non_closure_17->GetXaxis()->SetRangeUser(0, 600);
  non_closure_17->GetYaxis()->SetRangeUser(-5, 5);
  non_closure_17->GetXaxis()->SetNdivisions(505);
  non_closure_17->GetYaxis()->SetNdivisions(505);
  non_closure_17->GetXaxis()->SetTitleSize(0.05);
  non_closure_17->GetYaxis()->SetTitleSize(0.04);
  non_closure_17->GetXaxis()->SetTitleOffset(0.9);
  non_closure_17->GetYaxis()->SetTitleOffset(1.5);
  non_closure_17->GetXaxis()->SetTitle("p_{T}^{rec}");
  non_closure_17->GetYaxis()->SetTitle("non-closure uncertainty [%]");
  non_closure_17->SetTitle(" ");
  non_closure_17->SetMarkerStyle(20);
  non_closure_17->SetMarkerSize(0.8);

  resolution_18->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_18->GetXaxis()->SetNdivisions(505);
  resolution_18->GetYaxis()->SetNdivisions(505);
  resolution_18->GetXaxis()->SetTitleSize(0.05);
  resolution_18->GetYaxis()->SetTitleSize(0.04);
  resolution_18->GetXaxis()->SetTitleOffset(0.9);
  resolution_18->GetYaxis()->SetTitleOffset(1.5);
  resolution_18->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_18->GetYaxis()->SetTitle("mean #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_18->SetLineWidth(4);
  resolution_18->SetLineColor(kAzure+7);
  resolution_cor_18->SetLineWidth(4);
  resolution_cor_18->SetLineColor(kAzure+7);
  resolution_noJEC_18->SetLineWidth(4);
  resolution_noJEC_18->SetLineColor(kAzure+7);

  non_closure_18->GetXaxis()->SetRangeUser(0, 600);
  non_closure_18->GetYaxis()->SetRangeUser(-5, 5);
  non_closure_18->GetXaxis()->SetNdivisions(505);
  non_closure_18->GetYaxis()->SetNdivisions(505);
  non_closure_18->GetXaxis()->SetTitleSize(0.05);
  non_closure_18->GetYaxis()->SetTitleSize(0.04);
  non_closure_18->GetXaxis()->SetTitleOffset(0.9);
  non_closure_18->GetYaxis()->SetTitleOffset(1.5);
  non_closure_18->GetXaxis()->SetTitle("p_{T}^{rec}");
  non_closure_18->GetYaxis()->SetTitle("non-closure uncertainty [%]");
  non_closure_18->SetTitle(" ");
  non_closure_18->SetMarkerStyle(20);
  non_closure_18->SetMarkerSize(0.8);
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *c0 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_16->Draw("E0");
  zero_line->Draw("SAME");
  resolution_17->Draw("SAME E0");
  resolution_18->Draw("SAME E0");
  TLegend *leg0 = new TLegend(0.45,0.85,0.80,0.65);
  leg0->AddEntry(resolution_16,"Standard JEC applied for 16","l");
  leg0->AddEntry(resolution_17,"Standard JEC applied for 17","l");
  leg0->AddEntry(resolution_18,"Standard JEC applied for 18","l");
  leg0->SetTextSize(0.05);
  leg0->Draw();
  TLatex text0;
  text0.SetNDC(kTRUE);
  text0.SetTextFont(43);
  text0.SetTextSize(18);
  text0.DrawLatex(.2,.2, "lepton+jets");
  gPad->RedrawAxis();
  // first save without additional correction
  c0->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Resolution_Subjets/"+channel+"/pt_mean_compare_all_years_noAdditional.pdf");

  TCanvas *c1 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_cor_18->Draw("E1");
  zero_line->Draw("SAME");
  resolution_cor_17->Draw("SAME E1");
  resolution_cor_16->Draw("SAME E1");
  TLegend *leg1 = new TLegend(0.30,0.40,0.65,0.25);
  leg1->AddEntry(resolution_cor_16,"additional Corrections applied for 16","l");
  leg1->AddEntry(resolution_cor_17,"additional Corrections applied for 17","l");
  leg1->AddEntry(resolution_cor_18,"additional Corrections applied for 18","l");
  leg1->SetTextSize(0.05);
  leg1->Draw();
  TLatex text1;
  text1.SetNDC(kTRUE);
  text1.SetTextFont(43);
  text1.SetTextSize(18);
  text1.DrawLatex(.2,.2, "lepton+jets");
  gPad->RedrawAxis();
  // first save without additional correction
  c1->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Resolution_Subjets/"+channel+"/pt_mean_compare_all_years_Additional.pdf");

  TCanvas *c2 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_noJEC_18->Draw("E1");
  zero_line->Draw("SAME");
  resolution_noJEC_17->Draw("SAME E1");
  resolution_noJEC_16->Draw("SAME E1");
  TLegend *leg2 = new TLegend(0.45,0.85,0.80,0.65);
  leg2->AddEntry(resolution_noJEC_16,"no JEC applied for 16","l");
  leg2->AddEntry(resolution_noJEC_17,"no JEC applied for 17","l");
  leg2->AddEntry(resolution_noJEC_18,"no JEC applied for 18","l");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  TLatex text2;
  text2.SetNDC(kTRUE);
  text2.SetTextFont(43);
  text2.SetTextSize(18);
  text2.DrawLatex(.2,.2, "lepton+jets");
  gPad->RedrawAxis();
  // first save without additional correction
  c2->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Resolution_Subjets/"+channel+"/pt_mean_compare_all_years_noJEC.pdf");

  return 0;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
.██████  ███████ ████████     ██████  ███████ ███████  ██████
██       ██         ██        ██   ██ ██      ██      ██    ██
██   ███ █████      ██        ██████  █████   ███████ ██    ██
██    ██ ██         ██        ██   ██ ██           ██ ██    ██
.██████  ███████    ██        ██   ██ ███████ ███████  ██████
*/

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median){
  int n_bins = hists.size();
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TH1F * reso = new TH1F("reso"," ",n_bins, bins);
  vector<double> mean, error;

  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1 = new TF1("f1","gaus",-0.2,0.2);
    hists[i]->Fit("f1", "QR");
    double m = f1->GetParameter(1);
    double e = f1->GetParError(1);
    if(!use_median) mean.push_back(m);
    error.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  TAxis *xaxis = reso->GetXaxis();
  for(unsigned int i=0; i<n_bins; i++){
    reso->Fill(xaxis->GetBinCenter(i+1), mean[i]);
    reso->SetBinError(i+1, error[i]);
  }
  return reso;
}
//------------------------------------------------------------------------------
/*
.██████  ███████ ████████     ███    ███ ███████  █████  ███    ██ ███████
██       ██         ██        ████  ████ ██      ██   ██ ████   ██ ██
██   ███ █████      ██        ██ ████ ██ █████   ███████ ██ ██  ██ ███████
██    ██ ██         ██        ██  ██  ██ ██      ██   ██ ██  ██ ██      ██
.██████  ███████    ██        ██      ██ ███████ ██   ██ ██   ████ ███████
*/

vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error){
  int n_bins = hists.size();
  vector<double> mean, err;

  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1;
    if(do_ptrec) f1 = new TF1("f1","gaus");
    else         f1 = new TF1("f1","gaus", -0.2, 0.2);
    TString options = "QR";
    if(do_ptrec) options = "Q";
    hists[i]->Fit("f1", options);
    double m = f1->GetParameter(1);
    if(!use_median) mean.push_back(m);
    double e = f1->GetParError(1);
    err.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  if(error == "error") return err;
  else return mean;
}

//------------------------------------------------------------------------------
/*
███    ██  ██████  ███    ██        ██████ ██       ██████  ███████ ██    ██ ██████  ███████
████   ██ ██    ██ ████   ██       ██      ██      ██    ██ ██      ██    ██ ██   ██ ██
██ ██  ██ ██    ██ ██ ██  ██ █████ ██      ██      ██    ██ ███████ ██    ██ ██████  █████
██  ██ ██ ██    ██ ██  ██ ██       ██      ██      ██    ██      ██ ██    ██ ██   ██ ██
██   ████  ██████  ██   ████        ██████ ███████  ██████  ███████  ██████  ██   ██ ███████
*/
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown){
  int steps = 150;
  double stepsize = 700/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  // first generate new TGraph with all positive values points
  const int Npoints = uncert->GetN();
  double xPoint[Npoints];
  double yPoint[Npoints];
  for(unsigned int i=0; i<Npoints; i++){
    double y;
    uncert->GetPoint(i, xPoint[i], y);
    yPoint[i] = sqrt(y*y); // take abs for y-value
  }
  TGraph* uncert_abs = new TGraph( Npoints, xPoint, yPoint );

  //

  for(int i=0; i<steps; i++){
    x = stepsize * i;
    up[i] = sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    down[i] = - sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    xvalues[i] = x;
  }
  TGraph *uncertfit;
  if(updown == "area"){
    uncertfit = new TGraph(2*steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
      uncertfit->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
    }
  }
  else if(updown == "up"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
    }
  }
  else if(updown == "down"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],down[i]);
    }
  }
  return uncertfit;
}
