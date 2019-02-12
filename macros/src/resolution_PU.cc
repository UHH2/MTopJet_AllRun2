#include "../include/CentralInclude.h"

using namespace std;

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);

int main(int argc, char* argv[]){

  bool use_median = false;
  if(argc >1 && strcmp(argv[1], "median") == 0) use_median = true;

  // declare files
  TFile *f_tt = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  // set binning
  int n_ptbin = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TString bin_strings[] = {"0", "50", "80", "120", "170", "220", "270", "320", "370", "420", "600"};

  // this is for naming of pdf files
  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  // now get plots from file
  vector<vector<TH1F*>> reso_jec;
  vector<vector<TH1F*>> reso_raw;

  TString dir_jec = "RecGenHists_passedRec";
  TString dir_raw = "RecGenHists_passedRec_noJEC";
  vector<TString> PU_scenario = {"lowPU", "medPU", "highPU"};

  for(auto PU: PU_scenario){
    vector<TH1F*> r_jec;
    vector<TH1F*> r_raw;
    for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){
      TString name_jec = dir_jec + "_" + PU + "/PtResolution_" + std::to_string(ptbin);
      TString name_raw = dir_raw + "_" + PU + "/PtResolution_" + std::to_string(ptbin);
      r_jec.push_back( (TH1F*)f_tt->Get(name_jec) );
      r_raw.push_back( (TH1F*)f_tt->Get(name_raw) );
    }
    reso_jec.push_back(r_jec);
    reso_raw.push_back(r_raw);
  }

  // fit histograms and extract mean/median value
  vector<TH1F*> resolution_raw;
  vector<TH1F*> resolution_jec;

  for(unsigned int i=0; i<PU_scenario.size(); i++){
    resolution_jec.push_back(GetResoPlot(reso_jec[i], use_median));
    resolution_raw.push_back(GetResoPlot(reso_raw[i], use_median));
  }

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<PU_scenario.size(); i++){
    resolution_jec[i]->GetYaxis()->SetRangeUser(-0.1, 0.15);
    resolution_jec[i]->GetXaxis()->SetNdivisions(505);
    resolution_jec[i]->GetYaxis()->SetNdivisions(505);
    resolution_jec[i]->GetXaxis()->SetTitleSize(0.05);
    resolution_jec[i]->GetYaxis()->SetTitleSize(0.04);
    resolution_jec[i]->GetXaxis()->SetTitleOffset(0.9);
    resolution_jec[i]->GetYaxis()->SetTitleOffset(1.5);
    resolution_jec[i]->GetXaxis()->SetTitle("p_{T}^{gen}");
    resolution_jec[i]->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
    resolution_jec[i]->SetLineWidth(4);
    resolution_raw[i]->GetYaxis()->SetRangeUser(-0.1, 0.15);
    resolution_raw[i]->GetXaxis()->SetNdivisions(505);
    resolution_raw[i]->GetYaxis()->SetNdivisions(505);
    resolution_raw[i]->GetXaxis()->SetTitleSize(0.05);
    resolution_raw[i]->GetYaxis()->SetTitleSize(0.04);
    resolution_raw[i]->GetXaxis()->SetTitleOffset(0.9);
    resolution_raw[i]->GetYaxis()->SetTitleOffset(1.5);
    resolution_raw[i]->GetXaxis()->SetTitle("p_{T}^{gen}");
    resolution_raw[i]->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
    resolution_raw[i]->SetLineWidth(4);
  }
  resolution_jec[0]->SetLineColor(kRed+1);
  resolution_jec[1]->SetLineColor(kAzure+7);
  resolution_jec[2]->SetLineColor(kOrange+1);
  resolution_raw[0]->SetLineColor(kRed+1);
  resolution_raw[1]->SetLineColor(kAzure+7);
  resolution_raw[2]->SetLineColor(kOrange+1);

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  vector<TString> Leg_Text = {" 0 #leq N^{PV} < 10", "10 #leq N^{PV} < 20", "20 #leq N^{PV}"};

  TCanvas *c1 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_jec[0]->Draw("E1");
  zero_line->Draw("SAME");
  TLegend *leg1 = new TLegend(0.45,0.85,0.80,0.65);
  for(unsigned int i=0; i<PU_scenario.size(); i++){
    resolution_jec[i]->Draw("SAME E1");
    leg1->AddEntry(resolution_jec[i],Leg_Text[i],"l");
  }
  leg1->SetTextSize(0.05);
  leg1->Draw();
  TLatex text1;
  text1.SetNDC(kTRUE);
  text1.SetTextFont(43);
  text1.SetTextSize(18);
  text1.DrawLatex(.2,.2, "lepton+jets, Standard AK4 JEC applied");
  gPad->RedrawAxis();
  c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_"+mean_median+"_jec_PU.pdf");

  TCanvas *c2 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_raw[0]->Draw("E1");
  zero_line->Draw("SAME");
  TLegend *leg2 = new TLegend(0.45,0.85,0.80,0.65);
  for(unsigned int i=0; i<PU_scenario.size(); i++){
    resolution_raw[i]->Draw("SAME E1");
    leg2->AddEntry(resolution_raw[i],Leg_Text[i],"l");
  }
  leg2->SetTextSize(0.05);
  leg2->Draw();
  TLatex text2;
  text2.SetNDC(kTRUE);
  text2.SetTextFont(43);
  text2.SetTextSize(18);
  text2.DrawLatex(.2,.2, "lepton+jets, no JEC applied");
  gPad->RedrawAxis();
  c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_"+mean_median+"_raw_PU.pdf");


  return 0;
}
//------------------------------------------------------------------------------
/*
██████  ███████ ████████     ██████  ███████ ███████  ██████
██       ██         ██        ██   ██ ██      ██      ██    ██
██   ███ █████      ██        ██████  █████   ███████ ██    ██
██    ██ ██         ██        ██   ██ ██           ██ ██    ██
██████  ███████    ██        ██   ██ ███████ ███████  ██████
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
