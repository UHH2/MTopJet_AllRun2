#include "../include/CentralInclude.h"

using namespace std;
TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);


int main(int argc, char* argv[]){

  bool use_median = false;
  if(argc >1 && strcmp(argv[1], "median") == 0) use_median = true;

  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  const int n_bins = sizeof(bins)/sizeof(bins[0]) - 1;

  // declare files
  TFile *All = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_allHad.root");
  vector<string> dir = {"RecGenHists_rec0_gen400/", "RecGenHists_rec200_gen0/", "RecGenHists_rec300_gen0/", "RecGenHists_rec400_gen0/"};
  const int n_dir = dir.size();

  vector<TH1F*> resolution;

  for(unsigned int i=0; i<dir.size(); i++){
    vector<TH1F*> reso;
    for(int ptbin = 1; ptbin <= n_bins; ptbin++){
      TString name = dir[i] + "PtResolution_" + std::to_string(ptbin);
      reso.push_back( (TH1F*)All->Get(name) );
    }
    resolution.push_back( GetResoPlot(reso, use_median) );
  }

  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);

  TCanvas *c = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution[0]->SetTitle(" ");
  resolution[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution[0]->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution[0]->GetXaxis()->SetRangeUser(0, 600);
  resolution[0]->GetYaxis()->SetRangeUser(-0.05, 0.1);
  resolution[0]->GetYaxis()->SetNdivisions(505);
  resolution[0]->GetXaxis()->SetTitleSize(0.05);
  resolution[0]->GetYaxis()->SetTitleSize(0.04);
  resolution[0]->GetXaxis()->SetTitleOffset(0.9);
  resolution[0]->GetYaxis()->SetTitleOffset(1.5);

  for(unsigned int i=0; i<dir.size(); i++) resolution[i]->SetLineWidth(4);
  resolution[0]->SetLineColor(kRed+1);
  resolution[1]->SetLineColor(kOrange+1);
  resolution[2]->SetLineColor(kAzure+7);
  resolution[3]->SetLineColor(kGreen);


  resolution[0]->Draw("E1");

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);
  zero_line->Draw("SAME");

  for(unsigned int i=0; i<dir.size(); i++) resolution[i]->Draw("SAME E1");

  TLegend* leg = new TLegend(0.45,0.55,0.80,0.85);
  leg->AddEntry(resolution[0],"p_{T}^{rec} > 0, p_{T}^{gen} > 400","l");
  leg->AddEntry(resolution[1],"p_{T}^{rec} > 200, p_{T}^{gen} > 0","l");
  leg->AddEntry(resolution[2],"p_{T}^{rec} > 300, p_{T}^{gen} > 0","l");
  leg->AddEntry(resolution[3],"p_{T}^{rec} > 400, p_{T}^{gen} > 0","l");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextFont(43);
  text.SetTextSize(18);
  text.DrawLatex(.2,.2, "all hadronic");
  gPad->RedrawAxis();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/AllHadronic/pt_"+mean_median+"_ptcuts.pdf");

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
