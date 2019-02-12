#include "../include/CentralInclude.h"


using namespace std;

// compile with:
// g++ -o resolution_subjets_flavorJEC_updown resolution_subjets_flavorJEC_updown.cc `root-config --cflags --glibs`

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown);
vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error);
TGraph* GetNonClosure(TH1F* central, TH1F* up, TH1F* down, vector<double> ptrec);

int main(int argc, char* argv[]){

  bool use_median = false;
  if(argc >1 && strcmp(argv[1], "median") == 0) use_median = true;

  // declare files
  TFile *f_central = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *f_up = new TFile(dir+"JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *f_down = new TFile(dir+"JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  // set binning
  int n_ptbin = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};

  // this is for naming of pdf files
  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  // now get plots from file
  vector<TH1F*> reso_central, reso_up, reso_down, h_ptrec;

  std::string dir_jec = "RecGenHists_subjets/";

  std::string name_jec;
  std::string name_ptrec;

  for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){
    name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
    name_ptrec  = dir_jec + "PtRec_" + std::to_string(ptbin);
    reso_central.push_back( (TH1F*)f_central->Get(name_jec.c_str()) );
    reso_up.push_back( (TH1F*)f_up->Get(name_jec.c_str()) );
    reso_down.push_back( (TH1F*)f_down->Get(name_jec.c_str()) );
    h_ptrec.push_back( (TH1F*)f_central->Get(name_ptrec.c_str()) );
  }

  // fit histograms and extract mean/median value
  TH1F * resolution_central  = GetResoPlot(reso_central, use_median);
  TH1F * resolution_up  = GetResoPlot(reso_up, use_median);
  TH1F * resolution_down  = GetResoPlot(reso_down, use_median);

  vector<double> ptrec = GetMeans(h_ptrec, use_median, true, "mean");

  TGraph* nonclosure = GetNonClosure(resolution_central, resolution_up, resolution_down, ptrec);

  // non-closure in root file
  TString rootname = "files/NonClosure.root";
  TFile * outfile = new TFile(rootname,"RECREATE");
  nonclosure->Write("nonclosure");
  outfile->Close();

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // start plotting

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  resolution_central->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_central->GetXaxis()->SetNdivisions(505);
  resolution_central->GetYaxis()->SetNdivisions(505);
  resolution_central->GetXaxis()->SetTitleSize(0.05);
  resolution_central->GetYaxis()->SetTitleSize(0.04);
  resolution_central->GetXaxis()->SetTitleOffset(0.9);
  resolution_central->GetYaxis()->SetTitleOffset(1.5);
  resolution_central->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_central->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_central->SetLineWidth(4);
  resolution_central->SetLineColor(kAzure+7);
  resolution_up->SetLineWidth(4);
  resolution_up->SetLineColor(kOrange+1);
  resolution_down->SetLineWidth(4);
  resolution_down->SetLineColor(kRed+1);

  nonclosure->GetYaxis()->SetRangeUser(-0.1, 0.1);
  nonclosure->GetXaxis()->SetNdivisions(505);
  nonclosure->GetYaxis()->SetNdivisions(505);
  nonclosure->GetXaxis()->SetTitleSize(0.05);
  nonclosure->GetYaxis()->SetTitleSize(0.04);
  nonclosure->GetXaxis()->SetTitleOffset(0.9);
  nonclosure->GetYaxis()->SetTitleOffset(1.5);
  nonclosure->GetXaxis()->SetTitle("p_{T}^{rec}");
  nonclosure->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  nonclosure->SetMarkerSize(1);
  nonclosure->SetMarkerStyle(8);
  nonclosure->SetMarkerColor(kBlack);
  nonclosure->SetLineColor(kBlack);
  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------

  TCanvas *c0 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_central->Draw("E1");
  zero_line->Draw("SAME");
  resolution_central->Draw("SAME E1");
  resolution_up->Draw("SAME E1");
  resolution_down->Draw("SAME E1");
  TLegend *leg0 = new TLegend(0.45,0.85,0.80,0.65);
  leg0->AddEntry(resolution_up,"Up Variation","l");
  leg0->AddEntry(resolution_central,"Central","l");
  leg0->AddEntry(resolution_down,"Down Variation","l");
  leg0->SetTextSize(0.05);
  leg0->Draw();
  // TLatex text0;
  // text0.SetNDC(kTRUE);
  // text0.SetTextFont(43);
  // text0.SetTextSize(18);
  // text0.DrawLatex(.2,.2, "standard JEC");
  gPad->RedrawAxis();
  c0->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/UpDownVariation.pdf");

  TCanvas *c1 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  nonclosure->Draw("AP");
  zero_line->Draw("SAME");
  nonclosure->Draw("SAME P");
  gPad->RedrawAxis();
  c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/NonClosure.pdf");

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
//------------------------------------------------------------------------------
/*
 ██████  ███████ ████████     ███    ███ ███████  █████  ███    ██ ███████
██       ██         ██        ████  ████ ██      ██   ██ ████   ██ ██
██   ███ █████      ██        ██ ████ ██ █████   ███████ ██ ██  ██ ███████
██    ██ ██         ██        ██  ██  ██ ██      ██   ██ ██  ██ ██      ██
 ██████  ███████    ██        ██      ██ ███████ ██   ██ ██   ████ ███████
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


/*
███    ██  ██████  ███    ██        ██████ ██       ██████  ███████ ██    ██ ██████  ███████
████   ██ ██    ██ ████   ██       ██      ██      ██    ██ ██      ██    ██ ██   ██ ██
██ ██  ██ ██    ██ ██ ██  ██ █████ ██      ██      ██    ██ ███████ ██    ██ ██████  █████
██  ██ ██ ██    ██ ██  ██ ██       ██      ██      ██    ██      ██ ██    ██ ██   ██ ██
██   ████  ██████  ██   ████        ██████ ███████  ██████  ███████  ██████  ██   ██ ███████
*/

TGraph* GetNonClosure(TH1F* central, TH1F* up, TH1F* down, vector<double> ptrec){

  unsigned int n_bins = central->GetSize() - 2;

  // now find smalles deviation from 0 in every bin
  // if uncertainties include 0, set non-closure to zero
  vector<double> values;
  for(unsigned int i=1; i<=n_bins; i++){
    double c,u,d,content;
    double pt = ptrec[i-1];
    c = central->GetBinContent(i);
    u = up->GetBinContent(i);
    d = down->GetBinContent(i);
    if(u < 0 && d > 0) values.push_back(0);
    else if (d < 0 && u > 0) values.push_back(0);
    else{ // if uncertainties do not cover 0, get non closure
      if(c < u && c < d) content = c;
      else if(d < u && d < c) content = d;
      else if(u < d && u < c) content = u;
      else cout << "[ERROR] NO VALUE FOUND FOR NON-CLOSURE!" << endl;
      values.push_back(content);
    }
  }
  int npoints = n_bins + 2;
  TGraph* nonclosure = new TGraph();
  nonclosure->SetPoint(0, 0, values[0]); // to get a starting point at pt = 0
  for(unsigned int i=0; i<values.size(); i++){
    nonclosure->SetPoint(i+1, ptrec[i], values[i]);
  }
  nonclosure->SetPoint(n_bins+1, 600, values[n_bins-1]); // to get an ending point at pt = 600

  return nonclosure;
}
