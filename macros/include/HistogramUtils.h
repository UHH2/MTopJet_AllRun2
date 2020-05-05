#pragma once
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <iostream>
#include <tuple>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <TColor.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>

using namespace std;

/*
██████  ███████ ██████  ██ ███    ██
██   ██ ██      ██   ██ ██ ████   ██
██████  █████   ██████  ██ ██ ██  ██
██   ██ ██      ██   ██ ██ ██  ██ ██
██   ██ ███████ ██████  ██ ██   ████
*/


// -------------------------------------------------------------------------------------------
TH1F *rebin(TH1F *hist, int bin){
  TH1F *hist_bin = (TH1F*) hist->Clone();
  hist_bin->Rebin(bin);
  return hist_bin;
}

// -------------------------------------------------------------------------------------------
vector<TH1F*> rebin(vector<TH1F*> hist, int bin){
  vector<TH1F*> hist_bin_v;
  for(unsigned int i=0; i<hist.size(); i++) hist_bin_v.push_back(rebin(hist[i], bin));
  return hist_bin_v;
}

// -------------------------------------------------------------------------------------------
vector<vector<TH1F*>> rebin(vector<vector<TH1F*>> hist, int bin){
  vector<vector<TH1F*>> hist_bin_vv;
  for(unsigned int i=0; i<hist.size(); i++) hist_bin_vv.push_back(rebin(hist[i], bin));
  return hist_bin_vv;
}

/*
███    ██  ██████  ██████  ███    ███  █████  ██      ██ ███████ ███████
████   ██ ██    ██ ██   ██ ████  ████ ██   ██ ██      ██    ███  ██
██ ██  ██ ██    ██ ██████  ██ ████ ██ ███████ ██      ██   ███   █████
██  ██ ██ ██    ██ ██   ██ ██  ██  ██ ██   ██ ██      ██  ███    ██
██   ████  ██████  ██   ██ ██      ██ ██   ██ ███████ ██ ███████ ███████
*/


// -------------------------------------------------------------------------------------------
TH1F *normalize(TH1F *hist){
  TH1F *hist_norm = (TH1F*) hist->Clone();
  hist_norm->Scale(1/hist->Integral());
  return hist_norm;
}

// -------------------------------------------------------------------------------------------
vector<TH1F*> normalize(vector<TH1F*> hist){
  vector<TH1F*> hist_norm_v;
  for(unsigned int i=0; i<hist.size(); i++) hist_norm_v.push_back(normalize(hist[i]));
  return hist_norm_v;
}

// -------------------------------------------------------------------------------------------
vector<vector<TH1F*>> normalize(vector<vector<TH1F*>> hist){
  vector<vector<TH1F*>> hist_norm_vv;
  for(unsigned int i=0; i<hist.size(); i++) hist_norm_vv.push_back(normalize(hist[i]));
  return hist_norm_vv;
}


/*
███████ ██████  ██████   ██████  ██████
██      ██   ██ ██   ██ ██    ██ ██   ██
█████   ██████  ██████  ██    ██ ██████
██      ██   ██ ██   ██ ██    ██ ██   ██
███████ ██   ██ ██   ██  ██████  ██   ██
*/

// -------------------------------------------------------------------------------------------
vector<float> relative_error(TH1F* hist){
  TH1F* hist_ = (TH1F*) hist->Clone();
  int number_bins = hist_->GetNbinsX();

  vector<float> bin_rel_error_v;
  for(unsigned int bin=0; bin<number_bins; bin++) bin_rel_error_v.push_back(hist_->GetBinError(bin+1)/hist_->GetBinContent(bin+1));
  return bin_rel_error_v;
}

// -------------------------------------------------------------------------------------------
vector<float> relative_error(TH1F* hist, vector<float> error_v){
  TH1F* hist_ = (TH1F*) hist->Clone();
  int number_bins = hist_->GetNbinsX();

  vector<float> bin_rel_error_v;
  for(unsigned int bin=0; bin<number_bins; bin++) bin_rel_error_v.push_back(error_v[bin]/hist_->GetBinContent(bin+1));
  return bin_rel_error_v;
}

// -------------------------------------------------------------------------------------------
vector<float> relative_error(vector<float> content_v, vector<float> error_v){
  vector<float> bin_rel_error_v;
  for(unsigned int bin=0; bin<content_v.size(); bin++) bin_rel_error_v.push_back(error_v[bin]/content_v[bin]);
  return bin_rel_error_v;
}

// -------------------------------------------------------------------------------------------
/* Propagation of Uncertainties: many variables --------------------------------------------
dy = dely/delx_1 * dx_1 + dely/delx_2 * dx_2 + ... + dely/delx_n * dx_n
Gauss: dy^2 = (dely/delx_1 * dx_1)^2 + (dely/delx_2 * dx_2)^2 + ... + (dely/delx_n * dx_n)^2
...... variables are independent  of each other

PoU for normalized hists: ------------------------------------------------------------------
I = Integral of Histogram
bin x_1 and x_1,n(orm)
x_1,n = x_1/(Sum(all bins))
two different terms for PoU of dx_1,n:
delx_1,n/delx_1 * dx_1 = (I-x_1)/(I)^2 * dx_1
delx_1,n/delx_y * dx_y = -(x_1)/(I)^2 * dx_y (y = every bin except 1 in this example) */

vector<float> normalize_error(TH1F* hist){
  // error of each bin normalized and into vector
  TH1F* hist_ = (TH1F*) hist->Clone();
  float integral = hist_->Integral();
  float integral2 = integral*integral;
  int number_bins = hist_->GetNbinsX();

  vector<float> bin_error_norm, bin_content, bin_error;
  for(unsigned int bin=1; bin<number_bins+1; bin++){
    bin_content.push_back(hist_->GetBinContent(bin));
    bin_error.push_back(hist_->GetBinError(bin));
  }

  float error_term = 0;
  for(unsigned int bin=0; bin<number_bins; bin++){ // Calculate the error for each bin
    for(unsigned int i=0; i<number_bins; i++){ // Include contribution of each bin for PoU
      if(bin == i) error_term =+ pow(((integral-bin_content[bin])/integral2)*bin_error[i], 2);
      else         error_term =+ pow(((-1)*(bin_content[bin])/integral2)*bin_error[i], 2);
    }
    bin_error_norm.push_back(sqrt(error_term));
    error_term=0;
  }

  return bin_error_norm;
}

/*
██   ██ ██  ██████  ██   ██ ███████ ███████ ████████     ██████  ███████  █████  ██   ██
██   ██ ██ ██       ██   ██ ██      ██         ██        ██   ██ ██      ██   ██ ██  ██
███████ ██ ██   ███ ███████ █████   ███████    ██        ██████  █████   ███████ █████
██   ██ ██ ██    ██ ██   ██ ██           ██    ██        ██      ██      ██   ██ ██  ██
██   ██ ██  ██████  ██   ██ ███████ ███████    ██        ██      ███████ ██   ██ ██   ██
*/


// -------------------------------------------------------------------------------------------
double get_highest_peak(vector<TH1F*> hists){
  vector<double> max;
  double max_value;
  for(unsigned int i = 0; i<hists.size(); i++) max.push_back(hists[i]->GetMaximum());
  sort(max.begin(), max.end());
  double top_index = max.size()-1;
  return max[top_index];
}

/*
.█████  ██████  ██████  ██   ██ ██ ███████ ████████
██   ██ ██   ██ ██   ██ ██   ██ ██ ██         ██
███████ ██   ██ ██   ██ ███████ ██ ███████    ██
██   ██ ██   ██ ██   ██ ██   ██ ██      ██    ██
██   ██ ██████  ██████  ██   ██ ██ ███████    ██
*/


TH1F* AddHists(vector<TH1F*> hists_, int factor){
  TH1F* hist_added = (TH1F*) hists_[0]->Clone();
  for(unsigned int i=1; i<hists_.size();i++) hist_added->Add(hists_[i], factor);
  return hist_added;
}

// -------------------------------------------------------------------------------------------------------
TH1F* AddHists(TH1F* h1, TH1F* h2, int factor){
  TH1F* hist_added = (TH1F*) h1->Clone();
  hist_added->Add(h2, factor);
  return hist_added;
}
