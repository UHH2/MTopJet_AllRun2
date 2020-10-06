#pragma once
#include <TCanvas.h>
#include <TColor.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TVectorD.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
// #include <tuple>
#include <vector>

using namespace std;

/*
██████  ██ ███    ██
██   ██ ██ ████   ██
██████  ██ ██ ██  ██
██   ██ ██ ██  ██ ██
██████  ██ ██   ████
*/

// -------------------------------------------------------------------------------------------
vector<double> get_bin_content(TH1F* hist, int skip_bins=0){
  vector<double> bin_content_v;
  for(unsigned int bin=1+skip_bins; bin<hist->GetNbinsX()+1; bin++) bin_content_v.push_back(hist->GetBinContent(bin));
  return bin_content_v;
}

// -------------------------------------------------------------------------------------------
vector<double> get_bin_error(TH1F* hist, int skip_bins=0){
  vector<double> bin_error_v;
  for(unsigned int bin=1+skip_bins; bin<hist->GetNbinsX()+1; bin++) bin_error_v.push_back(hist->GetBinError(bin));
  return bin_error_v;
}

// -------------------------------------------------------------------------------------------
tuple<vector<double>, vector<double>> get_bin_content_and_error(TH1F* hist){
  vector<double> bin_error_v, bin_content_v;
  for(unsigned int bin=1; bin<hist->GetNbinsX()+1; bin++){
    bin_error_v.push_back(hist->GetBinError(bin));
    bin_content_v.push_back(hist->GetBinContent(bin));
  }

  return {bin_content_v, bin_error_v};
}

// -------------------------------------------------------------------------------------------
vector<int> bins_upper_limit(TH1F* hist, double limit){
  vector<int> bin_pass_v;
  for(int bin=0; bin < hist->GetNbinsX(); bin++){
    if((hist->GetBinContent(bin+1)>limit)) bin_pass_v.push_back(bin+1);
  }
  return bin_pass_v;
}

// -------------------------------------------------------------------------------------------
vector<int> bins_empty(TH1F* hist){
  vector<int> bin_empty_v;
  for(int bin=0; bin < hist->GetNbinsX(); bin++){
    if(!(abs(hist->GetBinContent(bin+1))>0)) bin_empty_v.push_back(bin+1);
  }
  return bin_empty_v;
}

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

/*
**********************************************************************************************
**************************** RELATIVE ERROR **************************************************
**********************************************************************************************
*/

// -------------------------------------------------------------------------------------------
vector<double> relative_error(TH1F* hist){
  TH1F* hist_ = (TH1F*) hist->Clone();
  int number_bins = hist_->GetNbinsX();

  vector<double> bin_rel_error_v;
  for(unsigned int bin=0; bin<number_bins; bin++) bin_rel_error_v.push_back(hist_->GetBinError(bin+1)/hist_->GetBinContent(bin+1));
  return bin_rel_error_v;
}

// -------------------------------------------------------------------------------------------
vector<double> relative_error(TH1F* hist, vector<double> error_v){
  TH1F* hist_ = (TH1F*) hist->Clone();
  int number_bins = hist_->GetNbinsX();

  vector<double> bin_rel_error_v;
  for(unsigned int bin=0; bin<number_bins; bin++) bin_rel_error_v.push_back(error_v[bin]/hist_->GetBinContent(bin+1));
  return bin_rel_error_v;
}

// -------------------------------------------------------------------------------------------
vector<double> relative_error(vector<double> content_v, vector<double> error_v){
  vector<double> bin_rel_error_v;
  for(unsigned int bin=0; bin<content_v.size(); bin++) bin_rel_error_v.push_back(error_v[bin]/content_v[bin]);
  return bin_rel_error_v;
}

/*
**********************************************************************************************
********************* Propagation of Uncertainties *******************************************
**********************************************************************************************
*/

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

vector<double> normalize_error(TH1F* hist){
  // error of each bin normalized and into vector
  TH1F* hist_ = (TH1F*) hist->Clone();
  double integral = hist_->Integral();
  double integral2 = integral*integral;
  int number_bins = hist_->GetNbinsX();

  vector<double> bin_error_norm, bin_content, bin_error;
  for(unsigned int bin=1; bin<number_bins+1; bin++){
    bin_content.push_back(hist_->GetBinContent(bin));
    bin_error.push_back(hist_->GetBinError(bin));
  }

  double error_term = 0;
  for(unsigned int bin=0; bin<number_bins; bin++){ // Calculate the error for each bin
    for(unsigned int i=0; i<number_bins; i++){ // Include contribution of each bin for PoU
      if(bin == i) error_term += pow(((integral-bin_content[bin])/integral2)*bin_error[i], 2.0);
      else         error_term += pow(((-1)*(bin_content[bin])/integral2)*bin_error[i], 2.0);
    }
    bin_error_norm.push_back(sqrt(error_term));
    error_term=0;
  }

  return bin_error_norm;
}


vector<vector<double>> normalize_error(vector<TH1F*> hist_v){
  // error of each bin normalized and into vector
  vector<vector<double>> all_errors;
  for(unsigned int i=0; i<hist_v.size(); i++){
    TH1F* hist_ = (TH1F*) hist_v[i]->Clone();
    double integral = hist_->Integral();
    double integral2 = integral*integral;
    int number_bins = hist_->GetNbinsX();

    vector<double> bin_error_norm, bin_content, bin_error;
    for(unsigned int bin=1; bin<number_bins+1; bin++){
      bin_content.push_back(hist_->GetBinContent(bin));
      bin_error.push_back(hist_->GetBinError(bin));
    }

    double error_term = 0;
    for(unsigned int bin=0; bin<number_bins; bin++){ // Calculate the error for each bin
      for(unsigned int i=0; i<number_bins; i++){ // Include contribution of each bin for PoU
        if(bin == i) error_term += pow(((integral-bin_content[bin])/integral2)*bin_error[i], 2.0);
        else         error_term += pow(((-1)*(bin_content[bin])/integral2)*bin_error[i], 2.0);
      }
      bin_error_norm.push_back(sqrt(error_term));
      error_term=0;
    }
    all_errors.push_back(bin_error_norm);
  }
  return all_errors;
}

/*
**********************************************************************************************
********************* Error from Hist difference *********************************************
**********************************************************************************************
*/

vector<double> bin_error_from_two_hists(TH1F *h1, TH1F *h2, int skip_first_bins=0){
  // vector<double> bin_content_h1, bin_content_h2;
  double bin_content_h1, bin_content_h2;
  vector<double> bin_error_v;
  int number_bins_h1 = h1->GetNbinsX();
  int number_bins_h2 = h2->GetNbinsX();
  if(number_bins_h1 != number_bins_h2) throw runtime_error("The two hists in bin_error_from_two_hists do not have the same binning");

  for(int bin=1+skip_first_bins; bin<number_bins_h1; bin++) bin_error_v.push_back(abs(bin_content_h1-bin_content_h2));
  return bin_error_v;
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
  for(unsigned int i=1; i<hists_.size();i++){
    hist_added->Add(hists_[i], factor);
  }
  return hist_added;
}

// -------------------------------------------------------------------------------------------------------
TH2F* AddHists(vector<TH2F*> hists_, int factor){
  TH2F* hist_added = (TH2F*) hists_[0]->Clone();
  for(unsigned int i=1; i<hists_.size();i++) hist_added->Add(hists_[i], factor);
  return hist_added;
}

// -------------------------------------------------------------------------------------------------------
TH1F* AddHists(TH1F* h1, TH1F* h2, int factor){
  TH1F* hist_added = (TH1F*) h1->Clone();
  hist_added->Add(h2, factor);
  return hist_added;
}

// -------------------------------------------------------------------------------------------------------
TH1F* add_second_year(TString year, TH1F* hist, TString dir, TString path, TString hist_name){
  TFile* file = new TFile(dir+year+"/muon/"+path);
  TH1F* new_hist = (TH1F*) file->Get(hist_name);
  new_hist->Add(hist, 1);
  return new_hist;
}

// -------------------------------------------------------------------------------------------------------
TH1F* add_second_year(TString year, TH1F* hist, TString dir, vector<TString> path_v, TString hist_name){
  vector<TFile*> file_v;
  vector<TH1F*> hist_v;
  for(unsigned int i=0; i<path_v.size(); i++) file_v.push_back(new TFile(dir+year+"/muon/"+path_v[i]));
  for(unsigned int i=0; i<file_v.size(); i++) hist_v.push_back((TH1F*)file_v[i]->Get(hist_name));
  TH1F* new_hist = AddHists(hist_v, 1);
  new_hist->Add(hist, 1);
  return new_hist;
}

/*
██████  ██       ██████  ████████
██   ██ ██      ██    ██    ██
██████  ██      ██    ██    ██
██      ██      ██    ██    ██
██      ███████  ██████     ██
*/

// -------------------------------------------------------------------------------------------------------
void plot_single_histogram(TH1F* hist, TString title, TString xAxis, int x_max, int color, TString save_path){
  hist->SetTitle(title);
  hist->GetXaxis()->SetRangeUser(0, x_max);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitle(xAxis);
  hist->GetYaxis()->SetTitle("");
  hist->SetLineWidth(2);
  hist->SetLineColor(color);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hist->Draw("HIST");
  A->SaveAs(save_path);
  delete A;
}
