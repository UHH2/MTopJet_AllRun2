#pragma once
#include "CentralInclude.h"
#include "HistogramUtils.h"
#include "Utils.h"

////////////////////////////////////////////////////////////////////////////////
// TMatrix m.Invert() does not work if full coloums and rows are empty.
// -----------------------------------------------------------------------------
// Difference between full Covarianz Matrix and
// a Covariance Matrix for only the considered bins?
// The functions in this header aim for the full covariance matrix considering
// all bins from the histograms. In case only a few bins are used,
// e.g. only the peak bins, the matrix is cut (TrimMatrix()) after
// normaliziations and inversion.

// void DrawCov(TH2D* h, TString save, TString axis, double zoff = 0.0);
void DrawCov(TH2* h, TString save, TString axis, double zoff = 0.0);
// void DrawCov(TH2F* h, TString save, TString axis, double zoff = 0.0);
void printCovMatrix(TMatrixD m, TString name, double factor=1, int prec=2, bool diag=true);
TH2D* TMatrixDtoTH2D(TMatrixD m, double nbins, double xmin, double xmax, TString name, bool debug=false);
TMatrixD NormCovMatrix(TH1* hist_, TMatrixD old_cov, bool width_, bool onlyDiag=false);
TMatrixD GetCovMatrixDiag(TH1F* hist, bool debug);
TMatrixD GetSystCovMatrix(TH1F* hist, bool debug);
TMatrixD GetCovMatrix(TH1F* hist, bool debug);
TMatrixD TH2ToTMatrix(TH2D* hist, double mDim);
TMatrixD TrimMatrix(TMatrixD matrix, vector<bool> skip, bool debug = false);
vector<bool> TrimSkipBin(TMatrixD matrix, vector<bool> vec, bool onlyZero, bool debug = false);

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TH2D* TMatrixDtoTH2D(TMatrixD m, double nbins, double xmin, double xmax, TString name, bool debug){
  if(debug) cout << "Convert TMatrixD to TH2D for " << name << endl;
  TH2D* h = new TH2D(name, "", nbins, xmin, xmax, nbins, xmin, xmax);
  for(unsigned int c=0; c<m.GetNcols(); c++){
    for(unsigned int r=0; r<m.GetNcols(); r++){
      h->SetBinContent(c+1, r+1, m[c][r]);
    }
  }
  return h;
}

TH2F* TMatrixDtoTH2F(TMatrixD m, double nbins, double xmin, double xmax, TString name, bool debug){
  if(debug) cout << "Convert TMatrixD to TH2F for " << name << endl;
  TH2F* h = new TH2F(name, "", nbins, xmin, xmax, nbins, xmin, xmax);
  for(unsigned int c=0; c<m.GetNcols(); c++){
    for(unsigned int r=0; r<m.GetNcols(); r++){
      h->SetBinContent(c+1, r+1, m[c][r]);
    }
  }
  return h;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void printCovMatrix(TMatrixD m, TString name, double factor, int prec, bool diag){

  Int_t space = prec+8;
  cout << endl << name << " - Number of bins x: " << setw(space) << m.GetNcols() << " y: " << setw(space) << m.GetNrows() << " | factor " << setw(space) << factor << endl << endl;
  cout << setw(space) << "y\\x" << " | ";
  if(m.GetNcols()!=m.GetNrows()) cerr << "Not a symmetric Matrix (" << name << ")" << endl;
  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << setw(space) << x << " | ";
  }
  cout << endl;

  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << setw(space) << x << " | ";
    for(unsigned int y=0; y<m.GetNrows(); y++){
      if(x==y){
        cout << RED<< setw(space) << Round(m[x][y]*factor, prec) << RESET << " | ";
      }
      else cout << setw(space) << Round(m[x][y]*factor, prec) << " | ";
    }
    cout << endl;
  }
  cout << endl;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TMatrixD GetCovMatrixDiag(TH1F* hist, bool debug){

  Int_t nbins = hist->GetNbinsX();
  Int_t mDim = 0;
  std::vector<int> indexEmptyBins;
  for(unsigned int i=1; i<=nbins; i++){
    if(hist->GetBinContent(i)>0) mDim++;
    else indexEmptyBins.push_back(i);
  }

  if(debug){
    cout << "Empty bins: ";
    for(unsigned int i=0; i<indexEmptyBins.size(); i++){
      cout << setw(4) << i;
    }
    cout << endl;
  }

  TMatrixD matrix = TMatrixD(nbins, nbins);
  for(int x=1; x<=nbins; x++){
    matrix[x-1][x-1] = (hist->GetBinContent(x)>0)?pow(hist->GetBinError(x),2):0.;
  }
  return matrix;
}

TMatrixD GetCovMatrix(TH1F* hist, bool debug){

  Int_t nbins = hist->GetNbinsX();
  if(debug){
    for(unsigned int i=1; i<=nbins; i++)
    {
      if(hist->GetBinContent(i)==0) cout << setw(4) << i;
    }
    cout << endl;
  }

  TMatrixD matrix = TMatrixD(nbins, nbins);
  for(int x=1; x<=nbins; x++)
  {
    matrix[x-1][x-1] = (hist->GetBinContent(x)>0)?pow(hist->GetBinError(x),2):0.;
  }
  return matrix;
}

TMatrixD GetSystCovMatrix(TH1F* hist, bool debug){

  Int_t nbins = hist->GetNbinsX();
  if(debug){
    for(unsigned int i=1; i<=nbins; i++)
    {
      if(hist->GetBinContent(i)==0) cout << setw(4) << i;
    }
    cout << endl;
  }

  TMatrixD matrix = TMatrixD(nbins, nbins);
  for(int x=1; x<=nbins; x++)
  {
    for(int y=1; y<=nbins; y++)
    {
      matrix[x-1][y-1] = hist->GetBinContent(x)*hist->GetBinContent(y);
    }
  }
  return matrix;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TMatrixD NormCovMatrix(TH1* hist_, TMatrixD old_cov, bool width_, bool onlyDiag){

  int mDim = old_cov.GetNrows();
  if(mDim!=old_cov.GetNcols()) throw runtime_error("WRONG DIM");

  TMatrixD new_cov = TMatrix(mDim, mDim);

  double integral = hist_->Integral();

  for(int i=1; i <= mDim; i++){
    for(int j=1; j <= mDim; j++){
      if(i!=j && onlyDiag) continue;
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= mDim; k++){
        for(int l=1; l <= mDim; l++){
          old_entry = old_cov[k-1][l-1];
          double binwidth_i = hist_->GetBinWidth(i);
          double binwidth_j = hist_->GetBinWidth(j);
          if(!width_){
            binwidth_i = 1.0;
            binwidth_j = 1.0;
          }
          if(i==k) derivation_i = (integral - hist_->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (hist_->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - hist_->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (hist_->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      new_cov[i-1][j-1] = sum;
    }
  }
  return new_cov;
}


TMatrixD NormCovMatrixAlt(TH1* hist_, TMatrixD old_cov, bool width_, bool onlyDiag){

  int mDim = old_cov.GetNrows();
  if(mDim!=old_cov.GetNcols()) throw runtime_error("WRONG DIM");

  TMatrixD new_cov = TMatrix(mDim, mDim);

  double integral = hist_->Integral();
  double integral2 = pow(integral,2);

  for(int i=1; i <= mDim; i++){

    double binwidth_i = 1.0;
    if(width_) binwidth_i = hist_->GetBinWidth(i);

    double counter_i2 = hist_->GetBinContent(i);
    double counter_i1 = integral - counter_i2;
    double denominator_i = (integral2) * (binwidth_i);

    double derivation_i_eq = (counter_i1) / denominator_i;
    double derivation_i_neq = - (counter_i2) / denominator_i;

    for(int j=1; j <= mDim; j++){
      if(i!=j && onlyDiag) continue;
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;

      double counter_j2 = hist_->GetBinContent(j);
      double counter_j1 = integral - counter_j2;

      double binwidth_j = 1.0;
      if(width_) binwidth_j = hist_->GetBinWidth(j);
      double denominator_j = (integral2) * (binwidth_j);

      double derivation_j_eq = (counter_j1) / denominator_j;
      double derivation_j_neq = - (counter_j2)  / denominator_j;


      for(int k=1; k <= mDim; k++){
        for(int l=1; l <= mDim; l++){
          old_entry = old_cov[k-1][l-1];
          derivation_i = i==k?derivation_i_eq:derivation_i_neq;
          derivation_j = j==l?derivation_j_eq:derivation_j_neq;
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      new_cov[i-1][j-1] = sum;
    }
  }
  return new_cov;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

// Cut covariance Matrix to the bins which are considered in analysis, e.g. the chi2.

TMatrixD TrimMatrix(TMatrixD matrix, vector<bool> skip, bool debug){

  Int_t nbins = matrix.GetNrows();
  vector<Int_t> usedEntries;
  for(int x=0; x<nbins; x++){
    if(matrix[x][x]!=0&&!skip[x]) usedEntries.push_back(x);
  }

  Int_t dim = usedEntries.size();
  if(debug) cout << "New Dim is " << dim << endl;
  TMatrixD trimmed = TMatrixD(dim, dim);
  for(int x=0; x<dim; x++){
    Int_t i = usedEntries[x];
    for(int y=0; y<dim; y++){
      Int_t j = usedEntries[y];
      trimmed[x][y] = matrix[i][j];
    }
  }

  return trimmed;
}

// Matrix is not invertable if rows or coloumns are 0. Sometimes, these cases needs to be trimmed
// before the matrix is inverted. TrimMatrix then creats the Matrix for the considered bins.

vector<bool> TrimSkipBin(TMatrixD matrix, vector<bool> vec, bool onlyZero, bool debug){

  Int_t nbins = matrix.GetNrows();
  vector<bool> vecTrim = vec;
  for(int x=nbins-1; x>=0; x--){
    if(onlyZero){
      if(matrix[x][x]==0) vecTrim.erase(vecTrim.begin()+x);
    }
    else{
      if(matrix[x][x]==0||vec[x]) vecTrim.erase(vecTrim.begin()+x);
    }
  }

  if(debug){
    for(int i: vec) cout << i << "\t"; cout << endl;
    for(int i: vecTrim) cout << i << "\t"; cout << endl;

    vector<int> test;
    for(int x=0; x<nbins; x++) test.push_back(x);
    for(int i: test) cout << i << "\t"; cout << endl;
    for(int x=nbins-1; x>=0; x--){
      if(onlyZero){
        if(matrix[x][x]==0) test.erase(test.begin()+x);
      }
      else{
        if(matrix[x][x]==0||vec[x]) test.erase(test.begin()+x);
      }
    }
    for(int i: test) cout << i << "\t"; cout << endl;
  }

  return vecTrim;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void DrawCov(TH2* h, TString save, TString axis, double zoff){

  TCanvas *B = new TCanvas(save, "", 800, 800); // name used to avoid Warning

  B->SetRightMargin(0.09);
  B->SetLeftMargin(0.15);
  B->SetRightMargin(0.20);

  h->GetXaxis()->SetTitle(axis);
  h->GetYaxis()->SetTitle(axis);
  TString norm = (save.Contains("norm")?"normalized ":"");
  TString ztitle = norm + "uncertainty";
  h->GetZaxis()->SetTitle(ztitle);
  h->GetZaxis()->SetTitleOffset(1.0+zoff);

  h->Draw("COLZ");
  B->SaveAs(save+".pdf");
  delete B;
  // h->Reset();
}

// void DrawCov(TH2F* h, TString save, TString axis, double zoff){
//
//   TCanvas *B = new TCanvas(save, "", 800, 800); // name used to avoid Warning
//
//   B->SetRightMargin(0.09);
//   B->SetLeftMargin(0.15);
//   B->SetRightMargin(0.20);
//
//   h->GetXaxis()->SetTitle(axis);
//   h->GetYaxis()->SetTitle(axis);
//   TString norm = (save.Contains("norm")?"normalized ":"");
//   TString ztitle = norm + "uncertainty";
//   h->GetZaxis()->SetTitle(ztitle);
//   h->GetZaxis()->SetTitleOffset(1.0+zoff);
//
//   h->Draw("COLZ");
//   B->SaveAs(save+".pdf");
//
//   h->Reset();
// }
