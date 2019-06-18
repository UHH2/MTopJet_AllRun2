#include "normalise.h"

Normalise::Normalise(TH1* hist, TH2* oldMatrix, double lower, double upper, bool width){
  width_ = width;
  SetHist(hist);
  SetRange(lower, upper);
  TH2D* newMatrix = MatrixErrorPropagation(oldMatrix);
  SetMatrix(newMatrix);
  Norm();
  UpdateHist();
  delete newMatrix;
}

Normalise::Normalise(TH1* hist, double lower, double upper, bool width){
  width_ = width;
  SetHist(hist);
  SetRange(lower, upper);
  TH2D* oldMatrix = CalcMatrix();
  TH2D* newMatrix = MatrixErrorPropagation(oldMatrix);
  SetMatrix(newMatrix);
  Norm();
  UpdateHist();
  delete oldMatrix;
  delete newMatrix;
}

Normalise::~Normalise(){
  delete hist_;
  delete matrix_;
}

TH1D* Normalise::GetHist(){
  return hist_;
}

TH2D* Normalise::GetMatrix(){
  return matrix_;
}

void Normalise::SetHist(TH1* hist){
  hist_ = (TH1D*)hist->Clone("hist_");
  return;
}

void Normalise::SetMatrix(TH2* matrix){
  matrix_ = (TH2D*)matrix->Clone("matrix_");
  return;
}

void Normalise::SetRange(double lower, double upper){
  lower_ = hist_->GetXaxis()->FindBin(lower);
  upper_ = hist_->GetXaxis()->FindBin(upper);
  return;
}

void Normalise::Norm(){
  if(width_) hist_->Scale(1/hist_->Integral(lower_, upper_),"width");
  else      hist_->Scale(1/hist_->Integral(lower_, upper_));
}

TH2D* Normalise::CalcMatrix(){
  int nbins = upper_;
  TH2D* cov = new TH2D("cov", "cov", nbins, 1, upper_+1, nbins, 1, upper_+1);
  for(int i=1; i <= upper_; i++){
    double error = hist_->GetBinError(i);
    cov->Fill(i,i, error * error);
  }
  return cov;
}


TH2D* Normalise::MatrixErrorPropagation(TH2*old_cov){
  int nbins = upper_;
  TH2D* new_cov = (TH2D*) old_cov->Clone();
  new_cov->Reset();
  double integral = hist_->Integral(lower_, upper_);
  // std::cout << "lower bin: "<< lower_ << std::endl;
  // std::cout << "upper bin: "<< upper_ << std::endl;
  for(int i=1; i <= upper_; i++){
    for(int j=1; j <= upper_; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= upper_; k++){
        for(int l=1; l <= upper_; l++){
          old_entry = old_cov->GetBinContent(k,l);
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
      // std::cout << sum << std::endl;
      new_cov->SetBinContent(i, j, sum);
    }
  }
  return new_cov;
}

void Normalise::UpdateHist(){
  int nbins = hist_->GetXaxis()->GetNbins();
  for(int i=0; i<=nbins; i++){
    if(i < lower_ || i > upper_){
      hist_->SetBinContent(i, 0.0);
      hist_->SetBinError(i, 0.0);
    }
    else{
      double error = sqrt(matrix_->GetBinContent(i,i));
      hist_->SetBinError(i, error);
    }
  }
  return;
}
