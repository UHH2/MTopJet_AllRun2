#include "modelsys.h"

ModelSys::ModelSys(TH1* unfolded, TH1D* truth){
  TH1D* diff = GetDiff(unfolded, truth);
  CalculateCov(diff);
}

TH2* ModelSys::GetCov(){
  return cov;
}

void ModelSys::CalculateCov(TH1D* hist){
  int nbins = hist->GetSize();
  cov = new TH2D("cov", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  for(int i=1; i <= nbins; i++){
    for(int j=1; i <= nbins; j++){
      double error = hist->GetBinError(i) * hist->GetBinError(j);
      cov->SetBinContent(i,i,error);
    }
  }
  return;
}

TH1D* ModelSys::GetDiff(TH1* unfolded, TH1D* truth){
  TH1D* diff = (TH1D*) truth->Clone();
  diff->Reset();
  int nbins = truth->GetSize() - 2;
  for(int i=1; i<=nbins; i++){
    double entry = unfolded->GetBinContent(i) - truth->GetBinContent(i);
    diff->SetBinContent(i, entry);
  }
  return diff;
}
