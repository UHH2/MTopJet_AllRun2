#include "smearing.h"


Smearing::Smearing(std::vector<TH1*> outputs){
  result = (TH1*) outputs[0]->Clone();
  result->Reset();
  int nbins = outputs[0]->FindBin(1000);
  for(int i=1; i<nbins; i++){
    TH1D* dummy = new TH1D("dummy", "dummy", 20, 0, 40000);
    for(int j=0; j<outputs.size(); j++){
      double value = outputs[j]->GetBinContent(i);
      dummy->Fill(value, 1);
    }
    TF1*f1 = new TF1("f1","gaus");
    dummy->Fit("f1");
    TF1* fitresult = dummy->GetFunction("f1");
    double mean = fitresult->GetParameter(1);
    double rms = fitresult->GetParameter(2);
    variations.push_back(dummy);
    fits.push_back(fitresult);
    result->SetBinContent(i, mean);
    result->SetBinError(i, rms);
  }
  return;
}

TH1* Smearing::GetResult(){
  return result;
}

std::vector<TF1*> Smearing::GetFits(){
  return fits;
}

std::vector<TH1D*> Smearing::GetVariations(){
  return variations;
}
