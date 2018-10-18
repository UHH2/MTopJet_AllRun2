#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include <time.h>
#include <vector>

using namespace std;

vector<double> RandomNumbers(int N, double mean, double width);
vector<TH1D*> GetVariations(int N, TH1D* original);

int main(int argc, char* argv[]){

  TFile *input = new TFile("Histograms_fine.root", "update");
  TH1D *h_original;
  input->GetObject("Generator_sig", h_original);

  vector<TH1D*> h_variations = GetVariations(200, h_original);

  for(int i=0; i<h_variations.size(); i++){
    TString name = "Generator_sig_";
    name += i;
    h_variations[i]->Write(name);
  }
  delete h_original;
  input->Close();
  return 0;
}


// this generates a gaussian distributed random number
vector<double> RandomNumbers(int N, double mean, double width){
  std::default_random_engine generator(time(0));
  std::normal_distribution<double> distribution(mean, width);
  vector<double> numbers;
  for(int i=0; i<N; i++){
    double number = distribution(generator);
    numbers.push_back(number);
  }
  return numbers;
}

vector<TH1D*> GetVariations(int N, TH1D* original){
  TH1D* h_original = (TH1D*) original->Clone();
  vector<TH1D*> h_variations;
  for(int i=0; i<N; i++){
    TString name = "Generator_sig_";
    name += i;
    TH1D* newHist = (TH1D*) original->Clone(name);
    h_variations.push_back(newHist);
  }

  int Nbins = h_original->GetSize();
  for(int bin=0; bin<=Nbins; bin++){
    double content = h_original->GetBinContent(bin);
    double error = h_original->GetBinError(bin);
    vector<double> newContent = RandomNumbers(N, content, error);
    for(int i=0; i<N; i++){
      h_variations[i]->SetBinContent(bin, newContent[i]);
    }
  }
  delete h_original;
  return h_variations;
}
