#include "../include/CentralInclude.h"


// compile with:
// g++ -o tt_pdf tt_pdf.cc `root-config --cflags --glibs`

std::vector< double > GetRMS(std::vector< TH1F* > hists);
std::vector< double > GetRMS(TH1F* central, std::vector< TH1F* > hists);
TH1F* GetVariationUp(TH1F* central, std::vector< double > rms, TString name);
TH1F* GetVariationDown(TH1F* central, std::vector< double > rms, TString name);

using namespace std;

int main(int argc, char* argv[]){


  TFile* TT_f = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TH1F* CentralA = (TH1F*)TT_f->Get("XCone_cor/M_jet1_");  // 20 GeV bins
  TH1F* CentralB = (TH1F*)TT_f->Get("XCone_cor/M_jet1");   // 10 GeV bins
  TH1F* CentralC = (TH1F*)TT_f->Get("XCone_cor/M_jet1_B"); //  5 GeV bins


  TString histdir, histnameA, histnameB, histnameC;
  histdir = "PDFHists/";
  histnameA = "M_jet1_A_PDF_";  // 20 GeV bins
  histnameB = "M_jet1_B_PDF_";  // 10 GeV bins
  histnameC = "M_jet1_C_PDF_";  //  5 GeV bins

  std::vector<TH1F*> HistsA, HistsB, HistsC;

  for(unsigned int i=1; i<= 100; i++){
    std::string num = std::to_string(i);
    HistsA.push_back((TH1F*)TT_f->Get(histdir + histnameA + num));
    HistsB.push_back((TH1F*)TT_f->Get(histdir + histnameB + num));
    HistsC.push_back((TH1F*)TT_f->Get(histdir + histnameC + num));
  }

  vector<double> RMS_A = GetRMS(CentralA, HistsA);
  vector<double> RMS_B = GetRMS(CentralB, HistsB);
  vector<double> RMS_C = GetRMS(CentralC, HistsC);


  TH1F* VariationA_up = GetVariationUp(CentralA, RMS_A, "M_jet1_");
  TH1F* VariationA_down = GetVariationDown(CentralA, RMS_A, "M_jet1_");
  TH1F* VariationB_up = GetVariationUp(CentralB, RMS_B, "M_jet1");
  TH1F* VariationB_down = GetVariationDown(CentralB, RMS_B, "M_jet1");
  TH1F* VariationC_up = GetVariationUp(CentralC, RMS_C, "M_jet1_B");
  TH1F* VariationC_down = GetVariationDown(CentralC, RMS_C, "M_jet1_B");

  TFile* PDF_up = new TFile(dir+"PDF_up/PDF_Variations.root","RECREATE");
  VariationA_up->Write();
  VariationB_up->Write();
  VariationC_up->Write();
  PDF_up->Close();

  TFile* PDF_down = new TFile(dir+"PDF_down/PDF_Variations.root","RECREATE");
  VariationA_down->Write();
  VariationB_down->Write();
  VariationC_down->Write();
  PDF_down->Close();

  return 0;
}

// there are two ways to calculate the RMS for the 100 Histograms per bin:
// 1) Fill the bin content of the 100 Hists in another histogram and use GetRMS2
// 2) Use formula for standard deviation and set mean value of this distribution
//    to bin content of central Hist (). This is written in the instruction in UHH2 wiki.

// Method 1
std::vector< double > GetRMS(std::vector< TH1F* > hists){
  std::vector< double > v_RMS;
  int Nbins = hists[0]->GetSize() - 2;
  for(int bin=1; bin<=Nbins; bin++){
    TH1F* newhist = new TH1F("newhist", "newhist", 1000, 0, 10000);
    for(int j=0; j<hists.size(); j++){
      double entry = hists[j]->GetBinContent(bin);
      newhist->Fill(entry);
    }
    v_RMS.push_back(newhist->GetRMS());
    delete newhist;
  }
  return v_RMS;
}

// Method 2
std::vector< double > GetRMS(TH1F* central, std::vector< TH1F* > hists){
  std::vector< double > v_RMS;
  int Nbins = central->GetSize() - 2;
  for(int bin=1; bin<=Nbins; bin++){
    double sum = 0;  //sum[(x_i - <x>)^2]
    for(int j=0; j<hists.size(); j++){
      double diff = hists[j]->GetBinContent(bin) - central->GetBinContent(bin);
      sum += diff*diff;
    }
    double rms = sqrt( sum/(hists.size()) );
    v_RMS.push_back(rms);
  }
  return v_RMS;
}

TH1F* GetVariationUp(TH1F* central, std::vector< double > rms, TString name){
  int Nbins = rms.size();
  TH1F* hist = (TH1F*) central->Clone(name);
  hist->Reset();
  for(unsigned int i=1; i<=Nbins; i++){
    double ctr = central->GetBinContent(i);
    hist->SetBinContent(i, ctr + rms[i-1]);
  }
  return hist;
}

TH1F* GetVariationDown(TH1F* central, std::vector< double > rms, TString name){
  int Nbins = rms.size();
  TH1F* hist = (TH1F*) central->Clone(name);
  hist->Reset();
  for(unsigned int i=1; i<=Nbins; i++){
    double ctr = central->GetBinContent(i);
    hist->SetBinContent(i, ctr - rms[i-1]);
  }
  return hist;
}
