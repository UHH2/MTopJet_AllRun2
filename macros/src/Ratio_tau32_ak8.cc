#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){
  vector<TString> year = {"2017", "2018"};

  TString ttbar_string = "uhh2.AnalysisModuleRunner.MC.TTbar.root";

  vector<TString> masscut = {"no_masscut ", "masscut_120", "masscut_130", "masscut_140", "masscut_150"};
  vector<TString> classes = {"AK8_tau_pass_rec", "AK8_tau_pass_rec_masscut_120", "AK8_tau_pass_rec_masscut_130", "AK8_tau_pass_rec_masscut_140", "AK8_tau_pass_rec_masscut_150"};
  vector<TH1F*> hists;

  TFile *file_ttbar;

  // At first only ttbar
  for(unsigned int i=0; i<year.size(); i++){

    cout<< "" << endl;
    cout<< "####################################################################" << endl;
    cout<< "############################# "+year[i]+" #################################" << endl;
    cout<< "####################################################################" << endl;
    cout<< "" << endl;

    file_ttbar = new TFile(dir+year[i]+"/muon/"+ttbar_string); // file for each year
    for(unsigned int j=0; j<classes.size(); j++) hists.push_back((TH1F*)file_ttbar->Get(classes[j]+"/events_matched_all")); // push_back Histogramm for N-subjettiness for each masscut

    for(unsigned int j=0; j<hists.size(); j++){
      double ratio = hists[j]->GetBinContent(2)/(hists[j]->GetBinContent(1)+hists[j]->GetBinContent(2));
      cout<< "--------------------- "+masscut[j]+" --------------------------------" << endl;
      cout<< "matched:     " << hists[j]->GetBinContent(2) <<"     Ratio to total #events: " << ratio<<endl;
    }

    hists = {}; // Empty Hist at the end of the loop to avoid Segmentation Violation
  }
}
