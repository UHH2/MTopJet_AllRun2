
#include "../include/CentralInclude.h"


using namespace std;

int main(int argc, char* argv[]){

  TString histname = "PreSel04_Muon/pt_1";
  TString histname_for_file = "pt_1";

  vector<TString> processes = {"DATA", "TTbar", "WJets", "SingleTop", "other"};

  TString name_out = "PtMuonSYS";
  TFile * file_out = new TFile(dir+name_out+".root","RECREATE");

  // first write nominal hist (all Processes)
  for(unsigned int i=0; i<processes.size(); i++){
    TString filename = dir;
    filename += "uhh2.AnalysisModuleRunner.";
    if(processes[i] == "DATA"){
      filename += "DATA.DATA.root";
    }
    else{
      filename += "MC.";
      filename += processes[i];
      filename += ".root";
    }
    TFile *file = new TFile(filename);
    TH1F* hist = (TH1F*)file->Get(histname);
    file_out->cd();
    hist->Write(histname_for_file+"__"+processes[i]);
  }

  TFile* btagup = new TFile(dir+"BTAG_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* btagdown = new TFile(dir+"BTAG_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  TH1F* histup = (TH1F*)btagup->Get(histname);
  TH1F* histdown = (TH1F*)btagdown->Get(histname);
  file_out->cd();

  histup->Write(histname_for_file+"__"+"TTbar"+"__BTAG__plus");
  histdown->Write(histname_for_file+"__"+"TTbar"+"__BTAG__minus");

  file_out->Close();

  return 0;
}
