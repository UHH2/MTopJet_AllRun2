#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  vector<TString> years = {"2016", "2017", "2018"};

  TString path       = "";
  TString data_file  = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TString mc_file    = "uhh2.AnalysisModuleRunner.MC.";

  TString hist

  vector<TString> names = {"TTbar", "SingleTop", "WJets", "QCD", "DYJets", "Diboson"};
  vector<TFile*>  files;
  vector<TH1*>    hists;



  for(TString year: years){

    for(TString name: names){
      TString file_folder = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/"+year+"/muon/"
      TFile* file = new TFile(file_folder+mc_file+name+".root");
    }

    TString other_file  = mc_file+"other.root";

    THStack all_hists   = new THStack("all", "");
    THStack other_hists = new THStack("other", "");

  }
}
