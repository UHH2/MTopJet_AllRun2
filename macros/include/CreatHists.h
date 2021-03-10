#include "CentralInclude.h"

using namespace std;

// channels
TString muon = "muon";
TString elec = "elec";

// Years
TString y16 = "_2016v3";
TString y17 = "_2017v2";
TString y18 = "_2018";

// Hist files
TString data_f  = "uhh2.AnalysisModuleRunner.DATA.DATA";
TString ttbar_f = "uhh2.AnalysisModuleRunner.MC.TTbar";
TString wjets_f = "uhh2.AnalysisModuleRunner.MC.WJets";
TString st_f    = "uhh2.AnalysisModuleRunner.MC.SingleTop";
TString other_f = "uhh2.AnalysisModuleRunner.MC.other";

// Hist paths
TString data_muon_16 = dir+muon+"/"+data_f+y16+".root";
TString data_muon_17 = dir+muon+"/"+data_f+y17+".root";
TString data_muon_18 = dir+muon+"/"+data_f+y18+".root";
TString data_elec_16 = dir+elec+"/"+data_f+y16+".root";
TString data_elec_17 = dir+elec+"/"+data_f+y17+".root";
TString data_elec_18 = dir+elec+"/"+data_f+y18+".root";
vector<TString> data_muon = {data_muon_16, data_muon_17, data_muon_18};
vector<TString> data_elec = {data_elec_16, data_elec_17, data_elec_18};

TString ttbar_muon_16 = dir+muon+"/"+ttbar_f+y16+".root";
TString ttbar_muon_17 = dir+muon+"/"+ttbar_f+y17+".root";
TString ttbar_muon_18 = dir+muon+"/"+ttbar_f+y18+".root";
TString ttbar_elec_16 = dir+elec+"/"+ttbar_f+y16+".root";
TString ttbar_elec_17 = dir+elec+"/"+ttbar_f+y17+".root";
TString ttbar_elec_18 = dir+elec+"/"+ttbar_f+y18+".root";
vector<TString> ttbar_muon = {ttbar_muon_16, ttbar_muon_17, ttbar_muon_18};
vector<TString> ttbar_elec = {ttbar_elec_16, ttbar_elec_17, ttbar_elec_18};

TString wjets_muon_16 = dir+muon+"/"+wjets_f+y16+".root";
TString wjets_muon_17 = dir+muon+"/"+wjets_f+y17+".root";
TString wjets_muon_18 = dir+muon+"/"+wjets_f+y18+".root";
TString wjets_elec_16 = dir+elec+"/"+wjets_f+y16+".root";
TString wjets_elec_17 = dir+elec+"/"+wjets_f+y17+".root";
TString wjets_elec_18 = dir+elec+"/"+wjets_f+y18+".root";
vector<TString> wjets_muon = {wjets_muon_16, wjets_muon_17, wjets_muon_18};
vector<TString> wjets_elec = {wjets_elec_16, wjets_elec_17, wjets_elec_18};

TString st_muon_16 = dir+muon+"/"+st_f+y16+".root";
TString st_muon_17 = dir+muon+"/"+st_f+y17+".root";
TString st_muon_18 = dir+muon+"/"+st_f+y18+".root";
TString st_elec_16 = dir+elec+"/"+st_f+y16+".root";
TString st_elec_17 = dir+elec+"/"+st_f+y17+".root";
TString st_elec_18 = dir+elec+"/"+st_f+y18+".root";
vector<TString> st_muon = {st_muon_16, st_muon_17, st_muon_18};
vector<TString> st_elec = {st_elec_16, st_elec_17, st_elec_18};

TString other_muon_16 = dir+muon+"/"+other_f+y16+".root";
TString other_muon_17 = dir+muon+"/"+other_f+y17+".root";
TString other_muon_18 = dir+muon+"/"+other_f+y18+".root";
TString other_elec_16 = dir+elec+"/"+other_f+y16+".root";
TString other_elec_17 = dir+elec+"/"+other_f+y17+".root";
TString other_elec_18 = dir+elec+"/"+other_f+y18+".root";
vector<TString> other_muon = {other_muon_16, other_muon_17, other_muon_18};
vector<TString> other_elec = {other_elec_16, other_elec_17, other_elec_18};

// Creating Hists
vector<TH1F*> get_hists(vector<TString> names, TString hist){
  vector<TFile*> files;
  vector<TH1F*> hists;
  for(TString file: names){files.push_back(new TFile(file));}
  for(TFile* file: files){hists.push_back((TH1F*)file->Get(hist));}
  hists.push_back(AddHists(hists, 1)); // Add combine (HistogramUtils)
  return(hists);
}
