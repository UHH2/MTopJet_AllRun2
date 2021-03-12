#include "../include/CentralInclude.h"

using namespace std;

TH1F* get_hist(TFile* file, TH1F* dummy, TString obs_name, TString sel_name, TString weightname);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  vector<TString> obsnames = {"Mass_Rec", "Pt_Rec", "mW", "ptsub1", "ptsub2", "ptsub3", "ptsub3_subsel"};
  vector<TH1F*> dummyhists;
  dummyhists.push_back(new TH1F("mjet", "mjet", 25, 0, 500));
  dummyhists.push_back(new TH1F("pt", "pt", 30, 400, 1000));
  dummyhists.push_back(new TH1F("mW", "mW", 36, 0, 180));
  dummyhists.push_back(new TH1F("ptsub1", "ptsub1", 35, 0, 700));
  dummyhists.push_back(new TH1F("ptsub2", "ptsub2", 25, 0, 500));
  dummyhists.push_back(new TH1F("ptsub3", "ptsub3", 30, 0, 300));
  dummyhists.push_back(new TH1F("ptsub3_subsel", "ptsub3_subsel", 10, 0, 30));

  TString directory = "/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/";
  TString prefix_mc = "uhh2.AnalysisModuleRunner.MC.";
  TString prefix_data = "uhh2.AnalysisModuleRunner.DATA.";

  vector<TString> channels = {"elec", "muon"};
  vector<TString> processes = {"DATA", "TTbar", "SingleTop", "WJets", "other"};
  vector<TString> years = {"2016v3", "2017v2", "2018"};
  double sf16 = 0.778;
  double sf17 = 0.870;
  double sf18 = 0.881;

  vector<TString> systematics = {"muid_up", "muid_down", "mutr_up", "mutr_down", "elid_up", "elid_down", "eltr_up", "eltr_down", "elreco_up", "elreco_down", "pu_up", "pu_down", "btag_up", "btag_down", "JEC_up", "JEC_down", "JER_up", "JER_down", "COR_up", "COR_down", "JMS_up", "JMS_down"};
  // vector<TString> systematics = {"JEC_up", "JEC_down", "JER_up", "JER_down", "COR_up", "COR_down", "JMS_up", "JMS_down"};

  TFile * outfile = new TFile("RecoLevelPlots.root","RECREATE");
  for(int i=0; i<obsnames.size(); i++){
    // Fill nominal hists
    for(auto process: processes){
      bool isfirsthist = true;
      TH1F* hist;
      for(auto year: years){
        double sf;
        if(year == "2016v3")      sf = sf16;
        else if(year == "2017v2") sf = sf17;
        else if(year == "2018")   sf = sf18;
        for(auto channel: channels){
          TString prefix = prefix_mc;
          if(process == "DATA") prefix = prefix_data;
          TFile* file = new TFile(dir+"/"+channel+"/"+prefix+process+"_"+year+".root");
          if(isfirsthist){
            hist = get_hist(file, dummyhists[i], obsnames[i], "passed_measurement_rec", "none");
            if(obsnames[i] == "ptsub3_subsel") hist = get_hist(file, dummyhists[i], "ptsub3", "passed_subptmigration_rec", "none");
            if(process == "TTbar") hist->Scale(sf);
            isfirsthist = false;
          }
          else{
            TH1F* h = get_hist(file, dummyhists[i], obsnames[i], "passed_measurement_rec", "none");
            if(obsnames[i] == "ptsub3_subsel") h = get_hist(file, dummyhists[i], "ptsub3", "passed_subptmigration_rec", "none");
            if(process == "TTbar") h->Scale(sf);
            hist->Add(h);
          }
        }
      }
      outfile->cd();
      hist->Write(obsnames[i]+"__"+process);
    }

    for(auto sys: systematics){
      bool isfirsthist = true;
      TH1F* hist;
      for(auto year: years){
        double sf;
        if(year == "2016v3")      sf = sf16;
        else if(year == "2017v2") sf = sf17;
        else if(year == "2018")   sf = sf18;
        for(auto channel: channels){
          TFile* file;
          TString weightname = "sf_"+sys;
          if(sys.Contains("JEC") || sys.Contains("JER")  || sys.Contains("COR")  || sys.Contains("JMS")){
            TString subdir = sys;
            if(subdir == "JMS_up") subdir = "JMS_upup";
            if(subdir == "JMS_down") subdir = "JMS_downdown";

            file = new TFile(dir+"/"+channel+"/"+subdir+"/"+prefix_mc+"TTbar"+"_"+year+".root");
            weightname = "none";
          }
          else{
            file = new TFile(dir+"/"+channel+"/"+prefix_mc+"TTbar"+"_"+year+".root");
          }

          if(isfirsthist){
            hist = get_hist(file, dummyhists[i], obsnames[i], "passed_measurement_rec", weightname);
            if(obsnames[i] == "ptsub3_subsel") hist = get_hist(file, dummyhists[i], "ptsub3", "passed_subptmigration_rec", "none");
            hist->Scale(sf);
            isfirsthist = false;
          }
          else{
            TH1F* h = get_hist(file, dummyhists[i], obsnames[i], "passed_measurement_rec", weightname);
            if(obsnames[i] == "ptsub3_subsel") h = get_hist(file, dummyhists[i], "ptsub3", "passed_subptmigration_rec", "none");
            h->Scale(sf);
            hist->Add(h);
          }
        }
      }

      TString sysstring;
      if     (sys.Contains("up"))   sysstring = sys.ReplaceAll("_up", "__plus");
      else if(sys.Contains("down")) sysstring = sys.ReplaceAll("_down", "__minus");

      cout << sysstring << endl;
      outfile->cd();
      hist->Write(obsnames[i]+"__"+"TTbar"+"__"+sysstring);
    }

  }

  outfile->Close();

return 0;
}

TH1F* get_hist(TFile* file, TH1F* dummyhist, TString obs_name, TString sel_name, TString weightname){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  TH1F* h_hist = (TH1F*) dummyhist->Clone();

  Double_t obs;
  Bool_t passed_measurement_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress(obs_name, &obs);
  tree->SetBranchAddress(sel_name, &passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  if(weightname != "none"){
    cout << "  - using weight" << endl;
    tree->SetBranchAddress(weightname,&additional_factor);
  }
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(weightname != "none") rec_weight *= additional_factor;
    weight = rec_weight * gen_weight;
    if(passed_measurement_rec) h_hist->Fill(obs, weight);
  }
  return h_hist;
}
