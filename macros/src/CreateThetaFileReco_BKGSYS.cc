#include "../include/CentralInclude.h"

using namespace std;

vector<TH1F*> get_hists(TFile* file, vector<TH1F*> dummy, vector<TString> obs_name, TString sel_name, TString weightname);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  vector<TString> obsnames = {"Mass_Rec", "Pt_Rec", "mW", "ptsub1", "ptsub2", "ptsub3", "ptsub3_subpt"};
  vector<TH1F*> dummyhists;
  dummyhists.push_back(new TH1F("mjet", "mjet", 25, 0, 500));
  dummyhists.push_back(new TH1F("pt", "pt", 30, 400, 1000));
  dummyhists.push_back(new TH1F("mW", "mW", 36, 0, 180));
  dummyhists.push_back(new TH1F("ptsub1", "ptsub1", 35, 0, 700));
  dummyhists.push_back(new TH1F("ptsub2", "ptsub2", 25, 0, 500));
  dummyhists.push_back(new TH1F("ptsub3", "ptsub3", 30, 0, 300));
  dummyhists.push_back(new TH1F("ptsub3_subpt", "ptsub3_subpt", 30, 0, 60));

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

  TFile * outfile = new TFile("RecoLevelPlots_SYSBKG.root","RECREATE");
  // Fill nominal hists
  cout << "Fill nominal hists" << endl;
  vector<TH1F*> h_tt, h_st, h_wj, h_ot;
  bool firstbkg=true;
  for(auto process: processes){
    vector<vector<TH1F*>> h_all_years;
    for(auto year: years){
      double sf;
      if(year == "2016v3")      sf = sf16;
      else if(year == "2017v2") sf = sf17;
      else if(year == "2018")   sf = sf18;
      TString prefix = prefix_mc;
      if(process == "DATA") prefix = prefix_data;
      TFile* file_el = new TFile(dir+"/elec/"+prefix+process+"_"+year+".root");
      TFile* file_mu = new TFile(dir+"/muon/"+prefix+process+"_"+year+".root");
      cout << " - fill " << process << " " << year << endl;
      vector<TH1F*> hists_el = get_hists(file_el, dummyhists, obsnames, "passed_measurement_rec", "none");
      vector<TH1F*> hists_mu = get_hists(file_mu, dummyhists, obsnames, "passed_measurement_rec", "none");
      vector<TH1F*> hists_combined;
      for(unsigned int i=0; i<hists_el.size(); i++){
        TH1F * hcombime = (TH1F*) hists_el[i]->Clone();
        hcombime->Add(hists_mu[i]);
        if(process == "TTbar") hcombime->Scale(sf);
        hists_combined.push_back(hcombime);
      }
      h_all_years.push_back(hists_combined);
    }
    for(unsigned int i=0; i<obsnames.size();i++){
      TH1F * h_16 = h_all_years[0][i];
      TH1F * h_17 = h_all_years[1][i];
      TH1F * h_18 = h_all_years[2][i];
      TH1F * hist = (TH1F*) h_16->Clone();
      hist->Add(h_17);
      hist->Add(h_18);
      outfile->cd();
      hist->Write(obsnames[i]+"__"+process);
      if(process == "TTbar")          h_tt.push_back((TH1F*) hist->Clone());
      else if(process == "SingleTop") h_st.push_back((TH1F*) hist->Clone());
      else if(process == "WJets")     h_wj.push_back((TH1F*) hist->Clone());
      else if(process == "other")     h_ot.push_back((TH1F*) hist->Clone());
    }
  }
  cout << "Calculate BKG rate uncerts" << endl;
  for(unsigned int i=0; i<obsnames.size();i++){
    TH1F* h_st_plus = (TH1F*) h_tt[i]->Clone();
    TH1F* h_st_minus = (TH1F*) h_tt[i]->Clone();
    TH1F* h_wj_plus = (TH1F*) h_tt[i]->Clone();
    TH1F* h_wj_minus = (TH1F*) h_tt[i]->Clone();
    TH1F* h_ot_plus = (TH1F*) h_tt[i]->Clone();
    TH1F* h_ot_minus = (TH1F*) h_tt[i]->Clone();
    h_st_plus ->Add(h_st[i],  0.23);
    h_st_minus->Add(h_st[i], -0.23);
    h_wj_plus ->Add(h_wj[i],  0.19);
    h_wj_minus->Add(h_wj[i], -0.19);
    h_ot_plus ->Add(h_ot[i],  1.0);
    h_ot_minus->Add(h_ot[i], -1.0);
    h_st_plus->Write(obsnames[i]+"__"+"TTbar"+"__SingleTopRate__plus");
    h_st_minus->Write(obsnames[i]+"__"+"TTbar"+"__SingleTopRate__minus");
    h_wj_plus->Write(obsnames[i]+"__"+"TTbar"+"__WJetsRate__plus");
    h_wj_minus->Write(obsnames[i]+"__"+"TTbar"+"__WJetsRate__minus");
    h_ot_plus->Write(obsnames[i]+"__"+"TTbar"+"__OtherRate__plus");
    h_ot_minus->Write(obsnames[i]+"__"+"TTbar"+"__OtherRate__minus");
  }


  cout << "Fill systematics" << endl;
  for(auto sys: systematics){
    for(auto process: processes){
      if(process == "DATA") continue;
      vector<vector<TH1F*>> h_all_years;
      for(auto year: years){
        double sf;
        if(year == "2016v3")      sf = sf16;
        else if(year == "2017v2") sf = sf17;
        else if(year == "2018")   sf = sf18;
        TFile *file_mu, *file_el;
        TString weightname = "sf_"+sys;
        if(sys.Contains("JEC") || sys.Contains("JER")  || sys.Contains("COR")  || sys.Contains("JMS")){
          TString subdir = sys;
          if(subdir == "JMS_up") subdir = "JMS_upup";
          if(subdir == "JMS_down") subdir = "JMS_downdown";
          file_el = new TFile(dir+"/elec/"+subdir+"/"+prefix_mc+process+"_"+year+".root");
          file_mu = new TFile(dir+"/muon/"+subdir+"/"+prefix_mc+process+"_"+year+".root");
          weightname = "none";
        }
        else{
          file_el = new TFile(dir+"/elec/"+prefix_mc+process+"_"+year+".root");
          file_mu = new TFile(dir+"/muon/"+prefix_mc+process+"_"+year+".root");
        }
        cout << " - fill " << sys << " " << process << " "<< year << endl;
        vector<TH1F*> hists_el = get_hists(file_el, dummyhists, obsnames, "passed_measurement_rec", weightname);
        vector<TH1F*> hists_mu = get_hists(file_mu, dummyhists, obsnames, "passed_measurement_rec", weightname);
        vector<TH1F*> hists_combined;
        for(unsigned int i=0; i<hists_el.size(); i++){
          TH1F * hcombime = (TH1F*) hists_el[i]->Clone();
          hcombime->Add(hists_mu[i]);
          if(process == "TTbar") hcombime->Scale(sf);
          hists_combined.push_back(hcombime);
        }
        h_all_years.push_back(hists_combined);
      }
      for(unsigned int i=0; i<obsnames.size();i++){
        TH1F * h_16 = h_all_years[0][i];
        TH1F * h_17 = h_all_years[1][i];
        TH1F * h_18 = h_all_years[2][i];
        TH1F * hist = (TH1F*) h_16->Clone();
        hist->Add(h_17);
        hist->Add(h_18);
        TString sysstring;
        TString sysdummy = sys;
        if     (sysdummy.Contains("up"))   sysstring = sysdummy.ReplaceAll("_up", "__plus");
        else if(sysdummy.Contains("down")) sysstring = sysdummy.ReplaceAll("_down", "__minus");
        // cout << sysstring << endl;
        outfile->cd();
        hist->Write(obsnames[i]+"__"+process+"__"+sysstring);
      }
    }
  }



  outfile->Close();

  return 0;
}

vector<TH1F*> get_hists(TFile* file, vector<TH1F*> dummys, vector<TString> obs_name, TString sel_name, TString weightname){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  vector<TH1F*> hists;
  for(auto h: dummys) hists.push_back((TH1F*) h->Clone());

  vector<Double_t> obs;
  for(unsigned int i=0; i<obs_name.size(); i++) obs.push_back(-1);

  Bool_t passed_measurement_rec, passed_subpt_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;
  int index_pt3=-1;
  int index_ptsub=-1;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress(sel_name, &passed_measurement_rec);
  tree->SetBranchAddress("passed_subptmigration_rec", &passed_subpt_rec);

  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  for(unsigned int i=0; i<obs_name.size(); i++){
    if(obs_name[i]=="ptsub3") index_pt3 = i;
    if(obs_name[i]=="ptsub3_subpt") index_ptsub = i;
    else tree->SetBranchAddress(obs_name[i], &obs[i]);
  }

  if(weightname != "none"){
    cout << "  - using weight" << endl;
    tree->SetBranchAddress(weightname,&additional_factor);
  }
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(weightname != "none") rec_weight *= additional_factor;
    weight = rec_weight * gen_weight;
    if(passed_measurement_rec){
      for(unsigned int i=0; i<obs_name.size(); i++){
        if(obs_name[i]!="ptsub3_subpt") hists[i]->Fill(obs[i], weight);
      }
    }
    if(passed_subpt_rec && index_ptsub!=-1 && index_pt3!=-1) hists[index_ptsub]->Fill(obs[index_pt3], weight);
  }
  return hists;
}
