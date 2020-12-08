#include "../include/CentralInclude.h"

using namespace std;

TH1F* get_hist(TFile* file, TString weightname);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){
  TFile* f_16_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
  TFile* f_16_mu_fsrup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
  TFile* f_16_mu_fsrdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
  TFile* f_16_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
  TFile* f_16_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
  TFile* f_16_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
  TFile* f_16_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");

  TFile* f_16_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
  TFile* f_16_el_fsrup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
  TFile* f_16_el_fsrdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
  TFile* f_16_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
  TFile* f_16_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
  TFile* f_16_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
  TFile* f_16_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");

  TFile* f_17_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
  TFile* f_17_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
  TFile* f_17_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
  TFile* f_17_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
  TFile* f_17_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");

  TFile* f_17_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
  TFile* f_17_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
  TFile* f_17_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
  TFile* f_17_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
  TFile* f_17_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");


  TFile* f_18_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
  TFile* f_18_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
  TFile* f_18_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
  TFile* f_18_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2018.root");
  TFile* f_18_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");


  TFile* f_18_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
  TFile* f_18_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
  TFile* f_18_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
  TFile* f_18_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2018.root");
  TFile* f_18_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");

  // 2016
  TH1F* h_16_mu_nominal = get_hist(f_16_mu, "none");
  TH1F* h_16_mu_fsrup2 = get_hist(f_16_mu_fsrup, "none");
  TH1F* h_16_mu_fsrdown2 = get_hist(f_16_mu_fsrdown, "none");
  TH1F* h_16_mu_st = get_hist(f_16_mu_st, "none");
  TH1F* h_16_mu_wj = get_hist(f_16_mu_wj, "none");
  TH1F* h_16_mu_ot = get_hist(f_16_mu_ot, "none");
  TH1F* h_16_mu_data = get_hist(f_16_mu_data, "none");
  TH1F* h_16_el_nominal = get_hist(f_16_el, "none");
  TH1F* h_16_el_fsrup2 = get_hist(f_16_el_fsrup, "none");
  TH1F* h_16_el_fsrdown2 = get_hist(f_16_el_fsrdown, "none");
  TH1F* h_16_el_st = get_hist(f_16_el_st, "none");
  TH1F* h_16_el_wj = get_hist(f_16_el_wj, "none");
  TH1F* h_16_el_ot = get_hist(f_16_el_ot, "none");
  TH1F* h_16_el_data = get_hist(f_16_el_data, "none");
  // add channels
  TH1F* h_16_nominal = (TH1F*) h_16_mu_nominal->Clone();
  h_16_nominal->Add(h_16_el_nominal);
  TH1F* h_16_fsrup2 = (TH1F*) h_16_mu_fsrup2->Clone();
  h_16_fsrup2->Add(h_16_el_fsrup2);
  TH1F* h_16_fsrdown2 = (TH1F*) h_16_mu_fsrdown2->Clone();
  h_16_fsrdown2->Add(h_16_el_fsrdown2);
  TH1F* h_16_bgr = (TH1F*) h_16_mu_st->Clone();
  h_16_bgr->Add(h_16_el_st);
  h_16_bgr->Add(h_16_mu_wj);
  h_16_bgr->Add(h_16_el_wj);
  h_16_bgr->Add(h_16_mu_ot);
  h_16_bgr->Add(h_16_el_ot);
  TH1F* h_16_data = (TH1F*) h_16_mu_data->Clone();
  h_16_data->Add(h_16_el_data);

  // 2017
  TH1F* h_17_mu_nominal = get_hist(f_17_mu, "none");
  TH1F* h_17_mu_st = get_hist(f_17_mu_st, "none");
  TH1F* h_17_mu_wj = get_hist(f_17_mu_wj, "none");
  TH1F* h_17_mu_ot = get_hist(f_17_mu_ot, "none");
  TH1F* h_17_mu_data = get_hist(f_17_mu_data, "none");
  TH1F* h_17_mu_fsrupsqrt2 = get_hist(f_17_mu, "sf_fsr_upsqrt2");
  TH1F* h_17_mu_fsrup2 = get_hist(f_17_mu, "sf_fsr_up2");
  TH1F* h_17_mu_fsrup4 = get_hist(f_17_mu, "sf_fsr_up4");
  TH1F* h_17_mu_fsrdownsqrt2 = get_hist(f_17_mu, "sf_fsr_downsqrt2");
  TH1F* h_17_mu_fsrdown2 = get_hist(f_17_mu, "sf_fsr_down2");
  TH1F* h_17_mu_fsrdown4 = get_hist(f_17_mu, "sf_fsr_down4");
  TH1F* h_17_el_nominal = get_hist(f_17_el, "none");
  TH1F* h_17_el_st = get_hist(f_17_el_st, "none");
  TH1F* h_17_el_wj = get_hist(f_17_el_wj, "none");
  TH1F* h_17_el_ot = get_hist(f_17_el_ot, "none");
  TH1F* h_17_el_data = get_hist(f_17_el_data, "none");
  TH1F* h_17_el_fsrupsqrt2 = get_hist(f_17_el, "sf_fsr_upsqrt2");
  TH1F* h_17_el_fsrup2 = get_hist(f_17_el, "sf_fsr_up2");
  TH1F* h_17_el_fsrup4 = get_hist(f_17_el, "sf_fsr_up4");
  TH1F* h_17_el_fsrdownsqrt2 = get_hist(f_17_el, "sf_fsr_downsqrt2");
  TH1F* h_17_el_fsrdown2 = get_hist(f_17_el, "sf_fsr_down2");
  TH1F* h_17_el_fsrdown4 = get_hist(f_17_el, "sf_fsr_down4");
  // add channels
  TH1F* h_17_nominal = (TH1F*) h_17_mu_nominal->Clone();
  h_17_nominal->Add(h_17_el_nominal);
  TH1F* h_17_fsrupsqrt2 = (TH1F*) h_17_mu_fsrupsqrt2->Clone();
  h_17_fsrupsqrt2->Add(h_17_el_fsrupsqrt2);
  TH1F* h_17_fsrup2 = (TH1F*) h_17_mu_fsrup2->Clone();
  h_17_fsrup2->Add(h_17_el_fsrup2);
  TH1F* h_17_fsrup4 = (TH1F*) h_17_mu_fsrup4->Clone();
  h_17_fsrup4->Add(h_17_el_fsrup4);
  TH1F* h_17_fsrdownsqrt2 = (TH1F*) h_17_mu_fsrdownsqrt2->Clone();
  h_17_fsrdownsqrt2->Add(h_17_el_fsrdownsqrt2);
  TH1F* h_17_fsrdown2 = (TH1F*) h_17_mu_fsrdown2->Clone();
  h_17_fsrdown2->Add(h_17_el_fsrdown2);
  TH1F* h_17_fsrdown4 = (TH1F*) h_17_mu_fsrdown4->Clone();
  h_17_fsrdown4->Add(h_17_el_fsrdown4);
  TH1F* h_17_bgr = (TH1F*) h_17_mu_st->Clone();
  h_17_bgr->Add(h_17_el_st);
  h_17_bgr->Add(h_17_mu_wj);
  h_17_bgr->Add(h_17_el_wj);
  h_17_bgr->Add(h_17_mu_ot);
  h_17_bgr->Add(h_17_el_ot);
  TH1F* h_17_data = (TH1F*) h_17_mu_data->Clone();
  h_17_data->Add(h_17_el_data);

  // 2018
  TH1F* h_18_mu_nominal = get_hist(f_18_mu, "none");
  TH1F* h_18_mu_st = get_hist(f_18_mu_st, "none");
  TH1F* h_18_mu_wj = get_hist(f_18_mu_wj, "none");
  TH1F* h_18_mu_ot = get_hist(f_18_mu_ot, "none");
  TH1F* h_18_mu_data = get_hist(f_18_mu_data, "none");
  TH1F* h_18_mu_fsrupsqrt2 = get_hist(f_18_mu, "sf_fsr_upsqrt2");
  TH1F* h_18_mu_fsrup2 = get_hist(f_18_mu, "sf_fsr_up2");
  TH1F* h_18_mu_fsrup4 = get_hist(f_18_mu, "sf_fsr_up4");
  TH1F* h_18_mu_fsrdownsqrt2 = get_hist(f_18_mu, "sf_fsr_downsqrt2");
  TH1F* h_18_mu_fsrdown2 = get_hist(f_18_mu, "sf_fsr_down2");
  TH1F* h_18_mu_fsrdown4 = get_hist(f_18_mu, "sf_fsr_down4");
  TH1F* h_18_el_nominal = get_hist(f_18_el, "none");
  TH1F* h_18_el_st = get_hist(f_18_el_st, "none");
  TH1F* h_18_el_wj = get_hist(f_18_el_wj, "none");
  TH1F* h_18_el_ot = get_hist(f_18_el_ot, "none");
  TH1F* h_18_el_data = get_hist(f_18_el_data, "none");
  TH1F* h_18_el_fsrupsqrt2 = get_hist(f_18_el, "sf_fsr_upsqrt2");
  TH1F* h_18_el_fsrup2 = get_hist(f_18_el, "sf_fsr_up2");
  TH1F* h_18_el_fsrup4 = get_hist(f_18_el, "sf_fsr_up4");
  TH1F* h_18_el_fsrdownsqrt2 = get_hist(f_18_el, "sf_fsr_downsqrt2");
  TH1F* h_18_el_fsrdown2 = get_hist(f_18_el, "sf_fsr_down2");
  TH1F* h_18_el_fsrdown4 = get_hist(f_18_el, "sf_fsr_down4");
  // add channels
  TH1F* h_18_nominal = (TH1F*) h_18_mu_nominal->Clone();
  h_18_nominal->Add(h_18_el_nominal);
  TH1F* h_18_fsrupsqrt2 = (TH1F*) h_18_mu_fsrupsqrt2->Clone();
  h_18_fsrupsqrt2->Add(h_18_el_fsrupsqrt2);
  TH1F* h_18_fsrup2 = (TH1F*) h_18_mu_fsrup2->Clone();
  h_18_fsrup2->Add(h_18_el_fsrup2);
  TH1F* h_18_fsrup4 = (TH1F*) h_18_mu_fsrup4->Clone();
  h_18_fsrup4->Add(h_18_el_fsrup4);
  TH1F* h_18_fsrdownsqrt2 = (TH1F*) h_18_mu_fsrdownsqrt2->Clone();
  h_18_fsrdownsqrt2->Add(h_18_el_fsrdownsqrt2);
  TH1F* h_18_fsrdown2 = (TH1F*) h_18_mu_fsrdown2->Clone();
  h_18_fsrdown2->Add(h_18_el_fsrdown2);
  TH1F* h_18_fsrdown4 = (TH1F*) h_18_mu_fsrdown4->Clone();
  h_18_fsrdown4->Add(h_18_el_fsrdown4);
  TH1F* h_18_bgr = (TH1F*) h_18_mu_st->Clone();
  h_18_bgr->Add(h_18_el_st);
  h_18_bgr->Add(h_18_mu_wj);
  h_18_bgr->Add(h_18_el_wj);
  h_18_bgr->Add(h_18_mu_ot);
  h_18_bgr->Add(h_18_el_ot);
  TH1F* h_18_data = (TH1F*) h_18_mu_data->Clone();
  h_18_data->Add(h_18_el_data);

  TFile* f_out = new TFile("FSR_hists.root","recreate");
  f_out->cd();
  h_16_nominal->Write("nominal_16");
  h_16_fsrup2->Write("fsrup2_16");
  h_16_fsrdown2->Write("fsrdown2_16");
  h_16_bgr->Write("bgr_16");
  h_16_data->Write("data_16");

  h_17_nominal->Write("nominal_17");
  h_17_fsrupsqrt2->Write("fsrupsqrt2_17");
  h_17_fsrup2->Write("fsrup2_17");
  h_17_fsrup4->Write("fsrup4_17");
  h_17_fsrdownsqrt2->Write("fsrdownsqrt2_17");
  h_17_fsrdown2->Write("fsrdown2_17");
  h_17_fsrdown4->Write("fsrdown4_17");
  h_17_bgr->Write("bgr_17");
  h_17_data->Write("data_17");

  h_18_nominal->Write("nominal_18");
  h_18_fsrupsqrt2->Write("fsrupsqrt2_18");
  h_18_fsrup2->Write("fsrup2_18");
  h_18_fsrup4->Write("fsrup4_18");
  h_18_fsrdownsqrt2->Write("fsrdownsqrt2_18");
  h_18_fsrdown2->Write("fsrdown2_18");
  h_18_fsrdown4->Write("fsrdown4_18");
  h_18_bgr->Write("bgr_18");
  h_18_data->Write("data_18");

  f_out->Close();

  return 0;
}

TH1F* get_hist(TFile* file, TString weightname){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  TH1F* h_hist = new TH1F("hist", "hist", 10, 0, 1);

  Float_t tau32;
  Bool_t passed_measurement_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("tau32_ak8_had",&tau32);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  if(weightname != "none"){
    cout << "  - using weight" << endl;
    tree->SetBranchAddress(weightname,&additional_factor);
  }
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;

    // Here, only rec weight is changed!
    if(weightname != "none") rec_weight *= additional_factor;


    // get weights for migration matrix
    weight = rec_weight * gen_weight;
    if(tau32 < 0) continue;
    else h_hist->Fill(tau32, weight);
  }
  return h_hist;
}
