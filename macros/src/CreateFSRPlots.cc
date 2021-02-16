#include "../include/CentralInclude.h"

using namespace std;

TH1F* get_hist(TFile* file, vector<double> mcut, TString weightname);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  vector<TString> masswindows = {"140toInf", "140to200", "140to160", "160to170", "170to180", "180to190", "190to200", "200toInf"};
  // vector<TString> masswindows = {"140to200"};

  vector<vector<double>> mcutvalues = {
    {140., 10000.},
    {140., 200.},
    {140., 160.},
    {160., 170.},
    {170., 180.},
    {180., 190.},
    {190., 200.},
    {200., 10000.}
  };
  // vector<vector<double>> mcutvalues = { {140., 200.} };

  for(unsigned int i=0; i<masswindows.size(); i++){

    TString mass_string = masswindows[i];
    vector<double> mcut = mcutvalues[i];
    cout << "mjet: " << mass_string << endl;


    TFile* f_16_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_fsrup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
    TFile* f_16_mu_fsrdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
    TFile* f_16_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
    TFile* f_16_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
    TFile* f_16_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
    TFile* f_16_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
    TFile* f_16_mu_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");

    TFile* f_16_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_fsrup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
    TFile* f_16_el_fsrdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
    TFile* f_16_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
    TFile* f_16_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
    TFile* f_16_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
    TFile* f_16_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
    TFile* f_16_el_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");

    TFile* f_17_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
    TFile* f_17_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
    TFile* f_17_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
    TFile* f_17_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
    TFile* f_17_mu_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");

    TFile* f_17_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
    TFile* f_17_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
    TFile* f_17_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
    TFile* f_17_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
    TFile* f_17_el_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");

    TFile* f_18_mu = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
    TFile* f_18_mu_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
    TFile* f_18_mu_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.other_2018.root");
    TFile* f_18_mu_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
    TFile* f_18_mu_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");

    TFile* f_18_el = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_st = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
    TFile* f_18_el_wj = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
    TFile* f_18_el_ot = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.other_2018.root");
    TFile* f_18_el_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
    TFile* f_18_el_jecup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jecdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_corup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_cordown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jmsup = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jmsdown = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");


    // 2016
    TH1F* h_16_mu_nominal = get_hist(f_16_mu, mcut, "none");
    TH1F* h_16_mu_fsrup2 = get_hist(f_16_mu_fsrup, mcut, "none");
    TH1F* h_16_mu_fsrdown2 = get_hist(f_16_mu_fsrdown, mcut, "none");
    TH1F* h_16_mu_st = get_hist(f_16_mu_st, mcut, "none");
    TH1F* h_16_mu_wj = get_hist(f_16_mu_wj, mcut, "none");
    TH1F* h_16_mu_ot = get_hist(f_16_mu_ot, mcut, "none");
    TH1F* h_16_mu_data = get_hist(f_16_mu_data, mcut, "none");
    TH1F* h_16_mu_jecup = get_hist(f_16_mu_jecup, mcut, "none");
    TH1F* h_16_mu_jecdown = get_hist(f_16_mu_jecdown, mcut, "none");
    TH1F* h_16_mu_corup = get_hist(f_16_mu_corup, mcut, "none");
    TH1F* h_16_mu_cordown = get_hist(f_16_mu_cordown, mcut, "none");
    TH1F* h_16_mu_jmsup = get_hist(f_16_mu_jmsup, mcut, "none");
    TH1F* h_16_mu_jmsdown = get_hist(f_16_mu_jmsdown, mcut, "none");
    TH1F* h_16_el_nominal = get_hist(f_16_el, mcut, "none");
    TH1F* h_16_el_fsrup2 = get_hist(f_16_el_fsrup, mcut, "none");
    TH1F* h_16_el_fsrdown2 = get_hist(f_16_el_fsrdown, mcut, "none");
    TH1F* h_16_el_st = get_hist(f_16_el_st, mcut, "none");
    TH1F* h_16_el_wj = get_hist(f_16_el_wj, mcut, "none");
    TH1F* h_16_el_ot = get_hist(f_16_el_ot, mcut, "none");
    TH1F* h_16_el_data = get_hist(f_16_el_data, mcut, "none");
    TH1F* h_16_el_jecup = get_hist(f_16_el_jecup, mcut, "none");
    TH1F* h_16_el_jecdown = get_hist(f_16_el_jecdown, mcut, "none");
    TH1F* h_16_el_corup = get_hist(f_16_el_corup, mcut, "none");
    TH1F* h_16_el_cordown = get_hist(f_16_el_cordown, mcut, "none");
    TH1F* h_16_el_jmsup = get_hist(f_16_el_jmsup, mcut, "none");
    TH1F* h_16_el_jmsdown = get_hist(f_16_el_jmsdown, mcut, "none");
    // add channels
    TH1F* h_16_nominal = (TH1F*) h_16_mu_nominal->Clone();
    h_16_nominal->Add(h_16_el_nominal);
    TH1F* h_16_fsrup2 = (TH1F*) h_16_mu_fsrup2->Clone();
    h_16_fsrup2->Add(h_16_el_fsrup2);
    TH1F* h_16_fsrdown2 = (TH1F*) h_16_mu_fsrdown2->Clone();
    h_16_fsrdown2->Add(h_16_el_fsrdown2);

    TH1F* h_16_jecup = (TH1F*) h_16_mu_jecup->Clone();
    h_16_jecup->Add(h_16_el_jecup);
    TH1F* h_16_jecdown = (TH1F*) h_16_mu_jecdown->Clone();
    h_16_jecdown->Add(h_16_el_jecdown);

    TH1F* h_16_corup = (TH1F*) h_16_mu_corup->Clone();
    h_16_corup->Add(h_16_el_corup);
    TH1F* h_16_cordown = (TH1F*) h_16_mu_cordown->Clone();
    h_16_cordown->Add(h_16_el_cordown);

    TH1F* h_16_jmsup = (TH1F*) h_16_mu_jmsup->Clone();
    h_16_jmsup->Add(h_16_el_jmsup);
    TH1F* h_16_jmsdown = (TH1F*) h_16_mu_jmsdown->Clone();
    h_16_jmsdown->Add(h_16_el_jmsdown);

    TH1F* h_16_st = (TH1F*) h_16_mu_st->Clone();
    h_16_st->Add(h_16_el_st);
    TH1F* h_16_wj = (TH1F*) h_16_mu_wj->Clone();
    h_16_wj->Add(h_16_el_wj);
    TH1F* h_16_ot = (TH1F*) h_16_mu_ot->Clone();
    h_16_ot->Add(h_16_el_ot);
    TH1F* h_16_data = (TH1F*) h_16_mu_data->Clone();
    h_16_data->Add(h_16_el_data);

    // 2017
    TH1F* h_17_mu_nominal = get_hist(f_17_mu, mcut, "none");
    TH1F* h_17_mu_st = get_hist(f_17_mu_st, mcut, "none");
    TH1F* h_17_mu_wj = get_hist(f_17_mu_wj, mcut, "none");
    TH1F* h_17_mu_ot = get_hist(f_17_mu_ot, mcut, "none");
    TH1F* h_17_mu_data = get_hist(f_17_mu_data, mcut, "none");
    TH1F* h_17_mu_jecup = get_hist(f_17_mu_jecup, mcut, "none");
    TH1F* h_17_mu_jecdown = get_hist(f_17_mu_jecdown, mcut, "none");
    TH1F* h_17_mu_corup = get_hist(f_17_mu_corup, mcut, "none");
    TH1F* h_17_mu_cordown = get_hist(f_17_mu_cordown, mcut, "none");
    TH1F* h_17_mu_jmsup = get_hist(f_17_mu_jmsup, mcut, "none");
    TH1F* h_17_mu_jmsdown = get_hist(f_17_mu_jmsdown, mcut, "none");
    TH1F* h_17_mu_fsrupsqrt2 = get_hist(f_17_mu, mcut, "sf_fsr_upsqrt2");
    TH1F* h_17_mu_fsrup2 = get_hist(f_17_mu, mcut, "sf_fsr_up2");
    TH1F* h_17_mu_fsrup4 = get_hist(f_17_mu, mcut, "sf_fsr_up4");
    TH1F* h_17_mu_fsrdownsqrt2 = get_hist(f_17_mu, mcut, "sf_fsr_downsqrt2");
    TH1F* h_17_mu_fsrdown2 = get_hist(f_17_mu, mcut, "sf_fsr_down2");
    TH1F* h_17_mu_fsrdown4 = get_hist(f_17_mu, mcut, "sf_fsr_down4");
    TH1F* h_17_el_nominal = get_hist(f_17_el, mcut, "none");
    TH1F* h_17_el_st = get_hist(f_17_el_st, mcut, "none");
    TH1F* h_17_el_wj = get_hist(f_17_el_wj, mcut, "none");
    TH1F* h_17_el_ot = get_hist(f_17_el_ot, mcut, "none");
    TH1F* h_17_el_data = get_hist(f_17_el_data, mcut, "none");
    TH1F* h_17_el_fsrupsqrt2 = get_hist(f_17_el, mcut, "sf_fsr_upsqrt2");
    TH1F* h_17_el_fsrup2 = get_hist(f_17_el, mcut, "sf_fsr_up2");
    TH1F* h_17_el_fsrup4 = get_hist(f_17_el, mcut, "sf_fsr_up4");
    TH1F* h_17_el_fsrdownsqrt2 = get_hist(f_17_el, mcut, "sf_fsr_downsqrt2");
    TH1F* h_17_el_fsrdown2 = get_hist(f_17_el, mcut, "sf_fsr_down2");
    TH1F* h_17_el_fsrdown4 = get_hist(f_17_el, mcut, "sf_fsr_down4");
    TH1F* h_17_el_jecup = get_hist(f_17_el_jecup, mcut, "none");
    TH1F* h_17_el_jecdown = get_hist(f_17_el_jecdown, mcut, "none");
    TH1F* h_17_el_corup = get_hist(f_17_el_corup, mcut, "none");
    TH1F* h_17_el_cordown = get_hist(f_17_el_cordown, mcut, "none");
    TH1F* h_17_el_jmsup = get_hist(f_17_el_jmsup, mcut, "none");
    TH1F* h_17_el_jmsdown = get_hist(f_17_el_jmsdown, mcut, "none");
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
    TH1F* h_17_jecup = (TH1F*) h_17_mu_jecup->Clone();
    h_17_jecup->Add(h_17_el_jecup);
    TH1F* h_17_jecdown = (TH1F*) h_17_mu_jecdown->Clone();
    h_17_jecdown->Add(h_17_el_jecdown);
    TH1F* h_17_corup = (TH1F*) h_17_mu_corup->Clone();
    h_17_corup->Add(h_17_el_corup);
    TH1F* h_17_cordown = (TH1F*) h_17_mu_cordown->Clone();
    h_17_cordown->Add(h_17_el_cordown);
    TH1F* h_17_jmsup = (TH1F*) h_17_mu_jmsup->Clone();
    h_17_jmsup->Add(h_17_el_jmsup);
    TH1F* h_17_jmsdown = (TH1F*) h_17_mu_jmsdown->Clone();
    h_17_jmsdown->Add(h_17_el_jmsdown);
    TH1F* h_17_st = (TH1F*) h_17_mu_st->Clone();
    h_17_st->Add(h_17_el_st);
    TH1F* h_17_wj = (TH1F*) h_17_mu_wj->Clone();
    h_17_wj->Add(h_17_el_wj);
    TH1F* h_17_ot = (TH1F*) h_17_mu_ot->Clone();
    h_17_ot->Add(h_17_el_ot);
    TH1F* h_17_data = (TH1F*) h_17_mu_data->Clone();
    h_17_data->Add(h_17_el_data);

    // 2018
    TH1F* h_18_mu_nominal = get_hist(f_18_mu, mcut, "none");
    TH1F* h_18_mu_st = get_hist(f_18_mu_st, mcut, "none");
    TH1F* h_18_mu_wj = get_hist(f_18_mu_wj, mcut, "none");
    TH1F* h_18_mu_ot = get_hist(f_18_mu_ot, mcut, "none");
    TH1F* h_18_mu_data = get_hist(f_18_mu_data, mcut, "none");
    TH1F* h_18_mu_jecup = get_hist(f_18_mu_jecup, mcut, "none");
    TH1F* h_18_mu_jecdown = get_hist(f_18_mu_jecdown, mcut, "none");
    TH1F* h_18_mu_corup = get_hist(f_18_mu_corup, mcut, "none");
    TH1F* h_18_mu_cordown = get_hist(f_18_mu_cordown, mcut, "none");
    TH1F* h_18_mu_jmsup = get_hist(f_18_mu_jmsup, mcut, "none");
    TH1F* h_18_mu_jmsdown = get_hist(f_18_mu_jmsdown, mcut, "none");
    TH1F* h_18_mu_fsrupsqrt2 = get_hist(f_18_mu, mcut, "sf_fsr_upsqrt2");
    TH1F* h_18_mu_fsrup2 = get_hist(f_18_mu, mcut, "sf_fsr_up2");
    TH1F* h_18_mu_fsrup4 = get_hist(f_18_mu, mcut, "sf_fsr_up4");
    TH1F* h_18_mu_fsrdownsqrt2 = get_hist(f_18_mu, mcut, "sf_fsr_downsqrt2");
    TH1F* h_18_mu_fsrdown2 = get_hist(f_18_mu, mcut, "sf_fsr_down2");
    TH1F* h_18_mu_fsrdown4 = get_hist(f_18_mu, mcut, "sf_fsr_down4");
    TH1F* h_18_el_nominal = get_hist(f_18_el, mcut, "none");
    TH1F* h_18_el_st = get_hist(f_18_el_st, mcut, "none");
    TH1F* h_18_el_wj = get_hist(f_18_el_wj, mcut, "none");
    TH1F* h_18_el_ot = get_hist(f_18_el_ot, mcut, "none");
    TH1F* h_18_el_data = get_hist(f_18_el_data, mcut, "none");
    TH1F* h_18_el_fsrupsqrt2 = get_hist(f_18_el, mcut, "sf_fsr_upsqrt2");
    TH1F* h_18_el_fsrup2 = get_hist(f_18_el, mcut, "sf_fsr_up2");
    TH1F* h_18_el_fsrup4 = get_hist(f_18_el, mcut, "sf_fsr_up4");
    TH1F* h_18_el_fsrdownsqrt2 = get_hist(f_18_el, mcut, "sf_fsr_downsqrt2");
    TH1F* h_18_el_fsrdown2 = get_hist(f_18_el, mcut, "sf_fsr_down2");
    TH1F* h_18_el_fsrdown4 = get_hist(f_18_el, mcut, "sf_fsr_down4");
    TH1F* h_18_el_jecup = get_hist(f_18_el_jecup, mcut, "none");
    TH1F* h_18_el_jecdown = get_hist(f_18_el_jecdown, mcut, "none");
    TH1F* h_18_el_corup = get_hist(f_18_el_corup, mcut, "none");
    TH1F* h_18_el_cordown = get_hist(f_18_el_cordown, mcut, "none");
    TH1F* h_18_el_jmsup = get_hist(f_18_el_jmsup, mcut, "none");
    TH1F* h_18_el_jmsdown = get_hist(f_18_el_jmsdown, mcut, "none");
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
    TH1F* h_18_jecup = (TH1F*) h_18_mu_jecup->Clone();
    h_18_jecup->Add(h_18_el_jecup);
    TH1F* h_18_jecdown = (TH1F*) h_18_mu_jecdown->Clone();
    h_18_jecdown->Add(h_18_el_jecdown);
    TH1F* h_18_corup = (TH1F*) h_18_mu_corup->Clone();
    h_18_corup->Add(h_18_el_corup);
    TH1F* h_18_cordown = (TH1F*) h_18_mu_cordown->Clone();
    h_18_cordown->Add(h_18_el_cordown);
    TH1F* h_18_jmsup = (TH1F*) h_18_mu_jmsup->Clone();
    h_18_jmsup->Add(h_18_el_jmsup);
    TH1F* h_18_jmsdown = (TH1F*) h_18_mu_jmsdown->Clone();
    h_18_jmsdown->Add(h_18_el_jmsdown);
    TH1F* h_18_st = (TH1F*) h_18_mu_st->Clone();
    h_18_st->Add(h_18_el_st);
    TH1F* h_18_wj = (TH1F*) h_18_mu_wj->Clone();
    h_18_wj->Add(h_18_el_wj);
    TH1F* h_18_ot = (TH1F*) h_18_mu_ot->Clone();
    h_18_ot->Add(h_18_el_ot);
    TH1F* h_18_data = (TH1F*) h_18_mu_data->Clone();
    h_18_data->Add(h_18_el_data);

    // combine 17 + 18
    TH1F* h_1718_nominal = (TH1F*) h_18_nominal->Clone();
    h_1718_nominal->Add(h_17_nominal);
    TH1F* h_1718_fsrupsqrt2 = (TH1F*) h_18_fsrupsqrt2->Clone();
    h_1718_fsrupsqrt2->Add(h_17_fsrupsqrt2);
    TH1F* h_1718_fsrup2 = (TH1F*) h_18_fsrup2->Clone();
    h_1718_fsrup2->Add(h_17_fsrup2);
    TH1F* h_1718_fsrup4 = (TH1F*) h_18_fsrup4->Clone();
    h_1718_fsrup4->Add(h_17_fsrup4);
    TH1F* h_1718_fsrdownsqrt2 = (TH1F*) h_18_fsrdownsqrt2->Clone();
    h_1718_fsrdownsqrt2->Add(h_17_fsrdownsqrt2);
    TH1F* h_1718_fsrdown2 = (TH1F*) h_18_fsrdown2->Clone();
    h_1718_fsrdown2->Add(h_17_fsrdown2);
    TH1F* h_1718_fsrdown4 = (TH1F*) h_18_fsrdown4->Clone();
    h_1718_fsrdown4->Add(h_17_fsrdown4);
    TH1F* h_1718_jecup = (TH1F*) h_18_jecup->Clone();
    h_1718_jecup->Add(h_17_jecup);
    TH1F* h_1718_jecdown = (TH1F*) h_18_jecdown->Clone();
    h_1718_jecdown->Add(h_17_jecdown);
    TH1F* h_1718_corup = (TH1F*) h_18_corup->Clone();
    h_1718_corup->Add(h_17_corup);
    TH1F* h_1718_cordown = (TH1F*) h_18_cordown->Clone();
    h_1718_cordown->Add(h_17_cordown);
    TH1F* h_1718_jmsup = (TH1F*) h_18_jmsup->Clone();
    h_1718_jmsup->Add(h_17_jmsup);
    TH1F* h_1718_jmsdown = (TH1F*) h_18_jmsdown->Clone();
    h_1718_jmsdown->Add(h_17_jmsdown);
    TH1F* h_1718_st = (TH1F*) h_18_st->Clone();
    h_1718_st->Add(h_17_st);
    TH1F* h_1718_wj = (TH1F*) h_18_wj->Clone();
    h_1718_wj->Add(h_17_wj);
    TH1F* h_1718_ot = (TH1F*) h_18_ot->Clone();
    h_1718_ot->Add(h_17_ot);
    TH1F* h_1718_data = (TH1F*) h_18_data->Clone();
    h_1718_data->Add(h_17_data);

    TFile* f_out = new TFile("FSR_hists_mjet"+mass_string+".root","recreate");
    f_out->cd();
    h_16_nominal->Write("nominal_2016");
    h_16_fsrup2->Write("fsrup2_2016");
    h_16_fsrdown2->Write("fsrdown2_2016");
    h_16_jecup->Write("jecup_2016");
    h_16_jecdown->Write("jecdown_2016");
    h_16_corup->Write("corup_2016");
    h_16_cordown->Write("cordown_2016");
    h_16_jmsup->Write("jmsup_2016");
    h_16_jmsdown->Write("jmsdown_2016");
    h_16_st->Write("bgr_st_2016");
    h_16_wj->Write("bgr_wj_2016");
    h_16_ot->Write("bgr_ot_2016");
    h_16_data->Write("data_2016");

    h_17_nominal->Write("nominal_2017");
    h_17_fsrupsqrt2->Write("fsrupsqrt2_2017");
    h_17_fsrup2->Write("fsrup2_2017");
    h_17_fsrup4->Write("fsrup4_2017");
    h_17_fsrdownsqrt2->Write("fsrdownsqrt2_2017");
    h_17_fsrdown2->Write("fsrdown2_2017");
    h_17_fsrdown4->Write("fsrdown4_2017");
    h_17_jecup->Write("jecup_2017");
    h_17_jecdown->Write("jecdown_2017");
    h_17_corup->Write("corup_2017");
    h_17_cordown->Write("cordown_2017");
    h_17_jmsup->Write("jmsup_2017");
    h_17_jmsdown->Write("jmsdown_2017");
    h_17_st->Write("bgr_st_2017");
    h_17_wj->Write("bgr_wj_2017");
    h_17_ot->Write("bgr_ot_2017");
    h_17_data->Write("data_2017");

    h_18_nominal->Write("nominal_2018");
    h_18_fsrupsqrt2->Write("fsrupsqrt2_2018");
    h_18_fsrup2->Write("fsrup2_2018");
    h_18_fsrup4->Write("fsrup4_2018");
    h_18_fsrdownsqrt2->Write("fsrdownsqrt2_2018");
    h_18_fsrdown2->Write("fsrdown2_2018");
    h_18_fsrdown4->Write("fsrdown4_2018");
    h_18_jecup->Write("jecup_2018");
    h_18_jecdown->Write("jecdown_2018");
    h_18_corup->Write("corup_2018");
    h_18_cordown->Write("cordown_2018");
    h_18_jmsup->Write("jmsup_2018");
    h_18_jmsdown->Write("jmsdown_2018");
    h_18_st->Write("bgr_st_2018");
    h_18_wj->Write("bgr_wj_2018");
    h_18_ot->Write("bgr_ot_2018");
    h_18_data->Write("data_2018");

    h_1718_nominal->Write("nominal_combine");
    h_1718_fsrupsqrt2->Write("fsrupsqrt2_combine");
    h_1718_fsrup2->Write("fsrup2_combine");
    h_1718_fsrup4->Write("fsrup4_combine");
    h_1718_fsrdownsqrt2->Write("fsrdownsqrt2_combine");
    h_1718_fsrdown2->Write("fsrdown2_combine");
    h_1718_fsrdown4->Write("fsrdown4_combine");
    h_1718_jecup->Write("jecup_combine");
    h_1718_jecdown->Write("jecdown_combine");
    h_1718_corup->Write("corup_combine");
    h_1718_cordown->Write("cordown_combine");
    h_1718_jmsup->Write("jmsup_combine");
    h_1718_jmsdown->Write("jmsdown_combine");
    h_1718_st->Write("bgr_st_combine");
    h_1718_wj->Write("bgr_wj_combine");
    h_1718_ot->Write("bgr_ot_combine");
    h_1718_data->Write("data_combine");

    f_out->Close();
  }
  return 0;
}

TH1F* get_hist(TFile* file, vector<double> mcut, TString weightname){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  TH1F* h_hist = new TH1F("hist", "hist", 100, 0, 1);

  Float_t tau32, mjet;
  Bool_t passed_measurement_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("tau32_ak8_had",&tau32);
  tree->SetBranchAddress("mass_ak8_had",&mjet);
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

    // if(weight > 50){
    //   cout << "  - Weight too large, skip event" << "(weight = " << weight << ")"<< endl;
    //   continue;
    // }

    if(tau32 < 0 || mjet < mcut[0] || mjet > mcut[1]) continue;
    else h_hist->Fill(tau32, weight);
  }
  return h_hist;
}
