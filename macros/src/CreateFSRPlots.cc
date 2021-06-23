#include "../include/CentralInclude.h"
#include "../include/Utils.h"

using namespace std;

TH1F* get_hist(TFile* file, vector<double> mcut, TString weightname);
TH1F* AddHists(TH1F* h1, TH1F* h2);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  TH1::AddDirectory(kFALSE);

  // vector<TString> masswindows = {"140toInf", "140to200", "140to160", "160to170", "170to180", "180to190", "190to200", "200toInf"};
  vector<TString> masswindows = {"140toInf"};

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
    auto start = high_resolution_clock::now(); // Calculation time - start

    TString mass_string = masswindows[i];
    vector<double> mcut = mcutvalues[i];
    cout << "mjet: " << mass_string << endl;

    // TString dir = "/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel";
    TString dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel";
    TFile* f_16_mu = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_fsrup = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
    TFile* f_16_mu_fsrdown = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
    TFile* f_16_mu_st = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
    TFile* f_16_mu_wj = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
    TFile* f_16_mu_ot = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
    TFile* f_16_mu_data = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
    TFile* f_16_mu_jecup = new TFile(dir+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jecdown = new TFile(dir+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_corup = new TFile(dir+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_cordown = new TFile(dir+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jmsup = new TFile(dir+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_mu_jmsdown = new TFile(dir+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");

    TFile* f_16_el = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_fsrup = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
    TFile* f_16_el_fsrdown = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
    TFile* f_16_el_st = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
    TFile* f_16_el_wj = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
    TFile* f_16_el_ot = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
    TFile* f_16_el_data = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
    TFile* f_16_el_jecup = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jecdown = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_corup = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_cordown = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jmsup = new TFile(dir+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
    TFile* f_16_el_jmsdown = new TFile(dir+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");

    TFile* f_17_mu = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_st = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
    TFile* f_17_mu_wj = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
    TFile* f_17_mu_ot = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
    TFile* f_17_mu_data = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
    TFile* f_17_mu_jecup = new TFile(dir+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jecdown = new TFile(dir+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_corup = new TFile(dir+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_cordown = new TFile(dir+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jmsup = new TFile(dir+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_mu_jmsdown = new TFile(dir+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");

    TFile* f_17_el = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_st = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
    TFile* f_17_el_wj = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
    TFile* f_17_el_ot = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
    TFile* f_17_el_data = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
    TFile* f_17_el_jecup = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jecdown = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_corup = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_cordown = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jmsup = new TFile(dir+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
    TFile* f_17_el_jmsdown = new TFile(dir+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");

    TFile* f_18_mu = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_st = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
    TFile* f_18_mu_wj = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
    TFile* f_18_mu_ot = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.MC.other_2018.root");
    TFile* f_18_mu_data = new TFile(dir+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
    TFile* f_18_mu_jecup = new TFile(dir+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jecdown = new TFile(dir+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_corup = new TFile(dir+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_cordown = new TFile(dir+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jmsup = new TFile(dir+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_mu_jmsdown = new TFile(dir+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");

    TFile* f_18_el = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_st = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
    TFile* f_18_el_wj = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
    TFile* f_18_el_ot = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.MC.other_2018.root");
    TFile* f_18_el_data = new TFile(dir+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
    TFile* f_18_el_jecup = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jecdown = new TFile(dir+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_corup = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_cordown = new TFile(dir+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jmsup = new TFile(dir+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
    TFile* f_18_el_jmsdown = new TFile(dir+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");


    // 2016
    TH1F* h_16_mu_nominal = get_hist(f_16_mu, mcut, "none");
    TH1F* h_16_mu_fsrup2 = get_hist(f_16_mu_fsrup, mcut, "none");
    TH1F* h_16_mu_fsrdown2 = get_hist(f_16_mu_fsrdown, mcut, "none");
    TH1F* h_16_mu_st = get_hist(f_16_mu_st, mcut, "none");
    TH1F* h_16_mu_wj = get_hist(f_16_mu_wj, mcut, "none");
    TH1F* h_16_mu_ot = get_hist(f_16_mu_ot, mcut, "none");
    TH1F* h_16_mu_data = get_hist(f_16_mu_data, mcut, "none");
    // TH1F* h_16_mu_jecup = get_hist(f_16_mu_jecup, mcut, "none");
    // TH1F* h_16_mu_jecdown = get_hist(f_16_mu_jecdown, mcut, "none");
    // TH1F* h_16_mu_corup = get_hist(f_16_mu_corup, mcut, "none");
    // TH1F* h_16_mu_cordown = get_hist(f_16_mu_cordown, mcut, "none");
    TH1F* h_16_mu_jmsup = get_hist(f_16_mu_jmsup, mcut, "none");
    TH1F* h_16_mu_jmsdown = get_hist(f_16_mu_jmsdown, mcut, "none");
    TH1F* h_16_el_nominal = get_hist(f_16_el, mcut, "none");
    TH1F* h_16_el_fsrup2 = get_hist(f_16_el_fsrup, mcut, "none");
    TH1F* h_16_el_fsrdown2 = get_hist(f_16_el_fsrdown, mcut, "none");
    TH1F* h_16_el_st = get_hist(f_16_el_st, mcut, "none");
    TH1F* h_16_el_wj = get_hist(f_16_el_wj, mcut, "none");
    TH1F* h_16_el_ot = get_hist(f_16_el_ot, mcut, "none");
    TH1F* h_16_el_data = get_hist(f_16_el_data, mcut, "none");
    // TH1F* h_16_el_jecup = get_hist(f_16_el_jecup, mcut, "none");
    // TH1F* h_16_el_jecdown = get_hist(f_16_el_jecdown, mcut, "none");
    // TH1F* h_16_el_corup = get_hist(f_16_el_corup, mcut, "none");
    // TH1F* h_16_el_cordown = get_hist(f_16_el_cordown, mcut, "none");
    TH1F* h_16_el_jmsup = get_hist(f_16_el_jmsup, mcut, "none");
    TH1F* h_16_el_jmsdown = get_hist(f_16_el_jmsdown, mcut, "none");
    // add channels
    TH1F* h_16_nominal = AddHists(h_16_mu_nominal, h_16_el_nominal);
    TH1F* h_16_fsrup2 = AddHists(h_16_mu_fsrup2, h_16_el_fsrup2);
    TH1F* h_16_fsrdown2 = AddHists(h_16_mu_fsrdown2, h_16_el_fsrdown2);
    // TH1F* h_16_jecup = AddHists(h_16_mu_jecup, h_16_el_jecup);
    // TH1F* h_16_jecdown = AddHists(h_16_mu_jecdown, h_16_el_jecdown);
    // TH1F* h_16_corup = AddHists(h_16_mu_corup, h_16_el_corup);
    // TH1F* h_16_cordown = AddHists(h_16_mu_cordown, h_16_el_cordown);
    TH1F* h_16_jmsup = AddHists(h_16_mu_jmsup, h_16_el_jmsup);
    TH1F* h_16_jmsdown = AddHists(h_16_mu_jmsdown, h_16_el_jmsdown);
    TH1F* h_16_st = AddHists(h_16_mu_st, h_16_el_st);
    TH1F* h_16_wj = AddHists(h_16_mu_wj, h_16_el_wj);
    TH1F* h_16_ot = AddHists(h_16_mu_ot, h_16_el_ot);
    TH1F* h_16_data = AddHists(h_16_mu_data, h_16_el_data);

    // 2017
    TH1F* h_17_mu_nominal = get_hist(f_17_mu, mcut, "none");
    TH1F* h_17_mu_st = get_hist(f_17_mu_st, mcut, "none");
    TH1F* h_17_mu_wj = get_hist(f_17_mu_wj, mcut, "none");
    TH1F* h_17_mu_ot = get_hist(f_17_mu_ot, mcut, "none");
    TH1F* h_17_mu_data = get_hist(f_17_mu_data, mcut, "none");
    // TH1F* h_17_mu_jecup = get_hist(f_17_mu_jecup, mcut, "none");
    // TH1F* h_17_mu_jecdown = get_hist(f_17_mu_jecdown, mcut, "none");
    // TH1F* h_17_mu_corup = get_hist(f_17_mu_corup, mcut, "none");
    // TH1F* h_17_mu_cordown = get_hist(f_17_mu_cordown, mcut, "none");
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
    // TH1F* h_17_el_jecup = get_hist(f_17_el_jecup, mcut, "none");
    // TH1F* h_17_el_jecdown = get_hist(f_17_el_jecdown, mcut, "none");
    // TH1F* h_17_el_corup = get_hist(f_17_el_corup, mcut, "none");
    // TH1F* h_17_el_cordown = get_hist(f_17_el_cordown, mcut, "none");
    TH1F* h_17_el_jmsup = get_hist(f_17_el_jmsup, mcut, "none");
    TH1F* h_17_el_jmsdown = get_hist(f_17_el_jmsdown, mcut, "none");
    // add channels
    TH1F* h_17_nominal = AddHists(h_17_mu_nominal, h_17_el_nominal);
    TH1F* h_17_fsrupsqrt2 = AddHists(h_17_mu_fsrupsqrt2, h_17_el_fsrupsqrt2);
    TH1F* h_17_fsrdownsqrt2 = AddHists(h_17_mu_fsrdownsqrt2, h_17_el_fsrdownsqrt2);
    TH1F* h_17_fsrup2 = AddHists(h_17_mu_fsrup2, h_17_el_fsrup2);
    TH1F* h_17_fsrdown2 = AddHists(h_17_mu_fsrdown2, h_17_el_fsrdown2);
    TH1F* h_17_fsrup4 = AddHists(h_17_mu_fsrup4, h_17_el_fsrup4);
    TH1F* h_17_fsrdown4 = AddHists(h_17_mu_fsrdown4, h_17_el_fsrdown4);
    // TH1F* h_17_jecup = AddHists(h_17_mu_jecup, h_17_el_jecup);
    // TH1F* h_17_jecdown = AddHists(h_17_mu_jecdown, h_17_el_jecdown);
    // TH1F* h_17_corup = AddHists(h_17_mu_corup, h_17_el_corup);
    // TH1F* h_17_cordown = AddHists(h_17_mu_cordown, h_17_el_cordown);
    TH1F* h_17_jmsup = AddHists(h_17_mu_jmsup, h_17_el_jmsup);
    TH1F* h_17_jmsdown = AddHists(h_17_mu_jmsdown, h_17_el_jmsdown);
    TH1F* h_17_st = AddHists(h_17_mu_st, h_17_el_st);
    TH1F* h_17_wj = AddHists(h_17_mu_wj, h_17_el_wj);
    TH1F* h_17_ot = AddHists(h_17_mu_ot, h_17_el_ot);
    TH1F* h_17_data = AddHists(h_17_mu_data, h_17_el_data);

    // 2018
    TH1F* h_18_mu_nominal = get_hist(f_18_mu, mcut, "none");
    TH1F* h_18_mu_st = get_hist(f_18_mu_st, mcut, "none");
    TH1F* h_18_mu_wj = get_hist(f_18_mu_wj, mcut, "none");
    TH1F* h_18_mu_ot = get_hist(f_18_mu_ot, mcut, "none");
    TH1F* h_18_mu_data = get_hist(f_18_mu_data, mcut, "none");
    // TH1F* h_18_mu_jecup = get_hist(f_18_mu_jecup, mcut, "none");
    // TH1F* h_18_mu_jecdown = get_hist(f_18_mu_jecdown, mcut, "none");
    // TH1F* h_18_mu_corup = get_hist(f_18_mu_corup, mcut, "none");
    // TH1F* h_18_mu_cordown = get_hist(f_18_mu_cordown, mcut, "none");
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
    // TH1F* h_18_el_jecup = get_hist(f_18_el_jecup, mcut, "none");
    // TH1F* h_18_el_jecdown = get_hist(f_18_el_jecdown, mcut, "none");
    // TH1F* h_18_el_corup = get_hist(f_18_el_corup, mcut, "none");
    // TH1F* h_18_el_cordown = get_hist(f_18_el_cordown, mcut, "none");
    TH1F* h_18_el_jmsup = get_hist(f_18_el_jmsup, mcut, "none");
    TH1F* h_18_el_jmsdown = get_hist(f_18_el_jmsdown, mcut, "none");
    // add channels
    TH1F* h_18_nominal = AddHists(h_18_mu_nominal, h_18_el_nominal);
    TH1F* h_18_fsrupsqrt2 = AddHists(h_18_mu_fsrupsqrt2, h_18_el_fsrupsqrt2);
    TH1F* h_18_fsrdownsqrt2 = AddHists(h_18_mu_fsrdownsqrt2, h_18_el_fsrdownsqrt2);
    TH1F* h_18_fsrup2 = AddHists(h_18_mu_fsrup2, h_18_el_fsrup2);
    TH1F* h_18_fsrdown2 = AddHists(h_18_mu_fsrdown2, h_18_el_fsrdown2);
    TH1F* h_18_fsrup4 = AddHists(h_18_mu_fsrup4, h_18_el_fsrup4);
    TH1F* h_18_fsrdown4 = AddHists(h_18_mu_fsrdown4, h_18_el_fsrdown4);
    // TH1F* h_18_jecup = AddHists(h_18_mu_jecup, h_18_el_jecup);
    // TH1F* h_18_jecdown = AddHists(h_18_mu_jecdown, h_18_el_jecdown);
    // TH1F* h_18_corup = AddHists(h_18_mu_corup, h_18_el_corup);
    // TH1F* h_18_cordown = AddHists(h_18_mu_cordown, h_18_el_cordown);
    TH1F* h_18_jmsup = AddHists(h_18_mu_jmsup, h_18_el_jmsup);
    TH1F* h_18_jmsdown = AddHists(h_18_mu_jmsdown, h_18_el_jmsdown);
    TH1F* h_18_st = AddHists(h_18_mu_st, h_18_el_st);
    TH1F* h_18_wj = AddHists(h_18_mu_wj, h_18_el_wj);
    TH1F* h_18_ot = AddHists(h_18_mu_ot, h_18_el_ot);
    TH1F* h_18_data = AddHists(h_18_mu_data, h_18_el_data);

    // combine 17 + 18
    TH1F* h_1718_mu_nominal = AddHists(h_17_mu_nominal, h_18_mu_nominal);
    TH1F* h_1718_mu_fsrupsqrt2 = AddHists(h_17_mu_fsrupsqrt2, h_18_mu_fsrupsqrt2);
    TH1F* h_1718_mu_fsrdownsqrt2 = AddHists(h_17_mu_fsrdownsqrt2, h_18_mu_fsrdownsqrt2);
    TH1F* h_1718_mu_fsrup2 = AddHists(h_17_mu_fsrup2, h_18_mu_fsrup2);
    TH1F* h_1718_mu_fsrdown2 = AddHists(h_17_mu_fsrdown2, h_18_mu_fsrdown2);
    TH1F* h_1718_mu_fsrup4 = AddHists(h_17_mu_fsrup4, h_18_mu_fsrup4);
    TH1F* h_1718_mu_fsrdown4 = AddHists(h_17_mu_fsrdown4, h_18_mu_fsrdown4);
    // TH1F* h_1718_mu_jecup = AddHists(h_17_mu_jecup, h_18_mu_jecup);
    // TH1F* h_1718_mu_jecdown = AddHists(h_17_mu_jecdown, h_18_mu_jecdown);
    // TH1F* h_1718_mu_corup = AddHists(h_17_mu_corup, h_18_mu_corup);
    // TH1F* h_1718_mu_cordown = AddHists(h_17_mu_cordown, h_18_mu_cordown);
    TH1F* h_1718_mu_jmsup = AddHists(h_17_mu_jmsup, h_18_mu_jmsup);
    TH1F* h_1718_mu_jmsdown = AddHists(h_17_mu_jmsdown, h_18_mu_jmsdown);
    TH1F* h_1718_mu_st = AddHists(h_17_mu_st, h_18_mu_st);
    TH1F* h_1718_mu_wj = AddHists(h_17_mu_wj, h_18_mu_wj);
    TH1F* h_1718_mu_ot = AddHists(h_17_mu_ot, h_18_mu_ot);
    TH1F* h_1718_mu_data = AddHists(h_17_mu_data, h_18_mu_data);

    TH1F* h_1718_el_nominal = AddHists(h_17_el_nominal, h_18_el_nominal);
    TH1F* h_1718_el_fsrupsqrt2 = AddHists(h_17_el_fsrupsqrt2, h_18_el_fsrupsqrt2);
    TH1F* h_1718_el_fsrdownsqrt2 = AddHists(h_17_el_fsrdownsqrt2, h_18_el_fsrdownsqrt2);
    TH1F* h_1718_el_fsrup2 = AddHists(h_17_el_fsrup2, h_18_el_fsrup2);
    TH1F* h_1718_el_fsrdown2 = AddHists(h_17_el_fsrdown2, h_18_el_fsrdown2);
    TH1F* h_1718_el_fsrup4 = AddHists(h_17_el_fsrup4, h_18_el_fsrup4);
    TH1F* h_1718_el_fsrdown4 = AddHists(h_17_el_fsrdown4, h_18_el_fsrdown4);
    // TH1F* h_1718_el_jecup = AddHists(h_17_el_jecup, h_18_el_jecup);
    // TH1F* h_1718_el_jecdown = AddHists(h_17_el_jecdown, h_18_el_jecdown);
    // TH1F* h_1718_el_corup = AddHists(h_17_el_corup, h_18_el_corup);
    // TH1F* h_1718_el_cordown = AddHists(h_17_el_cordown, h_18_el_cordown);
    TH1F* h_1718_el_jmsup = AddHists(h_17_el_jmsup, h_18_el_jmsup);
    TH1F* h_1718_el_jmsdown = AddHists(h_17_el_jmsdown, h_18_el_jmsdown);
    TH1F* h_1718_el_st = AddHists(h_17_el_st, h_18_el_st);
    TH1F* h_1718_el_wj = AddHists(h_17_el_wj, h_18_el_wj);
    TH1F* h_1718_el_ot = AddHists(h_17_el_ot, h_18_el_ot);
    TH1F* h_1718_el_data = AddHists(h_17_el_data, h_18_el_data);

    TH1F* h_1718_nominal = AddHists(h_17_nominal, h_18_nominal);
    TH1F* h_1718_fsrupsqrt2 = AddHists(h_17_fsrupsqrt2, h_18_fsrupsqrt2);
    TH1F* h_1718_fsrdownsqrt2 = AddHists(h_17_fsrdownsqrt2, h_18_fsrdownsqrt2);
    TH1F* h_1718_fsrup2 = AddHists(h_17_fsrup2, h_18_fsrup2);
    TH1F* h_1718_fsrdown2 = AddHists(h_17_fsrdown2, h_18_fsrdown2);
    TH1F* h_1718_fsrup4 = AddHists(h_17_fsrup4, h_18_fsrup4);
    TH1F* h_1718_fsrdown4 = AddHists(h_17_fsrdown4, h_18_fsrdown4);
    // TH1F* h_1718_jecup = AddHists(h_17_jecup, h_18_jecup);
    // TH1F* h_1718_jecdown = AddHists(h_17_jecdown, h_18_jecdown);
    // TH1F* h_1718_corup = AddHists(h_17_corup, h_18_corup);
    // TH1F* h_1718_cordown = AddHists(h_17_cordown, h_18_cordown);
    TH1F* h_1718_jmsup = AddHists(h_17_jmsup, h_18_jmsup);
    TH1F* h_1718_jmsdown = AddHists(h_17_jmsdown, h_18_jmsdown);
    TH1F* h_1718_st = AddHists(h_17_st, h_18_st);
    TH1F* h_1718_wj = AddHists(h_17_wj, h_18_wj);
    TH1F* h_1718_ot = AddHists(h_17_ot, h_18_ot);
    TH1F* h_1718_data = AddHists(h_17_data, h_18_data);

    TFile* f_out = new TFile("FSR_hists_mjet"+mass_string+".root","recreate");
    f_out->cd();

    ///////////////////////////////////////////////////
    ///                    2016                     ///
    ///////////////////////////////////////////////////
    h_16_mu_nominal->Write("nominal_muon_2016");
    h_16_mu_fsrup2->Write("fsrup2_muon_2016");
    h_16_mu_fsrdown2->Write("fsrdown2_muon_2016");
    // h_16_mu_jecup->Write("jecup_2016");
    // h_16_mu_jecdown->Write("jecdown_muon_2016");
    // h_16_mu_corup->Write("corup_muon_2016");
    // h_16_mu_cordown->Write("cordown_muon_2016");
    h_16_mu_jmsup->Write("jmsup_muon_2016");
    h_16_mu_jmsdown->Write("jmsdown_muon_2016");
    h_16_mu_st->Write("bgr_st_muon_2016");
    h_16_mu_wj->Write("bgr_wj_muon_2016");
    h_16_mu_ot->Write("bgr_ot_muon_2016");
    h_16_mu_data->Write("data_muon_2016");

    h_16_el_nominal->Write("nominal_elec_2016");
    h_16_el_fsrup2->Write("fsrup2_elec_2016");
    h_16_el_fsrdown2->Write("fsrdown2_elec_2016");
    // h_16_el_jecup->Write("jecup_2016");
    // h_16_el_jecdown->Write("jecdown_elec_2016");
    // h_16_el_corup->Write("corup_elec_2016");
    // h_16_el_cordown->Write("cordown_elec_2016");
    h_16_el_jmsup->Write("jmsup_elec_2016");
    h_16_el_jmsdown->Write("jmsdown_elec_2016");
    h_16_el_st->Write("bgr_st_elec_2016");
    h_16_el_wj->Write("bgr_wj_elec_2016");
    h_16_el_ot->Write("bgr_ot_elec_2016");
    h_16_el_data->Write("data_elec_2016");

    h_16_nominal->Write("nominal_2016");
    h_16_fsrup2->Write("fsrup2_2016");
    h_16_fsrdown2->Write("fsrdown2_2016");
    // h_16_jecup->Write("jecup_2016");
    // h_16_jecdown->Write("jecdown_2016");
    // h_16_corup->Write("corup_2016");
    // h_16_cordown->Write("cordown_2016");
    h_16_jmsup->Write("jmsup_2016");
    h_16_jmsdown->Write("jmsdown_2016");
    h_16_st->Write("bgr_st_2016");
    h_16_wj->Write("bgr_wj_2016");
    h_16_ot->Write("bgr_ot_2016");
    h_16_data->Write("data_2016");

    ///////////////////////////////////////////////////
    ///                    2017                     ///
    ///////////////////////////////////////////////////
    h_17_mu_nominal->Write("nominal_muon_2017");
    h_17_mu_fsrupsqrt2->Write("fsrupsqrt2_muon_2017");
    h_17_mu_fsrup2->Write("fsrup2_muon_2017");
    h_17_mu_fsrup4->Write("fsrup4_muon_2017");
    h_17_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon_2017");
    h_17_mu_fsrdown2->Write("fsrdown2_muon_2017");
    h_17_mu_fsrdown4->Write("fsrdown4_muon_2017");
    // h_17_mu_jecup->Write("jecup_muon_2017");
    // h_17_mu_jecdown->Write("jecdown_muon_2017");
    // h_17_mu_corup->Write("corup_muon_2017");
    // h_17_mu_cordown->Write("cordown_muon_2017");
    h_17_mu_jmsup->Write("jmsup_muon_2017");
    h_17_mu_jmsdown->Write("jmsdown_muon_2017");
    h_17_mu_st->Write("bgr_st_muon_2017");
    h_17_mu_wj->Write("bgr_wj_muon_2017");
    h_17_mu_ot->Write("bgr_ot_muon_2017");
    h_17_mu_data->Write("data_muon_2017");

    h_17_el_nominal->Write("nominal_elec_2017");
    h_17_el_fsrupsqrt2->Write("fsrupsqrt2_elec_2017");
    h_17_el_fsrup2->Write("fsrup2_elec_2017");
    h_17_el_fsrup4->Write("fsrup4_elec_2017");
    h_17_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec_2017");
    h_17_el_fsrdown2->Write("fsrdown2_elec_2017");
    h_17_el_fsrdown4->Write("fsrdown4_elec_2017");
    // h_17_el_jecup->Write("jecup_elec_2017");
    // h_17_el_jecdown->Write("jecdown_elec_2017");
    // h_17_el_corup->Write("corup_elec_2017");
    // h_17_el_cordown->Write("cordown_elec_2017");
    h_17_el_jmsup->Write("jmsup_elec_2017");
    h_17_el_jmsdown->Write("jmsdown_elec_2017");
    h_17_el_st->Write("bgr_st_elec_2017");
    h_17_el_wj->Write("bgr_wj_elec_2017");
    h_17_el_ot->Write("bgr_ot_elec_2017");
    h_17_el_data->Write("data_elec_2017");

    h_17_nominal->Write("nominal_2017");
    h_17_fsrupsqrt2->Write("fsrupsqrt2_2017");
    h_17_fsrup2->Write("fsrup2_2017");
    h_17_fsrup4->Write("fsrup4_2017");
    h_17_fsrdownsqrt2->Write("fsrdownsqrt2_2017");
    h_17_fsrdown2->Write("fsrdown2_2017");
    h_17_fsrdown4->Write("fsrdown4_2017");
    // h_17_jecup->Write("jecup_2017");
    // h_17_jecdown->Write("jecdown_2017");
    // h_17_corup->Write("corup_2017");
    // h_17_cordown->Write("cordown_2017");
    h_17_jmsup->Write("jmsup_2017");
    h_17_jmsdown->Write("jmsdown_2017");
    h_17_st->Write("bgr_st_2017");
    h_17_wj->Write("bgr_wj_2017");
    h_17_ot->Write("bgr_ot_2017");
    h_17_data->Write("data_2017");

    ///////////////////////////////////////////////////
    ///                    2018                     ///
    ///////////////////////////////////////////////////
    h_18_mu_nominal->Write("nominal_muon_2018");
    h_18_mu_fsrupsqrt2->Write("fsrupsqrt2_muon_2018");
    h_18_mu_fsrup2->Write("fsrup2_muon_2018");
    h_18_mu_fsrup4->Write("fsrup4_muon_2018");
    h_18_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon_2018");
    h_18_mu_fsrdown2->Write("fsrdown2_muon_2018");
    h_18_mu_fsrdown4->Write("fsrdown4_muon_2018");
    // h_18_mu_jecup->Write("jecup_muon_2018");
    // h_18_mu_jecdown->Write("jecdown_muon_2018");
    // h_18_mu_corup->Write("corup_muon_2018");
    // h_18_mu_cordown->Write("cordown_muon_2018");
    h_18_mu_jmsup->Write("jmsup_muon_2018");
    h_18_mu_jmsdown->Write("jmsdown_muon_2018");
    h_18_mu_st->Write("bgr_st_muon_2018");
    h_18_mu_wj->Write("bgr_wj_muon_2018");
    h_18_mu_ot->Write("bgr_ot_muon_2018");
    h_18_mu_data->Write("data_muon_2018");

    h_18_el_nominal->Write("nominal_elec_2018");
    h_18_el_fsrupsqrt2->Write("fsrupsqrt2_elec_2018");
    h_18_el_fsrup2->Write("fsrup2_elec_2018");
    h_18_el_fsrup4->Write("fsrup4_elec_2018");
    h_18_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec_2018");
    h_18_el_fsrdown2->Write("fsrdown2_elec_2018");
    h_18_el_fsrdown4->Write("fsrdown4_elec_2018");
    // h_18_el_jecup->Write("jecup_elec_2018");
    // h_18_el_jecdown->Write("jecdown_elec_2018");
    // h_18_el_corup->Write("corup_elec_2018");
    // h_18_el_cordown->Write("cordown_elec_2018");
    h_18_el_jmsup->Write("jmsup_elec_2018");
    h_18_el_jmsdown->Write("jmsdown_elec_2018");
    h_18_el_st->Write("bgr_st_elec_2018");
    h_18_el_wj->Write("bgr_wj_elec_2018");
    h_18_el_ot->Write("bgr_ot_elec_2018");
    h_18_el_data->Write("data_elec_2018");

    h_18_nominal->Write("nominal_2018");
    h_18_fsrupsqrt2->Write("fsrupsqrt2_2018");
    h_18_fsrup2->Write("fsrup2_2018");
    h_18_fsrup4->Write("fsrup4_2018");
    h_18_fsrdownsqrt2->Write("fsrdownsqrt2_2018");
    h_18_fsrdown2->Write("fsrdown2_2018");
    h_18_fsrdown4->Write("fsrdown4_2018");
    // h_18_jecup->Write("jecup_2018");
    // h_18_jecdown->Write("jecdown_2018");
    // h_18_corup->Write("corup_2018");
    // h_18_cordown->Write("cordown_2018");
    h_18_jmsup->Write("jmsup_2018");
    h_18_jmsdown->Write("jmsdown_2018");
    h_18_st->Write("bgr_st_2018");
    h_18_wj->Write("bgr_wj_2018");
    h_18_ot->Write("bgr_ot_2018");
    h_18_data->Write("data_2018");

    ///////////////////////////////////////////////////
    ///                    1718                     ///
    ///////////////////////////////////////////////////
    h_1718_mu_nominal->Write("nominal_muon_combine");
    h_1718_mu_fsrupsqrt2->Write("fsrupsqrt2_muon_combine");
    h_1718_mu_fsrup2->Write("fsrup2_muon_combine");
    h_1718_mu_fsrup4->Write("fsrup4_muon_combine");
    h_1718_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon_combine");
    h_1718_mu_fsrdown2->Write("fsrdown2_muon_combine");
    h_1718_mu_fsrdown4->Write("fsrdown4_muon_combine");
    // h_1718_mu_jecup->Write("jecup_muon_combine");
    // h_1718_mu_jecdown->Write("jecdown_muon_combine");
    // h_1718_mu_corup->Write("corup_muon_combine");
    // h_1718_mu_cordown->Write("cordown_muon_combine");
    h_1718_mu_jmsup->Write("jmsup_muon_combine");
    h_1718_mu_jmsdown->Write("jmsdown_muon_combine");
    h_1718_mu_st->Write("bgr_st_muon_combine");
    h_1718_mu_wj->Write("bgr_wj_muon_combine");
    h_1718_mu_ot->Write("bgr_ot_muon_combine");
    h_1718_mu_data->Write("data_muon_combine");

    h_1718_el_nominal->Write("nominal_elec_combine");
    h_1718_el_fsrupsqrt2->Write("fsrupsqrt2_elec_combine");
    h_1718_el_fsrup2->Write("fsrup2_elec_combine");
    h_1718_el_fsrup4->Write("fsrup4_elec_combine");
    h_1718_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec_combine");
    h_1718_el_fsrdown2->Write("fsrdown2_elec_combine");
    h_1718_el_fsrdown4->Write("fsrdown4_elec_combine");
    // h_1718_el_jecup->Write("jecup_elec_combine");
    // h_1718_el_jecdown->Write("jecdown_elec_combine");
    // h_1718_el_corup->Write("corup_elec_combine");
    // h_1718_el_cordown->Write("cordown_elec_combine");
    h_1718_el_jmsup->Write("jmsup_elec_combine");
    h_1718_el_jmsdown->Write("jmsdown_elec_combine");
    h_1718_el_st->Write("bgr_st_elec_combine");
    h_1718_el_wj->Write("bgr_wj_elec_combine");
    h_1718_el_ot->Write("bgr_ot_elec_combine");
    h_1718_el_data->Write("data_elec_combine");

    h_1718_nominal->Write("nominal_combine");
    h_1718_fsrupsqrt2->Write("fsrupsqrt2_combine");
    h_1718_fsrup2->Write("fsrup2_combine");
    h_1718_fsrup4->Write("fsrup4_combine");
    h_1718_fsrdownsqrt2->Write("fsrdownsqrt2_combine");
    h_1718_fsrdown2->Write("fsrdown2_combine");
    h_1718_fsrdown4->Write("fsrdown4_combine");
    // h_1718_jecup->Write("jecup_combine");
    // h_1718_jecdown->Write("jecdown_combine");
    // h_1718_corup->Write("corup_combine");
    // h_1718_cordown->Write("cordown_combine");
    h_1718_jmsup->Write("jmsup_combine");
    h_1718_jmsdown->Write("jmsdown_combine");
    h_1718_st->Write("bgr_st_combine");
    h_1718_wj->Write("bgr_wj_combine");
    h_1718_ot->Write("bgr_ot_combine");
    h_1718_data->Write("data_combine");

    f_out->Close();

    auto stop  = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "  - time needed: " << GREEN << duration.count()/1000 << "s" << RESET << endl;
  }
  return 0;
}

// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
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

// -------------------------------------------------------------------------------------------------------
TH1F* AddHists(TH1F* h1, TH1F* h2){
  TH1F* hist_added = (TH1F*) h1->Clone();
  hist_added->Add(h2);
  return hist_added;
}
