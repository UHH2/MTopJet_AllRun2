#include "../include/CentralInclude.h"
#include "../include/Utils.h"
#include "TSystem.h"


// #include "TString.h"

using namespace std;

TString study, bin, year;

bool studyPt = false;
bool studyMjet = false;
bool debug = false;

void FillAllHists(vector<TString> window, vector<vector<double>> cutvalues, TString var);
TH1F* get_hist(TFile* file, vector<double> cut, TString weightname, TString var, bool weight_cut=false);
TH1F* AddHists(TH1F* h1, TH1F* h2);
VecDD CutValuesFromWindow(vector<TString> windows);

// TString dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel";
TString directory = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel";
TFile* f_16_mu           = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_fsrup     = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
TFile* f_16_mu_fsrdown   = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
TFile* f_16_mu_st        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
TFile* f_16_mu_wj        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
TFile* f_16_mu_ot        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
TFile* f_16_mu_data      = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
TFile* f_16_mu_jecup     = new TFile(directory+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_jecdown   = new TFile(directory+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_corup     = new TFile(directory+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_cordown   = new TFile(directory+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_jmsup     = new TFile(directory+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_jmsdown   = new TFile(directory+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_mu_hdampup   = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2016v3.root");
TFile* f_16_mu_hdampdown = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2016v3.root");
TFile* f_16_mu_isrup     = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_isrup_2016v3.root");
TFile* f_16_mu_isrdown   = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_isrdown_2016v3.root");
TFile* f_16_mu_tuneup    = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2016v3.root");
TFile* f_16_mu_tunedown  = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2016v3.root");
TFile* f_16_mu_gluon     = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2016v3.root");
TFile* f_16_mu_qcd       = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2016v3.root");

TFile* f_16_el           = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_fsrup     = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
TFile* f_16_el_fsrdown   = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");
TFile* f_16_el_st        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2016v3.root");
TFile* f_16_el_wj        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2016v3.root");
TFile* f_16_el_ot        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.other_2016v3.root");
TFile* f_16_el_data      = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2016v3.root");
TFile* f_16_el_jecup     = new TFile(directory+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_jecdown   = new TFile(directory+"/elec/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_corup     = new TFile(directory+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_cordown   = new TFile(directory+"/elec/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_jmsup     = new TFile(directory+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_jmsdown   = new TFile(directory+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
TFile* f_16_el_hdampup   = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2016v3.root");
TFile* f_16_el_hdampdown = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2016v3.root");
TFile* f_16_el_isrup     = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_isrup_2016v3.root");
TFile* f_16_el_isrdown   = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_isrdown_2016v3.root");
TFile* f_16_el_tuneup    = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2016v3.root");
TFile* f_16_el_tunedown  = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2016v3.root");
TFile* f_16_el_gluon     = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2016v3.root");
TFile* f_16_el_qcd       = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2016v3.root");

TFile* f_17_mu           = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
// TFile* f_17_mu           = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2L2Nu_Mtt0000to0700_2017v2.root");
// TFile* f_17_mu_fsrup     = new TFile(directory+"/muon/FSRup_4/uhh2.AnalysisModuleRunner.MC.TTbar_2L2Nu_Mtt0000to0700_2017v2.root");
// TFile* f_17_mu_fsrup_test= new TFile(directory+"/muon/FSRup_4_test/uhh2.AnalysisModuleRunner.MC.TTbar_2L2Nu_Mtt0000to0700_2017v2.root");
TFile* f_17_mu_st        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
TFile* f_17_mu_wj        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
TFile* f_17_mu_ot        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
TFile* f_17_mu_data      = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
TFile* f_17_mu_jecup     = new TFile(directory+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_jecdown   = new TFile(directory+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_corup     = new TFile(directory+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_cordown   = new TFile(directory+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_jmsup     = new TFile(directory+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_jmsdown   = new TFile(directory+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_mu_hdampup   = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2017v2.root");
TFile* f_17_mu_hdampdown = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2017v2.root");
TFile* f_17_mu_tuneup    = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2017v2.root");
TFile* f_17_mu_tunedown  = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2017v2.root");
TFile* f_17_mu_gluon     = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2017v2.root");
TFile* f_17_mu_qcd       = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2017v2.root");

TFile* f_17_el           = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_st        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2017v2.root");
TFile* f_17_el_wj        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2017v2.root");
TFile* f_17_el_ot        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.other_2017v2.root");
TFile* f_17_el_data      = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2017v2.root");
TFile* f_17_el_jecup     = new TFile(directory+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_jecdown   = new TFile(directory+"/elec/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_corup     = new TFile(directory+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_cordown   = new TFile(directory+"/elec/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_jmsup     = new TFile(directory+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_jmsdown   = new TFile(directory+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
TFile* f_17_el_hdampup   = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2017v2.root");
TFile* f_17_el_hdampdown = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2017v2.root");
TFile* f_17_el_tuneup    = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2017v2.root");
TFile* f_17_el_tunedown  = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2017v2.root");
TFile* f_17_el_gluon     = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2017v2.root");
TFile* f_17_el_qcd       = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2017v2.root");

TFile* f_18_mu           = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_st        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
TFile* f_18_mu_wj        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
TFile* f_18_mu_ot        = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.other_2018.root");
TFile* f_18_mu_data      = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
TFile* f_18_mu_jecup     = new TFile(directory+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_jecdown   = new TFile(directory+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_corup     = new TFile(directory+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_cordown   = new TFile(directory+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_jmsup     = new TFile(directory+"/muon/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_jmsdown   = new TFile(directory+"/muon/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_mu_hdampup   = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2018.root");
TFile* f_18_mu_hdampdown = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2018.root");
TFile* f_18_mu_tuneup    = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2018.root");
TFile* f_18_mu_tunedown  = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2018.root");
TFile* f_18_mu_gluon     = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2018.root");
TFile* f_18_mu_qcd       = new TFile(directory+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2018.root");

TFile* f_18_el           = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_st        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.SingleTop_2018.root");
TFile* f_18_el_wj        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.WJets_2018.root");
TFile* f_18_el_ot        = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.other_2018.root");
TFile* f_18_el_data      = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.DATA.DATA_2018.root");
TFile* f_18_el_jecup     = new TFile(directory+"/elec/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_jecdown   = new TFile(directory+"/elec/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_corup     = new TFile(directory+"/elec/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_cordown   = new TFile(directory+"/elec/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_jmsup     = new TFile(directory+"/elec/JMS_upup/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_jmsdown   = new TFile(directory+"/elec/JMS_downdown/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");
TFile* f_18_el_hdampup   = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampup_2018.root");
TFile* f_18_el_hdampdown = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_hdampdown_2018.root");
TFile* f_18_el_tuneup    = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneUp_2018.root");
TFile* f_18_el_tunedown  = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_TuneDown_2018.root");
TFile* f_18_el_gluon     = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_GluonMove_2018.root");
TFile* f_18_el_qcd       = new TFile(directory+"/elec/uhh2.AnalysisModuleRunner.MC.TTbar_QCDbased_2018.root");

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  cout << "-------------------------------------------" << endl;
  cout << "Usage: ./CreatFSRPlots <binstudy> <bin>" << endl;
  cout << "binstudy: pt, mjet" << endl;
  cout << "bin: e.g. 140toInf" << endl;
  cout << "-------------------------------------------" << endl;

  TH1::AddDirectory(kFALSE);
  study         = argv[1];
  bin           = argv[2];

  studyPt = study.EqualTo("pt");
  studyMjet = study.EqualTo("mjet");

  if(!(studyPt||studyMjet)){cerr << "Bin study please! Check above." << endl; return 0;}

  // ---------------------------------------------------------------------------
  //                            Bin studies
  // ---------------------------------------------------------------------------

  vector<TString> window = {(TString) argv[2]};
  VecDD cutvalues = CutValuesFromWindow(window);

  // ---------------------------------------------------------------------------
  //                            Read Files
  // ---------------------------------------------------------------------------

  FillAllHists(window, cutvalues, "tau32");
  if(studyPt)   FillAllHists(window, cutvalues, "pt");
  if(studyMjet) FillAllHists(window, cutvalues,  "mjet");

}

// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
void FillAllHists(vector<TString> window, vector<vector<double>> cutvalues, TString var){

  bool isTau = var=="tau32";
  bool isMjet = var=="mjet";
  bool isPt = var=="pt";

  for(unsigned int i=0; i<window.size(); i++){

    auto start = high_resolution_clock::now(); // Calculation time - start

    TString bin_string = window[i];
    vector<double> cut = cutvalues[i];
    cout << study << " Bin: " << bin_string << endl;

    // ------------------------------------------------------------------------- 2016
    if(debug) cout << "Start with 2016 ... " << endl;
    if(debug) cout << "\t ... muon channel" << endl;
    TH1F* h_16_mu_nominal = get_hist(f_16_mu, cut, "none", var);
    TH1F* h_16_mu_fsrup2 = get_hist(f_16_mu_fsrup, cut, "none", var);
    TH1F* h_16_mu_fsrdown2 = get_hist(f_16_mu_fsrdown, cut, "none", var);
    TH1F* h_16_mu_st = get_hist(f_16_mu_st, cut, "none", var);
    TH1F* h_16_mu_wj = get_hist(f_16_mu_wj, cut, "none", var);
    TH1F* h_16_mu_ot = get_hist(f_16_mu_ot, cut, "none", var);
    TH1F* h_16_mu_data = get_hist(f_16_mu_data, cut, "none", var);
    TH1F* h_16_mu_jecup = get_hist(f_16_mu_jecup, cut, "none", var);
    TH1F* h_16_mu_jecdown = get_hist(f_16_mu_jecdown, cut, "none", var);
    TH1F* h_16_mu_corup = get_hist(f_16_mu_corup, cut, "none", var);
    TH1F* h_16_mu_cordown = get_hist(f_16_mu_cordown, cut, "none", var);
    TH1F* h_16_mu_jmsup = get_hist(f_16_mu_jmsup, cut, "none", var);
    TH1F* h_16_mu_jmsdown = get_hist(f_16_mu_jmsdown, cut, "none", var);
    TH1F* h_16_mu_hdampup = get_hist(f_16_mu_hdampup, cut, "none", var);
    TH1F* h_16_mu_hdampdown = get_hist(f_16_mu_hdampdown, cut, "none", var);
    TH1F* h_16_mu_isrup = get_hist(f_16_mu_isrup, cut, "none", var);
    TH1F* h_16_mu_isrdown = get_hist(f_16_mu_isrdown, cut, "none", var);
    TH1F* h_16_mu_tuneup = get_hist(f_16_mu_tuneup, cut, "none", var);
    TH1F* h_16_mu_tunedown = get_hist(f_16_mu_tunedown, cut, "none", var);
    TH1F* h_16_mu_gluon = get_hist(f_16_mu_gluon, cut, "none", var);
    TH1F* h_16_mu_qcd = get_hist(f_16_mu_qcd, cut, "none", var);

    if(debug) cout << "\t ... elec channel" << endl;
    TH1F* h_16_el_nominal = get_hist(f_16_el, cut, "none", var);
    TH1F* h_16_el_fsrup2 = get_hist(f_16_el_fsrup, cut, "none", var);
    TH1F* h_16_el_fsrdown2 = get_hist(f_16_el_fsrdown, cut, "none", var);
    TH1F* h_16_el_st = get_hist(f_16_el_st, cut, "none", var);
    TH1F* h_16_el_wj = get_hist(f_16_el_wj, cut, "none", var);
    TH1F* h_16_el_ot = get_hist(f_16_el_ot, cut, "none", var);
    TH1F* h_16_el_data = get_hist(f_16_el_data, cut, "none", var);
    TH1F* h_16_el_jecup = get_hist(f_16_el_jecup, cut, "none", var);
    TH1F* h_16_el_jecdown = get_hist(f_16_el_jecdown, cut, "none", var);
    TH1F* h_16_el_corup = get_hist(f_16_el_corup, cut, "none", var);
    TH1F* h_16_el_cordown = get_hist(f_16_el_cordown, cut, "none", var);
    TH1F* h_16_el_jmsup = get_hist(f_16_el_jmsup, cut, "none", var);
    TH1F* h_16_el_jmsdown = get_hist(f_16_el_jmsdown, cut, "none", var);
    TH1F* h_16_el_hdampup = get_hist(f_16_el_hdampup, cut, "none", var);
    TH1F* h_16_el_hdampdown = get_hist(f_16_el_hdampdown, cut, "none", var);
    TH1F* h_16_el_isrup = get_hist(f_16_el_isrup, cut, "none", var);
    TH1F* h_16_el_isrdown = get_hist(f_16_el_isrdown, cut, "none", var);
    TH1F* h_16_el_tuneup = get_hist(f_16_el_tuneup, cut, "none", var);
    TH1F* h_16_el_tunedown = get_hist(f_16_el_tunedown, cut, "none", var);
    TH1F* h_16_el_gluon = get_hist(f_16_el_gluon, cut, "none", var);
    TH1F* h_16_el_qcd = get_hist(f_16_el_qcd, cut, "none", var);

    // ------------------------------------------------------------ add channels

    if(debug) cout << "\t ... combine channels" << endl;
    TH1F* h_16_nominal = AddHists(h_16_mu_nominal, h_16_el_nominal);
    TH1F* h_16_fsrup2 = AddHists(h_16_mu_fsrup2, h_16_el_fsrup2);
    TH1F* h_16_fsrdown2 = AddHists(h_16_mu_fsrdown2, h_16_el_fsrdown2);
    TH1F* h_16_jecup = AddHists(h_16_mu_jecup, h_16_el_jecup);
    TH1F* h_16_jecdown = AddHists(h_16_mu_jecdown, h_16_el_jecdown);
    TH1F* h_16_corup = AddHists(h_16_mu_corup, h_16_el_corup);
    TH1F* h_16_cordown = AddHists(h_16_mu_cordown, h_16_el_cordown);
    TH1F* h_16_jmsup = AddHists(h_16_mu_jmsup, h_16_el_jmsup);
    TH1F* h_16_jmsdown = AddHists(h_16_mu_jmsdown, h_16_el_jmsdown);
    TH1F* h_16_hdampup = AddHists(h_16_mu_hdampup, h_16_el_hdampup);
    TH1F* h_16_hdampdown = AddHists(h_16_mu_hdampdown, h_16_el_hdampdown);
    TH1F* h_16_isrup = AddHists(h_16_mu_isrup, h_16_el_isrup);
    TH1F* h_16_isrdown = AddHists(h_16_mu_isrdown, h_16_el_isrdown);
    TH1F* h_16_tuneup = AddHists(h_16_mu_tuneup, h_16_el_tuneup);
    TH1F* h_16_tunedown = AddHists(h_16_mu_tunedown, h_16_el_tunedown);
    TH1F* h_16_st = AddHists(h_16_mu_st, h_16_el_st);
    TH1F* h_16_wj = AddHists(h_16_mu_wj, h_16_el_wj);
    TH1F* h_16_ot = AddHists(h_16_mu_ot, h_16_el_ot);
    TH1F* h_16_data = AddHists(h_16_mu_data, h_16_el_data);
    TH1F* h_16_gluon = AddHists(h_16_mu_gluon, h_16_el_gluon);
    TH1F* h_16_qcd = AddHists(h_16_mu_qcd, h_16_el_qcd);

    // ------------------------------------------------------------------------- 2017
    if(debug) cout << "Start with 2017 ... " << endl;
    if(debug) cout << "\t ... muon channel" << endl;
    TH1F* h_17_mu_nominal = get_hist(f_17_mu, cut, "none", var);
    TH1F* h_17_mu_st = get_hist(f_17_mu_st, cut, "none", var);
    TH1F* h_17_mu_wj = get_hist(f_17_mu_wj, cut, "none", var);
    TH1F* h_17_mu_ot = get_hist(f_17_mu_ot, cut, "none", var);
    TH1F* h_17_mu_data = get_hist(f_17_mu_data, cut, "none", var);
    TH1F* h_17_mu_jecup = get_hist(f_17_mu_jecup, cut, "none", var);
    TH1F* h_17_mu_jecdown = get_hist(f_17_mu_jecdown, cut, "none", var);
    TH1F* h_17_mu_corup = get_hist(f_17_mu_corup, cut, "none", var);
    TH1F* h_17_mu_cordown = get_hist(f_17_mu_cordown, cut, "none", var);
    TH1F* h_17_mu_jmsup = get_hist(f_17_mu_jmsup, cut, "none", var);
    TH1F* h_17_mu_jmsdown = get_hist(f_17_mu_jmsdown, cut, "none", var);
    TH1F* h_17_mu_fsrupsqrt2 = get_hist(f_17_mu, cut, "sf_fsr_upsqrt2", var);
    TH1F* h_17_mu_fsrup2 = get_hist(f_17_mu, cut, "sf_fsr_up2", var);
    TH1F* h_17_mu_fsrup4 = get_hist(f_17_mu, cut, "sf_fsr_up4", var);
    // TH1F* h_17_mu_fsrup4_cut = get_hist(f_17_mu, cut, "sf_fsr_up4", var, true);
    TH1F* h_17_mu_fsrdownsqrt2 = get_hist(f_17_mu, cut, "sf_fsr_downsqrt2", var);
    TH1F* h_17_mu_fsrdown2 = get_hist(f_17_mu, cut, "sf_fsr_down2", var);
    TH1F* h_17_mu_fsrdown4 = get_hist(f_17_mu, cut, "sf_fsr_down4", var);
    TH1F* h_17_mu_hdampup = get_hist(f_17_mu_hdampup, cut, "none", var);
    TH1F* h_17_mu_hdampdown = get_hist(f_17_mu_hdampdown, cut, "none", var);
    TH1F* h_17_mu_tuneup = get_hist(f_17_mu_tuneup, cut, "none", var);
    TH1F* h_17_mu_tunedown = get_hist(f_17_mu_tunedown, cut, "none", var);
    TH1F* h_17_mu_isrup = get_hist(f_17_mu, cut, "sf_isr_up2", var);
    TH1F* h_17_mu_isrdown = get_hist(f_17_mu, cut, "sf_isr_down2", var);
    TH1F* h_17_mu_gluon = get_hist(f_17_mu_gluon, cut, "none", var);
    TH1F* h_17_mu_qcd = get_hist(f_17_mu_qcd, cut, "none", var);

    if(debug) cout << "\t ... elec channel" << endl;
    TH1F* h_17_el_nominal = get_hist(f_17_el, cut, "none", var);
    TH1F* h_17_el_st = get_hist(f_17_el_st, cut, "none", var);
    TH1F* h_17_el_wj = get_hist(f_17_el_wj, cut, "none", var);
    TH1F* h_17_el_ot = get_hist(f_17_el_ot, cut, "none", var);
    TH1F* h_17_el_data = get_hist(f_17_el_data, cut, "none", var);
    TH1F* h_17_el_fsrupsqrt2 = get_hist(f_17_el, cut, "sf_fsr_upsqrt2", var);
    TH1F* h_17_el_fsrup2 = get_hist(f_17_el, cut, "sf_fsr_up2", var);
    TH1F* h_17_el_fsrup4 = get_hist(f_17_el, cut, "sf_fsr_up4", var);
    // TH1F* h_17_el_fsrup4_cut = get_hist(f_17_el, cut, "sf_fsr_up4", var, true);
    TH1F* h_17_el_fsrdownsqrt2 = get_hist(f_17_el, cut, "sf_fsr_downsqrt2", var);
    TH1F* h_17_el_fsrdown2 = get_hist(f_17_el, cut, "sf_fsr_down2", var);
    TH1F* h_17_el_fsrdown4 = get_hist(f_17_el, cut, "sf_fsr_down4", var);
    TH1F* h_17_el_jecup = get_hist(f_17_el_jecup, cut, "none", var);
    TH1F* h_17_el_jecdown = get_hist(f_17_el_jecdown, cut, "none", var);
    TH1F* h_17_el_corup = get_hist(f_17_el_corup, cut, "none", var);
    TH1F* h_17_el_cordown = get_hist(f_17_el_cordown, cut, "none", var);
    TH1F* h_17_el_jmsup = get_hist(f_17_el_jmsup, cut, "none", var);
    TH1F* h_17_el_jmsdown = get_hist(f_17_el_jmsdown, cut, "none", var);
    TH1F* h_17_el_hdampup = get_hist(f_17_el_hdampup, cut, "none", var);
    TH1F* h_17_el_hdampdown = get_hist(f_17_el_hdampdown, cut, "none", var);
    TH1F* h_17_el_tuneup = get_hist(f_17_el_tuneup, cut, "none", var);
    TH1F* h_17_el_tunedown = get_hist(f_17_el_tunedown, cut, "none", var);
    TH1F* h_17_el_isrup = get_hist(f_17_el, cut, "sf_isr_up2", var);
    TH1F* h_17_el_isrdown = get_hist(f_17_el, cut, "sf_isr_down2", var);
    TH1F* h_17_el_gluon = get_hist(f_17_el_gluon, cut, "none", var);
    TH1F* h_17_el_qcd = get_hist(f_17_el_qcd, cut, "none", var);

    // ------------------------------------------------------------ add channels

    if(debug) cout << "\t ... combine channels" << endl;
    TH1F* h_17_nominal = AddHists(h_17_mu_nominal, h_17_el_nominal);
    TH1F* h_17_fsrupsqrt2 = AddHists(h_17_mu_fsrupsqrt2, h_17_el_fsrupsqrt2);
    TH1F* h_17_fsrdownsqrt2 = AddHists(h_17_mu_fsrdownsqrt2, h_17_el_fsrdownsqrt2);
    TH1F* h_17_fsrup2 = AddHists(h_17_mu_fsrup2, h_17_el_fsrup2);
    TH1F* h_17_fsrdown2 = AddHists(h_17_mu_fsrdown2, h_17_el_fsrdown2);
    TH1F* h_17_fsrup4 = AddHists(h_17_mu_fsrup4, h_17_el_fsrup4);
    // TH1F* h_17_fsrup4_cut = AddHists(h_17_mu_fsrup4_cut, h_17_el_fsrup4_cut);
    TH1F* h_17_fsrdown4 = AddHists(h_17_mu_fsrdown4, h_17_el_fsrdown4);
    TH1F* h_17_jecup = AddHists(h_17_mu_jecup, h_17_el_jecup);
    TH1F* h_17_jecdown = AddHists(h_17_mu_jecdown, h_17_el_jecdown);
    TH1F* h_17_corup = AddHists(h_17_mu_corup, h_17_el_corup);
    TH1F* h_17_cordown = AddHists(h_17_mu_cordown, h_17_el_cordown);
    TH1F* h_17_jmsup = AddHists(h_17_mu_jmsup, h_17_el_jmsup);
    TH1F* h_17_jmsdown = AddHists(h_17_mu_jmsdown, h_17_el_jmsdown);
    TH1F* h_17_isrup = AddHists(h_17_mu_isrup, h_17_el_isrup);
    TH1F* h_17_isrdown = AddHists(h_17_mu_isrdown, h_17_el_isrdown);
    TH1F* h_17_tuneup = AddHists(h_17_mu_tuneup, h_17_el_tuneup);
    TH1F* h_17_tunedown = AddHists(h_17_mu_tunedown, h_17_el_tunedown);
    TH1F* h_17_hdampup = AddHists(h_17_mu_hdampup, h_17_el_hdampup);
    TH1F* h_17_hdampdown = AddHists(h_17_mu_hdampdown, h_17_el_hdampdown);
    TH1F* h_17_st = AddHists(h_17_mu_st, h_17_el_st);
    TH1F* h_17_wj = AddHists(h_17_mu_wj, h_17_el_wj);
    TH1F* h_17_ot = AddHists(h_17_mu_ot, h_17_el_ot);
    TH1F* h_17_data = AddHists(h_17_mu_data, h_17_el_data);
    TH1F* h_17_gluon = AddHists(h_17_mu_gluon, h_17_el_gluon);
    TH1F* h_17_qcd = AddHists(h_17_mu_qcd, h_17_el_qcd);

    // ------------------------------------------------------------------------- 2018
    if(debug) cout << "Start with 2016 ... " << endl;
    if(debug) cout << "\t ... muon channel" << endl;
    TH1F* h_18_mu_nominal = get_hist(f_18_mu, cut, "none", var);
    TH1F* h_18_mu_st = get_hist(f_18_mu_st, cut, "none", var);
    TH1F* h_18_mu_wj = get_hist(f_18_mu_wj, cut, "none", var);
    TH1F* h_18_mu_ot = get_hist(f_18_mu_ot, cut, "none", var);
    TH1F* h_18_mu_data = get_hist(f_18_mu_data, cut, "none", var);
    TH1F* h_18_mu_jecup = get_hist(f_18_mu_jecup, cut, "none", var);
    TH1F* h_18_mu_jecdown = get_hist(f_18_mu_jecdown, cut, "none", var);
    TH1F* h_18_mu_corup = get_hist(f_18_mu_corup, cut, "none", var);
    TH1F* h_18_mu_cordown = get_hist(f_18_mu_cordown, cut, "none", var);
    TH1F* h_18_mu_jmsup = get_hist(f_18_mu_jmsup, cut, "none", var);
    TH1F* h_18_mu_jmsdown = get_hist(f_18_mu_jmsdown, cut, "none", var);
    TH1F* h_18_mu_fsrupsqrt2 = get_hist(f_18_mu, cut, "sf_fsr_upsqrt2", var);
    TH1F* h_18_mu_fsrup2 = get_hist(f_18_mu, cut, "sf_fsr_up2", var);
    TH1F* h_18_mu_fsrup4 = get_hist(f_18_mu, cut, "sf_fsr_up4", var);
    // TH1F* h_18_mu_fsrup4_cut = get_hist(f_18_mu, cut, "sf_fsr_up4", var, true);
    TH1F* h_18_mu_fsrdownsqrt2 = get_hist(f_18_mu, cut, "sf_fsr_downsqrt2", var);
    TH1F* h_18_mu_fsrdown2 = get_hist(f_18_mu, cut, "sf_fsr_down2", var);
    TH1F* h_18_mu_fsrdown4 = get_hist(f_18_mu, cut, "sf_fsr_down4", var);
    TH1F* h_18_mu_hdampup = get_hist(f_18_mu_hdampup, cut, "none", var);
    TH1F* h_18_mu_hdampdown = get_hist(f_18_mu_hdampdown, cut, "none", var);
    TH1F* h_18_mu_tuneup = get_hist(f_18_mu_tuneup, cut, "none", var);
    TH1F* h_18_mu_tunedown = get_hist(f_18_mu_tunedown, cut, "none", var);
    TH1F* h_18_mu_isrup = get_hist(f_18_mu, cut, "sf_isr_up2", var);
    TH1F* h_18_mu_isrdown = get_hist(f_18_mu, cut, "sf_isr_down2", var);
    TH1F* h_18_mu_gluon = get_hist(f_18_mu_gluon, cut, "none", var);
    TH1F* h_18_mu_qcd = get_hist(f_18_mu_qcd, cut, "none", var);

    if(debug) cout << "\t ... elec channel" << endl;
    TH1F* h_18_el_nominal = get_hist(f_18_el, cut, "none", var);
    TH1F* h_18_el_st = get_hist(f_18_el_st, cut, "none", var);
    TH1F* h_18_el_wj = get_hist(f_18_el_wj, cut, "none", var);
    TH1F* h_18_el_ot = get_hist(f_18_el_ot, cut, "none", var);
    TH1F* h_18_el_data = get_hist(f_18_el_data, cut, "none", var);
    TH1F* h_18_el_fsrupsqrt2 = get_hist(f_18_el, cut, "sf_fsr_upsqrt2", var);
    TH1F* h_18_el_fsrup2 = get_hist(f_18_el, cut, "sf_fsr_up2", var);
    TH1F* h_18_el_fsrup4 = get_hist(f_18_el, cut, "sf_fsr_up4", var);
    // TH1F* h_18_el_fsrup4_cut = get_hist(f_18_el, cut, "sf_fsr_up4", var, true);
    TH1F* h_18_el_fsrdownsqrt2 = get_hist(f_18_el, cut, "sf_fsr_downsqrt2", var);
    TH1F* h_18_el_fsrdown2 = get_hist(f_18_el, cut, "sf_fsr_down2", var);
    TH1F* h_18_el_fsrdown4 = get_hist(f_18_el, cut, "sf_fsr_down4", var);
    TH1F* h_18_el_jecup = get_hist(f_18_el_jecup, cut, "none", var);
    TH1F* h_18_el_jecdown = get_hist(f_18_el_jecdown, cut, "none", var);
    TH1F* h_18_el_corup = get_hist(f_18_el_corup, cut, "none", var);
    TH1F* h_18_el_cordown = get_hist(f_18_el_cordown, cut, "none", var);
    TH1F* h_18_el_jmsup = get_hist(f_18_el_jmsup, cut, "none", var);
    TH1F* h_18_el_jmsdown = get_hist(f_18_el_jmsdown, cut, "none", var);
    TH1F* h_18_el_hdampup = get_hist(f_18_el_hdampup, cut, "none", var);
    TH1F* h_18_el_hdampdown = get_hist(f_18_el_hdampdown, cut, "none", var);
    TH1F* h_18_el_tuneup = get_hist(f_18_el_tuneup, cut, "none", var);
    TH1F* h_18_el_tunedown = get_hist(f_18_el_tunedown, cut, "none", var);
    TH1F* h_18_el_isrup = get_hist(f_18_el, cut, "sf_isr_up2", var);
    TH1F* h_18_el_isrdown = get_hist(f_18_el, cut, "sf_isr_down2", var);
    TH1F* h_18_el_gluon = get_hist(f_18_el_gluon, cut, "none", var);
    TH1F* h_18_el_qcd = get_hist(f_18_el_qcd, cut, "none", var);

    // ------------------------------------------------------------ add channels

    if(debug) cout << "\t ... combine channels" << endl;
    TH1F* h_18_nominal = AddHists(h_18_mu_nominal, h_18_el_nominal);
    TH1F* h_18_fsrupsqrt2 = AddHists(h_18_mu_fsrupsqrt2, h_18_el_fsrupsqrt2);
    TH1F* h_18_fsrdownsqrt2 = AddHists(h_18_mu_fsrdownsqrt2, h_18_el_fsrdownsqrt2);
    TH1F* h_18_fsrup2 = AddHists(h_18_mu_fsrup2, h_18_el_fsrup2);
    TH1F* h_18_fsrdown2 = AddHists(h_18_mu_fsrdown2, h_18_el_fsrdown2);
    TH1F* h_18_fsrup4 = AddHists(h_18_mu_fsrup4, h_18_el_fsrup4);
    // TH1F* h_18_fsrup4_cut = AddHists(h_18_mu_fsrup4_cut, h_18_el_fsrup4_cut);
    TH1F* h_18_fsrdown4 = AddHists(h_18_mu_fsrdown4, h_18_el_fsrdown4);
    TH1F* h_18_jecup = AddHists(h_18_mu_jecup, h_18_el_jecup);
    TH1F* h_18_jecdown = AddHists(h_18_mu_jecdown, h_18_el_jecdown);
    TH1F* h_18_corup = AddHists(h_18_mu_corup, h_18_el_corup);
    TH1F* h_18_cordown = AddHists(h_18_mu_cordown, h_18_el_cordown);
    TH1F* h_18_jmsup = AddHists(h_18_mu_jmsup, h_18_el_jmsup);
    TH1F* h_18_jmsdown = AddHists(h_18_mu_jmsdown, h_18_el_jmsdown);
    TH1F* h_18_isrup = AddHists(h_18_mu_isrup, h_18_el_isrup);
    TH1F* h_18_isrdown = AddHists(h_18_mu_isrdown, h_18_el_isrdown);
    TH1F* h_18_tuneup = AddHists(h_18_mu_tuneup, h_18_el_tuneup);
    TH1F* h_18_tunedown = AddHists(h_18_mu_tunedown, h_18_el_tunedown);
    TH1F* h_18_hdampup = AddHists(h_18_mu_hdampup, h_18_el_hdampup);
    TH1F* h_18_hdampdown = AddHists(h_18_mu_hdampdown, h_18_el_hdampdown);
    TH1F* h_18_st = AddHists(h_18_mu_st, h_18_el_st);
    TH1F* h_18_wj = AddHists(h_18_mu_wj, h_18_el_wj);
    TH1F* h_18_ot = AddHists(h_18_mu_ot, h_18_el_ot);
    TH1F* h_18_data = AddHists(h_18_mu_data, h_18_el_data);
    TH1F* h_18_gluon = AddHists(h_18_mu_gluon, h_18_el_gluon);
    TH1F* h_18_qcd = AddHists(h_18_mu_qcd, h_18_el_qcd);

    // ------------------------------------------------------------------------- combine 17 + 18

    TH1F* h_1718_mu_nominal = AddHists(h_17_mu_nominal, h_18_mu_nominal);
    TH1F* h_1718_mu_fsrupsqrt2 = AddHists(h_17_mu_fsrupsqrt2, h_18_mu_fsrupsqrt2);
    TH1F* h_1718_mu_fsrdownsqrt2 = AddHists(h_17_mu_fsrdownsqrt2, h_18_mu_fsrdownsqrt2);
    TH1F* h_1718_mu_fsrup2 = AddHists(h_17_mu_fsrup2, h_18_mu_fsrup2);
    TH1F* h_1718_mu_fsrdown2 = AddHists(h_17_mu_fsrdown2, h_18_mu_fsrdown2);
    TH1F* h_1718_mu_fsrup4 = AddHists(h_17_mu_fsrup4, h_18_mu_fsrup4);
    // TH1F* h_1718_mu_fsrup4_cut = AddHists(h_17_mu_fsrup4_cut, h_18_mu_fsrup4_cut);
    TH1F* h_1718_mu_fsrdown4 = AddHists(h_17_mu_fsrdown4, h_18_mu_fsrdown4);
    TH1F* h_1718_mu_jecup = AddHists(h_17_mu_jecup, h_18_mu_jecup);
    TH1F* h_1718_mu_jecdown = AddHists(h_17_mu_jecdown, h_18_mu_jecdown);
    TH1F* h_1718_mu_corup = AddHists(h_17_mu_corup, h_18_mu_corup);
    TH1F* h_1718_mu_cordown = AddHists(h_17_mu_cordown, h_18_mu_cordown);
    TH1F* h_1718_mu_jmsup = AddHists(h_17_mu_jmsup, h_18_mu_jmsup);
    TH1F* h_1718_mu_jmsdown = AddHists(h_17_mu_jmsdown, h_18_mu_jmsdown);
    TH1F* h_1718_mu_st = AddHists(h_17_mu_st, h_18_mu_st);
    TH1F* h_1718_mu_wj = AddHists(h_17_mu_wj, h_18_mu_wj);
    TH1F* h_1718_mu_ot = AddHists(h_17_mu_ot, h_18_mu_ot);
    TH1F* h_1718_mu_data = AddHists(h_17_mu_data, h_18_mu_data);
    TH1F* h_1718_mu_hdampup = AddHists(h_17_mu_hdampup, h_18_mu_hdampup);
    TH1F* h_1718_mu_hdampdown = AddHists(h_17_mu_hdampdown, h_18_mu_hdampdown);
    TH1F* h_1718_mu_tuneup = AddHists(h_17_mu_tuneup, h_18_mu_tuneup);
    TH1F* h_1718_mu_tunedown = AddHists(h_17_mu_tunedown, h_18_mu_tunedown);
    TH1F* h_1718_mu_isrup = AddHists(h_17_mu_isrup, h_18_mu_isrup);
    TH1F* h_1718_mu_isrdown = AddHists(h_17_mu_isrdown, h_18_mu_isrdown);
    TH1F* h_1718_mu_gluon = AddHists(h_17_mu_gluon, h_18_mu_gluon);
    TH1F* h_1718_mu_qcd = AddHists(h_17_mu_qcd, h_18_mu_qcd);

    TH1F* h_1718_el_nominal = AddHists(h_17_el_nominal, h_18_el_nominal);
    TH1F* h_1718_el_fsrupsqrt2 = AddHists(h_17_el_fsrupsqrt2, h_18_el_fsrupsqrt2);
    TH1F* h_1718_el_fsrdownsqrt2 = AddHists(h_17_el_fsrdownsqrt2, h_18_el_fsrdownsqrt2);
    TH1F* h_1718_el_fsrup2 = AddHists(h_17_el_fsrup2, h_18_el_fsrup2);
    TH1F* h_1718_el_fsrdown2 = AddHists(h_17_el_fsrdown2, h_18_el_fsrdown2);
    TH1F* h_1718_el_fsrup4 = AddHists(h_17_el_fsrup4, h_18_el_fsrup4);
    // TH1F* h_1718_el_fsrup4_cut = AddHists(h_17_el_fsrup4_cut, h_18_el_fsrup4_cut);
    TH1F* h_1718_el_fsrdown4 = AddHists(h_17_el_fsrdown4, h_18_el_fsrdown4);
    TH1F* h_1718_el_jecup = AddHists(h_17_el_jecup, h_18_el_jecup);
    TH1F* h_1718_el_jecdown = AddHists(h_17_el_jecdown, h_18_el_jecdown);
    TH1F* h_1718_el_corup = AddHists(h_17_el_corup, h_18_el_corup);
    TH1F* h_1718_el_cordown = AddHists(h_17_el_cordown, h_18_el_cordown);
    TH1F* h_1718_el_jmsup = AddHists(h_17_el_jmsup, h_18_el_jmsup);
    TH1F* h_1718_el_jmsdown = AddHists(h_17_el_jmsdown, h_18_el_jmsdown);
    TH1F* h_1718_el_st = AddHists(h_17_el_st, h_18_el_st);
    TH1F* h_1718_el_wj = AddHists(h_17_el_wj, h_18_el_wj);
    TH1F* h_1718_el_ot = AddHists(h_17_el_ot, h_18_el_ot);
    TH1F* h_1718_el_data = AddHists(h_17_el_data, h_18_el_data);
    TH1F* h_1718_el_hdampup = AddHists(h_17_el_hdampup, h_18_el_hdampup);
    TH1F* h_1718_el_hdampdown = AddHists(h_17_el_hdampdown, h_18_el_hdampdown);
    TH1F* h_1718_el_tuneup = AddHists(h_17_el_tuneup, h_18_el_tuneup);
    TH1F* h_1718_el_tunedown = AddHists(h_17_el_tunedown, h_18_el_tunedown);
    TH1F* h_1718_el_isrup = AddHists(h_17_el_isrup, h_18_el_isrup);
    TH1F* h_1718_el_isrdown = AddHists(h_17_el_isrdown, h_18_el_isrdown);
    TH1F* h_1718_el_gluon = AddHists(h_17_el_gluon, h_18_el_gluon);
    TH1F* h_1718_el_qcd = AddHists(h_17_el_qcd, h_18_el_qcd);

    TH1F* h_1718_nominal = AddHists(h_17_nominal, h_18_nominal);
    TH1F* h_1718_fsrupsqrt2 = AddHists(h_17_fsrupsqrt2, h_18_fsrupsqrt2);
    TH1F* h_1718_fsrdownsqrt2 = AddHists(h_17_fsrdownsqrt2, h_18_fsrdownsqrt2);
    TH1F* h_1718_fsrup2 = AddHists(h_17_fsrup2, h_18_fsrup2);
    TH1F* h_1718_fsrdown2 = AddHists(h_17_fsrdown2, h_18_fsrdown2);
    TH1F* h_1718_fsrup4 = AddHists(h_17_fsrup4, h_18_fsrup4);
    // TH1F* h_1718_fsrup4_cut = AddHists(h_17_fsrup4_cut, h_18_fsrup4_cut);
    TH1F* h_1718_fsrdown4 = AddHists(h_17_fsrdown4, h_18_fsrdown4);
    TH1F* h_1718_jecup = AddHists(h_17_jecup, h_18_jecup);
    TH1F* h_1718_jecdown = AddHists(h_17_jecdown, h_18_jecdown);
    TH1F* h_1718_corup = AddHists(h_17_corup, h_18_corup);
    TH1F* h_1718_cordown = AddHists(h_17_cordown, h_18_cordown);
    TH1F* h_1718_jmsup = AddHists(h_17_jmsup, h_18_jmsup);
    TH1F* h_1718_jmsdown = AddHists(h_17_jmsdown, h_18_jmsdown);
    TH1F* h_1718_st = AddHists(h_17_st, h_18_st);
    TH1F* h_1718_wj = AddHists(h_17_wj, h_18_wj);
    TH1F* h_1718_ot = AddHists(h_17_ot, h_18_ot);
    TH1F* h_1718_data = AddHists(h_17_data, h_18_data);
    TH1F* h_1718_hdampup = AddHists(h_17_hdampup, h_18_hdampup);
    TH1F* h_1718_hdampdown = AddHists(h_17_hdampdown, h_18_hdampdown);
    TH1F* h_1718_tuneup = AddHists(h_17_tuneup, h_18_tuneup);
    TH1F* h_1718_tunedown = AddHists(h_17_tunedown, h_18_tunedown);
    TH1F* h_1718_isrup = AddHists(h_17_isrup, h_18_isrup);
    TH1F* h_1718_isrdown = AddHists(h_17_isrdown, h_18_isrdown);
    TH1F* h_1718_gluon = AddHists(h_17_gluon, h_18_gluon);
    TH1F* h_1718_qcd = AddHists(h_17_qcd, h_18_qcd);

    if(debug) continue; // Dont update TFile while debugging
    TString file = "files/FSR_hists_"+study+bin_string+".root";
    // TString file = removeWeights?"files/FSR_hists_"+study+bin_string+"_noweights.root":"files/FSR_hists_"+study+bin_string+"_weights.root";
    TString option = gSystem->AccessPathName(file)?"recreate":"update";
    TFile* f_out = new TFile(file,option);
    f_out->cd();
    TString add = isMjet?"_mjet":isPt?"_pt":"";
    //
    // ///////////////////////////////////////////////////
    // ///                    2016                     ///
    // ///////////////////////////////////////////////////
    h_16_mu_nominal->Write("nominal_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_fsrup2->Write("fsrup2_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_fsrdown2->Write("fsrdown2_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_jecup->Write("jecup_2016", TObject::kOverwrite);
    h_16_mu_jecdown->Write("jecdown_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_corup->Write("corup_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_cordown->Write("cordown_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_jmsup->Write("jmsup_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_jmsdown->Write("jmsdown_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_hdampup->Write("hdampup_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_hdampdown->Write("hdampdown_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_isrup->Write("isrup_muon2"+add+"_2016", TObject::kOverwrite);
    h_16_mu_isrdown->Write("isrdown_muon2"+add+"_2016", TObject::kOverwrite);
    h_16_mu_tuneup->Write("tuneup_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_tunedown->Write("tunedown_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_st->Write("bgr_st_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_wj->Write("bgr_wj_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_ot->Write("bgr_ot_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_data->Write("data_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_gluon->Write("gluonmove_muon"+add+"_2016", TObject::kOverwrite);
    h_16_mu_qcd->Write("qcdbased_muon"+add+"_2016", TObject::kOverwrite);

    h_16_el_nominal->Write("nominal_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_fsrup2->Write("fsrup2_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_fsrdown2->Write("fsrdown2_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_jecup->Write("jecup_2016", TObject::kOverwrite);
    h_16_el_jecdown->Write("jecdown_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_corup->Write("corup_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_cordown->Write("cordown_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_jmsup->Write("jmsup_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_jmsdown->Write("jmsdown_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_hdampup->Write("hdampup_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_hdampdown->Write("hdampdown_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_isrup->Write("isrup_elec2"+add+"_2016", TObject::kOverwrite);
    h_16_el_isrdown->Write("isrdown_elec2"+add+"_2016", TObject::kOverwrite);
    h_16_el_tuneup->Write("tuneup_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_tunedown->Write("tunedown_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_st->Write("bgr_st_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_wj->Write("bgr_wj_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_ot->Write("bgr_ot_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_data->Write("data_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_gluon->Write("gluonmove_elec"+add+"_2016", TObject::kOverwrite);
    h_16_el_qcd->Write("qcdbased_elec"+add+"_2016", TObject::kOverwrite);

    h_16_nominal->Write("nominal"+add+"_2016", TObject::kOverwrite);
    h_16_fsrup2->Write("fsrup2"+add+"_2016", TObject::kOverwrite);
    h_16_fsrdown2->Write("fsrdown2"+add+"_2016", TObject::kOverwrite);
    h_16_jecup->Write("jecup"+add+"_2016", TObject::kOverwrite);
    h_16_jecdown->Write("jecdown"+add+"_2016", TObject::kOverwrite);
    h_16_corup->Write("corup"+add+"_2016", TObject::kOverwrite);
    h_16_cordown->Write("cordown"+add+"_2016", TObject::kOverwrite);
    h_16_jmsup->Write("jmsup"+add+"_2016", TObject::kOverwrite);
    h_16_jmsdown->Write("jmsdown"+add+"_2016", TObject::kOverwrite);
    h_16_hdampup->Write("hdampup"+add+"_2016", TObject::kOverwrite);
    h_16_hdampdown->Write("hdampdown"+add+"_2016", TObject::kOverwrite);
    h_16_isrup->Write("isrup2"+add+"_2016", TObject::kOverwrite);
    h_16_isrdown->Write("isrdown2"+add+"_2016", TObject::kOverwrite);
    h_16_tuneup->Write("tuneup"+add+"_2016", TObject::kOverwrite);
    h_16_tunedown->Write("tunedown"+add+"_2016", TObject::kOverwrite);
    h_16_st->Write("bgr_st"+add+"_2016", TObject::kOverwrite);
    h_16_wj->Write("bgr_wj"+add+"_2016", TObject::kOverwrite);
    h_16_ot->Write("bgr_ot"+add+"_2016", TObject::kOverwrite);
    h_16_data->Write("data"+add+"_2016", TObject::kOverwrite);
    h_16_gluon->Write("gluonmove"+add+"_2016", TObject::kOverwrite);
    h_16_qcd->Write("qcdbased"+add+"_2016", TObject::kOverwrite);

    ///////////////////////////////////////////////////
    ///                    2017                     ///
    ///////////////////////////////////////////////////
    h_17_mu_nominal->Write("nominal_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_fsrupsqrt2->Write("fsrupsqrt2_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_fsrup2->Write("fsrup2_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_fsrup4->Write("fsrup4_muon"+add+"_2017", TObject::kOverwrite);
    // h_17_mu_fsrup4_cut->Write("fsrup4_muon"+add+"_cut_2017", TObject::kOverwrite);
    h_17_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_fsrdown2->Write("fsrdown2_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_fsrdown4->Write("fsrdown4_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_jecup->Write("jecup_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_jecdown->Write("jecdown_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_corup->Write("corup_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_cordown->Write("cordown_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_jmsup->Write("jmsup_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_jmsdown->Write("jmsdown_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_hdampup->Write("hdampup_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_hdampdown->Write("hdampdown_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_isrup->Write("isrup_muon2"+add+"_2017", TObject::kOverwrite);
    h_17_mu_isrdown->Write("isrdown_muon2"+add+"_2017", TObject::kOverwrite);
    h_17_mu_tuneup->Write("tuneup_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_tunedown->Write("tunedown_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_st->Write("bgr_st_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_wj->Write("bgr_wj_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_ot->Write("bgr_ot_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_data->Write("data_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_gluon->Write("gluonmove_muon"+add+"_2017", TObject::kOverwrite);
    h_17_mu_qcd->Write("qcdbased_muon"+add+"_2017", TObject::kOverwrite);

    h_17_el_nominal->Write("nominal_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_fsrupsqrt2->Write("fsrupsqrt2_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_fsrup2->Write("fsrup2_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_fsrup4->Write("fsrup4_elec"+add+"_2017", TObject::kOverwrite);
    // h_17_el_fsrup4_cut->Write("fsrup4_elec"+add+"_cut_2017", TObject::kOverwrite);
    h_17_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_fsrdown2->Write("fsrdown2_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_fsrdown4->Write("fsrdown4_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_jecup->Write("jecup_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_jecdown->Write("jecdown_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_corup->Write("corup_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_cordown->Write("cordown_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_jmsup->Write("jmsup_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_jmsdown->Write("jmsdown_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_hdampup->Write("hdampup_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_hdampdown->Write("hdampdown_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_isrup->Write("isrup_elec2"+add+"_2017", TObject::kOverwrite);
    h_17_el_isrdown->Write("isrdown_elec2"+add+"_2017", TObject::kOverwrite);
    h_17_el_tuneup->Write("tuneup_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_tunedown->Write("tunedown_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_st->Write("bgr_st_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_wj->Write("bgr_wj_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_ot->Write("bgr_ot_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_data->Write("data_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_gluon->Write("gluonmove_elec"+add+"_2017", TObject::kOverwrite);
    h_17_el_qcd->Write("qcdbased_elec"+add+"_2017", TObject::kOverwrite);

    h_17_nominal->Write("nominal"+add+"_2017", TObject::kOverwrite);
    h_17_fsrupsqrt2->Write("fsrupsqrt2"+add+"_2017", TObject::kOverwrite);
    h_17_fsrup2->Write("fsrup2"+add+"_2017", TObject::kOverwrite);
    h_17_fsrup4->Write("fsrup4"+add+"_2017", TObject::kOverwrite);
    // h_17_fsrup4_cut->Write("fsrup4"+add+"_cut_2017", TObject::kOverwrite);
    h_17_fsrdownsqrt2->Write("fsrdownsqrt2"+add+"_2017", TObject::kOverwrite);
    h_17_fsrdown2->Write("fsrdown2"+add+"_2017", TObject::kOverwrite);
    h_17_fsrdown4->Write("fsrdown4"+add+"_2017", TObject::kOverwrite);
    h_17_jecup->Write("jecup"+add+"_2017", TObject::kOverwrite);
    h_17_jecdown->Write("jecdown"+add+"_2017", TObject::kOverwrite);
    h_17_corup->Write("corup"+add+"_2017", TObject::kOverwrite);
    h_17_cordown->Write("cordown"+add+"_2017", TObject::kOverwrite);
    h_17_jmsup->Write("jmsup"+add+"_2017", TObject::kOverwrite);
    h_17_jmsdown->Write("jmsdown"+add+"_2017", TObject::kOverwrite);
    h_17_hdampup->Write("hdampup"+add+"_2017", TObject::kOverwrite);
    h_17_hdampdown->Write("hdampdown"+add+"_2017", TObject::kOverwrite);
    h_17_isrup->Write("isrup2"+add+"_2017", TObject::kOverwrite);
    h_17_isrdown->Write("isrdown2"+add+"_2017", TObject::kOverwrite);
    h_17_tuneup->Write("tuneup"+add+"_2017", TObject::kOverwrite);
    h_17_tunedown->Write("tunedown"+add+"_2017", TObject::kOverwrite);
    h_17_st->Write("bgr_st"+add+"_2017", TObject::kOverwrite);
    h_17_wj->Write("bgr_wj"+add+"_2017", TObject::kOverwrite);
    h_17_ot->Write("bgr_ot"+add+"_2017", TObject::kOverwrite);
    h_17_data->Write("data"+add+"_2017", TObject::kOverwrite);
    h_17_gluon->Write("gluonmove"+add+"_2017", TObject::kOverwrite);
    h_17_qcd->Write("qcdbased"+add+"_2017", TObject::kOverwrite);

    ///////////////////////////////////////////////////
    ///                    2018                     ///
    ///////////////////////////////////////////////////
    h_18_mu_nominal->Write("nominal_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_fsrupsqrt2->Write("fsrupsqrt2_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_fsrup2->Write("fsrup2_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_fsrup4->Write("fsrup4_muon"+add+"_2018", TObject::kOverwrite);
    // h_18_mu_fsrup4_cut->Write("fsrup4_muon"+add+"_cut_2018", TObject::kOverwrite);
    h_18_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_fsrdown2->Write("fsrdown2_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_fsrdown4->Write("fsrdown4_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_jecup->Write("jecup_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_jecdown->Write("jecdown_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_corup->Write("corup_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_cordown->Write("cordown_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_jmsup->Write("jmsup_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_jmsdown->Write("jmsdown_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_hdampup->Write("hdampup_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_hdampdown->Write("hdampdown_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_isrup->Write("isrup_muon2"+add+"_2018", TObject::kOverwrite);
    h_18_mu_isrdown->Write("isrdown_muon2"+add+"_2018", TObject::kOverwrite);
    h_18_mu_tuneup->Write("tuneup_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_tunedown->Write("tunedown_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_st->Write("bgr_st_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_wj->Write("bgr_wj_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_ot->Write("bgr_ot_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_data->Write("data_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_gluon->Write("gluonmove_muon"+add+"_2018", TObject::kOverwrite);
    h_18_mu_qcd->Write("qcdbased_muon"+add+"_2018", TObject::kOverwrite);

    h_18_el_nominal->Write("nominal_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_fsrupsqrt2->Write("fsrupsqrt2_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_fsrup2->Write("fsrup2_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_fsrup4->Write("fsrup4_elec"+add+"_2018", TObject::kOverwrite);
    // h_18_el_fsrup4_cut->Write("fsrup4_elec"+add+"_cut_2018", TObject::kOverwrite);
    h_18_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_fsrdown2->Write("fsrdown2_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_fsrdown4->Write("fsrdown4_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_jecup->Write("jecup_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_jecdown->Write("jecdown_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_corup->Write("corup_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_cordown->Write("cordown_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_jmsup->Write("jmsup_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_jmsdown->Write("jmsdown_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_hdampup->Write("hdampup_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_hdampdown->Write("hdampdown_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_isrup->Write("isrup_elec2"+add+"_2018", TObject::kOverwrite);
    h_18_el_isrdown->Write("isrdown_elec2"+add+"_2018", TObject::kOverwrite);
    h_18_el_tuneup->Write("tuneup_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_tunedown->Write("tunedown_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_st->Write("bgr_st_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_wj->Write("bgr_wj_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_ot->Write("bgr_ot_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_data->Write("data_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_gluon->Write("gluonmove_elec"+add+"_2018", TObject::kOverwrite);
    h_18_el_qcd->Write("qcdbased_elec"+add+"_2018", TObject::kOverwrite);

    h_18_nominal->Write("nominal"+add+"_2018", TObject::kOverwrite);
    h_18_fsrupsqrt2->Write("fsrupsqrt2"+add+"_2018", TObject::kOverwrite);
    h_18_fsrup2->Write("fsrup2"+add+"_2018", TObject::kOverwrite);
    h_18_fsrup4->Write("fsrup4"+add+"_2018", TObject::kOverwrite);
    // h_18_fsrup4_cut->Write("fsrup4"+add+"_cut_2018", TObject::kOverwrite);
    h_18_fsrdownsqrt2->Write("fsrdownsqrt2"+add+"_2018", TObject::kOverwrite);
    h_18_fsrdown2->Write("fsrdown2"+add+"_2018", TObject::kOverwrite);
    h_18_fsrdown4->Write("fsrdown4"+add+"_2018", TObject::kOverwrite);
    h_18_jecup->Write("jecup"+add+"_2018", TObject::kOverwrite);
    h_18_jecdown->Write("jecdown"+add+"_2018", TObject::kOverwrite);
    h_18_corup->Write("corup"+add+"_2018", TObject::kOverwrite);
    h_18_cordown->Write("cordown"+add+"_2018", TObject::kOverwrite);
    h_18_jmsup->Write("jmsup"+add+"_2018", TObject::kOverwrite);
    h_18_jmsdown->Write("jmsdown"+add+"_2018", TObject::kOverwrite);
    h_18_hdampup->Write("hdampup"+add+"_2018", TObject::kOverwrite);
    h_18_hdampdown->Write("hdampdown"+add+"_2018", TObject::kOverwrite);
    h_18_isrup->Write("isrup2"+add+"_2018", TObject::kOverwrite);
    h_18_isrdown->Write("isrdown2"+add+"_2018", TObject::kOverwrite);
    h_18_tuneup->Write("tuneup"+add+"_2018", TObject::kOverwrite);
    h_18_tunedown->Write("tunedown"+add+"_2018", TObject::kOverwrite);
    h_18_st->Write("bgr_st"+add+"_2018", TObject::kOverwrite);
    h_18_wj->Write("bgr_wj"+add+"_2018", TObject::kOverwrite);
    h_18_ot->Write("bgr_ot"+add+"_2018", TObject::kOverwrite);
    h_18_data->Write("data"+add+"_2018", TObject::kOverwrite);
    h_18_gluon->Write("gluonmove"+add+"_2018", TObject::kOverwrite);
    h_18_qcd->Write("qcdbased"+add+"_2018", TObject::kOverwrite);

    ///////////////////////////////////////////////////
    ///                    1718                     ///
    ///////////////////////////////////////////////////
    h_1718_mu_nominal->Write("nominal_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_fsrupsqrt2->Write("fsrupsqrt2_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_fsrup2->Write("fsrup2_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_fsrup4->Write("fsrup4_muon"+add+"_combine", TObject::kOverwrite);
    // h_1718_mu_fsrup4_cut->Write("fsrup4_muon"+add+"_cut_combine", TObject::kOverwrite);
    h_1718_mu_fsrdownsqrt2->Write("fsrdownsqrt2_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_fsrdown2->Write("fsrdown2_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_fsrdown4->Write("fsrdown4_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_jecup->Write("jecup_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_jecdown->Write("jecdown_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_corup->Write("corup_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_cordown->Write("cordown_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_jmsup->Write("jmsup_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_jmsdown->Write("jmsdown_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_hdampup->Write("hdampup_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_hdampdown->Write("hdampdown_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_isrup->Write("isrup_muon2"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_isrdown->Write("isrdown_muon2"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_tuneup->Write("tuneup_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_tunedown->Write("tunedown_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_st->Write("bgr_st_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_wj->Write("bgr_wj_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_ot->Write("bgr_ot_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_data->Write("data_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_gluon->Write("gluonmove_muon"+add+"_combine", TObject::kOverwrite);
    h_1718_mu_qcd->Write("qcdbased_muon"+add+"_combine", TObject::kOverwrite);

    h_1718_el_nominal->Write("nominal_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_fsrupsqrt2->Write("fsrupsqrt2_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_fsrup2->Write("fsrup2_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_fsrup4->Write("fsrup4_elec"+add+"_combine", TObject::kOverwrite);
    // h_1718_el_fsrup4_cut->Write("fsrup4_elec"+add+"_cut_combine", TObject::kOverwrite);
    h_1718_el_fsrdownsqrt2->Write("fsrdownsqrt2_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_fsrdown2->Write("fsrdown2_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_fsrdown4->Write("fsrdown4_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_jecup->Write("jecup_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_jecdown->Write("jecdown_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_corup->Write("corup_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_cordown->Write("cordown_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_jmsup->Write("jmsup_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_jmsdown->Write("jmsdown_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_hdampup->Write("hdampup_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_hdampdown->Write("hdampdown_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_isrup->Write("isrup_elec2"+add+"_combine", TObject::kOverwrite);
    h_1718_el_isrdown->Write("isrdown_elec2"+add+"_combine", TObject::kOverwrite);
    h_1718_el_tuneup->Write("tuneup_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_tunedown->Write("tunedown_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_st->Write("bgr_st_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_wj->Write("bgr_wj_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_ot->Write("bgr_ot_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_data->Write("data_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_gluon->Write("gluonmove_elec"+add+"_combine", TObject::kOverwrite);
    h_1718_el_qcd->Write("qcdbased_elec"+add+"_combine", TObject::kOverwrite);

    h_1718_nominal->Write("nominal"+add+"_combine", TObject::kOverwrite);
    h_1718_fsrupsqrt2->Write("fsrupsqrt2"+add+"_combine", TObject::kOverwrite);
    h_1718_fsrup2->Write("fsrup2"+add+"_combine", TObject::kOverwrite);
    h_1718_fsrup4->Write("fsrup4"+add+"_combine", TObject::kOverwrite);
    // h_1718_fsrup4_cut->Write("fsrup4"+add+"_cut_combine", TObject::kOverwrite);
    h_1718_fsrdownsqrt2->Write("fsrdownsqrt2"+add+"_combine", TObject::kOverwrite);
    h_1718_fsrdown2->Write("fsrdown2"+add+"_combine", TObject::kOverwrite);
    h_1718_fsrdown4->Write("fsrdown4"+add+"_combine", TObject::kOverwrite);
    h_1718_jecup->Write("jecup"+add+"_combine", TObject::kOverwrite);
    h_1718_jecdown->Write("jecdown"+add+"_combine", TObject::kOverwrite);
    h_1718_corup->Write("corup"+add+"_combine", TObject::kOverwrite);
    h_1718_cordown->Write("cordown"+add+"_combine", TObject::kOverwrite);
    h_1718_jmsup->Write("jmsup"+add+"_combine", TObject::kOverwrite);
    h_1718_jmsdown->Write("jmsdown"+add+"_combine", TObject::kOverwrite);
    h_1718_hdampup->Write("hdampup"+add+"_combine", TObject::kOverwrite);
    h_1718_hdampdown->Write("hdampdown"+add+"_combine", TObject::kOverwrite);
    h_1718_isrup->Write("isrup2"+add+"_combine", TObject::kOverwrite);
    h_1718_isrdown->Write("isrdown2"+add+"_combine", TObject::kOverwrite);
    h_1718_tuneup->Write("tuneup"+add+"_combine", TObject::kOverwrite);
    h_1718_tunedown->Write("tunedown"+add+"_combine", TObject::kOverwrite);
    h_1718_st->Write("bgr_st"+add+"_combine", TObject::kOverwrite);
    h_1718_wj->Write("bgr_wj"+add+"_combine", TObject::kOverwrite);
    h_1718_ot->Write("bgr_ot"+add+"_combine", TObject::kOverwrite);
    h_1718_data->Write("data"+add+"_combine", TObject::kOverwrite);
    h_1718_gluon->Write("gluonmove"+add+"_combine", TObject::kOverwrite);
    h_1718_qcd->Write("qcdbased"+add+"_combine", TObject::kOverwrite);

    f_out->Close();

    auto stop  = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "  - time needed: " << GREEN << duration.count()/1000 << "s" << RESET << endl;
  }
  // return 0;
}

// -------------------------------------------------------------------------------------------------------
// ---                                                                                                 ---
// -------------------------------------------------------------------------------------------------------

TH1F* get_hist(TFile* file, vector<double> cut, TString weightname, TString var, bool weight_cut){

  bool plotTau = var=="tau32";
  bool plotMjet = var=="mjet";
  bool plotPt = var=="pt";

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  int xbins = plotTau?100:plotMjet?400:plotPt?1000:0;
  int xmin = plotTau?0:plotMjet?0:plotPt?0:0;
  int xmax = plotTau?1:plotMjet?500:plotPt?1000:0;
  TH1F* h_hist = new TH1F("hist", "hist", xbins, xmin, xmax);

  Float_t tau32, mjet, ptjet;
  Bool_t passed_measurement_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  map<TString, Float_t> cutvar, fillvar;
  cutvar["tau32"] = 0; cutvar["mjet"] = 0; cutvar["pt"] = 0;
  fillvar["tau32"] = 0; fillvar["mjet"] = 0; fillvar["pt"] = 0;
  if(debug && weightname == "none") return h_hist;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("tau32_ak8_had",&tau32);
  tree->SetBranchAddress("mass_ak8_had",&mjet);
  tree->SetBranchAddress("pt_ak8_had",&ptjet);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  if(weightname != "none"){
    cout << "  - using weight" << endl;
    tree->SetBranchAddress(weightname,&additional_factor);
  }
  tree->SetBranchStatus("*",1);

  int skipEvents = 0;
  // cout << "Number Entries " << tree->GetEntriesFast() << endl;
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(!passed_measurement_rec) continue;

    // Here, only rec weight is changed!
    if(weightname != "none") rec_weight *= additional_factor;

    // get weights for migration matrix
    weight = rec_weight * gen_weight;

    double ps_weight = gen_weight*additional_factor;
    if(ps_weight>300){
      cout << "gen weight greater 300 with " << ps_weight << "\t gen " << gen_weight  << "\t rec " << gen_weight << endl;
    }

    fillvar["mjet"] = mjet; fillvar["pt"] = ptjet; fillvar["tau32"] = tau32;
    cutvar["mjet"] = mjet; cutvar["pt"] = ptjet;

    if(studyPt && mjet<140) continue;
    if(fillvar[var] < 0 || cutvar[study] < cut[0] || cutvar[study] > cut[1]) continue;
    else {
      // if(weight_cut && weight>20){
      //   cout << "remove Bins with weight > 20" << endl;
      //   skipEvents++;
      //   if(debug) cout << "remove Bins with weight > 20 | " << weight << endl;
      //   continue;
      // }
      h_hist->Fill(fillvar[var], weight);
    }

  }
  // cout << h_hist->Integral() << endl;
  if(debug) cout << "\tSkip " << skipEvents << " Events" << endl;
  return h_hist;
}

// -------------------------------------------------------------------------------------------------------
// ---                                                                                                 ---
// -------------------------------------------------------------------------------------------------------

TH1F* AddHists(TH1F* h1, TH1F* h2){
  TH1F* hist_added = (TH1F*) h1->Clone();
  hist_added->Add(h2);
  return hist_added;
}

// -------------------------------------------------------------------------------------------------------
// ---                                                                                                 ---
// -------------------------------------------------------------------------------------------------------

VecDD CutValuesFromWindow(vector<TString> windows){
  VecDD cutvalues = {};
  string delimiter = "to";
  for(TString w: windows){
    string ws = (string) w;
    int pos = ws.find(delimiter);
    TString first = ws.substr(0, ws.find(delimiter));
    TString last = ws.erase(0, pos + delimiter.length());
    if(last=="Inf") last = "100000.";

    cutvalues.push_back({first.Atof(), last.Atof()});
  }
  return cutvalues;
}
