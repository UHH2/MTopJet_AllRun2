#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/CovarianzMatrix.h"

#include <sys/types.h>
#include <sys/stat.h>
#include "TSystem.h"

bool forTalk = false;
bool debug   = false;

using namespace std;

// -------------------------------------------------------------------------------------------------------
// ---------------------------------------- FUNCTIONS ----------------------------------------------------
// -------------------------------------------------------------------------------------------------------

void WriteHistogramToRoot(TH1F* hist, TString name, int xmin, int xmax){
  int nBins = xmax-xmin;
  // cout << name << setw(4) << nBins << setw(4) << xmin << setw(4) << xmax << endl;
  TH1F* wmass = new TH1F(name, name, nBins, xmin, xmax);
  for(unsigned int i=1; i<=wmass->GetNbinsX(); i++){
    double bin_center  = wmass->GetBinCenter(i);
    double new_content = hist->GetBinContent((int) bin_center+1);
    double new_error   = hist->GetBinError((int) bin_center+1);
    // if(i==1) cout << setw(4) << bin_center << setw(4) << (int) bin_center+1 << endl;
    wmass->SetBinContent(i, new_content);
    wmass->SetBinError(i, new_error);
  }
  wmass->Write(name, TObject::kOverwrite);
}

// -------------------------------------------------------------------------------------------------------
// ------------------------------------------ MAIN -------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]){

  // #################################################################################################
  // Settings ########################################################################################

  // Declare different variables used in the code ------------------------------
  print_seperater();

  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_SYS_MassPlot <rebin_number>\n";
    cout << "hist->Rebin(rebin_number)\n";
    return 0;
  }

  // Default -------------------------------------------------------------------
  cout.precision(6);
  Int_t   oldLevel     = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ...
  gErrorIgnoreLevel    = kError;          // suppress TCanvas output
  gStyle->SetOptTitle(0);
  SetupGlobalStyle();

  TString reconst      = "btag";            // match btag_cut btag_sel compare min_mass
  TString MC_uncert    = "central";         // central avg nominal
  bool    usePeak_in   = false;
  bool    one_Ptbin    = false;
  int bin_width = 1;

  // Input ---------------------------------------------------------------------
  // int bin_width   = stoi(argv[1]);

  // Directories ---------------------------------------------------------------

  TString save_path = get_save_path();
  save_path += "/JetCorrections";
  save_path = creat_folder_and_path(save_path, "MassPlots");
  // if reconst is not btag include another folder here
  save_path = creat_folder_and_path(save_path, "BinWidth_"+to_string(bin_width));
  cout << save_path << endl;
  creat_folder(save_path+"/muon");
  creat_folder(save_path+"/elec");
  creat_folder(save_path+"/combine");

  // Directories ---------------------------------------------------------------
  vector<int> index_year   = {0,1,2,3};
  vector<TString> years    = {"2016","2017","2018","combine"};
  vector<TString> channels = {"muon","elec","combine"};

  // -------------------------------------------------------------------------------------------------------
  // ------------------------------------------ MAIN -------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------

  map<TString, vector<int>> range;
  range["comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hh"] = {70, 105};
  range["comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hl"] = {75, 104};
  range["comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_lh"] = {62, 98};
  range["comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_ll"] = {63, 101};

  vector<double> means_nominal = {};
  vector<double> means_jec_up = {};
  vector<double> means_cor_up = {};
  vector<double> means_jec_down = {};
  vector<double> means_cor_down = {};

  for(int ptbin=0; ptbin<5; ptbin++){

    // Define ptbins -----------------------------------------------------------
    if((ptbin==4)) continue;
    TString addition="";
    addition="_hh";
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="";
    string add = (string) addition;

    cout << endl << "Start with bin " << addition << " ..." << endl;

    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
    // if(reconst=="btag")      w_mass = hist_class+"wmass_match";
    if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(reconst=="btag"&&ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(reconst=="btag"&&ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";
    if(reconst=="btag"&&ptbin==4) w_mass = hist_class+"wmass_match";
    // if(reconst=="btag")      w_mass = hist_class+"wjet_pt_match_S1divW_W";

    // #################################################################################################
    // Get Hists #######################################################################################
    cout << "Get Hists ..." << endl;

    vector<TH1F*> h_data_muon          = get_all_hists(data_muon, w_mass);
    vector<TH1F*> h_data_elec          = get_all_hists(data_elec, w_mass);
    vector<TH1F*> h_data_combine       = combine_channels(h_data_muon, h_data_elec);

    vector<TH1F*> h_ttbar_muon         = get_all_hists(ttbar_muon, w_mass);
    vector<TH1F*> h_ttbar_elec         = get_all_hists(ttbar_elec, w_mass);
    vector<TH1F*> h_ttbar_combine      = combine_channels(h_ttbar_muon, h_ttbar_elec);

    vector<TH1F*> h_st_muon            = get_all_hists(st_muon, w_mass);
    vector<TH1F*> h_st_elec            = get_all_hists(st_elec, w_mass);
    vector<TH1F*> h_st_combine         = combine_channels(h_st_muon, h_st_elec);

    vector<TH1F*> h_wjets_muon         = get_all_hists(wjets_muon, w_mass);
    vector<TH1F*> h_wjets_elec         = get_all_hists(wjets_elec, w_mass);
    vector<TH1F*> h_wjets_combine      = combine_channels(h_wjets_muon, h_wjets_elec);

    vector<TH1F*> h_other_muon         = get_all_hists(other_muon, w_mass);
    vector<TH1F*> h_other_elec         = get_all_hists(other_elec, w_mass);
    vector<TH1F*> h_other_combine      = combine_channels(h_other_muon, h_other_elec);

    // jec -----------------------------------------------------------------------
    vector<TH1F*> h_jec_up_muon        = get_all_hists(jec_up_muon, w_mass);
    vector<TH1F*> h_jec_up_elec        = get_all_hists(jec_up_elec, w_mass);
    vector<TH1F*> h_jec_up_combine     = combine_channels(h_jec_up_muon, h_jec_up_elec);

    vector<TH1F*> h_jec_down_muon      = get_all_hists(jec_down_muon, w_mass);
    vector<TH1F*> h_jec_down_elec      = get_all_hists(jec_down_elec, w_mass);
    vector<TH1F*> h_jec_down_combine   = combine_channels(h_jec_down_muon, h_jec_down_elec);

    // cor -----------------------------------------------------------------------
    vector<TH1F*> h_cor_up_muon        = get_all_hists(cor_up_muon, w_mass);
    vector<TH1F*> h_cor_up_elec        = get_all_hists(cor_up_elec, w_mass);
    vector<TH1F*> h_cor_up_combine     = combine_channels(h_cor_up_muon, h_cor_up_elec);

    vector<TH1F*> h_cor_down_muon      = get_all_hists(cor_down_muon, w_mass);
    vector<TH1F*> h_cor_down_elec      = get_all_hists(cor_down_elec, w_mass);
    vector<TH1F*> h_cor_down_combine   = combine_channels(h_cor_down_muon, h_cor_down_elec);

    // jer -----------------------------------------------------------------------
    vector<TH1F*> h_jer_up_muon        = get_all_hists(jer_up_muon, w_mass);
    vector<TH1F*> h_jer_up_elec        = get_all_hists(jer_up_elec, w_mass);
    vector<TH1F*> h_jer_up_combine     = combine_channels(h_jer_up_muon, h_jer_up_elec);

    vector<TH1F*> h_jer_down_muon      = get_all_hists(jer_down_muon, w_mass);
    vector<TH1F*> h_jer_down_elec      = get_all_hists(jer_down_elec, w_mass);
    vector<TH1F*> h_jer_down_combine   = combine_channels(h_jer_down_muon, h_jer_down_elec);

    // FSR -----------------------------------------------------------------------
    vector<TH1F*> h_fsr_up_sqrt2_muon    = get_all_hists(fsr_up_sqrt2_muon, w_mass);
    vector<TH1F*> h_fsr_up_sqrt2_elec    = get_all_hists(fsr_up_sqrt2_elec, w_mass);
    vector<TH1F*> h_fsr_up_sqrt2_combine = combine_channels(h_fsr_up_sqrt2_muon, h_fsr_up_sqrt2_elec);

    vector<TH1F*> h_fsr_down_sqrt2_muon      = get_all_hists(fsr_down_sqrt2_muon, w_mass);
    vector<TH1F*> h_fsr_down_sqrt2_elec      = get_all_hists(fsr_down_sqrt2_elec, w_mass);
    vector<TH1F*> h_fsr_down_sqrt2_combine   = combine_channels(h_fsr_down_sqrt2_muon, h_fsr_down_sqrt2_elec);

    vector<TH1F*> h_fsr_up_2_muon        = get_all_hists(fsr_up_2_muon, w_mass);
    vector<TH1F*> h_fsr_up_2_elec        = get_all_hists(fsr_up_2_elec, w_mass);
    vector<TH1F*> h_fsr_up_2_combine     = combine_channels(h_fsr_up_2_muon, h_fsr_up_2_elec);

    vector<TH1F*> h_fsr_down_2_muon      = get_all_hists(fsr_down_2_muon, w_mass);
    vector<TH1F*> h_fsr_down_2_elec      = get_all_hists(fsr_down_2_elec, w_mass);
    vector<TH1F*> h_fsr_down_2_combine   = combine_channels(h_fsr_down_2_muon, h_fsr_down_2_elec);

    vector<TH1F*> h_fsr_up_4_muon        = get_all_hists(fsr_up_4_muon, w_mass);
    vector<TH1F*> h_fsr_up_4_elec        = get_all_hists(fsr_up_4_elec, w_mass);
    vector<TH1F*> h_fsr_up_4_combine     = combine_channels(h_fsr_up_4_muon, h_fsr_up_4_elec);

    vector<TH1F*> h_fsr_down_4_muon      = get_all_hists(fsr_down_4_muon, w_mass);
    vector<TH1F*> h_fsr_down_4_elec      = get_all_hists(fsr_down_4_elec, w_mass);
    vector<TH1F*> h_fsr_down_4_combine   = combine_channels(h_fsr_down_4_muon, h_fsr_down_4_elec);

    // hdamp -----------------------------------------------------------------------
    vector<TH1F*> h_hdamp_up_muon        = get_all_hists(hdamp_up_muon, w_mass);
    vector<TH1F*> h_hdamp_up_elec        = get_all_hists(hdamp_up_elec, w_mass);
    vector<TH1F*> h_hdamp_up_combine     = combine_channels(h_hdamp_up_muon, h_hdamp_up_elec);

    vector<TH1F*> h_hdamp_down_muon      = get_all_hists(hdamp_down_muon, w_mass);
    vector<TH1F*> h_hdamp_down_elec      = get_all_hists(hdamp_down_elec, w_mass);
    vector<TH1F*> h_hdamp_down_combine   = combine_channels(h_hdamp_down_muon, h_hdamp_down_elec);

    // tune -----------------------------------------------------------------------
    vector<TH1F*> h_tune_up_muon        = get_all_hists(tune_up_muon, w_mass);
    vector<TH1F*> h_tune_up_elec        = get_all_hists(tune_up_elec, w_mass);
    vector<TH1F*> h_tune_up_combine     = combine_channels(h_tune_up_muon, h_tune_up_elec);

    vector<TH1F*> h_tune_down_muon      = get_all_hists(tune_down_muon, w_mass);
    vector<TH1F*> h_tune_down_elec      = get_all_hists(tune_down_elec, w_mass);
    vector<TH1F*> h_tune_down_combine   = combine_channels(h_tune_down_muon, h_tune_down_elec);

    // CR -------------------------------------------------------------------------
    vector<TH1F*> h_cr_qb_muon        = get_all_hists(cr_qb_muon, w_mass);
    vector<TH1F*> h_cr_qb_elec        = get_all_hists(cr_qb_elec, w_mass);
    vector<TH1F*> h_cr_qb_combine     = combine_channels(h_cr_qb_muon, h_cr_qb_elec);

    vector<TH1F*> h_cr_gm_muon      = get_all_hists(cr_gm_muon, w_mass);
    vector<TH1F*> h_cr_gm_elec      = get_all_hists(cr_gm_elec, w_mass);
    vector<TH1F*> h_cr_gm_combine   = combine_channels(h_cr_gm_muon, h_cr_gm_elec);

    // #################################################################################################
    // Combine Hists ###################################################################################
    cout << "Combine Hists ..." << endl;

    // Bkg -----------------------------------------------------------------------
    vector<TH1F*> h_bkg_muon    = AddHists(h_st_muon, h_wjets_muon, h_other_muon, 1);
    vector<TH1F*> h_bkg_elec    = AddHists(h_st_elec, h_wjets_elec, h_other_elec, 1);
    vector<TH1F*> h_bkg_combine = AddHists(h_st_combine, h_wjets_combine, h_other_combine, 1);

    // Data ----------------------------------------------------------------------
    if(debug) cout << "\t ... data - bkg" << endl;

    vector<TH1F*> all_data_muon              = AddHists(h_data_muon, h_bkg_muon, -1);
    vector<TH1F*> all_data_elec              = AddHists(h_data_elec, h_bkg_elec, -1);
    vector<TH1F*> all_data_combine           = AddHists(h_data_combine, h_bkg_combine, -1);

    vector<TH1F*> all_data_muon_rebin        = RebinVector(all_data_muon, bin_width);
    vector<TH1F*> all_data_elec_rebin        = RebinVector(all_data_elec, bin_width);
    vector<TH1F*> all_data_combine_rebin     = RebinVector(all_data_combine, bin_width);

    vector<TH1F*> all_data_muon_norm, all_data_elec_norm, all_data_combine_norm;
    for(auto hist: all_data_muon_rebin) all_data_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_data_elec_rebin) all_data_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_data_combine_rebin) all_data_combine_norm.push_back(Normalize(hist));



    // TTbar ---------------------------------------------------------------------
    if(debug) cout << "\t ... ttbar" << endl;

    vector<TH1F*> all_ttbar_muon             = h_ttbar_muon;
    vector<TH1F*> all_ttbar_elec             = h_ttbar_elec;
    vector<TH1F*> all_ttbar_combine          = h_ttbar_combine;

    vector<TH1F*> all_ttbar_muon_rebin       = RebinVector(all_ttbar_muon, bin_width);
    vector<TH1F*> all_ttbar_elec_rebin       = RebinVector(all_ttbar_elec, bin_width);
    vector<TH1F*> all_ttbar_combine_rebin    = RebinVector(all_ttbar_combine, bin_width);

    vector<TH1F*> all_ttbar_muon_norm, all_ttbar_elec_norm, all_ttbar_combine_norm;
    for(auto hist: all_ttbar_muon_rebin) all_ttbar_muon_norm.push_back(Normalize(hist)); // uncorrelated error
    for(auto hist: all_ttbar_elec_rebin) all_ttbar_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_ttbar_combine_rebin) all_ttbar_combine_norm.push_back(Normalize(hist));

    // JEC -----------------------------------------------------------------------
    if(debug) cout << "\t ... JEC" << endl;

    vector<TH1F*> all_jec_up_muon            = h_jec_up_muon;
    vector<TH1F*> all_jec_up_elec            = h_jec_up_elec;
    vector<TH1F*> all_jec_up_combine         = h_jec_up_combine;

    vector<TH1F*> all_jec_down_muon          = h_jec_down_muon;
    vector<TH1F*> all_jec_down_elec          = h_jec_down_elec;
    vector<TH1F*> all_jec_down_combine       = h_jec_down_combine;

    vector<TH1F*> all_jec_down_muon_rebin    = RebinVector(all_jec_down_muon, bin_width);
    vector<TH1F*> all_jec_down_elec_rebin    = RebinVector(all_jec_down_elec, bin_width);
    vector<TH1F*> all_jec_down_combine_rebin = RebinVector(all_jec_down_combine, bin_width);

    vector<TH1F*> all_jec_up_muon_rebin      = RebinVector(all_jec_up_muon, bin_width);
    vector<TH1F*> all_jec_up_elec_rebin      = RebinVector(all_jec_up_elec, bin_width);
    vector<TH1F*> all_jec_up_combine_rebin   = RebinVector(all_jec_up_combine, bin_width);

    vector<TH1F*> all_jec_down_muon_norm, all_jec_down_elec_norm, all_jec_down_combine_norm;
    for(auto hist: all_jec_down_muon_rebin) all_jec_down_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_jec_down_elec_rebin) all_jec_down_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_jec_down_combine_rebin) all_jec_down_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_jec_up_muon_norm, all_jec_up_elec_norm, all_jec_up_combine_norm;
    for(auto hist: all_jec_up_muon_rebin) all_jec_up_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_jec_up_elec_rebin) all_jec_up_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_jec_up_combine_rebin) all_jec_up_combine_norm.push_back(Normalize(hist));

    // COR -----------------------------------------------------------------------
    if(debug) cout << "\t ... XCone" << endl;

    vector<TH1F*> all_cor_up_muon            = h_cor_up_muon;
    vector<TH1F*> all_cor_up_elec            = h_cor_up_elec;
    vector<TH1F*> all_cor_up_combine         = h_cor_up_combine;

    vector<TH1F*> all_cor_down_muon          = h_cor_down_muon;
    vector<TH1F*> all_cor_down_elec          = h_cor_down_elec;
    vector<TH1F*> all_cor_down_combine       = h_cor_down_combine;

    vector<TH1F*> all_cor_up_muon_rebin      = RebinVector(all_cor_up_muon, bin_width);
    vector<TH1F*> all_cor_up_elec_rebin      = RebinVector(all_cor_up_elec, bin_width);
    vector<TH1F*> all_cor_up_combine_rebin   = RebinVector(all_cor_up_combine, bin_width);

    vector<TH1F*> all_cor_down_muon_rebin    = RebinVector(all_cor_down_muon, bin_width);
    vector<TH1F*> all_cor_down_elec_rebin    = RebinVector(all_cor_down_elec, bin_width);
    vector<TH1F*> all_cor_down_combine_rebin = RebinVector(all_cor_down_combine, bin_width);

    vector<TH1F*> all_cor_down_muon_norm, all_cor_down_elec_norm, all_cor_down_combine_norm;
    for(auto hist: all_cor_down_muon_rebin) all_cor_down_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_cor_down_elec_rebin) all_cor_down_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_cor_down_combine_rebin) all_cor_down_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_cor_up_muon_norm, all_cor_up_elec_norm, all_cor_up_combine_norm;
    for(auto hist: all_cor_up_muon_rebin) all_cor_up_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_cor_up_elec_rebin) all_cor_up_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_cor_up_combine_rebin) all_cor_up_combine_norm.push_back(Normalize(hist));

    // JEC -----------------------------------------------------------------------
    if(debug) cout << "\t ... JEC" << endl;

    vector<TH1F*> all_jer_up_muon            = h_jer_up_muon;
    vector<TH1F*> all_jer_up_elec            = h_jer_up_elec;
    vector<TH1F*> all_jer_up_combine         = h_jer_up_combine;

    vector<TH1F*> all_jer_down_muon          = h_jer_down_muon;
    vector<TH1F*> all_jer_down_elec          = h_jer_down_elec;
    vector<TH1F*> all_jer_down_combine       = h_jer_down_combine;

    vector<TH1F*> all_jer_down_muon_rebin    = RebinVector(all_jer_down_muon, bin_width);
    vector<TH1F*> all_jer_down_elec_rebin    = RebinVector(all_jer_down_elec, bin_width);
    vector<TH1F*> all_jer_down_combine_rebin = RebinVector(all_jer_down_combine, bin_width);

    vector<TH1F*> all_jer_up_muon_rebin      = RebinVector(all_jer_up_muon, bin_width);
    vector<TH1F*> all_jer_up_elec_rebin      = RebinVector(all_jer_up_elec, bin_width);
    vector<TH1F*> all_jer_up_combine_rebin   = RebinVector(all_jer_up_combine, bin_width);

    vector<TH1F*> all_jer_down_muon_norm, all_jer_down_elec_norm, all_jer_down_combine_norm;
    for(auto hist: all_jer_down_muon_rebin) all_jer_down_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_jer_down_elec_rebin) all_jer_down_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_jer_down_combine_rebin) all_jer_down_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_jer_up_muon_norm, all_jer_up_elec_norm, all_jer_up_combine_norm;
    for(auto hist: all_jer_up_muon_rebin) all_jer_up_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_jer_up_elec_rebin) all_jer_up_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_jer_up_combine_rebin) all_jer_up_combine_norm.push_back(Normalize(hist));

    // hdamp -----------------------------------------------------------------------
    if(debug) cout << "\t ... hdamp" << endl;

    vector<TH1F*> all_hdamp_up_muon            = h_hdamp_up_muon;
    vector<TH1F*> all_hdamp_up_elec            = h_hdamp_up_elec;
    vector<TH1F*> all_hdamp_up_combine         = h_hdamp_up_combine;

    vector<TH1F*> all_hdamp_down_muon          = h_hdamp_down_muon;
    vector<TH1F*> all_hdamp_down_elec          = h_hdamp_down_elec;
    vector<TH1F*> all_hdamp_down_combine       = h_hdamp_down_combine;

    vector<TH1F*> all_hdamp_down_muon_rebin    = RebinVector(all_hdamp_down_muon, bin_width);
    vector<TH1F*> all_hdamp_down_elec_rebin    = RebinVector(all_hdamp_down_elec, bin_width);
    vector<TH1F*> all_hdamp_down_combine_rebin = RebinVector(all_hdamp_down_combine, bin_width);

    vector<TH1F*> all_hdamp_up_muon_rebin      = RebinVector(all_hdamp_up_muon, bin_width);
    vector<TH1F*> all_hdamp_up_elec_rebin      = RebinVector(all_hdamp_up_elec, bin_width);
    vector<TH1F*> all_hdamp_up_combine_rebin   = RebinVector(all_hdamp_up_combine, bin_width);

    vector<TH1F*> all_hdamp_down_muon_norm, all_hdamp_down_elec_norm, all_hdamp_down_combine_norm;
    for(auto hist: all_hdamp_down_muon_rebin) all_hdamp_down_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_hdamp_down_elec_rebin) all_hdamp_down_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_hdamp_down_combine_rebin) all_hdamp_down_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_hdamp_up_muon_norm, all_hdamp_up_elec_norm, all_hdamp_up_combine_norm;
    for(auto hist: all_hdamp_up_muon_rebin) all_hdamp_up_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_hdamp_up_elec_rebin) all_hdamp_up_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_hdamp_up_combine_rebin) all_hdamp_up_combine_norm.push_back(Normalize(hist));

    // CR --------------------------------------------------------------------------
    if(debug) cout << "\t ... hdamp" << endl;

    vector<TH1F*> all_cr_qb_muon            = h_cr_qb_muon;
    vector<TH1F*> all_cr_qb_elec            = h_cr_qb_elec;
    vector<TH1F*> all_cr_qb_combine         = h_cr_qb_combine;

    vector<TH1F*> all_cr_gm_muon          = h_cr_gm_muon;
    vector<TH1F*> all_cr_gm_elec          = h_cr_gm_elec;
    vector<TH1F*> all_cr_gm_combine       = h_cr_gm_combine;

    vector<TH1F*> all_cr_gm_muon_rebin    = RebinVector(all_cr_gm_muon, bin_width);
    vector<TH1F*> all_cr_gm_elec_rebin    = RebinVector(all_cr_gm_elec, bin_width);
    vector<TH1F*> all_cr_gm_combine_rebin = RebinVector(all_cr_gm_combine, bin_width);

    vector<TH1F*> all_cr_qb_muon_rebin      = RebinVector(all_cr_qb_muon, bin_width);
    vector<TH1F*> all_cr_qb_elec_rebin      = RebinVector(all_cr_qb_elec, bin_width);
    vector<TH1F*> all_cr_qb_combine_rebin   = RebinVector(all_cr_qb_combine, bin_width);

    vector<TH1F*> all_cr_gm_muon_norm, all_cr_gm_elec_norm, all_cr_gm_combine_norm;
    for(auto hist: all_cr_gm_muon_rebin) all_cr_gm_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_cr_gm_elec_rebin) all_cr_gm_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_cr_gm_combine_rebin) all_cr_gm_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_cr_qb_muon_norm, all_cr_qb_elec_norm, all_cr_qb_combine_norm;
    for(auto hist: all_cr_qb_muon_rebin) all_cr_qb_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_cr_qb_elec_rebin) all_cr_qb_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_cr_qb_combine_rebin) all_cr_qb_combine_norm.push_back(Normalize(hist));

    // tune -----------------------------------------------------------------------
    if(debug) cout << "\t ... tune" << endl;

    vector<TH1F*> all_tune_up_muon            = h_tune_up_muon;
    vector<TH1F*> all_tune_up_elec            = h_tune_up_elec;
    vector<TH1F*> all_tune_up_combine         = h_tune_up_combine;

    vector<TH1F*> all_tune_down_muon          = h_tune_down_muon;
    vector<TH1F*> all_tune_down_elec          = h_tune_down_elec;
    vector<TH1F*> all_tune_down_combine       = h_tune_down_combine;

    vector<TH1F*> all_tune_down_muon_rebin    = RebinVector(all_tune_down_muon, bin_width);
    vector<TH1F*> all_tune_down_elec_rebin    = RebinVector(all_tune_down_elec, bin_width);
    vector<TH1F*> all_tune_down_combine_rebin = RebinVector(all_tune_down_combine, bin_width);

    vector<TH1F*> all_tune_up_muon_rebin      = RebinVector(all_tune_up_muon, bin_width);
    vector<TH1F*> all_tune_up_elec_rebin      = RebinVector(all_tune_up_elec, bin_width);
    vector<TH1F*> all_tune_up_combine_rebin   = RebinVector(all_tune_up_combine, bin_width);

    vector<TH1F*> all_tune_down_muon_norm, all_tune_down_elec_norm, all_tune_down_combine_norm;
    for(auto hist: all_tune_down_muon_rebin) all_tune_down_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_tune_down_elec_rebin) all_tune_down_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_tune_down_combine_rebin) all_tune_down_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_tune_up_muon_norm, all_tune_up_elec_norm, all_tune_up_combine_norm;
    for(auto hist: all_tune_up_muon_rebin) all_tune_up_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_tune_up_elec_rebin) all_tune_up_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_tune_up_combine_rebin) all_tune_up_combine_norm.push_back(Normalize(hist));

    // fsr -----------------------------------------------------------------------
    if(debug) cout << "\t ... fsr" << endl;

    vector<TH1F*> all_fsr_up_sqrt2_muon            = h_fsr_up_sqrt2_muon;
    vector<TH1F*> all_fsr_up_sqrt2_elec            = h_fsr_up_sqrt2_elec;
    vector<TH1F*> all_fsr_up_sqrt2_combine         = h_fsr_up_sqrt2_combine;

    vector<TH1F*> all_fsr_down_sqrt2_muon          = h_fsr_down_sqrt2_muon;
    vector<TH1F*> all_fsr_down_sqrt2_elec          = h_fsr_down_sqrt2_elec;
    vector<TH1F*> all_fsr_down_sqrt2_combine       = h_fsr_down_sqrt2_combine;

    vector<TH1F*> all_fsr_down_sqrt2_muon_rebin    = RebinVector(all_fsr_down_sqrt2_muon, bin_width);
    vector<TH1F*> all_fsr_down_sqrt2_elec_rebin    = RebinVector(all_fsr_down_sqrt2_elec, bin_width);
    vector<TH1F*> all_fsr_down_sqrt2_combine_rebin = RebinVector(all_fsr_down_sqrt2_combine, bin_width);

    vector<TH1F*> all_fsr_up_sqrt2_muon_rebin      = RebinVector(all_fsr_up_sqrt2_muon, bin_width);
    vector<TH1F*> all_fsr_up_sqrt2_elec_rebin      = RebinVector(all_fsr_up_sqrt2_elec, bin_width);
    vector<TH1F*> all_fsr_up_sqrt2_combine_rebin   = RebinVector(all_fsr_up_sqrt2_combine, bin_width);

    vector<TH1F*> all_fsr_down_sqrt2_muon_norm, all_fsr_down_sqrt2_elec_norm, all_fsr_down_sqrt2_combine_norm;
    for(auto hist: all_fsr_down_sqrt2_muon_rebin) all_fsr_down_sqrt2_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_sqrt2_elec_rebin) all_fsr_down_sqrt2_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_sqrt2_combine_rebin) all_fsr_down_sqrt2_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_fsr_up_sqrt2_muon_norm, all_fsr_up_sqrt2_elec_norm, all_fsr_up_sqrt2_combine_norm;
    for(auto hist: all_fsr_up_sqrt2_muon_rebin) all_fsr_up_sqrt2_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_sqrt2_elec_rebin) all_fsr_up_sqrt2_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_sqrt2_combine_rebin) all_fsr_up_sqrt2_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_fsr_up_2_muon            = h_fsr_up_2_muon;
    vector<TH1F*> all_fsr_up_2_elec            = h_fsr_up_2_elec;
    vector<TH1F*> all_fsr_up_2_combine         = h_fsr_up_2_combine;

    vector<TH1F*> all_fsr_down_2_muon          = h_fsr_down_2_muon;
    vector<TH1F*> all_fsr_down_2_elec          = h_fsr_down_2_elec;
    vector<TH1F*> all_fsr_down_2_combine       = h_fsr_down_2_combine;

    vector<TH1F*> all_fsr_down_2_muon_rebin    = RebinVector(all_fsr_down_2_muon, bin_width);
    vector<TH1F*> all_fsr_down_2_elec_rebin    = RebinVector(all_fsr_down_2_elec, bin_width);
    vector<TH1F*> all_fsr_down_2_combine_rebin = RebinVector(all_fsr_down_2_combine, bin_width);

    vector<TH1F*> all_fsr_up_2_muon_rebin      = RebinVector(all_fsr_up_2_muon, bin_width);
    vector<TH1F*> all_fsr_up_2_elec_rebin      = RebinVector(all_fsr_up_2_elec, bin_width);
    vector<TH1F*> all_fsr_up_2_combine_rebin   = RebinVector(all_fsr_up_2_combine, bin_width);

    vector<TH1F*> all_fsr_down_2_muon_norm, all_fsr_down_2_elec_norm, all_fsr_down_2_combine_norm;
    for(auto hist: all_fsr_down_2_muon_rebin) all_fsr_down_2_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_2_elec_rebin) all_fsr_down_2_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_2_combine_rebin) all_fsr_down_2_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_fsr_up_2_muon_norm, all_fsr_up_2_elec_norm, all_fsr_up_2_combine_norm;
    for(auto hist: all_fsr_up_2_muon_rebin) all_fsr_up_2_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_2_elec_rebin) all_fsr_up_2_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_2_combine_rebin) all_fsr_up_2_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_fsr_up_4_muon            = h_fsr_up_4_muon;
    vector<TH1F*> all_fsr_up_4_elec            = h_fsr_up_4_elec;
    vector<TH1F*> all_fsr_up_4_combine         = h_fsr_up_4_combine;

    vector<TH1F*> all_fsr_down_4_muon          = h_fsr_down_4_muon;
    vector<TH1F*> all_fsr_down_4_elec          = h_fsr_down_4_elec;
    vector<TH1F*> all_fsr_down_4_combine       = h_fsr_down_4_combine;

    vector<TH1F*> all_fsr_down_4_muon_rebin    = RebinVector(all_fsr_down_4_muon, bin_width);
    vector<TH1F*> all_fsr_down_4_elec_rebin    = RebinVector(all_fsr_down_4_elec, bin_width);
    vector<TH1F*> all_fsr_down_4_combine_rebin = RebinVector(all_fsr_down_4_combine, bin_width);

    vector<TH1F*> all_fsr_up_4_muon_rebin      = RebinVector(all_fsr_up_4_muon, bin_width);
    vector<TH1F*> all_fsr_up_4_elec_rebin      = RebinVector(all_fsr_up_4_elec, bin_width);
    vector<TH1F*> all_fsr_up_4_combine_rebin   = RebinVector(all_fsr_up_4_combine, bin_width);

    vector<TH1F*> all_fsr_down_4_muon_norm, all_fsr_down_4_elec_norm, all_fsr_down_4_combine_norm;
    for(auto hist: all_fsr_down_4_muon_rebin) all_fsr_down_4_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_4_elec_rebin) all_fsr_down_4_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_down_4_combine_rebin) all_fsr_down_4_combine_norm.push_back(Normalize(hist));

    vector<TH1F*> all_fsr_up_4_muon_norm, all_fsr_up_4_elec_norm, all_fsr_up_4_combine_norm;
    for(auto hist: all_fsr_up_4_muon_rebin) all_fsr_up_4_muon_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_4_elec_rebin) all_fsr_up_4_elec_norm.push_back(Normalize(hist));
    for(auto hist: all_fsr_up_4_combine_rebin) all_fsr_up_4_combine_norm.push_back(Normalize(hist));

    // Set Bin Error -------------------------------------------------------------
    TMatrixD m_cov_data    = GetCovMatrix(all_data_combine_norm[3], debug);
    TMatrixD m_cov_ttbar   = GetCovMatrix(all_ttbar_combine_norm[3], debug);
    TMatrixD m_cov_JECup   = GetCovMatrix(all_jec_up_combine_norm[3], debug);
    TMatrixD m_cov_JECdown = GetCovMatrix(all_jec_down_combine_norm[3], debug);
    TMatrixD m_cov_CORup   = GetCovMatrix(all_cor_up_combine_norm[3], debug);
    TMatrixD m_cov_CORdown = GetCovMatrix(all_cor_down_combine_norm[3], debug);
    TMatrixD m_cov_JERup   = GetCovMatrix(all_jer_up_combine_norm[3], debug);
    TMatrixD m_cov_JERdown = GetCovMatrix(all_jer_down_combine_norm[3], debug);

    // TMatrixD NormCovMatrix(TH1* hist_, TMatrixD old_cov, bool width_, bool onlyDiag)
    // For plot only diagonal important
    int mDim = all_data_combine_norm[3]->GetNbinsX();
    TMatrixD m_cov_data_norm    = NormCovMatrix(all_data_combine_norm[3],     m_cov_data,    false, true);
    TMatrixD m_cov_ttbar_norm   = NormCovMatrix(all_ttbar_combine_norm[3],    m_cov_ttbar,   false, true);
    TMatrixD m_cov_JECup_norm   = NormCovMatrix(all_jec_up_combine_norm[3],   m_cov_JECup,   false, true);
    TMatrixD m_cov_JECdown_norm = NormCovMatrix(all_jec_down_combine_norm[3], m_cov_JECdown, false, true);
    TMatrixD m_cov_CORup_norm   = NormCovMatrix(all_cor_up_combine_norm[3],   m_cov_CORup,   false, true);
    TMatrixD m_cov_CORdown_norm = NormCovMatrix(all_cor_down_combine_norm[3], m_cov_CORdown, false, true);
    TMatrixD m_cov_JERup_norm   = NormCovMatrix(all_jer_up_combine_norm[3],   m_cov_JERup,   false, true);
    TMatrixD m_cov_JERdown_norm = NormCovMatrix(all_jer_down_combine_norm[3], m_cov_JERdown, false, true);

    for(unsigned int i = 1; i<=mDim; i++){
      all_data_combine_norm[3]->SetBinError(i, sqrt(m_cov_data_norm[i-1][i-1]));
      all_ttbar_combine_norm[3]->SetBinError(i, sqrt(m_cov_ttbar_norm[i-1][i-1]));
      all_jec_up_combine_norm[3]->SetBinError(i, sqrt(m_cov_JECup_norm[i-1][i-1]));
      all_jec_down_combine_norm[3]->SetBinError(i, sqrt(m_cov_JECdown_norm[i-1][i-1]));
      all_cor_up_combine_norm[3]->SetBinError(i, sqrt(m_cov_CORup_norm[i-1][i-1]));
      all_cor_down_combine_norm[3]->SetBinError(i, sqrt(m_cov_CORdown_norm[i-1][i-1]));
      all_jer_up_combine_norm[3]->SetBinError(i, sqrt(m_cov_JERup_norm[i-1][i-1]));
      all_jer_down_combine_norm[3]->SetBinError(i, sqrt(m_cov_JERdown_norm[i-1][i-1]));
    }

    // -------------------------------------------------------------------------------------------------------
    // ------------------------------------------ Plots ------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------

    for(int year=0; year<4; year++){
      string sYear = (string) years[year];
      cout << "Start with year "+sYear+" ..." << endl;

      // Masspeack -------------------------------------------------------------
      cout << "\t ... Masspeak Bins" << endl;

      // Get bins with bin-content>Limit
      vector<double> PeakLimit;
      vector<int> peak_bins_muon, peak_bins_elec, peak_bins_combine;
      double limit = (year==3)?100:0;
      peak_bins_muon    = bins_upper_limit(all_data_muon_rebin[year], limit);
      peak_bins_elec    = bins_upper_limit(all_data_elec_rebin[year], limit);
      peak_bins_combine = bins_upper_limit(all_data_combine_rebin[year], limit);

      // Bin Center
      vector<double> bin_centers_muon    = bin_center_upper_limit(all_data_muon_rebin[year],    peak_bins_muon,    bin_width);
      vector<double> bin_centers_elec    = bin_center_upper_limit(all_data_elec_rebin[year],    peak_bins_elec,    bin_width);
      vector<double> bin_centers_combine = bin_center_upper_limit(all_data_combine_rebin[year], peak_bins_combine, bin_width);

      if(year==3){
        int nBins = (int) bin_centers_combine[1]-(int) bin_centers_combine[0];
        TString option = gSystem->AccessPathName("files/PaperPlots_Peak.root")?"recreate":"update";
        TFile *f_out = new TFile("files/PaperPlots_Peak.root", option);
        f_out->cd();

        // void WriteHistogramToRoot(TH1F* hist, TString name, int xmin, int xmax){
        WriteHistogramToRoot(all_data_combine_norm[year],     "wmass"+addition+"__DATA",         range[w_mass][0], range[w_mass][1]);
        WriteHistogramToRoot(all_ttbar_combine_norm[year],    "wmass"+addition+"__TTbar",        range[w_mass][0], range[w_mass][1]);
        WriteHistogramToRoot(all_jec_up_combine_norm[year],   "wmass"+addition+"__TTbarJECup",   range[w_mass][0], range[w_mass][1]);
        WriteHistogramToRoot(all_jec_down_combine_norm[year], "wmass"+addition+"__TTbarJECdown", range[w_mass][0], range[w_mass][1]);
        WriteHistogramToRoot(all_cor_up_combine_norm[year],   "wmass"+addition+"__TTbarCORup",   range[w_mass][0], range[w_mass][1]);
        WriteHistogramToRoot(all_cor_down_combine_norm[year], "wmass"+addition+"__TTbarCORdown", range[w_mass][0], range[w_mass][1]);
        f_out->Close();
      }

      TString option = gSystem->AccessPathName("files/WMassPlots.root")?"recreate":"update";
      TFile *f_out = new TFile("files/WMassPlots.root", option);
      f_out->cd();

      // void WriteHistogramToRoot(TH1F* hist, TString name, int xmin, int xmax){
      TString name = sYear+addition;
      map<TString, double> means;

      means[name+"_comb_data"] = trunc_mean(all_data_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_ttbar"] = trunc_mean(all_ttbar_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_jec_up"] = trunc_mean(all_jec_up_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_jec_down"] = trunc_mean(all_jec_down_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_cor_up"] = trunc_mean(all_cor_up_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_cor_down"] = trunc_mean(all_cor_down_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_jer_up"] = trunc_mean(all_jer_up_combine_norm[year], range[w_mass][0], range[w_mass][1]);
      means[name+"_comb_jer_down"] = trunc_mean(all_jer_down_combine_norm[year], range[w_mass][0], range[w_mass][1]);

      means_nominal.push_back(means[name+"_comb_ttbar"]);
      means_jec_up.push_back(means[name+"_comb_jec_up"]);
      means_jec_down.push_back(means[name+"_comb_jec_down"]);
      means_cor_up.push_back(means[name+"_comb_cor_up"]);
      means_cor_down.push_back(means[name+"_comb_cor_down"]);

      //       bin             JECup            COR up           Nominal          COR down          JEC down alt
      // ----------------------------------------------------------------------------------------------
      //  hh           85.9054           85.6499           85.5152             85.37           85.0989
      // ----------------------------------------------------------------------------------------------
      //  hl           86.1001           85.8945           85.7095           85.5485           85.3465
      // ----------------------------------------------------------------------------------------------
      //  lh           86.1654           85.8911           85.6818           85.4935            85.172
      // ----------------------------------------------------------------------------------------------
      //  ll             86.08           85.8316           85.6485           85.4794           85.2079
      //
      //  bin             JECup            COR up           Nominal          COR down          JEC down neu
      // ----------------------------------------------------------------------------------------------
      //  hh            87.246           86.9342           86.7664           86.6002           86.2651
      // ----------------------------------------------------------------------------------------------
      //  hl            87.393           87.1675           86.9716            86.772           86.5533
      // ----------------------------------------------------------------------------------------------
      //  lh           87.5371           87.1833           86.9299           86.7161            86.365
      // ----------------------------------------------------------------------------------------------
      //  ll             87.42           87.1159           86.9017           86.7041           86.3986
      //
      //
      // $\ptW<300\GeV$, $\ptratio<0.7$  &  ${+}0.65$    &   ${+}0.35$   &  ${+}0.14$  &  ${-}0.06$    &  ${-}0.37$   \\ neu berechnet
      // $\ptW<300\GeV$, $\ptratio>0.7$  &  ${+}0.77$    &   ${+}0.42$   &  ${+}0.16$  &  ${-}0.05$    &  ${-}0.40$   \\
      // $\ptW>300\GeV$, $\ptratio<0.7$  &  ${+}0.63$    &   ${+}0.40$   &  ${+}0.21$  &  ${+}0.01$    &  ${-}0.21$   \\
      // $\ptW>300\GeV$, $\ptratio>0.7$  &  ${+}0.48$    &   ${+}0.17$   &   0         &  ${-}0.17$    &  ${-}0.50$

      //   $\ptW<300\GeV$, $\ptratio<0.7$  &  ${+}0.67$    &   ${+}0.35$   &  ${+}0.13$  &  ${-}0.07$    &  ${-}0.37$   \\ alt mit gerundet gerechnet
      //   $\ptW<300\GeV$, $\ptratio>0.7$  &  ${+}0.77$    &   ${+}0.41$   &  ${+}0.16$  &  ${-}0.05$    &  ${-}0.40$   \\
      //   $\ptW>300\GeV$, $\ptratio<0.7$  &  ${+}0.62$    &   ${+}0.40$   &  ${+}0.20$  &  ${+}0.00$    &  ${-}0.22$   \\
      //   $\ptW>300\GeV$, $\ptratio>0.7$  &  ${+}0.48$    &   ${+}0.16$   &   0         &  ${-}0.17$    &  ${-}0.50$

      //      bin             JECup            COR up           Nominal          COR down          JEC down ausgeben
      // ll              0.65              0.35              0.14             -0.06             -0.37
      // lh              0.77              0.42              0.16             -0.05             -0.40
      // hl              0.63              0.40              0.21              0.01             -0.21
      // hh              0.48              0.17              0.00             -0.17             -0.50
      //
      //   $\ptW<300\GeV$, $\ptratio<0.7$  &  ${+}0.56$    &   ${+}0.32$   &  ${+}0.13$  &  ${-}0.04$    &  ${-}0.31$   \\
      //   $\ptW<300\GeV$, $\ptratio>0.7$  &  ${+}0.65$    &   ${+}0.38$   &  ${+}0.17$  &  ${-}0.02$    &  ${-}0.34$   \\
      //   $\ptW>300\GeV$, $\ptratio<0.7$  &  ${+}0.58$    &   ${+}0.38$   &  ${+}0.19$  &  ${+}0.03$    &  ${-}0.17$   \\
      //   $\ptW>300\GeV$, $\ptratio>0.7$  &  ${+}0.39$    &   ${+}0.13$   &   0         &  ${-}0.15$    &  ${-}0.42$

      if(sYear=="combine"){
        for(auto m: means){
          cout << setw(30) << m.first << setw(15) << m.second << endl;
        }

        cout << endl;
        cout << " --- Absolut difference ---" << endl;
        cout << endl;
        cout << setw(18) << addition << setw(18) << "Data" << setw(18) << "Nominal" << setw(18) << "JEC up" << setw(18) << "JEC down" << setw(18) << "COR up" << setw(18) << "COR down" << endl;
        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "Data";
        cout << RED << setw(18) << means[name+"_comb_data"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "Nominal";
        cout << setw(18) << means[name+"_comb_ttbar"]-means[name+"_comb_data"];
        cout << RED << setw(18) << means[name+"_comb_ttbar"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "JEC up";
        cout << setw(18) << means[name+"_comb_jec_up"]-means[name+"_comb_data"];
        cout << setw(18) << means[name+"_comb_jec_up"]-means[name+"_comb_ttbar"];
        cout << RED << setw(18) << means[name+"_comb_jec_up"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "JEC down";
        cout << setw(18) << means[name+"_comb_jec_down"]-means[name+"_comb_data"];
        cout << setw(18) << means[name+"_comb_jec_down"]-means[name+"_comb_ttbar"];
        cout << setw(18) << means[name+"_comb_jec_down"]-means[name+"_comb_jec_up"];
        cout << RED << setw(18) << means[name+"_comb_jec_down"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "COR up";
        cout << setw(18) << means[name+"_comb_cor_up"]-means[name+"_comb_data"];
        cout << setw(18) << means[name+"_comb_cor_up"]-means[name+"_comb_ttbar"];
        cout << setw(18) << means[name+"_comb_cor_up"]-means[name+"_comb_jec_up"];
        cout << setw(18) << means[name+"_comb_cor_up"]-means[name+"_comb_jec_down"];
        cout << RED << setw(18) << means[name+"_comb_cor_up"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "COR down";
        cout << setw(18) << means[name+"_comb_cor_down"]-means[name+"_comb_data"];
        cout << setw(18) << means[name+"_comb_cor_down"]-means[name+"_comb_ttbar"];
        cout << setw(18) << means[name+"_comb_cor_down"]-means[name+"_comb_jec_up"];
        cout << setw(18) << means[name+"_comb_cor_down"]-means[name+"_comb_jec_down"];
        cout << setw(18) << means[name+"_comb_cor_down"]-means[name+"_comb_cor_up"];
        cout << RED << setw(18) << means[name+"_comb_cor_down"] << RESET;
        cout << endl;

        cout << endl;
        cout << " --- Relative difference ---" << endl;
        cout << endl;
        cout << setw(18) << addition << setw(18) << "Data" << setw(18) << "Nominal" << setw(18) << "JEC up" << setw(18) << "JEC down" << setw(18) << "COR up" << setw(18) << "COR down" << endl;
        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "Data";
        cout << RED << setw(18) << means[name+"_comb_data"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "Nominal";
        cout << setw(18) << (fabs(means[name+"_comb_ttbar"]-means[name+"_comb_data"]))/means[name+"_comb_ttbar"]*100;
        cout << RED << setw(18) << means[name+"_comb_ttbar"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "JEC up";
        cout << setw(18) << (fabs(means[name+"_comb_jec_up"]-means[name+"_comb_data"]))/means[name+"_comb_jec_up"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_jec_up"]-means[name+"_comb_ttbar"]))/means[name+"_comb_jec_up"]*100;
        cout << RED << setw(18) << means[name+"_comb_jec_up"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "JEC down";
        cout << setw(18) << (fabs(means[name+"_comb_jec_down"]-means[name+"_comb_data"]))/means[name+"_comb_jec_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_jec_down"]-means[name+"_comb_ttbar"]))/means[name+"_comb_jec_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_jec_down"]-means[name+"_comb_jec_up"]))/means[name+"_comb_jec_down"]*100;
        cout << RED << setw(18) << means[name+"_comb_jec_down"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "COR up";
        cout << setw(18) << (fabs(means[name+"_comb_cor_up"]-means[name+"_comb_data"]))/means[name+"_comb_cor_up"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_up"]-means[name+"_comb_ttbar"]))/means[name+"_comb_cor_up"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_up"]-means[name+"_comb_jec_up"]))/means[name+"_comb_cor_up"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_up"]-means[name+"_comb_jec_down"]))/means[name+"_comb_cor_up"]*100;
        cout << RED << setw(18) << means[name+"_comb_cor_up"] << RESET;
        cout << endl;

        cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
        cout << setw(18) << "COR down";
        cout << setw(18) << (fabs(means[name+"_comb_cor_down"]-means[name+"_comb_data"]))/means[name+"_comb_cor_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_down"]-means[name+"_comb_ttbar"]))/means[name+"_comb_cor_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_down"]-means[name+"_comb_jec_up"]))/means[name+"_comb_cor_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_down"]-means[name+"_comb_jec_down"]))/means[name+"_comb_cor_down"]*100;
        cout << setw(18) << (fabs(means[name+"_comb_cor_down"]-means[name+"_comb_cor_up"]))/means[name+"_comb_cor_down"]*100;
        cout << RED << setw(18) << means[name+"_comb_cor_down"] << RESET;
        cout << endl;
      }

      WriteHistogramToRoot(all_data_combine_norm[year],     "wmass"+addition+"__DATA__"+sYear,         (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_ttbar_combine_norm[year],    "wmass"+addition+"__TTbar__"+sYear,        (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_jec_up_combine_norm[year],   "wmass"+addition+"__TTbarJECup__"+sYear,   (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_jec_down_combine_norm[year], "wmass"+addition+"__TTbarJECdown__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_cor_up_combine_norm[year],   "wmass"+addition+"__TTbarCORup__"+sYear,   (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_cor_down_combine_norm[year], "wmass"+addition+"__TTbarCORdown__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_jer_up_combine_norm[year],   "wmass"+addition+"__TTbarJERup__"+sYear,   (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_jer_down_combine_norm[year], "wmass"+addition+"__TTbarJERdown__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_hdamp_up_combine_norm[year],   "wmass"+addition+"__TTbarhdampup__"+sYear,   (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_hdamp_down_combine_norm[year], "wmass"+addition+"__TTbarhdampdown__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_tune_up_combine_norm[year],    "wmass"+addition+"__TTbartuneup__"+sYear,    (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_tune_down_combine_norm[year],  "wmass"+addition+"__TTbartunedown__"+sYear,  (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_cr_gm_combine_norm[year],      "wmass"+addition+"__TTbargluonmove__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_cr_qb_combine_norm[year],      "wmass"+addition+"__TTbarQCDbased__"+sYear,  (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_up_sqrt2_combine_norm[year],   "wmass"+addition+"__TTbarFSRup_sqrt2__"+sYear,   (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_down_sqrt2_combine_norm[year], "wmass"+addition+"__TTbarFSRdown_sqrt2__"+sYear, (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_up_2_combine_norm[year],       "wmass"+addition+"__TTbarFSRup_2__"+sYear,       (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_down_2_combine_norm[year],     "wmass"+addition+"__TTbarFSRdown_2__"+sYear,     (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_up_4_combine_norm[year],       "wmass"+addition+"__TTbarFSRup_4__"+sYear,       (int) bin_centers_combine[0], (int) bin_centers_combine[1]);
      WriteHistogramToRoot(all_fsr_down_4_combine_norm[year],     "wmass"+addition+"__TTbarFSRdown_4__"+sYear,     (int) bin_centers_combine[0], (int) bin_centers_combine[1]);

      f_out->Close();

    } // year
  } // pt bin

  vector<TString> bins = {"hh", "hl", "lh", "ll"};
  cout << endl;
  cout << setw(4) << "bin" << setw(18) << "JECup" << setw(18) << "COR up" << setw(18) << "Nominal" << setw(18) << "COR down" << setw(18) << "JEC down" << endl;
  for(int b=0; b<bins.size(); b++){
    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << setw(4) << bins[b] << setw(18) << means_jec_up[b] << setw(18) << means_cor_up[b] << setw(18) << means_nominal[b] << setw(18) << means_cor_down[b] << setw(18) << means_jec_down[b] << endl;
  }

  cout << endl;
  cout << "absolut" << endl;
  cout << setw(4) << "bin" << setw(18) << "JECup" << setw(18) << "COR up" << setw(18) << "Nominal" << setw(18) << "COR down" << setw(18) << "JEC down" << endl;
  for(int b=0; b<bins.size(); b++){
    cout << setw(4) << bins[b] << setw(18) << means_jec_up[b]-means_nominal[0] << setw(18) << means_cor_up[b]-means_nominal[0] << setw(18) << means_nominal[b]-means_nominal[0] << setw(18) << means_cor_down[b]-means_nominal[0] << setw(18) << means_jec_down[b]-means_nominal[0] << endl;
  }

  cout << endl;
  cout << "absolut rounded" << endl;
  cout << setw(4) << "bin" << setw(18) << "JECup" << setw(18) << "COR up" << setw(18) << "Nominal" << setw(18) << "COR down" << setw(18) << "JEC down" << endl;
  for(int b=0; b<bins.size(); b++){
    cout << setw(4) << bins[b] << setw(18) << dtos(means_jec_up[b]-means_nominal[0], 2) << setw(18) << dtos(means_cor_up[b]-means_nominal[0], 2) << setw(18) << dtos(means_nominal[b]-means_nominal[0], 2) << setw(18) << dtos(means_cor_down[b]-means_nominal[0], 2) << setw(18) << dtos(means_jec_down[b]-means_nominal[0], 2) << endl;
  }
} // main
