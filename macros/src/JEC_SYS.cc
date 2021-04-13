#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>

typedef std::vector<double> VecD;
typedef std::vector<vector<double>> VecDD;

typedef std::map<TString, double> MapD;
typedef std::map<TString, VecD> MapVD;

Double_t plane(Double_t *x, Double_t *par){
  return 0.;
}

VecD ExtractJMSValues(TF1* chi2function);
void SetFunctionParameters(TF1* chi2function, VecDD paramters);
void drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange, TString save_path);
void JMSinTXT(double nom, double up, double down, TString variation, TString save);
void printVector(VecD vector);


// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  bool debug = false;
  TString save_path = get_save_path();

  // Modification ####################################################################################
  // #################################################################################################
  // Changes of how many pt bins to use // CHANGE_PT
  // #################################################################################################
  // Modification ####################################################################################

  // #################################################################################################
  // Declare different variables used in the code ####################################################
  print_seperater();
  cout.precision(6);

  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_SYS <year>\n";
    return 0;
  }

  // #################################################################################################
  // Only one fit for all bins #######################################################################
  if(debug) cout << "String into bool" << endl;

  // Default -------------------------------------------------------------------
  Int_t   oldLevel     = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ... - functions_explain
  gErrorIgnoreLevel    = kWarning;          // suppress TCanvas output
  TString reconst      = "btag";            // match btag_cut btag_sel compare min_mass
  int     bin_width    = 1;
  bool    usePeak_in   = true;
  bool    into_latex   = false;
  bool    print_table  = false;
  bool    onlyData     = false;
  bool    useOnly_lin  = false;
  bool    pt_bins      = true;

  // Input ---------------------------------------------------------------------
  TString year         = argv[1];

  // #################################################################################################
  // cout settings ###################################################################################
  cout << "============================= General Settings\n";
  cout << "debug:     " << debug       << endl;
  cout << "Latex:     " << into_latex  << endl;
  cout << "Table:     " << print_table << endl;
  cout << "\n============================= Calculation Settings\n";
  cout << "Peak:      " << usePeak_in  << endl;
  cout << "Pt bins:   " << pt_bins     << endl;
  cout << "Only lin:  " << useOnly_lin << endl;
  cout << "Bin Width: " << bin_width   << endl;



  // #################################################################################################
  // Define variables from input #####################################################################

  // Rebin -------------------------------------------------------------------------------------------
  if(debug) cout << "Set Number Bins" << endl;
  int     number_bins     = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Set Year" << endl;

  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16  = true;
  else if(strcmp(year, "2017")==0)     is17  = true;
  else if(strcmp(year, "2018")==0)     is18  = true;
  else if(strcmp(year, "combined")==0) isAll = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or combined)");

  // same fit ----------------------------------------------------------------------------------------
  // fits are ordered in by the number of parameters. Lin fit is the first option -> index_lin =0
  // Used later in Code
  int same_fit;
  if(useOnly_lin) same_fit=0;

  /*
  ██████  ██ ██████  ███████  ██████ ████████  ██████  ██████  ██ ███████ ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██      ██
  ██   ██ ██ ██████  █████   ██         ██    ██    ██ ██████  ██ █████   ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██           ██
  ██████  ██ ██   ██ ███████  ██████    ██     ██████  ██   ██ ██ ███████ ███████
  */
  if(debug) cout << "Directiories" << endl;

  TString save_path_general = save_path+"/Plots/JEC_SYS";
  save_path_general = creat_folder_and_path(save_path_general, "chi2");
  if(pt_bins) save_path_general = creat_folder_and_path(save_path_general, "pt_bins"); // CHANGE_PT
  else        save_path_general = creat_folder_and_path(save_path_general, "no_bins");

  // #################################################################################################
  // creat subdirectories ############################################################################
  if(debug) cout << "Sub-Directiories" << endl;
  save_path_general = creat_folder_and_path(save_path_general, year);
  save_path_general = creat_folder_and_path(save_path_general, reconst);
  save_path_general = creat_folder_and_path(save_path_general, "rebin"+str_number_bins);
  if(useOnly_lin)       save_path_general = creat_folder_and_path(save_path_general, "linear");
  if(usePeak_in&&isAll) save_path_general = creat_folder_and_path(save_path_general, "masspeak");
  creat_folder(save_path_general, "single_bins");
  creat_folder(save_path_general, "projection");

  /*
  ██████  ████████     ██       ██████   ██████  ██████
  ██   ██    ██        ██      ██    ██ ██    ██ ██   ██
  ██████     ██        ██      ██    ██ ██    ██ ██████
  ██         ██        ██      ██    ██ ██    ██ ██
  ██         ██        ███████  ██████   ██████  ██
  */
  /* Goes down to LATEX */
  if(debug) cout << "Start: Pt-Loop" << endl;
  cout << "\nPath: "+save_path_general << endl;

  TString addition="";
  TString chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll;
  VecD chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll;

  double number_bins_total = 0;
  double number_not_lin    = 0;

  for(int ptbin=0; ptbin<5; ptbin++){
    /* Four ptbins: "W-pT" and ratio of "subjet_high-pT/W-pT" - each two bins || CHANGE_PT */

    if(pt_bins){
      if(ptbin==4)    continue;
    } else{
      if(!(ptbin==4)) continue;
    }
    TString ptbin_str = to_string(ptbin);
    if(ptbin==0) addition="_hh";
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="";
    creat_folder(save_path_general, "single_bins/"+addition);

    cout << "\n!----- "+addition+" -----!\n";

    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */
    cout << '\n';
    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
    /* Options for other pT bin in storage */
    /* CHANGE_PT */
    if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(reconst=="btag"&&ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(reconst=="btag"&&ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";
    if(reconst=="btag"&&ptbin==4) w_mass = hist_class+"wmass_match";

    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    vector<TFile*>  file_bkg_v;
    vector<TH1F*>   hists_bkg_v;
    TString prefix = "uhh2.AnalysisModuleRunner.MC.";
    vector<TString> channels = {"elec", "muon"};
    vector<TString> path_bkg_v = {"other", "SingleTop", "WJets"};
    vector<TString> years;
    if(is16) years = {"2016v3"};
    else if(is17) years = {"2017v2"};
    else if(is18) years = {"2018"};
    else if(isAll) years = {"2016v3", "2017v2", "2018"};

    for(auto y: years){
      for(auto channel: channels){
        for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+"/"+channel+"/"+prefix+path_bkg_v[i]+"_"+y+".root"));
      }
    }
    for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
    TH1F *bkg = AddHists(hists_bkg_v, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    vector<TFile*> data_files;
    vector<TH1F*> data_hists;
    TString data_path        = "uhh2.AnalysisModuleRunner.DATA.DATA";
    for(auto y: years){
      for(auto channel: channels){
        data_files.push_back(new TFile(dir+"/"+channel+"/"+data_path+"_"+y+".root"));
      }
    }

    for(unsigned int i=0; i<data_files.size(); i++) data_hists.push_back((TH1F*)data_files[i]->Get(w_mass));
    TH1F   *data             = AddHists(data_hists, 1);

    TH1F   *data_rebin       = rebin(data, bin_width);
    TH1F   *data_norm        = normalize(data);
    TH1F   *data_rebin_norm  = normalize(data_rebin);

    VecD data_rebin_norm_err = normalize_error(data_rebin);
    for(int ipar=0; ipar<data_rebin_norm->GetNbinsX(); ipar++)
    {
      data_rebin_norm->SetBinError(ipar+1, data_rebin_norm_err[ipar]);
    }

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;
    vector<TFile*> ttbar_files, ttbar_files_jecup, ttbar_files_jecdown, ttbar_files_corup, ttbar_files_cordown;
    vector<TH1F*> ttbar_hists, ttbar_hists_jecup, ttbar_hists_jecdown, ttbar_hists_corup, ttbar_hists_cordown;
    TString ttbar_path        = "uhh2.AnalysisModuleRunner.MC.TTbar";
    for(auto y: years){
      for(auto channel: channels){
        ttbar_files.push_back(new TFile(dir+"/"+channel+"/"+ttbar_path+"_"+y+".root"));
        ttbar_files_jecup.push_back(new TFile(dir+"/"+channel+"/JEC_up/"+ttbar_path+"_"+y+".root"));
        ttbar_files_jecdown.push_back(new TFile(dir+"/"+channel+"/JEC_down/"+ttbar_path+"_"+y+".root"));
        ttbar_files_corup.push_back(new TFile(dir+"/"+channel+"/COR_up/"+ttbar_path+"_"+y+".root"));
        ttbar_files_cordown.push_back(new TFile(dir+"/"+channel+"/COR_down/"+ttbar_path+"_"+y+".root"));
      }
    }

    for(unsigned int i=0; i<ttbar_files.size(); i++) ttbar_hists.push_back((TH1F*)ttbar_files[i]->Get(w_mass));
    for(unsigned int i=0; i<ttbar_files_jecup.size(); i++) ttbar_hists_jecup.push_back((TH1F*)ttbar_files_jecup[i]->Get(w_mass));
    for(unsigned int i=0; i<ttbar_files_jecdown.size(); i++) ttbar_hists_jecdown.push_back((TH1F*)ttbar_files_jecdown[i]->Get(w_mass));
    for(unsigned int i=0; i<ttbar_files_corup.size(); i++) ttbar_hists_corup.push_back((TH1F*)ttbar_files_corup[i]->Get(w_mass));
    for(unsigned int i=0; i<ttbar_files_cordown.size(); i++) ttbar_hists_cordown.push_back((TH1F*)ttbar_files_cordown[i]->Get(w_mass));

    TH1F   *ttbar             = AddHists(ttbar_hists, 1);
    TH1F   *JECup             = AddHists(ttbar_hists_jecup, 1);
    TH1F   *JECdown           = AddHists(ttbar_hists_jecdown, 1);
    TH1F   *XConeup             = AddHists(ttbar_hists_corup, 1);
    TH1F   *XConedown           = AddHists(ttbar_hists_cordown, 1);

    ttbar->Add(bkg, 1);
    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);
    XConeup->Add(bkg, 1);
    XConedown->Add(bkg, 1);

    TH1F* ttbar_rebin          = rebin(ttbar, bin_width);
    TH1F* ttbar_norm           = normalize(ttbar);
    TH1F* ttbar_rebin_norm     = normalize(ttbar_rebin);

    TH1F* JECup_norm           = normalize(JECup);
    TH1F* JECdown_norm         = normalize(JECdown);
    TH1F* JECup_rebin          = rebin(JECup, bin_width);
    TH1F* JECdown_rebin        = rebin(JECdown, bin_width);
    TH1F* JECup_rebin_norm     = rebin(JECup_norm, bin_width);
    TH1F* JECdown_rebin_norm   = rebin(JECdown_norm, bin_width);

    TH1F* XConeup_norm         = normalize(XConeup);
    TH1F* XConedown_norm       = normalize(XConedown);
    TH1F* XConeup_rebin        = rebin(XConeup, bin_width);
    TH1F* XConedown_rebin      = rebin(XConedown, bin_width);
    TH1F* XConeup_rebin_norm   = rebin(XConeup_norm, bin_width);
    TH1F* XConedown_rebin_norm = rebin(XConedown_norm, bin_width);

    VecD ttbar_rebin_norm_err     = normalize_error(ttbar_rebin);
    VecD JECup_rebin_norm_err     = normalize_error(JECup_rebin);
    VecD JECdown_rebin_norm_err   = normalize_error(JECdown_rebin);
    VecD XConeup_rebin_norm_err   = normalize_error(XConeup_rebin);
    VecD XConedown_rebin_norm_err = normalize_error(XConedown_rebin);

    /*
    .██    ██ ███████ ███████ ██████      ██████  ██ ███    ██ ███████
    .██    ██ ██      ██      ██   ██     ██   ██ ██ ████   ██ ██
    .██    ██ ███████ █████   ██   ██     ██████  ██ ██ ██  ██ ███████
    .██    ██      ██ ██      ██   ██     ██   ██ ██ ██  ██ ██      ██
    . ██████  ███████ ███████ ██████      ██████  ██ ██   ████ ███████
    */

    // #################################################################################################
    // Empty Data Bins #################################################################################
    if(debug) cout << "Empty Bins Data" << endl;

    vector<int> empty_bins_v = bins_empty(data_rebin_norm); // Get empty bins
    int number_empty_bins    = empty_bins_v.size();

    // #################################################################################################
    // Masspeack #######################################################################################
    if(debug) cout << "Masspeak Bins" << endl;

    bool usePeak =false;
    VecD PeakLimit;
    if(isAll&&usePeak_in) usePeak = true;

    double Limit;
    if(!pt_bins) Limit = 190;
    else         Limit = 100;
    vector<int> peak_bins_v;
    if(usePeak) peak_bins_v = bins_upper_limit(data_rebin, Limit); // Get bins withc bin-content>Limit

    // Debug Peak & Empty ######################################################
    if(debug){
      cout << "number peak bins:  " << peak_bins_v.size();
      cout << "number empty bins: " << number_empty_bins << endl;
      for(int bin=0; bin < peak_bins_v.size(); bin++) cout << "Peak bins: " << peak_bins_v[bin] << endl;
      cout << endl;
      for(int bin=0; bin < number_empty_bins; bin++)  cout << "Empty Bins:  " << empty_bins_v[bin] << "\n";
    }

    // #################################################################################################
    // Settings Plots ##################################################################################
    /*
    The Mass Plots are drawn in a seperated file - JEC_SYS_MassPlots.cc
    */
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetLegendBorderSize(0);

    // Legend ------------------------------------------------------------------------------------------
    TLegend *leg;


    /*
    ██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██ ███████
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██ ██
    ██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██ ███████
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██      ██
    ██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████ ███████
    */
    /*
    Looking at each bin seperatly. Getting the bin content of the nominal file and
    the variations and making a fit through the values. The variation UP and DOWN
    are assigned to the values +1 and -1.
    Afterward the fit parameters will be stored in vectors for later purpose.
    */

    if(debug) print_seperater();
    if(debug) cout << "Start: Dependency in WJet bins" << endl;

    // #################################################################################################
    // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
    if(debug) cout << "Title for bins" << endl;

    vector<TString> mass_bin_title;
    for(int bin=0; bin<number_bins; bin++) mass_bin_title.push_back(to_string(bin*bin_width)+" < m_{Wjet} < "+to_string((bin+1)*bin_width));

    vector<TString> number_bin;
    for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

    // Debug Mass Titles #######################################################
    if(debug){
      cout << "number_bins " << number_bins << " | number titles " << mass_bin_title.size() << " | bin width " << bin_width << "\n\n";
      for(unsigned int bin=0; bin<mass_bin_title.size(); bin++){
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())) continue;
        if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue; // cross-check for empty bins
        cout << mass_bin_title[bin] << '\n';
      }
    }

    // #################################################################################################
    // Cout Bin Error ##################################################################################
    if(debug||print_table){
      cout << "" << '\n' << fixed;
      cout << "|================================================================================================|" << endl;
      cout << "| ----------------------------------------- Bin Content -----------------------------------------|" << endl;
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << "|  Bin | Data Content |   JEC down   |  XCone down  |   nominal    |    JEC up    |   XCone up   |" << endl;
      cout << "| -----|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
      // In this loop no additional treatment is  necessary.
      for(unsigned int i=0; i<number_bins; i++){
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), i+1) != peak_bins_v.end())) continue;
        int bin;
        double data, Jd, tt, Ju, Xd, Xu;
        bin  = i+1;
        data = data_rebin_norm->GetBinContent(i+1);
        Jd   = JECdown_rebin_norm->GetBinContent(i+1);
        tt   = ttbar_rebin_norm->GetBinContent(i+1);
        Ju   = JECup_rebin_norm->GetBinContent(i+1);
        Xd   = XConedown_rebin_norm->GetBinContent(i+1);
        Xu   = XConeup_rebin_norm->GetBinContent(i+1);

        if( i % 2 == 0) cout << DGRAY;
        else            cout << RESET;
        cout << "|" << setw(6) << centered(to_string(bin)) << "|" << setw(14) << centered(to_string(data)) << "|";
        cout << setw(14) << centered(to_string(Jd)) << "|";
        cout << setw(14) << centered(to_string(Xd)) << "|" << setw(14) << centered(to_string(tt)) << "|";
        cout << setw(14) << centered(to_string(Ju)) << "|" << setw(14) << centered(to_string(Xu)) << "|\n";
      }
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << RESET;

      // #################################################################################################
      // Cout Bin Error ##################################################################################

      cout << "" << '\n' << fixed;
      cout << "|================================================================================================|" << endl;
      cout << "| ------------------------------------------ Bin Errors -----------------------------------------|" << endl;
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << "|  Bin | Data Content | JEC down Err | nominal(cen) |  JEC up Err  |   Average    |   cen/avg    |" << endl;
      cout << "| -----|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
      // In this loop no additional treatment is  necessary.
      for(unsigned int i=0; i<number_bins; i++){
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), i+1) != peak_bins_v.end())) continue;
        string bin, data_str, Jd_str, tt_str, Ju_str, avg_str, cenavg_str;
        double avg = (JECdown_rebin_norm_err[i]+ttbar_rebin_norm_err[i]+JECup_rebin_norm_err[i])/3;
        double cenavg;
        if(abs(avg) > 0) cenavg = ttbar_rebin_norm_err[i]/avg;                 // make sure bin error is not 0
        else             cenavg = avg;
        if(i+1<10)     bin = to_string(i+1);                       // i+1 to get bin number
        if(i+1>=10)    bin = to_string(i+1);                        // +start_bin to skip first bins
        data_str           = to_string(data_rebin_norm_err[i]);
        Jd_str             = to_string(JECdown_rebin_norm_err[i]);
        tt_str             = to_string(ttbar_rebin_norm_err[i]);
        Ju_str             = to_string(JECup_rebin_norm_err[i]);
        avg_str            = to_string(avg);
        cenavg_str         = to_string(cenavg);

        if( i % 2 == 0) cout << DGRAY;
        else            cout << RESET;
        cout << "|" << setw(6) << centered(bin) << "|" << setw(14) << centered(data_str) << "|";
        cout << setw(14) << centered(Jd_str) << "|";
        cout << setw(14) << centered(tt_str) << "|" << setw(14) << centered(Ju_str) << "|";
        cout << setw(14) << centered(avg_str) << "|" << setw(14) << centered(cenavg_str) << "|\n";
      }
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << RESET;

      // #################################################################################################
      // XCone Bin Content and error #####################################################################
      cout << "" << '\n' << fixed;
      cout << "|================================================================================================|" << endl;
      cout << "| ------------------------------------------ Bin Errors -----------------------------------------|" << endl;
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << "|  Bin | Data Content | XC down Err  | nominal(cen) |  XC up Err   |   Average    |   cen/avg    |" << endl;
      cout << "| -----|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
      // In this loop no additional treatment is  necessary.
      for(unsigned int i=0; i<number_bins; i++){
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), i+1) != peak_bins_v.end())) continue;
        string bin, data_str, Xd_str, tt_str, Xu_str, avg_str, cenavg_str;
        double avg = (XConedown_rebin_norm_err[i]+ttbar_rebin_norm_err[i]+XConeup_rebin_norm_err[i])/3;
        double cenavg;
        if(abs(avg) > 0) cenavg = ttbar_rebin_norm_err[i]/avg;      // make sure bin error is not 0
        else             cenavg = avg;
        if(i+1<10)     bin = to_string(i+1);                        // i+1 to get bin number
        if(i+1>=10)    bin = to_string(i+1);                        // +start_bin to skip first bins
        data_str   = to_string(data_rebin_norm_err[i]);
        Xd_str     = to_string(XConedown_rebin_norm_err[i]);
        tt_str     = to_string(ttbar_rebin_norm_err[i]);
        Xu_str     = to_string(XConeup_rebin_norm_err[i]);
        avg_str    = to_string(avg);
        cenavg_str = to_string(cenavg);

        if( i % 2 == 0) cout << DGRAY;
        else            cout << RESET;
        cout << "|" << setw(6) << centered(bin) << "|" << setw(14) << centered(data_str) << "|";
        cout << setw(14) << centered(Xd_str) << "|";
        cout << setw(14) << centered(tt_str) << "|" << setw(14) << centered(Xu_str) << "|";
        cout << setw(14) << centered(avg_str) << "|" << setw(14) << centered(cenavg_str) << "|\n";
      }
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << RESET << "\n\n" << endl;
    }

    // #################################################################################################
    // Bin Error one sigma estimation ##################################################################
    if(debug) cout << "Define BinError" << endl;

    VecD bin_error_one_sigma_avg, bin_error_one_sigma_central;
    for(int bin=0; bin<number_bins; bin++){
      bin_error_one_sigma_central.push_back(ttbar_rebin_norm_err[bin]);
      bin_error_one_sigma_avg.push_back((XConedown_rebin_norm_err[bin]+ttbar_rebin_norm_err[bin]+XConeup_rebin_norm_err[bin])/3);
    }

    /*
    ███████ ██ ████████     ██    ██  █████  ██      ██    ██ ███████ ███████
    ██      ██    ██        ██    ██ ██   ██ ██      ██    ██ ██      ██
    █████   ██    ██        ██    ██ ███████ ██      ██    ██ █████   ███████
    ██      ██    ██         ██  ██  ██   ██ ██      ██    ██ ██           ██
    ██      ██    ██          ████   ██   ██ ███████  ██████  ███████ ███████
    */

    // #################################################################################################
    // Some Declarations ###############################################################################
    if(debug) cout << "Define fit input" << endl;

    // Fit ---------------------------------------------------------------------------------------------
    /* For similiar functions (xy2 & x2y) or (xyy2 & xx2y) the probability for y2 functions is always higher.
    Therefore, functions with only x2 are excluded to make the code more clean.*/
    TGraph2DErrors *bin_fit;
    TString str_fit_lin  = "[0] + [1]*x + [2]*y";                     // 3
    TString str_fit_xy2  = "[0] + [1]*x + [2]*y*y";                   // 3
    TString str_fit_x2y  = "[0] + [1]*x*x + [2]*y";                   // 3
    TString str_fit_xyy2 = "[0] + [1]*x + [2]*y + [3]*y*y";           // 4
    TString str_fit_xx2y = "[0] + [1]*x + [2]*x*x + [3]*y";           // 4
    TString str_fit_quad = "[0] + [1]*x*x + [2]*y*y";                 // 3
    TString str_fit_poly = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y"; // 5
    vector<vector<TString>> str_fits  = {{str_fit_lin}, {str_fit_xy2, str_fit_x2y}, {str_fit_xyy2, str_fit_xx2y}, {str_fit_quad}, {str_fit_poly}};
    vector<vector<TString>> name_fits = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
    vector<vector<int>> ndfs  = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
    vector<vector<int>> n_fit_para  = {{3}, {3, 3}, {4, 4}, {3}, {5}};
    vector<int> all_n_para;                     // Fill later - Declare befor loop
    vector<TF2*> all_fits;                      // Fill later - Declare befor loop
    vector<TString> str_all_fits;               // Fill later - Declare befor loop
    VecD chi2_parameters;             // Fill later - Declare befor loop
    VecDD fits_parameters;     // Fill later - Declare befor loop - Debug purpose only

    double npar=0;
    int n_factors = 5; // JECdown, XConedown, ttbar, XConeup,JECup

    // Order important! - Used for TGraph2DErrors in Loop ----------------------
    vector<TH1F*> hists              = {ttbar_rebin_norm, JECdown_rebin_norm, JECup_rebin_norm, XConedown_rebin_norm, XConeup_rebin_norm};
    VecDD hists_err = {ttbar_rebin_norm_err, JECdown_rebin_norm_err, JECup_rebin_norm_err, XConedown_rebin_norm_err, XConeup_rebin_norm_err};
    VecDD hists_cont, hists_cont_err; // Used in Projection, Filled in Loop

    VecD factor_x          = {0.0, -1.0,  1.0,  0.0,  0.0};
    VecD factor_y          = {0.0,  0.0,  0.0, -1.0,  1.0};
    VecD error_x           = {0.0,  0.0,  0.0,  0.0,  0.0};
    VecD error_y           = {0.0,  0.0,  0.0,  0.0,  0.0};

    // #################################################################################################
    // Hists for Fits ##################################################################################
    TH1I used_fits   = TH1I("used_fits", "used_fits", number_bins, 1, number_bins+1);
    TH1I number_fits = TH1I("number_fits", "number_fits", 5, 0, 5);

    /*
    ███████ ██ ████████     ██       ██████   ██████  ██████      ███████ ████████  █████  ██████  ████████
    ██      ██    ██        ██      ██    ██ ██    ██ ██   ██     ██         ██    ██   ██ ██   ██    ██
    █████   ██    ██        ██      ██    ██ ██    ██ ██████      ███████    ██    ███████ ██████     ██
    ██      ██    ██        ██      ██    ██ ██    ██ ██               ██    ██    ██   ██ ██   ██    ██
    ██      ██    ██        ███████  ██████   ██████  ██          ███████    ██    ██   ██ ██   ██    ██
    */

    if(debug) cout << "Start: Bin Loop" << endl;

    for(int bin=0; bin < number_bins; bin++){

      // #################################################################################################
      // Skip bin ########################################################################################
      // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue;
      if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())) continue;
      cout << "\n";
      if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      if(debug) cout << mass_bin_title[bin] << endl;
      number_bins_total++;

      // #################################################################################################
      // Fill Vectors ####################################################################################
      if(debug) cout << "Get Bin Content" << endl;

      VecD bin_content, bin_content_err;
      for(int factor=0; factor<n_factors; factor++){
        bin_content.push_back(hists[factor]->GetBinContent(bin+1));
        bin_content_err.push_back(hists_err[factor][bin]);
      }

      // Later used in Projection
      hists_cont.push_back(bin_content);
      hists_cont_err.push_back(bin_content_err);

      // #################################################################################################
      // Set TGraph ######################################################################################
      if(debug) cout << "TGraphError - Single bin " << bin+1 << endl;
      bin_fit = new TGraph2DErrors(5, &factor_x[0], &factor_y[0],&bin_content[0],&error_x[0],&error_y[0],&bin_content_err[0]);
      bin_fit->SetName(number_bin[bin]); // To avoid Warning: Replacing existing TGraph2D

      // #################################################################################################
      // set fit function ################################################################################
      if(debug) cout << "\nDefine fit function" << endl;
      // TF2 *fit_con   = new TF2("fit_con",  "[0]");
      // Declared here, because they are reseted each bin
      TF2 *fit_lin  = new TF2("fit_lin",  str_fits[0][0]);
      TF2 *fit_xy2  = new TF2("fit_xy2",  str_fits[1][0]);
      TF2 *fit_x2y  = new TF2("fit_x2y",  str_fits[1][1]);
      TF2 *fit_xyy2 = new TF2("fit_xyy2", str_fits[2][0]);
      TF2 *fit_xx2y = new TF2("fit_xx2y", str_fits[2][1]);
      TF2 *fit_quad = new TF2("fit_quad", str_fits[3][0]);
      TF2 *fit_poly = new TF2("fit_poly", str_fits[4][0]);
      vector<vector<TF2*>> fits = {{fit_lin}, {fit_xy2, fit_x2y}, {fit_xyy2, fit_xx2y}, {fit_quad}, {fit_poly}};

      // #################################################################################################
      // Fit Bin #########################################################################################
      if(debug) cout << "Start: Fitting" << endl;

      /*
      ███████ ██ ████████
      ██      ██    ██
      █████   ██    ██
      ██      ██    ██
      ██      ██    ██
      */

      double chi2_fit, chi2_fit1, chi2_fit2;
      double prob, prob1, prob2;
      TF2 *used_fit;
      int order_index = 0; int index_fit = 0;
      for(unsigned int order=0; order<fits.size(); order++){
        if(order!=same_fit && useOnly_lin) continue;
        if(debug) cout << "\nfit a "+name_fits[order][0]+" function ------------\n";
        bool twoFits = fits[order].size()==2;

        bin_fit->Fit(fits[order][0], "Q");
        if(twoFits) bin_fit->Fit(fits[order][1], "Q");

        // -----------------------------------------------------------------------
        chi2_fit1 = fits[order][0]->GetChisquare();
        if(twoFits) chi2_fit2 = fits[order][1]->GetChisquare();

        // -----------------------------------------------------------------------
        prob1 = TMath::Prob(chi2_fit1,  ndfs[order][0]);
        if(twoFits) prob2 = TMath::Prob(chi2_fit2, ndfs[order][1]);

        prob = prob1;
        chi2_fit = chi2_fit1;
        index_fit = 0;
        if(twoFits && (prob1 < prob2)){
          prob = prob2;
          chi2_fit = chi2_fit2;
          index_fit = 1;
        }
        if(debug) cout << "ndfs: " << ndfs[order][index_fit] << " | chi2: " << chi2_fit <<" | prob: " << prob << endl;

        // -----------------------------------------------------------------------
        if(prob>0.05 || ndfs[order][0]==0 || useOnly_lin){ // prob > 0.05 is common in data science
          cout << "Bin " << number_bin[bin] << ": A " << GREEN << name_fits[order][index_fit] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;

          used_fit = new TF2(*fits[order][index_fit]);        // used in Loop
          order_index=order;                                  // used in Loop - SetParameters for Chi2

          all_fits.push_back(fits[order][index_fit]);         // TF2*    - used in Chi2
          str_all_fits.push_back(str_fits[order][index_fit]); // TString - used in Chi2
          all_n_para.push_back(n_fit_para[order][index_fit]); // int     - used in Projection

          number_fits.AddBinContent(order+1);                // Dummy Hist - total number
          used_fits.AddBinContent(bin+1, order+1);           // Dummy Hist - Each bin
          break;
        }
        else{
          if(order==0) number_not_lin++;
          cout << "Bin " << number_bin[bin] << ": A "<< RED << name_fits[order][index_fit] << RESET << " fit is rejected (" << prob << ")" << endl;
        }
      }

      // #######################################################################
      // Plot ##################################################################
      if(debug) cout << "Plot Fitting" << endl;

      // -----------------------------------------------------------------------
      // Plot Points -----------------------------------------------------------
      double max_bin_content = GetMaxValue(bin_content);
      double min_bin_content = GetMinValue(bin_content);
      used_fit->SetRange(-1, -1, min_bin_content*0.8, 1, 1, max_bin_content*1.2);

      // Additional Points to plot the z axis properly -------------------------
      /*
      The additional points are necessary to scale the z axis properly. One can not zoom out of the z axis without
      cutting of the errors above the highest- and below the lowest point. The trick is to add two additional points
      and stretch the z axis. Afterward one can zoom in (and not out).
      */
      bin_fit->SetPoint(5, 0., 0., max_bin_content*1.8);
      bin_fit->SetPointError(5, 0., 0., 0.);
      bin_fit->SetPoint(6, 0., 0., min_bin_content*0.2);
      bin_fit->SetPointError(6, 0., 0., 0.);

      bin_fit->SetTitle(mass_bin_title[bin]);
      bin_fit->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
      bin_fit->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
      bin_fit->GetHistogram()->GetZaxis()->SetTitle("#DeltaN/N");
      bin_fit->GetHistogram()->GetXaxis()->SetTitleOffset(2.5);
      bin_fit->GetHistogram()->GetYaxis()->SetTitleOffset(2.5);
      bin_fit->GetHistogram()->GetZaxis()->SetTitleOffset(2.6);
      bin_fit->GetHistogram()->GetXaxis()->SetTitleSize(0.026);
      bin_fit->GetHistogram()->GetYaxis()->SetTitleSize(0.026);
      bin_fit->GetHistogram()->GetZaxis()->SetTitleSize(0.026);
      bin_fit->SetMarkerStyle(20);

      // -----------------------------------------------------------------------
      // PolyLine to Mark Points in 2D Graph -----------------------------------
      TPolyLine3D *xup= new TPolyLine3D(2);
      xup->SetPoint(0, 0, 1, min_bin_content*0.8); // min_bin_content*0.8
      xup->SetPoint(1, 0, 1, bin_content[4]-bin_content_err[4]);
      xup->SetLineStyle(2);
      xup->SetLineColor(kRed);
      xup->SetLineWidth(1);

      TPolyLine3D *xdo= new TPolyLine3D(2);
      xdo->SetPoint(0, 0, -1, min_bin_content*0.8); // min_bin_content*0.8
      xdo->SetPoint(1, 0, -1, bin_content[3]-bin_content_err[3]);
      xdo->SetLineStyle(2);
      xdo->SetLineColor(kRed);
      xdo->SetLineWidth(1);

      TPolyLine3D *jup= new TPolyLine3D(2);
      jup->SetPoint(0, 1, 0, min_bin_content*0.8); // min_bin_content*0.8
      jup->SetPoint(1, 1, 0, bin_content[2]-bin_content_err[2]);
      jup->SetLineStyle(2);
      jup->SetLineColor(kRed);
      jup->SetLineWidth(1);

      TPolyLine3D *jdo= new TPolyLine3D(2);
      jdo->SetPoint(0, -1, 0, min_bin_content*0.8); // min_bin_content*0.8
      jdo->SetPoint(1, -1, 0, bin_content[1]-bin_content_err[1]);
      jdo->SetLineStyle(2);
      jdo->SetLineColor(kRed);
      jdo->SetLineWidth(1);

      TPolyLine3D *mid= new TPolyLine3D(2);
      mid->SetPoint(0, 0, 0, min_bin_content*0.8); // min_bin_content*0.8
      mid->SetPoint(1, 0, 0, bin_content[0]-bin_content_err[0]);
      mid->SetLineStyle(2);
      mid->SetLineColor(kRed);
      mid->SetLineWidth(1);

      TPolyLine3D *xax= new TPolyLine3D(2);
      xax->SetPoint(0, -1.05, 0, min_bin_content*0.8); // min_bin_content*0.8
      xax->SetPoint(1,  1.05, 0, min_bin_content*0.8);
      xax->SetLineStyle(3);
      xax->SetLineColor(kBlack);
      xax->SetLineWidth(1);
      TPolyLine3D *yax= new TPolyLine3D(2);
      yax->SetPoint(0, 0, -1.05, min_bin_content*0.8); // min_bin_content*0.8
      yax->SetPoint(1, 0,  1.05, min_bin_content*0.8);
      yax->SetLineStyle(3);
      yax->SetLineColor(kBlack);
      yax->SetLineWidth(1);

      // -----------------------------------------------------------------------

      TCanvas *B = new TCanvas(number_bin[bin]+"B"+ptbin_str, "B", 800, 800); // name used to avoid Warning
      bin_fit->GetHistogram()->GetXaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetYaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetZaxis()->SetRangeUser(min_bin_content*0.8, max_bin_content*1.2);

      B->SetRightMargin(0.09);
      B->SetLeftMargin(0.15);

      bin_fit->Draw("err P");
      jup->Draw("same");
      jdo->Draw("same");
      mid->Draw("same");
      xup->Draw("same");
      xdo->Draw("same");
      yax->Draw("same");
      xax->Draw("same");

      B->SaveAs(save_path_general+"/single_bins/"+addition+"/Bin"+number_bin[bin]+"_Points.pdf"); // Double / (//) does not affect the dir
      used_fit->Draw("same surf4");
      B->SaveAs(save_path_general+"/single_bins/"+addition+"/Bin"+number_bin[bin]+"_Combined.pdf");
      B->Clear();
      used_fit->SetRange(-2, -2, -1, 2, 2, 1);
      used_fit->Draw("surf4");
      B->SaveAs(save_path_general+"/single_bins/"+addition+"/Bin"+number_bin[bin]+"_Fit.pdf");

      B->Clear();

      // Remove Additional Points for fit --------------------------------------
      bin_fit->RemovePoint(6); // order is on purpose
      bin_fit->RemovePoint(5);

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      if(debug) cout << "HeatMap" << endl;

      TH2F* heat_map = new TH2F("", "", 3, -1.5, 1.5, 3, -1.5, 1.5);
      heat_map->SetBinContent(1,1,0);
      heat_map->SetBinContent(1,2,XConeup_rebin_norm->GetBinContent(bin));
      heat_map->SetBinContent(1,3,0);
      heat_map->SetBinContent(2,1,JECup_rebin_norm->GetBinContent(bin));
      heat_map->SetBinContent(2,2,ttbar_rebin_norm->GetBinContent(bin));
      heat_map->SetBinContent(2,3,JECdown_rebin_norm->GetBinContent(bin));
      heat_map->SetBinContent(3,1,0);
      heat_map->SetBinContent(3,2,XConeup_rebin_norm->GetBinContent(bin));
      heat_map->SetBinContent(3,3,0);

      TCanvas *B1 = new TCanvas("B1", "", 800, 800); // name used to avoid Warning
      bin_fit->GetHistogram()->GetXaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetYaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetZaxis()->SetRangeUser(min_bin_content*0.8, max_bin_content*1.2);

      B1->SetRightMargin(0.09);
      B1->SetLeftMargin(0.15);
      B1->SetRightMargin(0.20);

      double tt_bin = ttbar_rebin_norm->GetBinContent(bin);
      heat_map->GetZaxis()->SetRangeUser(tt_bin-tt_bin*0.15, tt_bin+tt_bin*0.15);

      gStyle->SetPalette(kCoffee);
      heat_map->Draw("COLZ TEXT ERROR");

      B1->SaveAs(save_path_general+"/single_bins/"+addition+"/HeatMap_"+number_bin[bin]+".pdf");
      B1->Clear();
      gStyle->SetPalette(kDeepSea);

      // delete to avoid same name
      delete heat_map;

      // #################################################################################################
      // Fill Fit Parameters #############################################################################
      if(debug) cout << "Fill Fit Parameters" << endl;

      // Data Bin Content in Vector ############################################
      if(debug) cout << "Data Bin Content in Vector" << endl;
      VecD data_rebin_norm_bin_content;
      for(int bin=1; bin < number_bins+1; bin++) data_rebin_norm_bin_content.push_back(data_rebin_norm->GetBinContent(bin));

      /*
      .███████ ███████ ████████     ██████   █████  ██████
      .██      ██         ██        ██   ██ ██   ██ ██   ██
      .███████ █████      ██        ██████  ███████ ██████
      .     ██ ██         ██        ██      ██   ██ ██   ██
      .███████ ███████    ██        ██      ██   ██ ██   ██
      */

      // Fill Parameter for Chi2 ###############################################
      VecD fit_para; // Ablage - Cannot use chi2_para, since it contains ALL parameters
      double sigma_tot_square;
      if(onlyData) sigma_tot_square = pow(data_rebin_norm_err[bin], 2);
      else         sigma_tot_square = pow(data_rebin_norm_err[bin], 2)+pow(ttbar_rebin_norm_err[bin],2);

      fit_para.push_back(data_rebin_norm_bin_content[bin]);
      chi2_parameters.push_back(data_rebin_norm_bin_content[bin]);
      for(int ipar=0; ipar<n_fit_para[order_index][index_fit]; ipar++){
        fit_para.push_back(used_fit->GetParameter(ipar));
        chi2_parameters.push_back(used_fit->GetParameter(ipar));
      }
      fit_para.push_back(sigma_tot_square);
      chi2_parameters.push_back(sigma_tot_square);

      fits_parameters.push_back(fit_para);
    }
    cout << endl;

    /*
    ███████ ██ ████████     ██       ██████   ██████  ██████      ███████ ███    ██ ██████
    ██      ██    ██        ██      ██    ██ ██    ██ ██   ██     ██      ████   ██ ██   ██
    █████   ██    ██        ██      ██    ██ ██    ██ ██████      █████   ██ ██  ██ ██   ██
    ██      ██    ██        ██      ██    ██ ██    ██ ██          ██      ██  ██ ██ ██   ██
    ██      ██    ██        ███████  ██████   ██████  ██          ███████ ██   ████ ██████
    */

    int npar2 = 0;
    for(unsigned int fit=0;fit<fits_parameters.size(); fit++){
      for(unsigned int ipar=0;ipar<fits_parameters[fit].size(); ipar++) npar2++;
    }
    if(debug) cout << "Number Parameter: " << npar2 << endl;
    if(debug) cout << "Chi2 Parameter:   " << chi2_parameters.size() << endl;


    // #################################################################################################
    // Hists info fits #################################################################################
    used_fits.SetMarkerStyle(8);  // data hist style
    used_fits.SetMarkerColor(kBlack);
    used_fits.SetTitle("");
    used_fits.GetXaxis()->SetTitle("Bin");
    used_fits.GetYaxis()->SetTitle("Fit");
    used_fits.GetYaxis()->SetRangeUser(0, 7);
    TCanvas *H = new TCanvas("H"+ptbin_str, "H", 600, 600);
    used_fits.Draw("P");
    leg = new TLegend(0.6,0.65,0.8,0.85);
    for(unsigned int i=0; i<str_fits.size();i++) leg->AddEntry((TObject*)0,to_string(i+1)+": "+str_fits[i][0],"");
    leg->SetTextSize(0.015);
    leg->Draw();
    gPad->RedrawAxis();
    H->SaveAs(save_path_general+"/used_fits_per_bin"+addition+".pdf");
    leg->Clear();

    number_fits.SetTitle("");
    number_fits.GetXaxis()->SetTitle("Fit");
    number_fits.GetYaxis()->SetTitle("Number");
    TCanvas *H1 = new TCanvas("H1"+ptbin_str, "H1", 600, 600);
    number_fits.Draw("Hist");
    leg = new TLegend(0.4,0.65,0.7,0.85);
    for(unsigned int i=0; i<str_fits.size();i++) leg->AddEntry((TObject*)0,to_string(i)+": "+str_fits[i][0],"");
    leg->SetTextSize(0.025);
    leg->Draw();
    gPad->RedrawAxis();
    H1->SaveAs(save_path_general+"/number_fits"+addition+".pdf");
    leg->Clear();

    /*
    .██████ ██   ██ ██ ██████
    ██      ██   ██ ██      ██
    ██      ███████ ██  █████
    ██      ██   ██ ██ ██
    .██████ ██   ██ ██ ███████
    */
    if(debug) cout << "\nChi2 - Parameters\n";

    TString str_chi2_function = ""; TString plus = " + ";
    TString par0 = "[0]"; TString par1 = "[1]"; TString par2 = "[2]"; TString par3 = "[3]"; TString par4 = "[4]"; // Max only 5 para.
    vector<TString> parameters = {"[0]", "[1]", "[2]", "[3]", "[4]"};

    // #################################################################################################
    // Replace Parameters ##############################################################################
    if(debug) cout << "Chi2 - Replace parameters and build string\n";
    if(!usePeak && str_all_fits.size()!= number_bins-number_empty_bins) throw runtime_error("Something is wrong with the number of fits");

    int count_para = 0;
    for(unsigned int fit=0; fit<str_all_fits.size(); fit++){
      TString str_fit_full = " + ((["+to_string(count_para)+"] - ";
      TString str_fit = str_all_fits[fit];
      /*
      Treat the first fit differently. Replacing the parametes leads to an overlap.
      [1]-[1]*x-[2]*y -> [2]-[2]*x-[2]*y -> [3]-[3]*x-[3]*y -> etc.
      in "Chi2 All" it is done backwards. IMPLEMENT HERE, TOO!
      */
      if(fit==0){
        if(!str_fit.Contains("y*y")){
          str_fit = "[1] - [2]*x - [3]*y";
          count_para = 3;
        }
        else if(str_fit.Contains("y*y") && !str_fit.Contains("]*y ") && !str_fit.Contains("x*x")){
          str_fit = "[1] - [2]*x - [3]*y*y";
          count_para = 3;
        }
        else if(str_fit.Contains("y*y") && str_fit.Contains("]*y ") && !str_fit.Contains("x*x")){
          str_fit = "[1] - [2]*x - [3]*y - [4]*y*y";
          count_para = 4;
        }
        else if(str_fit.Contains("y*y") && !str_fit.Contains("]*y ") && str_fit.Contains("x*x")){
          str_fit = "[1] - [2]*x*x - [3]*y*y";
          count_para = 3;
        }
        else if(str_fit.Contains("y*y") && str_fit.Contains("]*y ") && str_fit.Contains("x*x")){
          str_fit = "[1] - [2]*x - [3]*y - [4]*x*x - [5]*y*y";
          count_para = 5;
        }
      }
      else{
        for(unsigned int ipar=0; ipar<parameters.size();ipar++){
          if(str_fit.Contains(parameters[ipar])){
            if(fit!=0){
              count_para++;
              str_fit.ReplaceAll(parameters[ipar], "["+to_string(count_para)+"]");
            }
            str_fit.ReplaceAll("+", "-");
          }
          else break; // Each fit contains the first three parameters
        }
      }
      str_fit_full += str_fit;

      count_para++; // data bin err
      str_fit_full += "))^2/["+to_string(count_para)+"]";
      count_para++; // data bin content - for next iteration
      if(debug) cout << str_fit_full << "\n";
      str_chi2_function += str_fit_full;
    }

    // #################################################################################################
    // Set all parameters ##############################################################################
    if(debug) cout << "Chi2 - Define Funtion and Set Parameter\n";
    TF2 *chi2_function = new TF2("chi2_function"+addition, str_chi2_function, -10, 10, -10, 10);
    for(unsigned int ipar=0; ipar<chi2_parameters.size(); ipar++) chi2_function->SetParameter(ipar, chi2_parameters[ipar]);
    chi2_function->SetTitle("");
    chi2_function->GetXaxis()->SetTitle("JEC");
    chi2_function->GetYaxis()->SetTitle("Additional XCone correction");
    chi2_function->SetFillStyle(1000);
    chi2_function->SetLineWidth(1);

    // #################################################################################################
    // Add chi2 functions ##############################################################################
    if(debug) cout << "Chi2 - Safe Chi2 from single ptbins for Chi2-All\n";
    if(ptbin==0){
      chi2_parameters_hh = chi2_parameters;
      chi2_str_hh = str_chi2_function;
    }
    else if(ptbin==1){
      chi2_parameters_hl = chi2_parameters;
      chi2_str_hl = str_chi2_function;
    }
    else if(ptbin==2){
      chi2_parameters_lh = chi2_parameters;
      chi2_str_lh = str_chi2_function;
    }
    else if(ptbin==3){
      chi2_parameters_ll = chi2_parameters;
      chi2_str_ll = str_chi2_function;
    }

    // #################################################################################################
    // Minimum #########################################################################################
    if(debug) cout << "Chi2 - Get Minimum of Chi2\n";
    /* Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447. */
    double twoD_minX, twoD_minY;
    double twoD_minZ = chi2_function->GetMinimumXY(twoD_minX,twoD_minY);
    cout  << "z_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n' << endl;
    int number_para = chi2_function->GetNpar();

    fstream jec_txt;
    jec_txt.open(save_path_general+"/jec_factor"+addition+".txt", ios::out);
    jec_txt << number_bins << " bins [0, 180] GeV\n\n";
    jec_txt << "zmin: " << twoD_minZ << "\n(" << setw(7) << Round(twoD_minX, 4) << ", " << setw(7) << Round(twoD_minY, 4) << ")" << endl;
    jec_txt.close();

    // #################################################################################################
    // Evaluate 1sigma #################################################################################
    if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";

    // Get Points --------------------------------------------------------------
    auto start = high_resolution_clock::now(); // Calculation time - start
    VecDD points = FindXY(chi2_function, twoD_minZ+2.3, twoD_minX-2, twoD_minX+2, twoD_minY-2, twoD_minY+2, 1000, 0.001);
    auto stop  = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Numeric solution for 1\u03C3 estimation took " << GREEN << duration.count()/1000 << "s" << RESET << endl;
    if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;

    // Draw Points -------------------------------------------------------------
    if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";
    TGraph2D *zmin_point = new TGraph2D();
    zmin_point->SetPoint(0, twoD_minX, twoD_minY, twoD_minZ);
    zmin_point->SetMarkerColor(kBlack);
    zmin_point->SetMarkerStyle(kFullCircle);
    zmin_point->SetMarkerSize(0.5);

    TGraph2D *sigma_points = new TGraph2D();
    for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
    sigma_points->SetMarkerColor(kRed);
    sigma_points->SetMarkerStyle(kFullCircle);
    sigma_points->SetMarkerSize(0.1);

    chi2_function->SetTitle("");
    chi2_function->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
    chi2_function->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
    chi2_function->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
    chi2_function->GetHistogram()->GetZaxis()->CenterTitle();
    chi2_function->GetHistogram()->GetZaxis()->SetTitleOffset(0.5);
    chi2_function->SetFillStyle(1000);
    chi2_function->SetLineWidth(1);
    chi2_function->SetRange(-3, -3, 3, 3);

    // #################################################################################################
    // Plot ############################################################################################
    if(debug) cout<< "\n2D Chi2 - plot\n";

    TString option = "cont4z";
    TString option_add = "";

    chi2_function->SetContour(50);        // Contours
    gStyle->SetPalette(kDeepSea);         // kDeepSea kGreyScale kRust
    if(ptbin==0) TColor::InvertPalette(); // CHANGE_PT
    // if(ptbin==4) TColor::InvertPalette(); // CHANGE_PT

    TCanvas *D = new TCanvas("D"+ptbin_str,"D", 600, 600);
    D->SetRightMargin(0.12);
    D->SetLogz();
    chi2_function->Draw("cont4z");
    D->SaveAs(save_path_general+"/chi2_cont4z"+addition+".pdf");
    D->SetTheta(90);
    D->SetPhi(0);
    zmin_point->Draw("SAME P");
    sigma_points->Draw("SAME P");
    D->SaveAs(save_path_general+"/chi2_cont4z_points"+addition+".pdf");
    D->Clear();

    /*
    ██████  ██████   ██████       ██ ███████  ██████ ████████ ██  ██████  ███    ██
    ██   ██ ██   ██ ██    ██      ██ ██      ██         ██    ██ ██    ██ ████   ██
    ██████  ██████  ██    ██      ██ █████   ██         ██    ██ ██    ██ ██ ██  ██
    ██      ██   ██ ██    ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██
    ██      ██   ██  ██████   █████  ███████  ██████    ██    ██  ██████  ██   ████
    */
    cout << "\n####################################################################\n";
    cout << "######################### Projection ###############################\n";
    cout << "####################################################################\n";

    vector<TString> factor_name = {"JEC", "XC"};
    // #################################################################################################
    // Which correction ################################################################################
    for(int factor=0; factor<2; factor++){
      if(debug) cout << "------------------------------------------------------- " << factor_name[factor] << endl;
      VecD oneD_factor     = {-1,0,1};
      VecD oneD_factor_err = { 0,0,0};

      // #################################################################################################
      // Bins ############################################################################################
      int empty_bins = 0; // to skip empty bins
      for(int bin=0; bin<number_bins; bin++){
        // if(debug) cout << "------------------------------------------------------- " << bin+1 << endl;
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())){ // Search_peackmethod
          empty_bins++;
          continue;
        }
        if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()){
          empty_bins++;
          continue;
        }

        int used_bin = bin-empty_bins;
        VecD oneD_events, oneD_events_err;
        if(factor==0){
          oneD_events     = {hists_cont[used_bin][1], hists_cont[used_bin][0], hists_cont[used_bin][2]};
          oneD_events_err = {hists_cont_err[used_bin][1], hists_cont_err[used_bin][0], hists_cont_err[used_bin][2]};
        } else if(factor==1){
          oneD_events     = {hists_cont[used_bin][3], hists_cont[used_bin][0], hists_cont[used_bin][4]};
          oneD_events_err = {hists_cont_err[used_bin][3], hists_cont_err[used_bin][0], hists_cont_err[used_bin][4]};
        }
        // XCone -----------------------------------------------------------------------------------------
        TString str_fit = str_all_fits[used_bin];
        if(factor==0) str_fit.ReplaceAll("y", "0");
        else if(factor==1){
          str_fit.ReplaceAll("x", "0");
          str_fit.ReplaceAll("y", "x"); // 1D Graph only takes x as argument
        }

        TF1 *fit = new TF1("fit", str_fit, -2, 2);
        TGraphErrors *graph = new TGraphErrors(3, &oneD_factor[0], &oneD_events[0], &oneD_factor_err[0], &oneD_events_err[0]);
        for(int ipar=0; ipar<all_n_para[used_bin]; ipar++) fit->SetParameter(ipar, all_fits[used_bin]->GetParameter(ipar));

        // #################################################################################################
        // Plot ############################################################################################
        graph->SetTitle(mass_bin_title[bin]);
        TCanvas *E2 = new TCanvas(number_bin[bin]+"E2_"+factor_name[factor]+to_string(ptbin), "E2", 600, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        if(factor==0) graph->GetXaxis()->SetTitle("JEC factor");
        else          graph->GetXaxis()->SetTitle("Additional XCone correction factor");
        graph->GetYaxis()->SetTitle("#Delta N/N");
        graph->GetYaxis()->SetTitleOffset(2.2);

        graph->Draw("AP");
        fit->Draw("SAME");
        E2->SaveAs(save_path_general+"/projection/Bin"+number_bin[bin]+"_"+factor_name[factor]+addition+".pdf");
      }
    }

    /*
    ██       █████  ████████ ███████ ██   ██
    ██      ██   ██    ██    ██       ██ ██
    ██      ███████    ██    █████     ███
    ██      ██   ██    ██    ██       ██ ██
    ███████ ██   ██    ██    ███████ ██   ██
    */

    if(into_latex){
      cout << "\n####################################################################\n";
      cout << "############################ Latex #################################\n";
      cout << "####################################################################\n";
      // suppress output with &> /dev/null
      if(usePeak){
        if(debug) system("python python/single_bins_latex.py "+save_path_general+" "+addition+" projection");
        else      system("python python/single_bins_latex.py "+save_path_general+" "+addition+" projection &> /dev/null");
      } else{
        if(debug) system("python python/single_bins_latex.py "+save_path_general+" "+addition+" projection");
        else      system("python python/single_bins_latex.py "+save_path_general+" "+addition+" projection &> /dev/null");
      }
    }
  }

  cout << "\nNumber bins:" << setw(4) << number_bins_total << endl;
  cout << "Not linear: "   << setw(4) << number_not_lin    << " (" << (number_not_lin/number_bins_total)*100 << "%)" << endl;
  if(!pt_bins) return 0;
  /*
  .█████  ██████  ██████       ██████ ██   ██ ██ ██████
  ██   ██ ██   ██ ██   ██     ██      ██   ██ ██      ██
  ███████ ██   ██ ██   ██     ██      ███████ ██  █████
  ██   ██ ██   ██ ██   ██     ██      ██   ██ ██ ██
  ██   ██ ██████  ██████       ██████ ██   ██ ██ ███████
  */
  if(debug) cout<< "\nAll 2D Chi2 - Start";
  cout << endl;

  vector<TString>        all_string     = {chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll};
  VecDD all_parameters = {chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll};

  // #################################################################################################
  // Get String ######################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - Build string";

  TString full_chi2 = chi2_str_hh;
  int count_para = 0;
  for(int function=0; function<4; function++){
    if(function == 0) continue;

    count_para += all_parameters[function-1].size();
    for(int ipar=all_parameters[function].size()-1; ipar>=0; ipar--){
      TString ipar_string = "["+to_string(ipar)+"]";

      if(all_string[function].Contains(ipar_string)){
        all_string[function] = all_string[function].ReplaceAll(ipar_string, "["+to_string(ipar+count_para)+"]");
      }
    }
    full_chi2 += all_string[function];
  }
  TF2 *full_chi2_function = new TF2("chi2_function", full_chi2, -3, 3, -3, 3);

  // #################################################################################################
  // Set Parameter ###################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - Set Parameter";
  cout << endl;
  double ipar = 0;
  for(unsigned int function=0; function<all_parameters.size(); function++){
    for(unsigned int parameter=0; parameter<all_parameters[function].size(); parameter++){
      full_chi2_function->SetParameter(ipar, all_parameters[function][parameter]);
      ipar++;
    }
  }

  // #################################################################################################
  // Minimum #########################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - Get Minimum\n";
  /* Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447. */
  double twoD_minX, twoD_minY;
  double twoD_minZ = full_chi2_function->GetMinimumXY(twoD_minX,twoD_minY);
  cout  << "z_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n';
  int number_para = full_chi2_function->GetNpar();

  fstream jec_txt;
  jec_txt.open(save_path_general+"/jec_factor_all.txt", ios::out);
  jec_txt << number_bins << " bins [0, 180] GeV\n\n";
  jec_txt << "zmin: " << twoD_minZ << "\nx   :  " << twoD_minX << "\ny   :  " << twoD_minY << endl;
  jec_txt.close();

  // #################################################################################################
  // Evaluate 1sigma #################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - Ellipse\n";

  // Get Points --------------------------------------------------------------------------------------
  VecDD points;
  auto start = high_resolution_clock::now(); // Calculation time - start
  points = FindXY(full_chi2_function, twoD_minZ+2.3, twoD_minX-1.5, twoD_minX+1.5, twoD_minY-1.5, twoD_minY+1.5, 1000, 0.1, true);
  auto stop = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  cout << "Numeric solution for 1\u03C3 estimation took " << GREEN << duration.count()/1000 << "s" << RESET << endl;
  if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;

  // Draw Points -------------------------------------------------------------------------------------
  TGraph2D *zmin_point = new TGraph2D();
  zmin_point->SetPoint(0, twoD_minX, twoD_minY, twoD_minZ);
  zmin_point->SetMarkerColor(kBlack);
  zmin_point->SetMarkerStyle(kFullCircle);
  zmin_point->SetMarkerSize(0.5);

  TGraph2D *sigma_points = new TGraph2D();
  for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
  sigma_points->SetMarkerColor(kRed);
  sigma_points->SetMarkerStyle(kFullCircle);
  sigma_points->SetMarkerSize(0.1);

  // #################################################################################################
  // Plot ############################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - plot\n";

  TString option = "cont4z";
  TString option_add = "";
  full_chi2_function->SetTitle("");

  full_chi2_function->GetXaxis()->SetTitle("JEC factor");
  full_chi2_function->GetYaxis()->SetTitle("Additional XCone correction factor");
  full_chi2_function->GetZaxis()->SetTitle("#chi^{2}");
  full_chi2_function->GetZaxis()->CenterTitle();
  full_chi2_function->GetZaxis()->SetTitleOffset(0.5);

  full_chi2_function->SetContour(60);   // Contours
  gStyle->SetPalette(kDeepSea);         // kDeepSea kGreyScale kRust
  TColor::InvertPalette(); // CHANGE_PT

  TCanvas *Z = new TCanvas("Z","Z", 600, 600);
  Z->SetRightMargin(0.12);
  Z->SetLogz();
  full_chi2_function->Draw("cont4z");
  full_chi2_function->SetMaximum(1100); // One has to redefine the axis after drawing
  full_chi2_function->SetMinimum(200);
  Z->SaveAs(save_path_general+"/chi2_all.pdf");
  Z->SetTheta(90);
  Z->SetPhi(0);
  zmin_point->Draw("SAME P");
  sigma_points->Draw("SAME P");
  Z->SaveAs(save_path_general+"/chi2_all_points.pdf");
  Z->Clear();

  /*
  ██████  ██████   ██████       ██ ███████  ██████ ████████ ██  ██████  ███    ██
  ██   ██ ██   ██ ██    ██      ██ ██      ██         ██    ██ ██    ██ ████   ██
  ██████  ██████  ██    ██      ██ █████   ██         ██    ██ ██    ██ ██ ██  ██
  ██      ██   ██ ██    ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██
  ██      ██   ██  ██████   █████  ███████  ██████    ██    ██  ██████  ██   ████
  */

  TString full_chi2_jec = full_chi2;
  TString full_chi2_nan = full_chi2;
  TString full_chi2_cor = full_chi2;

  full_chi2_jec.ReplaceAll("y", "("+to_string(twoD_minY)+")");

  full_chi2_cor.ReplaceAll("x", to_string(twoD_minX));
  full_chi2_cor.ReplaceAll("y", "x");


  TF1 *full_chi2_function_jec = new TF1("chi2_function_jec", full_chi2_jec, -3, 3);
  TF1 *full_chi2_function_cor = new TF1("chi2_function_cor", full_chi2_cor, -3, 3);

  double ipar_jec = 0;
  for(unsigned int function=0; function<all_parameters.size(); function++){
    for(unsigned int parameter=0; parameter<all_parameters[function].size(); parameter++){
      full_chi2_function_jec->SetParameter(ipar_jec, all_parameters[function][parameter]);
      ipar_jec++;
    }
  }

  SetFunctionParameters(full_chi2_function_jec, all_parameters);
  SetFunctionParameters(full_chi2_function_cor, all_parameters);

  VecD JMS_JEC = ExtractJMSValues(full_chi2_function_jec);
  VecD JMS_COR = ExtractJMSValues(full_chi2_function_cor);

  // printVector(JMS_JEC);
  // printVector(JMS_COR);

  drawChi2Projection(full_chi2_function_jec, "JEC", {twoD_minX-0.4, twoD_minX+0.4}, {215, 235}, save_path_general);
  drawChi2Projection(full_chi2_function_cor, "COR", {twoD_minY-1, twoD_minY+1}, {215, 235}, save_path_general);

  JMSinTXT(JMS_JEC[0], JMS_JEC[1], JMS_JEC[2], "JEC", save_path_general);
  JMSinTXT(JMS_COR[0], JMS_COR[1], JMS_COR[2], "COR", save_path_general);

  double sigma_jec = (JMS_JEC[1]+JMS_JEC[2])/2;
  double sigma_cor = (JMS_COR[1]+JMS_COR[2])/2;
  double sigma_jms = sqrt(pow(sigma_jec/JMS_JEC[0],2)+pow(sigma_cor/JMS_COR[0],2));

  cout << pow(sigma_cor/JMS_COR[0],2) << endl;
  cout << pow(sigma_cor/JMS_JEC[0],2) << endl;
  cout << "Error f_{JMS}: " << sigma_jms << endl;

  /*
  ███████ ██      ██      ██ ██████  ███████ ███████
  ██      ██      ██      ██ ██   ██ ██      ██
  █████   ██      ██      ██ ██████  ███████ █████
  ██      ██      ██      ██ ██           ██ ██
  ███████ ███████ ███████ ██ ██      ███████ ███████
  */

  // #################################################################################################
  // Extreme Points ##################################################################################
  if(debug) cout<< "\nFind extrme points of ellipse";

  double dmax1, dmax2, dmin1, dmin2;
  VecDD Ps_max, Ps_min;
  unsigned int n_points = points.size();
  cout << "Looking at " << n_points << " points" << endl;

  VecD xValues, yValues, dist_right, dist_left, dist_all;
  // MapVD dist_right, dist_left, dist_all;
  for(unsigned int point=0; point<n_points; point++){
    double x=points[point][0];
    double y=points[point][1];
    dist_all.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    // dist_all["x"].push_back(x);
    // dist_all["y"].push_back(y);
    if(x>twoD_minX) dist_right.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    if(x<twoD_minX) dist_left.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    xValues.push_back(x);
    yValues.push_back(y);
  }

  // #################################################################################################
  // find index ######################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - find index";

  int index_max_l_s = ReturnIndex_High(dist_left);
  int index_max_r_s = ReturnIndex_High(dist_right);
  int index_min_l_s = ReturnIndex_Low(dist_left);
  int index_min_r_s = ReturnIndex_Low(dist_right);

  cout << setw(10) << dist_left[index_max_l_s] << setw(10) << dist_left[index_min_l_s] << setw(10) << dist_right[index_max_r_s] << setw(10) << dist_right[index_min_r_s] << endl;

  int index_max_l, index_max_r, index_min_l, index_min_r;

  // Compare ------------------------------------------------------------------------------------------
  for(unsigned int i=0; i<dist_all.size(); i++){
    if(dist_left[index_max_l_s]==dist_all[i]) { index_max_l = i; cout << i << xValues[i] << " " << yValues[i] << " max l" << endl;}
    if(dist_right[index_max_r_s]==dist_all[i]){ index_max_r = i; cout << i << xValues[i] << " " << yValues[i] << " max r" << endl;}
    if(dist_left[index_min_l_s]==dist_all[i]) { index_min_l = i; cout << i << xValues[i] << " " << yValues[i] << " min l" << endl;}
    if(dist_right[index_min_r_s]==dist_all[i]){ index_min_r = i; cout << i << xValues[i] << " " << yValues[i] << " min r" << endl;}
  }
  if(debug) cout << "l-max: " << index_max_l << " | r-max: " << index_max_r << " | l-min: " << index_min_l << " | r_min: " << index_min_r << endl;

  Ps_max.push_back(points[index_max_l]);
  Ps_max.push_back(points[index_max_r]);
  Ps_max.push_back(points[index_min_l]);
  Ps_max.push_back(points[index_min_r]);

  cout << "Number Points: " << Ps_max.size() << endl;
  for(unsigned int point=0; point<Ps_max.size(); point++) cout << "(" << setw(9) << Ps_max[point][0] << ", " << setw(9) << Ps_max[point][1] << ")\n";

  // #################################################################################################
  // in file #########################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - in file";

  jec_txt.open(save_path_general+"/jec_factor_all.txt", ios::out | ios::app);
  jec_txt << setw(4) << "\nuu: " << "(" << setw(7) << Round(Ps_max[3][0], 4) << ", " << setw(7) << Round(Ps_max[3][1], 4) << ")" <<  endl;
  jec_txt << setw(4) << "dd: "   << "(" << setw(7) << Round(Ps_max[2][0], 4) << ", " << setw(7) << Round(Ps_max[2][1], 4) << ")" <<  endl;
  jec_txt << setw(4) << "ud: "   << "(" << setw(7) << Round(Ps_max[1][0], 4) << ", " << setw(7) << Round(Ps_max[1][1], 4) << ")" <<  endl;
  jec_txt << setw(4) << "du: "   << "(" << setw(7) << Round(Ps_max[0][0], 4) << ", " << setw(7) << Round(Ps_max[0][1], 4) << ")" <<  endl;
  jec_txt.close();

  // #################################################################################################
  // Settings ########################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - settings";

  TGraph2D *extreme_points = new TGraph2D();
  for(unsigned int i=0; i<Ps_min.size(); i++) extreme_points->SetPoint(i, Ps_min[i][0], Ps_min[i][1], Ps_min[i][2]);
  for(unsigned int i=0; i<Ps_max.size(); i++) extreme_points->SetPoint(i, Ps_max[i][0], Ps_max[i][1], Ps_max[i][2]);

  extreme_points->SetMarkerColor(kBlue);
  extreme_points->SetMarkerStyle(kFullCircle);
  extreme_points->SetMarkerSize(0.4);


  // #################################################################################################
  // Plot ############################################################################################
  if(debug) cout<< "\nAll 2D Chi2 - plot";

  gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
  TCanvas *Z1 = new TCanvas("Z1","Z1", 600, 600);
  Z1->SetRightMargin(0.12);
  Z1->SetLogz();
  full_chi2_function->Draw("cont4z");
  Z1->SetTheta(90);
  Z1->SetPhi(0);
  zmin_point->Draw("SAME P");
  sigma_points->Draw("SAME P");
  extreme_points->Draw("SAME P");
  Z1->SaveAs(save_path_general+"/chi2_all_extreme_points.pdf");
  Z1->Clear();
}





// #################################################################################################
// #################################################################################################
// #################################################################################################

VecD ExtractJMSValues(TF1* chi2function){
  double minChi2   = chi2function->GetMinimum();
  double f_JMS     = chi2function->GetX(minChi2, -3, 3);
  double f_up      = chi2function->GetX(minChi2+1, f_JMS, 3);
  double f_down    = chi2function->GetX(minChi2+1, -3, f_JMS);

  double sigmaup   = f_up - f_JMS;
  double sigmadown = f_JMS - f_down;

  VecD values      = {f_JMS, sigmaup, sigmadown, f_up, f_down, minChi2};
  return values;
}

// ------------------------------------------------------------------------------------
void SetFunctionParameters(TF1* chi2function, VecDD parameters){
  double ipar = 0;
  for(unsigned int function=0; function<parameters.size(); function++){
    for(unsigned int parameter=0; parameter<parameters[function].size(); parameter++){
      chi2function->SetParameter(ipar, parameters[function][parameter]);
      ipar++;
    }
  }
}

// ------------------------------------------------------------------------------------
void drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange, TString save_path){
  TCanvas* ccor = new TCanvas(xaxis, xaxis, 600, 600);
  chi2function->SetTitle("");
  chi2function->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
  chi2function->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
  chi2function->GetXaxis()->SetTitle("f^{JMS}_{"+xaxis+"}");
  chi2function->GetYaxis()->SetTitle("#chi^{2}");
  chi2function->GetYaxis()->SetTitleOffset(1.2);
  chi2function->Draw();
  ccor->SaveAs(save_path+"/Chi2_projection_"+xaxis+".pdf");
}

// ------------------------------------------------------------------------------------
void JMSinTXT(double nom, double up, double down, TString variation, TString save){
  ofstream textfile;
  textfile.open(save+"/factor_projection_"+variation+".txt");
  textfile << right;
  textfile << setw(10) << "nom" << setw(10) << "up" << setw(10) << "down" << "\n";
  textfile << setw(10) <<  nom  << setw(10) <<  up  << setw(10) <<  down;
  textfile << "\n";
  textfile.close();
}

// ------------------------------------------------------------------------------------
void printVector(VecD vector){
  cout << endl;
  for(auto value: vector) cout << setw(10) << value << " ";
  cout << endl;
}
