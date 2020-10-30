#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  print_seperater();

  bool debug = false;
  cout.precision(6);

  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_comparison_uncert <year>";
    return 0;
  }

  // #################################################################################################
  // Print Options ###################################################################################
  Int_t oldLevel    = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;          // suppress TCanvas output

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;

  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16  = true;
  else if(strcmp(year, "2017")==0)     is17  = true;
  else if(strcmp(year, "2018")==0)     is18  = true;
  else if(strcmp(year, "combined")==0) isAll = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or combined)");

  /*
  ██████  ██ ██████  ███████  ██████ ████████  ██████  ██████  ██ ███████ ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██      ██
  ██   ██ ██ ██████  █████   ██         ██    ██    ██ ██████  ██ █████   ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██           ██
  ██████  ██ ██   ██ ███████  ██████    ██     ██████  ██   ██ ██ ███████ ███████
  */

  // #################################################################################################
  // creating subdirectories. Necessary, because different binning examined ##########################
  /*
  Explanation for mkdir(char **, mode_t mode) and mode_t: https://jameshfisher.com/2017/02/24/what-is-mode_t/
  example:    0002 is -------w-
  .           0003 is -------wx
  .           0004 is ------r--
  .           0005 is ------r-x
  .           0006 is ------rw-
  .           0007 is ------rwx
  positions:  xyzw - x: type (folder, etc.) not necessary here;- y: owner;- z: group;- w: other;
  */

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS";
  save_path_general = creat_folder_and_path(save_path_general, "uncert");
  save_path_general = creat_folder_and_path(save_path_general, "pt_bins");
  save_path_general = creat_folder_and_path(save_path_general, year);

  // #################################################################################################
  // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
  int number_bins = 180;
  if(debug) cout << "Title for bins" << endl;
  vector<TString> mass_bin;
  for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin)+" < m_{Wjet} < "+to_string(bin+1));

  if(debug){
    cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << 1 << "\n\n";
    for(unsigned int bin=0; bin<mass_bin.size(); bin++){
      cout << mass_bin[bin] << '\n';
    }
  }

  vector<TString> number_bin;
  for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

  // #################################################################################################
  // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
  TString addition="";
  TString chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll;
  vector<double> chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll;
  vector<vector<vector<double>>> all_bins_combi_err;
  vector<vector<TH1F*>>          all_bins_combi;
  vector<vector<int>>            all_peak_bins;

  // Used in FIT LOOP to count number of overlaping points in total ----------------------------------
  int overlap_n6=0    ; int overlap_n7=0    ; int overlap_67=0    ; int overlap_all=0;
  int overlap_sys_n6=0; int overlap_sys_n7=0; int overlap_sys_67=0; int overlap_sys_all=0;
  int total_points=0  ; int total_points_sys=0;

  /*
  .███████ ████████  █████  ██████  ████████     ██████  ████████ ██████  ██ ███    ██
  .██         ██    ██   ██ ██   ██    ██        ██   ██    ██    ██   ██ ██ ████   ██
  .███████    ██    ███████ ██████     ██        ██████     ██    ██████  ██ ██ ██  ██
  .     ██    ██    ██   ██ ██   ██    ██        ██         ██    ██   ██ ██ ██  ██ ██
  .███████    ██    ██   ██ ██   ██    ██        ██         ██    ██████  ██ ██   ████
  */

  for(int ptbin=0; ptbin<4; ptbin++){
    // print_seperater();
    // print_seperater();
    // cout << "------------------------------------ " << ptbin << endl;
    // print_seperater();
    // print_seperater();

    if(ptbin==0) addition="_hh";
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="_full"; // not included yet - Break Segmentation
    creat_folder(save_path_general, addition);

    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */

    cout << '\n';
    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";

    if(ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";
    if(ptbin==4) w_mass = hist_class+"wmass_match";


    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    vector<TFile*> file_bkg_v;
    vector<TH1F*> hists_bkg_v;
    vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
    for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+year+"/muon/"+path_bkg_v[i]));
    for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
    TH1F *bkg = AddHists(hists_bkg_v, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TString data_path        = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    TFile  *data_file        = new TFile(dir+year+"/muon/"+data_path);
    TH1F   *data             = (TH1F*)data_file->Get(w_mass);

    TH1F   *data_norm            = normalize(data);
    vector<double> data_norm_err = normalize_error(data);

    // Masspeack ---------------------------------------------------------------
    if(debug) cout << "Masspeak Bins" << endl;
    double Limit = 75;
    vector<int> peak_bins_v = bins_upper_limit(data, Limit); // Get bins with bin-content>Limit
    if(debug) cout << "Masspeakbins: " << peak_bins_v.size() << endl;
    all_peak_bins.push_back(peak_bins_v);

    /*
    .████████ ████████ ██████   █████  ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ███████ ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ██   ██ ██   ██
    */
    // #################################################################################################
    // TTbar ###########################################################################################

    if(debug) cout << "TTbar" << endl;

    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile *ttbar_file  = new TFile(dir+year+"/muon/"+ttbar_path);
    TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);
    ttbar->Add(bkg, 1);
    TH1F* ttbar_norm              = normalize(ttbar);
    vector<double> ttbar_norm_err = normalize_error(ttbar);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC" << endl;

    TFile *JECup_file   = new TFile(dir+year+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"); // CHANGE_NORMAL
    TFile *JECdown_file = new TFile(dir+year+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"); // CHANGE_NORMAL
    TH1F  *JECup        = (TH1F*)JECup_file->Get(w_mass);
    TH1F  *JECdown      = (TH1F*)JECdown_file->Get(w_mass);
    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);
    TH1F* JECup_norm                = normalize(JECup);
    TH1F* JECdown_norm              = normalize(JECdown);
    vector<double> JECup_norm_err   = normalize_error(JECup);
    vector<double> JECdown_norm_err = normalize_error(JECdown);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone" << endl;

    TFile *XConeup_file   = new TFile(dir+year+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"); // CHANGE_NORMAL
    TFile *XConedown_file = new TFile(dir+year+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"); // CHANGE_NORMAL
    TH1F  *XConeup        = (TH1F*)XConeup_file->Get(w_mass);
    TH1F  *XConedown      = (TH1F*)XConedown_file->Get(w_mass);
    XConeup->Add(bkg, 1);
    XConedown->Add(bkg, 1);
    TH1F* XConeup_norm                  = normalize(XConeup);
    TH1F* XConedown_norm                = normalize(XConedown);
    vector<double> XConeup_norm_err     = normalize_error(XConeup);
    vector<double> XConedown_norm_err   = normalize_error(XConedown);

    /*
    . ██  ██████  █████  ███████
    .███ ██      ██   ██ ██
    . ██ ███████  ██████ ███████
    . ██ ██    ██     ██      ██
    . ██  ██████  █████  ███████
    */
    if(debug) cout << "1695" << endl;

    TString ttbar_path_1695 = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root"; // CHANGE_MTOP
    TFile *ttbar_file_1695  = new TFile(dir+year+"/muon/"+ttbar_path_1695);
    TH1F  *ttbar_1695       = (TH1F*)ttbar_file_1695->Get(w_mass);
    ttbar_1695->Add(bkg, 1);
    TH1F          *ttbar_1695_norm     = normalize(ttbar_1695);
    vector<double> ttbar_1695_norm_err = normalize_error(ttbar_1695);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC 1695" << endl;

    TFile *JECup_1695_file   = new TFile(dir+year+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root");
    TFile *JECdown_1695_file = new TFile(dir+year+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root");
    TH1F  *JECup_1695        = (TH1F*)JECup_1695_file->Get(w_mass);
    TH1F  *JECdown_1695      = (TH1F*)JECdown_1695_file->Get(w_mass);
    JECup_1695->Add(bkg, 1);
    JECdown_1695->Add(bkg, 1);
    TH1F          *JECup_1695_norm       = normalize(JECup_1695);
    TH1F          *JECdown_1695_norm     = normalize(JECdown_1695);
    vector<double> JECup_1695_norm_err   = normalize_error(JECup_1695);
    vector<double> JECdown_1695_norm_err = normalize_error(JECdown_1695);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone 1695" << endl;

    TFile *XConeup_1695_file    = new TFile(dir+year+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root");
    TFile *XConeudown_1695_file = new TFile(dir+year+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root");
    TH1F  *XConeup_1695         = (TH1F*)XConeup_1695_file->Get(w_mass);
    TH1F  *XConedown_1695       = (TH1F*)XConeudown_1695_file->Get(w_mass);
    XConeup_1695->Add(bkg, 1);
    XConedown_1695->Add(bkg, 1);
    TH1F* XConeup_1695_norm                  = normalize(XConeup_1695);
    TH1F* XConedown_1695_norm                = normalize(XConedown_1695);
    vector<double> XConeup_1695_norm_err     = normalize_error(XConeup_1695);
    vector<double> XConedown_1695_norm_err   = normalize_error(XConedown_1695);

    /*
    . ██ ███████ ███████ ███████
    .███      ██ ██      ██
    . ██     ██  ███████ ███████
    . ██    ██        ██      ██
    . ██    ██   ███████ ███████
    */
    if(debug) cout << "1755" << endl;

    TString ttbar_path_1755 = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root"; // CHANGE_MTOP
    TFile *ttbar_file_1755  = new TFile(dir+year+"/muon/"+ttbar_path_1755);
    TH1F  *ttbar_1755       = (TH1F*)ttbar_file_1755->Get(w_mass);
    ttbar_1755->Add(bkg, 1);
    TH1F  *ttbar_1755_norm  = normalize(ttbar_1755);
    vector<double> ttbar_1755_norm_err = normalize_error(ttbar_1755);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC 1755" << endl;

    TFile *JECup_1755_file   = new TFile(dir+year+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
    TFile *JECdown_1755_file = new TFile(dir+year+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
    TH1F  *JECup_1755        = (TH1F*)JECup_1755_file->Get(w_mass);
    TH1F  *JECdown_1755      = (TH1F*)JECdown_1755_file->Get(w_mass);
    JECup_1755->Add(bkg, 1);
    JECdown_1755->Add(bkg, 1);
    TH1F* JECup_1755_norm                = normalize(JECup_1755);
    TH1F* JECdown_1755_norm              = normalize(JECdown_1755);
    vector<double> JECup_1755_norm_err   = normalize_error(JECup_1755);
    vector<double> JECdown_1755_norm_err = normalize_error(JECdown_1755);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone 1755" << endl;

    TFile *XConeup_1755_file   = new TFile(dir+year+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
    TFile *XConedown_1755_file = new TFile(dir+year+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
    TH1F  *XConeup_1755        = (TH1F*)XConeup_1755_file->Get(w_mass);
    TH1F  *XConedown_1755      = (TH1F*)XConedown_1755_file->Get(w_mass);
    XConeup_1755->Add(bkg, 1);
    XConedown_1755->Add(bkg, 1);
    TH1F* XConeup_1755_norm                  = normalize(XConeup_1755);
    TH1F* XConedown_1755_norm                = normalize(XConedown_1755);
    vector<double> XConeup_1755_norm_err     = normalize_error(XConeup_1755);
    vector<double> XConedown_1755_norm_err   = normalize_error(XConedown_1755);

    /*
    . █████  ██      ██          ██ ███    ██      ██████  ███    ██ ███████
    .██   ██ ██      ██          ██ ████   ██     ██    ██ ████   ██ ██
    .███████ ██      ██          ██ ██ ██  ██     ██    ██ ██ ██  ██ █████
    .██   ██ ██      ██          ██ ██  ██ ██     ██    ██ ██  ██ ██ ██
    .██   ██ ███████ ███████     ██ ██   ████      ██████  ██   ████ ███████
    */

    vector<TH1F*> hist_all_v = {ttbar, ttbar_1695, ttbar_1755};
    vector<TH1F*> hist_jecup_v = {JECup, JECup_1695, JECup_1755};
    vector<TH1F*> hist_jecdown_v = {JECdown, JECdown_1695, JECdown_1755};
    vector<TH1F*> hist_xcup_v = {XConeup, XConeup_1695, XConeup_1755};
    vector<TH1F*> hist_xcdown_v = {XConedown, XConedown_1695, XConedown_1755};

    TH1F *hist_all     = AddHists(hist_all_v, 1);
    TH1F *hist_jup     = AddHists(hist_jecup_v, 1);
    TH1F *hist_jdown   = AddHists(hist_jecdown_v, 1);
    TH1F *hist_xcup    = AddHists(hist_xcup_v, 1);
    TH1F *hist_xcdown  = AddHists(hist_xcdown_v, 1);

    TH1F *hist_all_norm     = normalize(hist_all);
    TH1F *hist_jup_norm     = normalize(hist_jup);
    TH1F *hist_jdown_norm   = normalize(hist_jdown);
    TH1F *hist_xcup_norm    = normalize(hist_xcup);
    TH1F *hist_xcdown_norm  = normalize(hist_xcdown);

    vector<double> err_all_norm     = normalize_error(hist_all_norm);
    vector<double> err_jup_norm     = normalize_error(hist_jup_norm);
    vector<double> err_jdown_norm   = normalize_error(hist_jdown_norm);
    vector<double> err_xcup_norm    = normalize_error(hist_xcup_norm);
    vector<double> err_xcdown_norm  = normalize_error(hist_xcdown_norm);

    all_bins_combi.push_back({hist_all_norm, hist_jup_norm, hist_jdown_norm, hist_xcup_norm, hist_xcdown_norm});
    all_bins_combi_err.push_back({err_all_norm, err_jup_norm, err_jdown_norm, err_xcup_norm, err_xcdown_norm});

    /*
    ██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██
    ██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
    ██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████
    */

    TGraph2DErrors *bin_fit_nom, *bin_fit_1695, *bin_fit_1755;
    TString str_fit_lin  = "[0] + [1]*x + [2]*y";                     // 3
    TString str_fit_xy2  = "[0] + [1]*x + [2]*y*y";                   // 3
    TString str_fit_x2y  = "[0] + [1]*x*x + [2]*y";                   // 3
    TString str_fit_xyy2 = "[0] + [1]*x + [2]*y + [3]*y*y";           // 4
    TString str_fit_xx2y = "[0] + [1]*x + [2]*x*x + [3]*y";           // 4
    TString str_fit_quad = "[0] + [1]*x*x + [2]*y*y";                 // 3
    TString str_fit_poly = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y"; // 5
    vector<vector<TString>> str_fits   = {{str_fit_lin}, {str_fit_xy2, str_fit_x2y}, {str_fit_xyy2, str_fit_xx2y}, {str_fit_quad}, {str_fit_poly}};
    vector<vector<TString>> name_fits  = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
    vector<vector<int>>     ndfs       = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
    vector<vector<int>>     n_fit_para = {{3}, {3, 3}, {4, 4}, {3}, {5}};
    // vector<TF2*>            all_fits;                                     // Fill later - Declare befor loop
    // vector<TString>         str_all_fits;                                 // Fill later - Declare befor loop
    vector<int>             all_n_para, all_n_para_1695, all_n_para_1755; // Fill later - Declare befor loop
    int n_factors = 5; // ttbar, JECdown, JECup, XConedown, XConeup

    vector<double> factor_x = {0.0, -1.0,  1.0,  0.0,  0.0};
    vector<double> factor_y = {0.0,  0.0,  0.0, -1.0,  1.0};
    vector<double> error_x  = {0.0,  0.0,  0.0,  0.0,  0.0};
    vector<double> error_y  = {0.0,  0.0,  0.0,  0.0,  0.0};

    vector<TH1F*> hists              = {ttbar_norm, JECdown_norm, JECup_norm, XConedown_norm, XConeup_norm};
    vector<vector<double>> hists_err = {ttbar_norm_err, JECdown_norm_err, JECup_norm_err, XConedown_norm_err, XConeup_norm_err};
    vector<vector<double>> hists_cont, hists_cont_err;           // Used in Projection, Filled in Loop

    vector<TH1F*> hists_1695              = {ttbar_1695_norm, JECdown_1695_norm, JECup_1695_norm, XConedown_1695_norm, XConeup_1695_norm};
    vector<vector<double>> hists_1695_err = {ttbar_1695_norm_err, JECdown_1695_norm_err, JECup_1695_norm_err, XConedown_1695_norm_err, XConeup_1695_norm_err};
    vector<vector<double>> hists_1695_cont, hists_1695_cont_err; // Used in Projection, Filled in Loop

    vector<TH1F*> hists_1755              = {ttbar_1755_norm, JECdown_1755_norm, JECup_1755_norm, XConedown_1755_norm, XConeup_1755_norm};
    vector<vector<double>> hists_1755_err = {ttbar_1755_norm_err, JECdown_1755_norm_err, JECup_1755_norm_err, XConedown_1755_norm_err, XConeup_1755_norm_err};
    vector<vector<double>> hists_1755_cont, hists_1755_cont_err; // Used in Projection, Filled in Loop

    vector<TGraphErrors*> bin_graph_jec_nom, bin_graph_jec_1695, bin_graph_jec_1755;
    vector<TGraphErrors*> bin_graph_xc_nom, bin_graph_xc_1695, bin_graph_xc_1755;

    vector<TF1*> bin_function_jec_nom, bin_function_jec_1695, bin_function_jec_1755;
    vector<TF1*> bin_function_xc_nom, bin_function_xc_1695, bin_function_xc_1755;

    /*
    .██████  ██████  ██    ██ ███    ██ ████████
    ██      ██    ██ ██    ██ ████   ██    ██
    ██      ██    ██ ██    ██ ██ ██  ██    ██
    ██      ██    ██ ██    ██ ██  ██ ██    ██
    .██████  ██████   ██████  ██   ████    ██
    */

    // #################################################################################################
    // Number overlap ##################################################################################

    // Used in FIT LOOP to count number of overlaping points in total ----------------------------------
    int overlap_ptbin_n6=0    ; int overlap_ptbin_n7=0    ; int overlap_ptbin_67=0    ; int overlap_ptbin_all=0;
    int overlap_ptbin_sys_n6=0; int overlap_ptbin_sys_n7=0; int overlap_ptbin_sys_67=0; int overlap_ptbin_sys_all=0;
    int total_points_ptbin=0  ; int total_points_ptbin_sys=0;

    for(int bin=0; bin < number_bins; bin++){
      if(!(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())) continue;
      if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      for(int factor=0; factor<5; factor++){ // First only nominal
        // if(factor==1||factor==2) continue;
        // if(factor==3||factor==4) continue;
        // if(factor==1) continue;

        double bin_content_nom  = hists[factor]->GetBinContent(bin+1);
        double bin_error_nom    = hists_err[factor][bin];
        double bin_content_1695 = hists_1695[factor]->GetBinContent(bin+1);
        double bin_error_1695   = hists_1695_err[factor][bin];
        double bin_content_1755 = hists_1755[factor]->GetBinContent(bin+1);
        double bin_error_1755   = hists_1755_err[factor][bin];

        double diff_n6 = fabs(bin_content_nom-bin_content_1695);
        double diff_n7 = fabs(bin_content_nom-bin_content_1755);
        double diff_67 = fabs(bin_content_1695-bin_content_1755);

        // double err_n6 = 1*sqrt(pow(bin_error_nom,2)+pow(bin_error_1695,2));
        // double err_n7 = 1*sqrt(pow(bin_error_nom,2)+pow(bin_error_1755,2));
        // double err_67 = 1*sqrt(pow(bin_error_1695,2)+pow(bin_error_1695,2));

        // double err_n6 = 2*sqrt(pow(bin_error_nom,2)+pow(bin_error_1695,2));
        // double err_n7 = 2*sqrt(pow(bin_error_nom,2)+pow(bin_error_1755,2));
        // double err_67 = 2*sqrt(pow(bin_error_1695,2)+pow(bin_error_1695,2));

        double err_n6 = 3*sqrt(pow(bin_error_nom,2)+pow(bin_error_1695,2));
        double err_n7 = 3*sqrt(pow(bin_error_nom,2)+pow(bin_error_1755,2));
        double err_67 = 3*sqrt(pow(bin_error_1695,2)+pow(bin_error_1695,2));

        // #################################################################################################
        // Start ###########################################################################################
        if((diff_n6-err_n6) < 0){
          overlap_ptbin_n6+=1;
          overlap_n6+=1;
        }
        if((diff_n7-err_n7) < 0){
          overlap_ptbin_n7+=1;
          overlap_n7+=1;
        }
        if((diff_67-err_67) < 0){
          overlap_ptbin_67+=1;
          overlap_67+=1;
        }
        total_points_ptbin+=1; // only looking at two smaples at a time
        total_points+=1;
      }
    }

    double ratio_n6 = (double)overlap_ptbin_n6/(double)total_points_ptbin;
    double ratio_n7 = (double)overlap_ptbin_n7/(double)total_points_ptbin;
    double ratio_67 = (double)overlap_ptbin_67/(double)total_points_ptbin;

    cout << "-------------------------- "+addition << endl;
    cout << "Total points: " << total_points_ptbin << " | #Bins:      " << peak_bins_v.size() << endl;
    cout << "Overlap n6:   " <<   overlap_ptbin_n6 << " | Percentage: " <<       ratio_n6*100 << endl;
    cout << "Overlap n7:   " <<   overlap_ptbin_n7 << " | Percentage: " <<       ratio_n7*100 << endl;
    cout << "Overlap 67:   " <<   overlap_ptbin_67 << " | Percentage: " <<       ratio_67*100 << endl;


    // // RESET -----------------------------------------------------------------------------------------------
    // overlap_ptbin_n6=0    ; overlap_ptbin_n7=0    ; overlap_ptbin_67=0    ; overlap_ptbin_all=0;
    // overlap_ptbin_sys_n6=0; overlap_ptbin_sys_n7=0; overlap_ptbin_sys_67=0; overlap_ptbin_sys_all=0;
    // total_points_ptbin=0  ; total_points_ptbin_sys=0;

    //   /*
    //   .███████ ████████  █████  ██████  ████████     ███    ███  █████  ███████ ███████
    //   .██         ██    ██   ██ ██   ██    ██        ████  ████ ██   ██ ██      ██
    //   .███████    ██    ███████ ██████     ██        ██ ████ ██ ███████ ███████ ███████
    //   .     ██    ██    ██   ██ ██   ██    ██        ██  ██  ██ ██   ██      ██      ██
    //   .███████    ██    ██   ██ ██   ██    ██        ██      ██ ██   ██ ███████ ███████
    //   */
    //
    //   for(unsigned int mass=0; mass<3; mass++){ // 0-nominal; 1-1695; 2-1755
    //     vector<TF2*>    all_fits;      // Fill later - Declare befor loop
    //     vector<TString> str_all_fits;  // Fill later - Declare befor loop
    //     /*
    //     ███████ ██ ████████     ██       ██████   ██████  ██████
    //     ██      ██    ██        ██      ██    ██ ██    ██ ██   ██
    //     █████   ██    ██        ██      ██    ██ ██    ██ ██████
    //     ██      ██    ██        ██      ██    ██ ██    ██ ██
    //     ██      ██    ██        ███████  ██████   ██████  ██
    //     */
    //     TLegend *leg;
    //     print_seperater();
    //     if(debug) cout << "Start: looping over bins - mass " << mass << endl;
    //     for(int bin=0; bin < number_bins; bin++){
    //       if(!(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())) continue;
    //       if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
    //       if(debug) cout << mass_bin[bin] << endl;
    //
    //       // #################################################################################################
    //       // Fill Vectors ####################################################################################
    //       if(debug) cout << "Get Bin Content" << endl;
    //
    //       int bin_h = bin+1;
    //       vector<double> bin_content, bin_content_err;
    //
    //       for(int factor=0; factor<n_factors; factor++){
    //         if(mass==0){
    //           bin_content.push_back(hists[factor]->GetBinContent(bin+1));
    //           bin_content_err.push_back(hists_err[factor][bin]);
    //         }
    //         else if(mass==1){
    //           bin_content.push_back(hists_1695[factor]->GetBinContent(bin+1));
    //           bin_content_err.push_back(hists_1695_err[factor][bin]);
    //         }
    //         else if(mass==2){
    //           bin_content.push_back(hists_1755[factor]->GetBinContent(bin+1));
    //           bin_content_err.push_back(hists_1755_err[factor][bin]);
    //         }
    //       }
    //
    //       if(mass==0){
    //         hists_cont.push_back(bin_content);
    //         hists_cont_err.push_back(bin_content_err);
    //       }
    //       else if(mass==1){
    //         hists_1695_cont.push_back(bin_content);
    //         hists_1695_cont_err.push_back(bin_content_err);
    //       }
    //       else if(mass==2){
    //         hists_1755_cont.push_back(bin_content);
    //         hists_1755_cont_err.push_back(bin_content_err);
    //       }
    //
    //       // #################################################################################################
    //       // Set TGraph ######################################################################################
    //       if(debug) cout << "TGraphError - Single bin " << bin_h << endl;
    //       TGraph2DErrors *bin_fit = new TGraph2DErrors(5, &factor_x[0], &factor_y[0],&bin_content[0],&error_x[0],&error_y[0],&bin_content_err[0]);
    //       bin_fit->SetTitle(mass_bin[bin]);
    //       bin_fit->SetName(number_bin[bin]);
    //       bin_fit->SetMarkerStyle(kFullCircle);
    //       bin_fit->SetMarkerColor(kBlack);
    //       bin_fit->SetMarkerSize(0.5);
    //       bin_fit->SetLineColor(kBlack);
    //
    //       // #################################################################################################
    //       // set fit function ################################################################################
    //       if(debug) cout << "\nDefine fit function" << endl;
    //       // TF2 *fit_con   = new TF2("fit_con",  "[0]");
    //       // Declared here, because they are reseted each bin
    //       TF2 *fit_lin  = new TF2("fit_lin",  str_fits[0][0]);
    //       TF2 *fit_xy2  = new TF2("fit_xy2",  str_fits[1][0]);
    //       TF2 *fit_x2y  = new TF2("fit_x2y",  str_fits[1][1]);
    //       TF2 *fit_xyy2 = new TF2("fit_xyy2", str_fits[2][0]);
    //       TF2 *fit_xx2y = new TF2("fit_xx2y", str_fits[2][1]);
    //       TF2 *fit_quad = new TF2("fit_quad", str_fits[3][0]);
    //       TF2 *fit_poly = new TF2("fit_poly", str_fits[4][0]);
    //       vector<vector<TF2*>> fits = {{fit_lin}, {fit_xy2, fit_x2y}, {fit_xyy2, fit_xx2y}, {fit_quad}, {fit_poly}};
    //
    //       // #################################################################################################
    //       // Fit Bin #########################################################################################
    //       if(debug) cout << "Start: Fitting" << endl;
    //
    //       /*
    //       ███████ ██ ████████
    //       ██      ██    ██
    //       █████   ██    ██
    //       ██      ██    ██
    //       ██      ██    ██
    //       */
    //
    //       double chi2_fit, chi2_fit1, chi2_fit2;
    //       double prob, prob1, prob2;
    //       TF2 *used_fit;
    //       int order_index = 0; int index_fit = 0;
    //       for(unsigned int order=0; order<fits.size(); order++){
    //         if(debug) cout << "\nfit a "+name_fits[order][0]+" function ------------\n";
    //         bool twoFits = fits[order].size()==2;
    //
    //         bin_fit->Fit(fits[order][0], "Q");
    //         if(twoFits) bin_fit->Fit(fits[order][1], "Q");
    //
    //         // -----------------------------------------------------------------------
    //         if(debug) cout << "Chi2" << endl;
    //         chi2_fit1 = fits[order][0]->GetChisquare();
    //         if(twoFits) chi2_fit2 = fits[order][1]->GetChisquare();
    //
    //         // -----------------------------------------------------------------------
    //         if(debug) cout << "prob" << endl;
    //
    //         prob1 = TMath::Prob(chi2_fit1,  ndfs[order][0]);
    //         if(twoFits) prob2 = TMath::Prob(chi2_fit2, ndfs[order][1]);
    //
    //         prob = prob1;
    //         chi2_fit = chi2_fit1;
    //         index_fit = 0;
    //         if(twoFits && (prob1 < prob2)){
    //           prob = prob2;
    //           chi2_fit = chi2_fit2;
    //           index_fit = 1;
    //         }
    //         if(debug) cout << "ndfs: " << ndfs[order][index_fit] << " | chi2: " << chi2_fit <<" | prob: " << prob << endl;
    //
    //         // -----------------------------------------------------------------------
    //         if(prob>0.05 || ndfs[order][0]==0){ // prob > 0.05 is common in data science
    //           cout << "Bin " << number_bin[bin] << ": A " << GREEN << name_fits[order][index_fit] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
    //           used_fit = new TF2(*fits[order][index_fit]);        // used in Loop
    //           all_fits.push_back(fits[order][index_fit]);         // TF2*    - used in Projection
    //           str_all_fits.push_back(str_fits[order][index_fit]); // TString - used in Projection
    //           all_n_para.push_back(n_fit_para[order][index_fit]); // int     - used in Projection
    //           break;
    //         }
    //       }
    //     }
    //
    //     /*
    //     ███████ ██ ████████     ███████ ███    ██ ██████
    //     ██      ██    ██        ██      ████   ██ ██   ██
    //     █████   ██    ██        █████   ██ ██  ██ ██   ██
    //     ██      ██    ██        ██      ██  ██ ██ ██   ██
    //     ██      ██    ██        ███████ ██   ████ ██████
    //     */
    //
    //     /*
    //     ██████  ██████   ██████       ██ ███████  ██████ ████████ ██  ██████  ███    ██
    //     ██   ██ ██   ██ ██    ██      ██ ██      ██         ██    ██ ██    ██ ████   ██
    //     ██████  ██████  ██    ██      ██ █████   ██         ██    ██ ██    ██ ██ ██  ██
    //     ██      ██   ██ ██    ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██
    //     ██      ██   ██  ██████   █████  ███████  ██████    ██    ██  ██████  ██   ████
    //     */
    //     cout << "\n####################################################################\n";
    //     cout << "######################### Projection ###############################\n";
    //     cout << "####################################################################\n";
    //
    //     vector<TString> factor_name = {"JEC", "XC"};
    //     // #################################################################################################
    //     // Which correction ################################################################################
    //     for(int factor=0; factor<2; factor++){
    //       if(debug) cout << "------------------------------------------------------- " << factor_name[factor] << endl;
    //       vector<double> oneD_factor     = {-1,0,1};
    //       vector<double> oneD_factor_err = { 0,0,0};
    //
    //       // #################################################################################################
    //       // Bins ############################################################################################
    //       int empty_bins = 0; // to skip empty bins
    //       for(int bin=0; bin<number_bins; bin++){
    //         if(!(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())){ // Search_peackmethod
    //           empty_bins++;
    //           continue;
    //         }
    //
    //         int used_bin = bin-empty_bins;
    //         vector<double> oneD_events, oneD_events_err;
    //         vector<vector<double>> hists_cont_gen, hists_cont_gen_err;
    //         if(mass==0){
    //           hists_cont_gen    =hists_cont;
    //           hists_cont_gen_err=hists_cont_err;
    //         } else if(mass==1){
    //           hists_cont_gen    =hists_1695_cont;
    //           hists_cont_gen_err=hists_1695_cont_err;
    //         } else if(mass==2){
    //           hists_cont_gen    =hists_1755_cont;
    //           hists_cont_gen_err=hists_1755_cont_err;
    //         }
    //         if(factor==0){
    //           oneD_events     = {hists_cont_gen[used_bin][1], hists_cont_gen[used_bin][0], hists_cont_gen[used_bin][2]};
    //           oneD_events_err = {hists_cont_gen_err[used_bin][1], hists_cont_gen_err[used_bin][0], hists_cont_gen_err[used_bin][2]};
    //         } else if(factor==1){
    //           oneD_events     = {hists_cont_gen[used_bin][3], hists_cont_gen[used_bin][0], hists_cont_gen[used_bin][4]};
    //           oneD_events_err = {hists_cont_gen_err[used_bin][3], hists_cont_gen_err[used_bin][0], hists_cont_gen_err[used_bin][4]};
    //         }
    //
    //         // XCone -----------------------------------------------------------------------------------------
    //         TString str_fit = str_all_fits[used_bin];
    //         if     (factor==0) str_fit.ReplaceAll("y", "0");
    //         else if(factor==1){
    //           str_fit.ReplaceAll("x", "0");
    //           str_fit.ReplaceAll("y", "x"); // 1D Graph only takes x as argument
    //         }
    //
    //         TF1 *fit = new TF1("fit", str_fit, -2, 2);
    //         TGraphErrors *graph = new TGraphErrors(3, &oneD_factor[0], &oneD_events[0], &oneD_factor_err[0], &oneD_events_err[0]);
    //         for(int ipar=0; ipar<all_n_para[used_bin]; ipar++) fit->SetParameter(ipar, all_fits[used_bin]->GetParameter(ipar));
    //
    //         if     (mass==0){
    //           if(factor==0) bin_graph_jec_nom.push_back(graph);
    //           if(factor==1) bin_graph_xc_nom.push_back(graph);
    //
    //           if(factor==0) bin_function_jec_nom.push_back(fit);
    //           if(factor==1) bin_function_xc_nom.push_back(fit);
    //         }
    //         else if(mass==1) {
    //           if(factor==0) bin_graph_jec_1695.push_back(graph);
    //           if(factor==1) bin_graph_xc_1695.push_back(graph);
    //
    //           if(factor==0) bin_function_jec_1695.push_back(fit);
    //           if(factor==1) bin_function_xc_1695.push_back(fit);
    //         }
    //         else if(mass==2) {
    //           if(factor==0) bin_graph_jec_1755.push_back(graph);
    //           if(factor==1) bin_graph_xc_1755.push_back(graph);
    //
    //           if(factor==0) bin_function_jec_1755.push_back(fit);
    //           if(factor==1) bin_function_xc_1755.push_back(fit);
    //         }
    //       }
    //     }
    //   }
    //
    //   /*
    //   ███    ███  █████  ███████ ███████     ███████ ███    ██ ██████
    //   ████  ████ ██   ██ ██      ██          ██      ████   ██ ██   ██
    //   ██ ████ ██ ███████ ███████ ███████     █████   ██ ██  ██ ██   ██
    //   ██  ██  ██ ██   ██      ██      ██     ██      ██  ██ ██ ██   ██
    //   ██      ██ ██   ██ ███████ ███████     ███████ ██   ████ ██████
    //   */
    //   cout << peak_bins_v.size() << endl;
    //   cout << bin_function_jec_nom.size()  << "     " << bin_graph_jec_nom.size() << endl;
    //   cout << bin_function_jec_1695.size() << "     " << bin_graph_jec_1695.size() << endl;
    //   cout << bin_function_jec_1755.size() << "     " << bin_graph_jec_1755.size() << endl;
    //
    //   int correct_bin = peak_bins_v[0];
    //
    //   for(unsigned int factor=0; factor<2; factor++){
    //     for(unsigned int i=0; i<peak_bins_v.size(); i++){
    //       TGraphErrors* bin_graph_nom, *bin_graph_1695, *bin_graph_1755;
    //       TF1* bin_function_nom, *bin_function_1695, *bin_function_1755;
    //       if(factor==0){
    //         bin_function_nom  =bin_function_jec_nom[i];
    //         bin_function_1695 =bin_function_jec_1695[i];
    //         bin_function_1755 =bin_function_jec_1755[i];
    //
    //         bin_graph_nom     =bin_graph_jec_nom[i];
    //         bin_graph_1695    =bin_graph_jec_1695[i];
    //         bin_graph_1755    =bin_graph_jec_1755[i];
    //       }
    //       else{
    //         bin_function_nom  =bin_function_xc_nom[i];
    //         bin_function_1695 =bin_function_xc_1695[i];
    //         bin_function_1755 =bin_function_xc_1755[i];
    //
    //         bin_graph_nom     =bin_graph_xc_nom[i];
    //         bin_graph_1695    =bin_graph_xc_1695[i];
    //         bin_graph_1755    =bin_graph_xc_1755[i];
    //       }
    //
    //
    //       bin_graph_nom->SetTitle(mass_bin[i+correct_bin]);
    //       if(factor==0) bin_graph_nom->GetXaxis()->SetTitle("JEC factor");
    //       else          bin_graph_nom->GetXaxis()->SetTitle("Additional XCone correction factor");
    //       bin_graph_nom->SetMarkerColor(kBlack);
    //       bin_graph_nom->SetMarkerStyle(kFullCircle);
    //       bin_graph_nom->SetMarkerSize(0.5);
    //       bin_graph_nom->SetLineColor(kBlack);
    //       bin_graph_nom->GetXaxis()->SetLimits(-1.1,1.1);
    //       bin_function_nom->SetLineColor(kBlack);
    //
    //       double xn, yn;
    //       bin_graph_nom->GetPoint(1, xn, yn);
    //       bin_graph_nom->GetYaxis()->SetRangeUser(yn/2, yn*1.5);
    //
    //       bin_graph_1695->SetMarkerColor(kRed);
    //       bin_graph_1695->SetMarkerStyle(kFullCircle);
    //       bin_graph_1695->SetMarkerSize(0.5);
    //       bin_graph_1695->SetLineColor(kRed);
    //       bin_function_1695->SetLineColor(kRed);
    //
    //       bin_graph_1755->SetMarkerColor(kBlue);
    //       bin_graph_1755->SetMarkerStyle(kFullCircle);
    //       bin_graph_1755->SetMarkerSize(0.5);
    //       bin_graph_1755->SetLineColor(kBlue);
    //       bin_function_1755->SetLineColor(kBlue);
    //
    //       TCanvas *E2;
    //       TLegend *leg;
    //       if(factor==0) E2 = new TCanvas(number_bin[i+correct_bin]+"E2_JEC_"+to_string(ptbin), "E2", 600, 600);
    //       else          E2 = new TCanvas(number_bin[i+correct_bin]+"E2_XC_"+to_string(ptbin), "E2", 600, 600);
    //       bin_graph_nom->Draw("AP");
    //       bin_graph_1695->Draw("SAME P");
    //       bin_graph_1755->Draw("SAME P");
    //
    //       bin_function_nom->Draw("SAME");
    //       bin_function_1695->Draw("SAME");
    //       bin_function_1755->Draw("SAME");
    //
    //       leg = new TLegend(0.20,0.15,0.40,0.3);
    //       leg->SetBorderSize(0);
    //       leg->AddEntry(bin_graph_nom ,"1725","pl");
    //       leg->AddEntry(bin_graph_1695,"1695","pl");
    //       leg->AddEntry(bin_graph_1755,"1755","pl");
    //       leg->Draw();
    //
    //       if(factor==0) E2->SaveAs(save_path_general+"/"+addition+"/Bin"+number_bin[i+correct_bin]+"_JEC.pdf");
    //       else          E2->SaveAs(save_path_general+"/"+addition+"/Bin"+number_bin[i+correct_bin]+"_XC.pdf");
    //       leg->Clear();
    //     }
    //   }
    //
    //   /*
    //   ██       █████  ████████ ███████ ██   ██
    //   ██      ██   ██    ██    ██       ██ ██
    //   ██      ███████    ██    █████     ███
    //   ██      ██   ██    ██    ██       ██ ██
    //   ███████ ██   ██    ██    ███████ ██   ██
    //   */
    //
    //   cout << "\n####################################################################\n";
    //   cout << "############################ Latex #################################\n";
    //   cout << "####################################################################\n";
    //   // suppress output with &> /dev/null
    //
    //   // system("python python/projection_comparison_all.py");
    //
  }
  /*
  ███████ ███    ██ ██████      ██████  ████████ ██████  ██ ███    ██
  ██      ████   ██ ██   ██     ██   ██    ██    ██   ██ ██ ████   ██
  █████   ██ ██  ██ ██   ██     ██████     ██    ██████  ██ ██ ██  ██
  ██      ██  ██ ██ ██   ██     ██         ██    ██   ██ ██ ██  ██ ██
  ███████ ██   ████ ██████      ██         ██    ██████  ██ ██   ████
  */

  double ratio_all_n6 = (double)overlap_n6/(double)total_points;
  double ratio_all_n7 = (double)overlap_n7/(double)total_points;
  double ratio_all_67 = (double)overlap_67/(double)total_points;

  cout << "-------------------------- all" << endl;
  cout << "Total points: " << total_points << endl;
  cout << "Overlap n6:   " <<   overlap_n6 << " | Percentage: " <<       ratio_all_n6*100 << endl;
  cout << "Overlap n7:   " <<   overlap_n7 << " | Percentage: " <<       ratio_all_n7*100 << endl;
  cout << "Overlap 67:   " <<   overlap_67 << " | Percentage: " <<       ratio_all_67*100 << endl;

  // #################################################################################################
  // #################################################################################################
  // #################################################################################################
  // #################################################################################################

  // /*
  // ███████ ██ ████████      █████  ██      ██
  // ██      ██    ██        ██   ██ ██      ██
  // █████   ██    ██        ███████ ██      ██
  // ██      ██    ██        ██   ██ ██      ██
  // ██      ██    ██        ██   ██ ███████ ███████
  // */
  // cout << "Bluuuuuuuub " << all_bins_combi.size() << all_bins_combi_err.size() << endl;
  // save_path_general = creat_folder_and_path(save_path_general, "add_all");
  // for(unsigned int ihist=0; ihist<all_bins_combi.size(); ihist++){
  //
  //   if(ihist==0) addition="_hh";
  //   if(ihist==1) addition="_hl";
  //   if(ihist==2) addition="_lh";
  //   if(ihist==3) addition="_ll";
  //   creat_folder(save_path_general, addition);
  //
  //   TH1F *ttbar_norm     = all_bins_combi[ihist][0];
  //   TH1F *JECup_norm     = all_bins_combi[ihist][1];
  //   TH1F *JECdown_norm   = all_bins_combi[ihist][2];
  //   TH1F *XConeup_norm   = all_bins_combi[ihist][3];
  //   TH1F *XConedown_norm = all_bins_combi[ihist][4];
  //
  //   vector<double> ttbar_norm_err     = all_bins_combi_err[ihist][0];
  //   vector<double> JECup_norm_err     = all_bins_combi_err[ihist][1];
  //   vector<double> JECdown_norm_err   = all_bins_combi_err[ihist][2];
  //   vector<double> XConeup_norm_err   = all_bins_combi_err[ihist][3];
  //   vector<double> XConedown_norm_err = all_bins_combi_err[ihist][4];
  //
  //   /*
  //   ██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██
  //   ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██
  //   ██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██
  //   ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
  //   ██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████
  //   */
  //
  //   TGraph2DErrors *bin_fit_nom, *bin_fit_1695, *bin_fit_1755;
  //   TString str_fit_lin  = "[0] + [1]*x + [2]*y";                     // 3
  //   TString str_fit_xy2  = "[0] + [1]*x + [2]*y*y";                   // 3
  //   TString str_fit_x2y  = "[0] + [1]*x*x + [2]*y";                   // 3
  //   TString str_fit_xyy2 = "[0] + [1]*x + [2]*y + [3]*y*y";           // 4
  //   TString str_fit_xx2y = "[0] + [1]*x + [2]*x*x + [3]*y";           // 4
  //   TString str_fit_quad = "[0] + [1]*x*x + [2]*y*y";                 // 3
  //   TString str_fit_poly = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y"; // 5
  //   vector<vector<TString>> str_fits   = {{str_fit_lin}, {str_fit_xy2, str_fit_x2y}, {str_fit_xyy2, str_fit_xx2y}, {str_fit_quad}, {str_fit_poly}};
  //   vector<vector<TString>> name_fits  = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
  //   vector<vector<int>>     ndfs       = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
  //   vector<vector<int>>     n_fit_para = {{3}, {3, 3}, {4, 4}, {3}, {5}};
  //   vector<TF2*>            all_fits;                                     // Fill later - Declare befor loop
  //   vector<TString>         str_all_fits;                                 // Fill later - Declare befor loop
  //   vector<int>             all_n_para, all_n_para_1695, all_n_para_1755; // Fill later - Declare befor loop
  //   int n_factors = 5; // ttbar, JECdown, JECup, XConedown, XConeup
  //
  //   vector<double> factor_x = {0.0, -1.0,  1.0,  0.0,  0.0};
  //   vector<double> factor_y = {0.0,  0.0,  0.0, -1.0,  1.0};
  //   vector<double> error_x  = {0.0,  0.0,  0.0,  0.0,  0.0};
  //   vector<double> error_y  = {0.0,  0.0,  0.0,  0.0,  0.0};
  //
  //   vector<TH1F*> hists              = {ttbar_norm, JECdown_norm, JECup_norm, XConedown_norm, XConeup_norm};
  //   vector<vector<double>> hists_err = {ttbar_norm_err, JECdown_norm_err, JECup_norm_err, XConedown_norm_err, XConeup_norm_err};
  //   vector<vector<double>> hists_cont, hists_cont_err;           // Used in Projection, Filled in Loop
  //
  //   /*
  //   ███████ ██ ████████     ██       ██████   ██████  ██████
  //   ██      ██    ██        ██      ██    ██ ██    ██ ██   ██
  //   █████   ██    ██        ██      ██    ██ ██    ██ ██████
  //   ██      ██    ██        ██      ██    ██ ██    ██ ██
  //   ██      ██    ██        ███████  ██████   ██████  ██
  //   */
  //
  //   for(int bin=0; bin < number_bins; bin++){
  //     if(!(find(all_peak_bins[ihist].begin(), all_peak_bins[ihist].end(), bin+1) != all_peak_bins[ihist].end())) continue;
  //     if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
  //     if(debug) cout << mass_bin[bin] << endl;
  //
  //     // #################################################################################################
  //     // Fill Vectors ####################################################################################
  //     if(debug) cout << "Get Bin Content" << endl;
  //
  //     int bin_h = bin+1;
  //     vector<double> bin_content, bin_content_err;
  //
  //     for(int factor=0; factor<n_factors; factor++){
  //       bin_content.push_back(hists[factor]->GetBinContent(bin+1));
  //       bin_content_err.push_back(hists_err[factor][bin]);
  //     }
  //
  //     hists_cont.push_back(bin_content);
  //     hists_cont_err.push_back(bin_content_err);
  //
  //     // #################################################################################################
  //     // Set TGraph ######################################################################################
  //     if(debug) cout << "TGraphError - Single bin " << bin_h << endl;
  //     TGraph2DErrors *bin_fit = new TGraph2DErrors(5, &factor_x[0], &factor_y[0],&bin_content[0],&error_x[0],&error_y[0],&bin_content_err[0]);
  //     bin_fit->SetTitle(mass_bin[bin]);
  //     bin_fit->SetName(number_bin[bin]);
  //     bin_fit->SetMarkerStyle(kFullCircle);
  //     bin_fit->SetMarkerColor(kBlack);
  //     bin_fit->SetMarkerSize(0.5);
  //     bin_fit->SetLineColor(kBlack);
  //
  //     // #################################################################################################
  //     // set fit function ################################################################################
  //     if(debug) cout << "\nDefine fit function" << endl;
  //     // TF2 *fit_con   = new TF2("fit_con",  "[0]");
  //     // Declared here, because they are reseted each bin
  //     TF2 *fit_lin  = new TF2("fit_lin",  str_fits[0][0]);
  //     TF2 *fit_xy2  = new TF2("fit_xy2",  str_fits[1][0]);
  //     TF2 *fit_x2y  = new TF2("fit_x2y",  str_fits[1][1]);
  //     TF2 *fit_xyy2 = new TF2("fit_xyy2", str_fits[2][0]);
  //     TF2 *fit_xx2y = new TF2("fit_xx2y", str_fits[2][1]);
  //     TF2 *fit_quad = new TF2("fit_quad", str_fits[3][0]);
  //     TF2 *fit_poly = new TF2("fit_poly", str_fits[4][0]);
  //     vector<vector<TF2*>> fits = {{fit_lin}, {fit_xy2, fit_x2y}, {fit_xyy2, fit_xx2y}, {fit_quad}, {fit_poly}};
  //
  //     // #################################################################################################
  //     // Fit Bin #########################################################################################
  //     if(debug) cout << "Start: Fitting" << endl;
  //
  //     /*
  //     ███████ ██ ████████
  //     ██      ██    ██
  //     █████   ██    ██
  //     ██      ██    ██
  //     ██      ██    ██
  //     */
  //
  //     double chi2_fit, chi2_fit1, chi2_fit2;
  //     double prob, prob1, prob2;
  //     TF2 *used_fit;
  //     int order_index = 0; int index_fit = 0;
  //     for(unsigned int order=0; order<fits.size(); order++){
  //       if(debug) cout << "\nfit a "+name_fits[order][0]+" function ------------\n";
  //       bool twoFits = fits[order].size()==2;
  //
  //       bin_fit->Fit(fits[order][0], "Q");
  //       if(twoFits) bin_fit->Fit(fits[order][1], "Q");
  //
  //       // -----------------------------------------------------------------------
  //       if(debug) cout << "Chi2" << endl;
  //       chi2_fit1 = fits[order][0]->GetChisquare();
  //       if(twoFits) chi2_fit2 = fits[order][1]->GetChisquare();
  //
  //       // -----------------------------------------------------------------------
  //       if(debug) cout << "prob" << endl;
  //
  //       prob1 = TMath::Prob(chi2_fit1,  ndfs[order][0]);
  //       if(twoFits) prob2 = TMath::Prob(chi2_fit2, ndfs[order][1]);
  //
  //       prob = prob1;
  //       chi2_fit = chi2_fit1;
  //       index_fit = 0;
  //       if(twoFits && (prob1 < prob2)){
  //         prob = prob2;
  //         chi2_fit = chi2_fit2;
  //         index_fit = 1;
  //       }
  //       if(debug) cout << "ndfs: " << ndfs[order][index_fit] << " | chi2: " << chi2_fit <<" | prob: " << prob << endl;
  //
  //       // -----------------------------------------------------------------------
  //       if(prob>0.05 || ndfs[order][0]==0){ // prob > 0.05 is common in data science
  //         cout << "Bin " << number_bin[bin] << ": A " << GREEN << name_fits[order][index_fit] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
  //         used_fit = new TF2(*fits[order][index_fit]);        // used in Loop
  //         all_fits.push_back(fits[order][index_fit]);         // TF2*    - used in Projection
  //         str_all_fits.push_back(str_fits[order][index_fit]); // TString - used in Projection
  //         all_n_para.push_back(n_fit_para[order][index_fit]); // int     - used in Projection
  //         break;
  //       }
  //     }
  //   }
  //   /*
  //   ███████ ██ ████████     ███████ ███    ██ ██████
  //   ██      ██    ██        ██      ████   ██ ██   ██
  //   █████   ██    ██        █████   ██ ██  ██ ██   ██
  //   ██      ██    ██        ██      ██  ██ ██ ██   ██
  //   ██      ██    ██        ███████ ██   ████ ██████
  //   */
  //
  //   /*
  //   ██████  ██████   ██████       ██ ███████  ██████ ████████ ██  ██████  ███    ██
  //   ██   ██ ██   ██ ██    ██      ██ ██      ██         ██    ██ ██    ██ ████   ██
  //   ██████  ██████  ██    ██      ██ █████   ██         ██    ██ ██    ██ ██ ██  ██
  //   ██      ██   ██ ██    ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██
  //   ██      ██   ██  ██████   █████  ███████  ██████    ██    ██  ██████  ██   ████
  //   */
  //   cout << "\n####################################################################\n";
  //   cout << "######################### Projection ###############################\n";
  //   cout << "####################################################################\n";
  //
  //   vector<TString> factor_name = {"JEC", "XC"};
  //   // #################################################################################################
  //   // Which correction ################################################################################
  //   for(int factor=0; factor<2; factor++){
  //     if(debug) cout << "------------------------------------------------------- " << factor_name[factor] << endl;
  //     vector<double> oneD_factor     = {-1,0,1};
  //     vector<double> oneD_factor_err = { 0,0,0};
  //
  //     // #################################################################################################
  //     // Bins ############################################################################################
  //     int empty_bins = 0; // to skip empty bins
  //     for(int bin=0; bin<number_bins; bin++){
  //       if(!(find(all_peak_bins[ihist].begin(), all_peak_bins[ihist].end(), bin+1) != all_peak_bins[ihist].end())){ // Search_peackmethod
  //         empty_bins++;
  //         continue;
  //       }
  //
  //       int used_bin = bin-empty_bins;
  //       vector<double> oneD_events, oneD_events_err;
  //       vector<vector<double>> hists_cont_gen, hists_cont_gen_err;
  //
  //       if(factor==0){
  //         oneD_events     = {hists_cont[used_bin][1], hists_cont[used_bin][0], hists_cont[used_bin][2]};
  //         oneD_events_err = {hists_cont_err[used_bin][1], hists_cont_err[used_bin][0], hists_cont_err[used_bin][2]};
  //       } else if(factor==1){
  //         oneD_events     = {hists_cont[used_bin][3], hists_cont[used_bin][0], hists_cont[used_bin][4]};
  //         oneD_events_err = {hists_cont_err[used_bin][3], hists_cont_err[used_bin][0], hists_cont_err[used_bin][4]};
  //       }
  //
  //       // XCone -----------------------------------------------------------------------------------------
  //       TString str_fit = str_all_fits[used_bin];
  //       if     (factor==0) str_fit.ReplaceAll("y", "0");
  //       else if(factor==1){
  //         str_fit.ReplaceAll("x", "0");
  //         str_fit.ReplaceAll("y", "x"); // 1D Graph only takes x as argument
  //       }
  //
  //       TF1 *fit = new TF1("fit", str_fit, -2, 2);
  //       TGraphErrors *graph = new TGraphErrors(3, &oneD_factor[0], &oneD_events[0], &oneD_factor_err[0], &oneD_events_err[0]);
  //       for(int ipar=0; ipar<all_n_para[used_bin]; ipar++) fit->SetParameter(ipar, all_fits[used_bin]->GetParameter(ipar));
  //
  //       // #################################################################################################
  //       // Plot ############################################################################################
  //       graph->SetTitle(mass_bin[bin]);
  //       TCanvas *E2 = new TCanvas(number_bin[bin]+"E2_"+factor_name[factor]+to_string(ihist), "E2", 600, 600);
  //       if(factor==0) graph->GetXaxis()->SetTitle("JEC");
  //       else          graph->GetXaxis()->SetTitle("XCone");
  //       graph->Draw("AP");
  //       fit->Draw("SAME");
  //       E2->SaveAs(save_path_general+"/"+addition+"/Bin"+number_bin[bin]+"_"+factor_name[factor]+addition+".pdf");
  //     }
  //   }
  // }

}
