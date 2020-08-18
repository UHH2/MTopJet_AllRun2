#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  /*
  Explanation:
  1) This script loops many times over bins. The bin number starts at 1 (e.g. for the function GetBinContent)
  On the other hand, the bin content or titles are stored in a vector. The vector starts at 0.
  Mostly, the loops start at bin=0. Therefor, using the function GetBinContent the bin needs to be bin+1
  */

  // Changes of how many pt bins to use - Search for CHANGE_PT
  // #################################################################################################
  // Declare different variables used in the code ####################################################
  TString chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll;
  vector<double> chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll;
  print_seperater();

  bool debug        = false;
  TString reconst   = "btag";    // match btag_cut btag_sel compare min_mass
  TString MC_uncert = "central"; // central avg nominal
  cout.precision(6);

  if(argc != 5){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number> <masspeak> <useOnly_lin>\n";
    cout << "rebin_number:        hist->Rebin(rebin_number)\n";
    cout << "Selection:           match, btag_cut, btag_sel, compare, min_mass\n";
    cout << "Use Masspeak:        0(false), 1(true)\n";
    return 0;
  }

  // #################################################################################################
  // Only one fit for all bins #######################################################################
  if(debug) cout << "String into bool" << endl;

  bool usePeak_in  = stob(argv[3]);
  bool useOnly_lin = stob(argv[4]);
  bool into_latex  = true;






  // #################################################################################################
  // cout setting ####################################################################################
  cout << "Peak:    " << usePeak_in  << endl;
  cout << "Only:    " << useOnly_lin << endl;
  cout << "Latex:   " << into_latex  << endl;
  cout << "debug:   " << debug       << endl;

  // #################################################################################################
  // Print Options ###################################################################################
  Int_t oldLevel = gErrorIgnoreLevel;
  // Set by: gErrorIgnoreLevel = ...

  // Rebin -------------------------------------------------------------------------------------------
  /*
  Right now every Data bin will be excluded, if the BinContent is zero. This leads to a problem,
  since in the first bins (at a binning of 1 GeV) the bin Content is (e.g.) 000001000.
  */
  int bin_width = atoi(argv[2]); // atoi() c-string representing integer into int
  int number_bins = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;

  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0) is16   = true;
  else if(strcmp(year, "2017")==0) is17   = true;
  else if(strcmp(year, "2018")==0) is18   = true;
  else if(strcmp(year,  "all")==0) isAll  = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or all)");

  vector<TString> all_years;
  if(is16 || is17 || is18) all_years = {year};
  else if(isAll) all_years = {"2016", "2017", "2018"};
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018, all or 1718)");

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

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/pt_bins"; // CHANGE_PT
  // TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS"; // CHANGE_PT

  // #################################################################################################
  // creat subdirectories ############################################################################

  save_path_general = creat_folder_and_path(save_path_general, year);
  save_path_general = creat_folder_and_path(save_path_general, reconst);
  save_path_general = creat_folder_and_path(save_path_general, "rebin"+str_number_bins);
  if(usePeak_in&&isAll) save_path_general = creat_folder_and_path(save_path_general, "masspeak");
  creat_folder(save_path_general, "single_bins");
  creat_folder(save_path_general, "projection");

  TString addition="";
  for(int ptbin=0; ptbin<4; ptbin++){ // CHANGE_PT
    TString ptbin_str = to_string(ptbin);
    if(ptbin==0) addition="_hh"; // CHANGE_PT
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    int same_fit;
    if(useOnly_lin){addition += "_lin"; same_fit=0;}
    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */

    vector<TH1F*> allyears_bkg,  allyears_data,  allyears_ttbar;
    vector<TH1F*> allyears_JECup, allyears_JECdown, allyears_XCup, allyears_XCdown;

    cout << '\n'; // CHANGE_PT
    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
    // if(reconst=="min_mass" ) w_mass = hist_class+"wmass_min";                      // minimal mass from subjet combination
    // if(reconst=="btag")      w_mass = hist_class+"wmass_match";                    // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
    // if(reconst=="btag_cut")  w_mass = hist_class+"wmass_btagcut";                  // same mass as above mit with btag_high > 0.7
    // if(reconst=="btag_sel")  w_mass = hist_class+"wmass_btagcut_one_btag_subjet";  // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
    // //                                                                                Selection: btag_high>0.7; ONE subjet with dr(ak4, subjet)<0.4; ONE high btag
    // if(reconst=="compare")   w_mass = hist_class+"wmass_compare";                  // mass closest to Wmass from subjet combination
    // if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptbin_low";
    // if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptbin_high";

    if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(reconst=="btag"&&ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(reconst=="btag"&&ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";


    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      vector<TFile*> file_bkg_v;
      vector<TH1F*> hists_bkg_v;
      vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
      for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+all_years[iyear]+"/muon/"+path_bkg_v[i]));
      for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
      TH1F *bkg = AddHists(hists_bkg_v, 1);
      allyears_bkg.push_back(bkg);
    }

    TH1F *bkg = AddHists(allyears_bkg, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *data_file      = new TFile(dir+all_years[iyear]+"/muon/"+data_path);
      TH1F  *data           = (TH1F*)data_file->Get(w_mass);
      allyears_data.push_back(data);
    }

    TH1F *data            = AddHists(allyears_data, 1);
    TH1F *data_rebin      = rebin(data, bin_width);
    TH1F *data_norm       = normalize(data);
    TH1F *data_rebin_norm = normalize(data_rebin);

    vector<double> data_rebin_norm_err = normalize_error(data_rebin);
    cout<< data_rebin_norm->GetNbinsX()<<"   "<< data_rebin_norm_err.size()<< endl;
    for(int ipar=0; ipar<data_rebin_norm->GetNbinsX(); ipar++) data_rebin_norm->SetBinError(ipar+1, data_rebin_norm_err[ipar]);

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *ttbar_file  = new TFile(dir+all_years[iyear]+"/muon/"+ttbar_path);
      TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);
      allyears_ttbar.push_back(ttbar);
    }

    TH1F  *ttbar           = AddHists(allyears_ttbar, 1);
    ttbar->Add(bkg, 1);
    TH1F* ttbar_rebin      = rebin(ttbar, bin_width);
    TH1F* ttbar_norm       = normalize(ttbar);
    TH1F* ttbar_rebin_norm = normalize(ttbar_rebin);

    vector<double> ttbar_rebin_norm_err = normalize_error(ttbar_rebin);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *JECup_file   = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
      TFile *JECdown_file = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");
      TH1F  *JECup        = (TH1F*)JECup_file->Get(w_mass);
      TH1F  *JECdown      = (TH1F*)JECdown_file->Get(w_mass);
      allyears_JECup.push_back(JECup);
      allyears_JECdown.push_back(JECdown);
    }

    TH1F *JECup   = AddHists(allyears_JECup, 1);
    TH1F *JECdown = AddHists(allyears_JECdown, 1);

    // Add hists from JEC and background ---------------------------------------------------------------
    if(debug) cout << "JEC+Background" << endl;

    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);

    TH1F* JECup_norm         = normalize(JECup);
    TH1F* JECdown_norm       = normalize(JECdown);
    TH1F* JECup_rebin        = rebin(JECup, bin_width);
    TH1F* JECdown_rebin      = rebin(JECdown, bin_width);
    TH1F* JECup_rebin_norm   = rebin(JECup_norm, bin_width);
    TH1F* JECdown_rebin_norm = rebin(JECdown_norm, bin_width);

    vector<double> JECup_rebin_norm_err   = normalize_error(JECup_rebin);
    vector<double> JECdown_rebin_norm_err = normalize_error(JECdown_rebin);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *XCup_file   = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConeup.root");
      TFile *XCdown_file = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConedown.root");
      TH1F  *XCup        = (TH1F*)XCup_file->Get(w_mass);
      TH1F  *XCdown      = (TH1F*)XCdown_file->Get(w_mass);
      allyears_XCup.push_back(XCup);
      allyears_XCdown.push_back(XCdown);
    }

    TH1F *XConeup   = AddHists(allyears_XCup, 1);
    TH1F *XConedown = AddHists(allyears_XCdown, 1);

    // Add hists from FSR and background ---------------------------------------------------------------
    if(debug) cout << "XCone+Background" << endl;

    XConeup->Add(bkg, 1);
    XConedown->Add(bkg, 1);

    TH1F* XConeup_norm         = normalize(XConeup);
    TH1F* XConedown_norm       = normalize(XConedown);
    TH1F* XConeup_rebin        = rebin(XConeup, bin_width);
    TH1F* XConedown_rebin      = rebin(XConedown, bin_width);
    TH1F* XConeup_rebin_norm   = rebin(XConeup_norm, bin_width);
    TH1F* XConedown_rebin_norm = rebin(XConedown_norm, bin_width);

    vector<double> XConeup_rebin_norm_err   = normalize_error(XConeup_rebin);
    vector<double> XConedown_rebin_norm_err = normalize_error(XConedown_rebin);

    vector<TH1F*> corrections_up_norm = {JECup_rebin_norm, XConeup_rebin_norm};
    vector<TH1F*> corrections_down_norm = {JECdown_rebin_norm, XConedown_rebin_norm};

    // #################################################################################################
    // Empty Data Bins #################################################################################
    if(debug) cout << "Empty Bins Data" << endl;

    vector<int> empty_bins_v;
    for(int bin=0; bin < number_bins; bin++){
      if(!(abs(data_rebin_norm->GetBinContent(bin+1))>0)){
        empty_bins_v.push_back(bin+1);
        if(debug) cout << "Empty Bins: " << empty_bins_v[bin] << "\n";
      }
    }
    int number_empty_bins = empty_bins_v.size();
    if(debug) cout << "Number Empty Bins: " << number_empty_bins << "\n";

    // #################################################################################################
    // Masspeack #######################################################################################
    if(debug) cout << "Masspeak Bins" << endl;

    bool usePeak;
    if(isAll&&usePeak_in) usePeak = true;
    vector<double> PeakLimit;
    double Limit = 0;
    if(usePeak){
      // PeakLimit = {200., 400., 600., 800., 1000.}; // CHANGE_PT
      // Limit     = PeakLimit[bin_width-1]; // PeakWidth is 1, 2, 3, 4, 5 -> Just take -1 for element
      Limit = 75;
    }

    vector<int> peak_bins_v; // Search_peackmethod
    if(usePeak){
      for(int bin=0; bin < number_bins; bin++){
        if((data_rebin->GetBinContent(bin+1)>Limit)){
          peak_bins_v.push_back(bin+1);
          if(debug) cout << "Peak Bins: " << peak_bins_v[bin] << "\n";
        }
      }
    }

    /*
    ██████  ██       ██████  ████████ ███████
    ██   ██ ██      ██    ██    ██    ██
    ██████  ██      ██    ██    ██    ███████
    ██      ██      ██    ██    ██         ██
    ██      ███████  ██████     ██    ███████
    */
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetLegendBorderSize(0);
    /*
    Plotting the histograms included above.
    */
    if(debug) cout << "Plots Sensitivity" << endl;

    // #################################################################################################
    // Settings ########################################################################################
    /* The vectors are only used here */
    vector<TH1F*> ttbar_ = {ttbar_rebin, ttbar_rebin_norm};
    vector<TH1F*> data_ = {data_rebin, data_rebin_norm};
    vector<vector<TH1F*>> corrections_up_ = {{JECup_rebin, XConeup_rebin}, corrections_up_norm};
    vector<vector<TH1F*>> corrections_down_ = {{JECdown_rebin, XConedown_rebin}, corrections_down_norm};

    for(unsigned int norm=0; norm<ttbar_.size(); norm++){
      ttbar_[norm]->SetTitle("");
      ttbar_[norm]->GetXaxis()->SetRangeUser(0, 180);
      ttbar_[norm]->GetYaxis()->SetRangeUser(0, ttbar_[norm]->GetMaximum()*1.2);
      ttbar_[norm]->GetXaxis()->SetNdivisions(505);
      ttbar_[norm]->GetYaxis()->SetNdivisions(505);
      ttbar_[norm]->GetXaxis()->SetTitleSize(0.05);
      ttbar_[norm]->GetYaxis()->SetTitleSize(0.04);
      ttbar_[norm]->GetXaxis()->SetTitleOffset(0.9);
      ttbar_[norm]->GetYaxis()->SetTitleOffset(1.5);
      ttbar_[norm]->GetXaxis()->SetTitle("m_{Wjet}");
      ttbar_[norm]->GetYaxis()->SetTitle("");
      ttbar_[norm]->SetLineWidth(2);  // ttbar hist style
      ttbar_[norm]->SetLineColor(kRed);

      // Data --------------------------------------------------------------------------------------------
      data_[norm]->SetMarkerStyle(8);  // data hist style
      data_[norm]->SetMarkerColor(kBlack);
    }

    // Legend ------------------------------------------------------------------------------------------
    TLegend *leg;

    /* The Mass Plots are drawn in a seperated file - JEC_SYS_MassPlots.cc */
    /*
    ██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██ ███████
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██ ██
    ██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██ ███████
    ██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██      ██
    ██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████ ███████
    */
    /*
    Looking at each bin seperatly. Getting the bin content of the nominal file and
    the variations and making a fit through the values. The variation UP and DOWN are
    assigned to the values +1 and -1.
    Afterward the fit parameters will be stored in vectors for later purpose.
    */

    if(debug) print_seperater();
    if(debug) cout << "Start: Dependency in WJet bins" << endl;

    // #################################################################################################
    // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
    if(debug) cout << "Title for bins" << endl;
    vector<TString> mass_bin;
    for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin*bin_width)+" < m_{Wjet} < "+to_string((bin+1)*bin_width));

    if(debug){
      cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << bin_width << "\n\n";
      for(unsigned int bin=0; bin<mass_bin.size(); bin++){
        if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue;
        cout << mass_bin[bin] << '\n';
      }
    }

    vector<TString> number_bin;
    for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

    // #################################################################################################
    // Bin Content and error ###########################################################################
    if(debug) print_seperater();
    if(debug) cout << "Vector for Content and Error" << endl;

    vector<TH1F*> XCone_ordered = {XConedown_rebin_norm, XConeup_rebin_norm};
    vector<TH1F*> JEC_ordered = {JECdown_rebin_norm, JECup_rebin_norm};

    // for the chi2 Method -----------------------------------------------------------------------------
    if(debug) cout << "Vectors for Chi2" << endl;
    vector<double> data_rebin_norm_bin_content;
    for(int bin=1; bin < number_bins+1; bin++) data_rebin_norm_bin_content.push_back(data_rebin_norm->GetBinContent(bin));

    // #################################################################################################
    // Cout Bin Error ##################################################################################
    if(debug){
      cout << "" << '\n' << fixed;
      cout << "|================================================================================================|" << endl;
      cout << "| ----------------------------------------- Bin Content -----------------------------------------|" << endl;
      cout << "| -----------------------------------------------------------------------------------------------|" << endl;
      cout << "|  Bin | Data Content |   JEC down   |  XCone down  |   nominal    |    JEC up    |   XCone up   |" << endl;
      cout << "| -----|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
      // In this loop no additional treatment is  necessary.
      for(unsigned int i=0; i<number_bins; i++){
        int bin;
        double data, Jd, tt, Ju, Xd, Xu;
        bin  = i+1;
        data = data_rebin_norm_bin_content[i];
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
        string bin, data_str, Jd_str, tt_str, Ju_str, avg_str, cenavg_str;
        double avg = (JECdown_rebin_norm_err[i]+ttbar_rebin_norm_err[i]+JECup_rebin_norm_err[i])/3;
        double cenavg;
        if(abs(avg) > 0) cenavg = ttbar_rebin_norm_err[i]/avg;                 // make sure bin error is not 0
        else             cenavg = avg;
        if(i+1<10)     bin = to_string(i+1);                       // i+1 to get bin number
        if(i+1>=10)    bin = to_string(i+1);                        // +start_bin to skip first bins
        data_str           = to_string(data_rebin_norm_bin_content[i]);
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
        string bin, data_str, Xd_str, tt_str, Xu_str, avg_str, cenavg_str;
        double avg = (XConedown_rebin_norm_err[i]+ttbar_rebin_norm_err[i]+XConeup_rebin_norm_err[i])/3;
        double cenavg;
        if(abs(avg) > 0) cenavg = ttbar_rebin_norm_err[i]/avg;                // make sure bin error is not 0
        else             cenavg = avg;
        if(i+1<10)     bin = to_string(i+1);                       // i+1 to get bin number
        if(i+1>=10)    bin = to_string(i+1);                        // +start_bin to skip first bins
        data_str   = to_string(data_rebin_norm_bin_content[i]);
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
    vector<double> bin_error_one_sigma_avg, bin_error_one_sigma_central;
    for(int bin=0; bin<number_bins; bin++){
      bin_error_one_sigma_central.push_back(ttbar_rebin_norm_err[bin]);
      bin_error_one_sigma_avg.push_back((XConedown_rebin_norm_err[bin]+ttbar_rebin_norm_err[bin]+XConeup_rebin_norm_err[bin])/3);
    }

    // #################################################################################################
    // Some Declarations ###############################################################################
    if(debug) cout << "Define fit input" << endl;

    // Fit ---------------------------------------------------------------------------------------------
    /* For similiar functions (xy2 & x2y) or (xyy2 & xx2y) the probability for y2 functions is always higher.
    Therefore, functions with only x2 are excluded to make the code more clean.*/
    TGraph2DErrors *bin_fit;
    TString str_fit_lin  = "[0] + [1]*x + [2]*y"; // 3
    TString str_fit_xy2  = "[0] + [1]*x + [2]*y*y"; // 3
    TString str_fit_x2y  = "[0] + [1]*x*x + [2]*y"; // 3
    TString str_fit_xyy2 = "[0] + [1]*x + [2]*y + [3]*y*y"; // 4
    TString str_fit_xx2y = "[0] + [1]*x + [2]*x*x + [3]*y"; // 4
    TString str_fit_quad = "[0] + [1]*x*x + [2]*y*y"; // 3
    TString str_fit_poly = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y"; // 5
    vector<vector<TString>> str_fits  = {{str_fit_lin}, {str_fit_xy2, str_fit_x2y}, {str_fit_xyy2, str_fit_xx2y}, {str_fit_quad}, {str_fit_poly}};
    vector<vector<TString>> name_fits = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
    double npar=0;
    vector<vector<int>> ndfs  = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
    vector<vector<int>> n_fit_para  = {{3}, {3, 3}, {4, 4}, {3}, {5}};
    vector<int> all_n_para;
    vector<TString> str_all_fits; vector<TF2*> all_fits; // Fill later
    vector<vector<double>> fits_parameters;     // filled in for-loop - structer purpose
    vector<vector<double>> fits_parameters_err; // filled in for-loop - structer purpose
    vector<double> chi2_parameters;

    int n_factors = 5; // JECdown, XConedown, ttbar, XConeup,JECup
    vector<TH1F*> hists              = {ttbar_rebin_norm, JECdown_rebin_norm, JECup_rebin_norm, XConedown_rebin_norm, XConeup_rebin_norm};
    vector<vector<double>> hists_err = {ttbar_rebin_norm_err, JECdown_rebin_norm_err, JECup_rebin_norm_err, XConedown_rebin_norm_err, XConeup_rebin_norm_err};
    vector<vector<double>> hists_cont, hists_cont_err;
    vector<double> factor_x          = {0.0, -1.0,  1.0,  0.0,  0.0};
    vector<double> factor_y          = {0.0,  0.0,  0.0, -1.0,  1.0};
    vector<double> error_x           = {0.0,  0.0,  0.0, 0.0,  0.0};
    vector<double> error_y           = {0.0,  0.0,  0.0, 0.0,  0.0};

    // #################################################################################################
    // Hists for Fits ##################################################################################
    TH1I used_fits = TH1I("used_fits", "used_fits", number_bins, 1, number_bins+1);
    TH1I number_fits = TH1I("number_fits", "number_fits", 5, 0, 5);

    /*
    ███████ ██ ████████     ██       ██████   ██████  ██████
    ██      ██    ██        ██      ██    ██ ██    ██ ██   ██
    █████   ██    ██        ██      ██    ██ ██    ██ ██████
    ██      ██    ██        ██      ██    ██ ██    ██ ██
    ██      ██    ██        ███████  ██████   ██████  ██
    */
    if(debug) cout << "Start: looping over bins" << endl;

    for(int bin=0; bin < number_bins; bin++){
      if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue;
      if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())) continue;  // Search_peackmethod
      cout << "\n";
      if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      if(debug) cout << mass_bin[bin] << endl;
      // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      // #################################################################################################
      // Fill Vectors ####################################################################################
      if(debug) cout << "Get Bin Content" << endl;

      vector<double> bin_content, bin_content_err;
      for(int factor=0; factor<n_factors; factor++){
        bin_content.push_back(hists[factor]->GetBinContent(bin+1));
        bin_content_err.push_back(hists_err[factor][bin]);
      }
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

        double prob = prob1;
        if(twoFits && (prob1 < prob2)){
          prob = prob2;
          index_fit = 1;
        }
        if(debug) cout << "ndfs: " << ndfs[order][index_fit] << " | chi2: " << chi2_fit <<" | prob: " << prob << endl;

        // -----------------------------------------------------------------------
        if(prob>0.05 || ndfs[order][0]==0 || useOnly_lin){
          cout << "Bin " << number_bin[bin] << ": A " << GREEN << name_fits[order][index_fit] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
          used_fit = new TF2(*fits[order][index_fit]); // in if clause
          str_all_fits.push_back(str_fits[order][index_fit]);
          all_fits.push_back(fits[order][index_fit]);
          npar += n_fit_para[order][index_fit];
          all_n_para.push_back(n_fit_para[order][index_fit]);
          order_index=order;
          number_fits.AddBinContent(order+1);
          if(twoFits){
            if(index_fit==0) used_fits.AddBinContent(bin+1, order+1-0.3);
            else             used_fits.AddBinContent(bin+1, order+1+0.3);
          } else             used_fits.AddBinContent(bin+1, order+1);
          break;
        }
        else cout << "Bin " << number_bin[bin] << ": "<< RED << name_fits[order][index_fit] << RESET << " fit is rejected (" << prob << ")" << endl;
      }

      // #################################################################################################
      // Plot  ###########################################################################################
      if(debug) cout << "Plot Fitting" << endl;

      // -------------------------------------------------------------------------------------------------
      // Plot Points -------------------------------------------------------------------------------------
      double max_bin_content = GetMaxValue(bin_content);
      double min_bin_content = GetMinValue(bin_content);
      used_fit->SetRange(-1, -1, min_bin_content*0.8, 1, 1, max_bin_content*1.2);

      // Additional Points to plot the z axis properly ---------------------------------------------------
      /* The additional points are necessary to scale the z axis properly. One can not zoom out of the z axis without
      cutting of the errors above the highest- and below the lowest point. The trick is to add two additional points
      to streth the z axis. Afterward one can zoom in (and not out).*/
      bin_fit->SetPoint(5, 0., 0., max_bin_content*1.8);
      bin_fit->SetPointError(5, 0., 0., 0.);
      bin_fit->SetPoint(6, 0., 0., min_bin_content*0.2);
      bin_fit->SetPointError(6, 0., 0., 0.);

      bin_fit->SetTitle(mass_bin[bin]);
      bin_fit->GetHistogram()->GetXaxis()->SetTitle("JEC");
      bin_fit->GetHistogram()->GetYaxis()->SetTitle("XCone");
      bin_fit->GetHistogram()->GetXaxis()->SetTitleOffset(-0.5);
      bin_fit->GetHistogram()->GetYaxis()->SetTitleOffset(-0.5);

      bin_fit->SetMarkerStyle(20);
      TCanvas *B = new TCanvas(number_bin[bin]+"B"+ptbin_str, "B", 600, 600);
      bin_fit->GetHistogram()->GetXaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetYaxis()->SetRangeUser(-2, 2);
      bin_fit->GetHistogram()->GetZaxis()->SetRangeUser(min_bin_content*0.8, max_bin_content*1.2);

      bin_fit->Draw("err P");
      B->SaveAs(save_path_general+"/single_bins/Bin"+number_bin[bin]+"_Points.pdf");
      used_fit->Draw("same surf4");
      B->SaveAs(save_path_general+"/single_bins/Bin"+number_bin[bin]+"_Combined"+addition+".pdf");
      B->Clear();
      used_fit->SetRange(-2, -2, -1, 2, 2, 1);
      used_fit->Draw("surf4");
      B->SaveAs(save_path_general+"/single_bins/Bin"+number_bin[bin]+"_Fit"+addition+".pdf");
      B->Clear();

      // Remove Additional Points for fit ----------------------------------------------------------------
      bin_fit->RemovePoint(6);
      bin_fit->RemovePoint(5);

      // #################################################################################################
      // Fill Fit Parameters #############################################################################
      if(debug) cout << "Fill Fit Parameters" << endl;

      vector<double> stacker, stacker_err; // Ablage
      stacker.push_back(data_rebin_norm_bin_content[bin]);
      chi2_parameters.push_back(data_rebin_norm_bin_content[bin]);
      for(int ipar=0; ipar<n_fit_para[order_index][index_fit]; ipar++){
        chi2_parameters.push_back(used_fit->GetParameter(ipar));
        stacker.push_back(used_fit->GetParameter(ipar));
        stacker_err.push_back(used_fit->GetParError(ipar));
      }
      stacker.push_back(data_rebin_norm_err[bin]);
      chi2_parameters.push_back(data_rebin_norm_err[bin]);

      fits_parameters.push_back(stacker);
      fits_parameters_err.push_back(stacker_err);
    }
    cout << endl;

    int npar2 = 0;
    for(unsigned int fit=0;fit<fits_parameters.size(); fit++){
      for(unsigned int ipar=0;ipar<fits_parameters[fit].size(); ipar++) npar2++;
    }
    if(debug) cout << "Number Parameter: " << npar2 << endl;
    if(debug) cout << "Chi2 Parameter: " << chi2_parameters.size() << endl;


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
    leg = new TLegend(0.4,0.65,0.7,0.85);
    for(unsigned int i=0; i<str_fits.size();i++) leg->AddEntry((TObject*)0,to_string(i+1)+": "+str_fits[i][0],"");
    leg->SetTextSize(0.025);
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
    int count_para = 0;

    // #################################################################################################
    // Replace Parameters ##############################################################################
    if(debug) cout << "Chi2 - Replace parameters\n";
    if(!usePeak && str_all_fits.size()!= number_bins-number_empty_bins) throw runtime_error("Something is wrong with the number of fits");
    // else if(usePeak && str_all_fits.size()!=end_bin-start_bin+1)        throw runtime_error("Something is wrong with the number of fits");

    for(unsigned int fit=0; fit<str_all_fits.size(); fit++){
      TString str_fit_full = " + ((["+to_string(count_para)+"] - ";
      TString str_fit = str_all_fits[fit];
      /* Treat the first fit differently. Replacing the parametes leads to an overlap.
      [1]-[1]*x-[2]*y -> [2]-[2]*x-[2]*y -> [3]-[3]*x-[3]*y -> etc. */
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
      str_fit_full += ")/["+to_string(count_para)+"])^2";
      count_para++; // data bin content - for next iteration
      if(debug) cout << str_fit_full << "\n";
      str_chi2_function += str_fit_full;
    }

    // #################################################################################################
    // Set all parameters ##############################################################################
    TF2 *chi2_function = new TF2("chi2_function"+addition, str_chi2_function, -10, 10, -10, 10);
    for(unsigned int ipar=0; ipar<chi2_parameters.size(); ipar++) chi2_function->SetParameter(ipar, chi2_parameters[ipar]);
    chi2_function->SetTitle("");
    chi2_function->GetXaxis()->SetTitle("JEC");
    chi2_function->GetYaxis()->SetTitle("XCone");
    chi2_function->SetFillStyle(1000);
    chi2_function->SetLineWidth(1);

    // #################################################################################################
    // Add chi2 functions ##############################################################################
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
    /* Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447. */
    double twoD_minX, twoD_minY;
    double twoD_minZ = chi2_function->GetMinimumXY(twoD_minX,twoD_minY);
    cout  << "z_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n' << endl;
    int number_para = chi2_function->GetNpar();

    fstream jec_txt;
    jec_txt.open(save_path_general+"/jec_factor"+addition+".txt", ios::out);
    jec_txt << number_bins << " bins [0, 180] GeV\n\n";
    jec_txt << "zmin: " << twoD_minZ << "\nx   :  " << twoD_minX << "\ny   :  " << twoD_minY << endl;
    jec_txt.close();

    // #################################################################################################
    // Evaluate 1sigma #################################################################################

    // Draw Points -------------------------------------------------------------------------------------
    vector<vector<double>> points;
    auto start = high_resolution_clock::now(); // Calculation time - start
    points = FindXY(chi2_function, twoD_minZ+2.3, twoD_minX-2, twoD_minX+2, twoD_minY-2, twoD_minY+2, 1000, 0.001);

    auto stop = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Numeric solution for 1sigma estimation took " << GREEN << duration.count()/1000 << "s" << RESET << endl;
    if(debug) cout << "Number Points at 1sigma: " << points.size() << endl;
    // if(debug) for(unsigned int i=0; i<points.size(); i++) cout << "x: " << points[i][0] << " | y: " << points[i][1] << " | z: " << points[i][2] << endl;

    TGraph2D *zmin_point = new TGraph2D();
    zmin_point->SetPoint(0, twoD_minX, twoD_minY, twoD_minZ);
    zmin_point->SetMarkerColor(kBlack);
    zmin_point->SetMarkerStyle(kFullCircle);
    zmin_point->SetMarkerSize(0.5);
    TGraph2D *sigma_points = new TGraph2D();
    for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
    // for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1]);

    sigma_points->SetMarkerColor(kRed);
    sigma_points->SetMarkerStyle(kFullCircle);
    sigma_points->SetMarkerSize(0.1);

    chi2_function->SetRange(-3, -3, 3, 3);
    // #################################################################################################
    // Plot ############################################################################################
    if(debug) cout<< "\n2D Chi2 - plot\n";
    gErrorIgnoreLevel = kWarning; // suppress TCanvas output
    TString option = "cont4z";
    TString option_add = "";
    // Contours ----------------------------------------------------------------------------------------
    chi2_function->SetContour(50);
    gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
    if(ptbin==0) TColor::InvertPalette(); // CHANGE_PT
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

    gErrorIgnoreLevel = oldLevel; // Reset TCanvas output

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
      vector<double> oneD_factor     = {-1,0,1};
      vector<double> oneD_factor_err = { 0,0,0};

      // #################################################################################################
      // Bins ############################################################################################
      int empty_bins = 0; // to skip empty bins
      for(int bin=0; bin<number_bins; bin++){

        if(debug) cout << "------------------------------------------------------- " << bin+1 << endl;
        if(usePeak && !(find(peak_bins_v.begin(), peak_bins_v.end(), bin+1) != peak_bins_v.end())){// Search_peackmethod
          empty_bins++;
          continue;
        }
        if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()){
          empty_bins++;
          continue;
        }

        int used_bin = bin-empty_bins;
        vector<double> oneD_events, oneD_events_err;
        if(factor==0){
          oneD_events     = {hists_cont[used_bin][1], hists_cont[used_bin][0], hists_cont[used_bin][2]};
          oneD_events_err = {hists_cont_err[used_bin][1], hists_cont_err[used_bin][0], hists_cont_err[used_bin][2]};
        } else if(factor==1){
          oneD_events     = {hists_cont[used_bin][3], hists_cont[used_bin][0], hists_cont[used_bin][4]};
          oneD_events_err = {hists_cont_err[used_bin][3], hists_cont_err[used_bin][0], hists_cont_err[used_bin][4]};
        }
        // if(debug) cout << setw(8) << oneD_events[0] << "   " << setw(8) << oneD_events[1] << "   " << setw(8) << oneD_events[2] << "   \n";
        // if(debug) cout << setw(8) << oneD_events_err[0] << "   " << setw(8) << oneD_events_err[1] << "   " << setw(8) << oneD_events_err[2] << "   \n";

        // XCone -----------------------------------------------------------------------------------------
        TString str_fit = str_all_fits[used_bin];
        if(factor==0) str_fit.ReplaceAll("y", "0");
        else if(factor==1){
          str_fit.ReplaceAll("x", "0");
          str_fit.ReplaceAll("y", "x");
        }

        TF1 *fit = new TF1("fit", str_fit, -2, 2);
        TGraphErrors *graph = new TGraphErrors(3, &oneD_factor[0], &oneD_events[0], &oneD_factor_err[0], &oneD_events_err[0]);
        for(int ipar=0; ipar<all_n_para[used_bin]; ipar++) fit->SetParameter(ipar, all_fits[used_bin]->GetParameter(ipar));

        // #################################################################################################
        // Plot ############################################################################################
        graph->SetTitle(mass_bin[bin]);
        gErrorIgnoreLevel = kWarning; // suppress TCanvas output
        TCanvas *E2 = new TCanvas(number_bin[bin]+"E2_"+factor_name[factor]+to_string(ptbin), "E2", 600, 600);
        if(factor==0) graph->GetXaxis()->SetTitle("JEC");
        else          graph->GetXaxis()->SetTitle("XCone");
        graph->Draw("AP");
        fit->Draw("SAME");
        E2->SaveAs(save_path_general+"/projection/Bin"+number_bin[bin]+"_"+factor_name[factor]+addition+".pdf");
        gErrorIgnoreLevel = oldLevel; // Reset TCanvas output
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
        if(debug) system("python python/single_bins_peak_latex.py "+year+" "+to_string(number_bins));
        else      system("python python/single_bins_peak_latex.py "+year+" "+to_string(number_bins)+" &> /dev/null");
      } else{
        if(debug) system("python python/single_bins_latex.py "+year+" "+to_string(number_bins)+" "+addition);
        else      system("python python/single_bins_latex.py "+year+" "+to_string(number_bins)+" "+addition+" &> /dev/null");
      }
    }
  }


  /*
  .█████  ██████  ██████       ██████ ██   ██ ██ ██████
  ██   ██ ██   ██ ██   ██     ██      ██   ██ ██      ██
  ███████ ██   ██ ██   ██     ██      ███████ ██  █████
  ██   ██ ██   ██ ██   ██     ██      ██   ██ ██ ██
  ██   ██ ██████  ██████       ██████ ██   ██ ██ ███████
  */

  cout << chi2_str_hh.Length() << "       " << chi2_parameters_hh.size() << endl;
  cout << chi2_str_hl.Length() << "       " << chi2_parameters_hl.size() << endl;
  cout << chi2_str_lh.Length() << "       " << chi2_parameters_lh.size() << endl;
  cout << chi2_str_ll.Length() << "       " << chi2_parameters_ll.size() << endl;
  cout << endl;

  vector<TString> all_string = {chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll};
  vector<vector<double>> all_parameters = {chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll};


  // #################################################################################################
  // Plot ############################################################################################
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
  // Plot ############################################################################################
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
  /* Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447. */
  double twoD_minX, twoD_minY;
  double twoD_minZ = full_chi2_function->GetMinimumXY(twoD_minX,twoD_minY);
  cout  << "z_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n' << endl;
  int number_para = full_chi2_function->GetNpar();

  fstream jec_txt;
  jec_txt.open(save_path_general+"/jec_factor_all.txt", ios::out);
  jec_txt << number_bins << " bins [0, 180] GeV\n\n";
  jec_txt << "zmin: " << twoD_minZ << "\nx   :  " << twoD_minX << "\ny   :  " << twoD_minY << endl;
  jec_txt.close();

  // #################################################################################################
  // Evaluate 1sigma #################################################################################

  // Draw Points -------------------------------------------------------------------------------------
  vector<vector<double>> points;
  auto start = high_resolution_clock::now(); // Calculation time - start
  points = FindXY(full_chi2_function, twoD_minZ+2.3, twoD_minX-1.5, twoD_minX+1.5, twoD_minY-1.5, twoD_minY+1.5, 60000, 0.0001);
  auto stop = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  cout << "Numeric solution for 1sigma estimation took " << GREEN << duration.count()/1000 << "s" << RESET << endl;
  if(debug) cout << "Number Points at 1sigma: " << points.size() << endl;

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
  if(debug) cout<< "\n2D Chi2 All - plot\n";

  gErrorIgnoreLevel = kWarning; // suppress TCanvas output
  TString option = "cont4z";
  TString option_add = "";
  // Contours ----------------------------------------------------------------------------------------
  full_chi2_function->SetContour(50);
  gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
  TCanvas *Z = new TCanvas("Z","Z", 600, 600);
  Z->SetRightMargin(0.12);
  Z->SetLogz();
  full_chi2_function->Draw("cont4z");
  Z->SaveAs(save_path_general+"/chi2_all.pdf");
  Z->SetTheta(90);
  Z->SetPhi(0);
  zmin_point->Draw("SAME P");
  sigma_points->Draw("SAME P");
  Z->SaveAs(save_path_general+"/chi2_all_points.pdf");
  Z->Clear();

  gErrorIgnoreLevel = oldLevel; // Reset TCanvas output

  /*
  ███████ ██      ██      ██ ██████  ███████ ███████
  ██      ██      ██      ██ ██   ██ ██      ██
  █████   ██      ██      ██ ██████  ███████ █████
  ██      ██      ██      ██ ██           ██ ██
  ███████ ███████ ███████ ██ ██      ███████ ███████
  */

  // #################################################################################################
  // Extreme Points ##################################################################################
  if(debug) cout<< "\n Find extrme points of ellipse\n";

  double dmax1, dmax2, dmin1, dmin2;
  vector<vector<double>> Ps_max, Ps_min;
  unsigned int n_points = points.size();
  cout << "Looking at " << n_points << " points" << endl;

  vector<double> dist_right, dist_left, dist_all, xValues, yValues;
  for(unsigned int point=0; point<n_points; point++){
    double x=points[point][0];
    double y=points[point][1];
    dist_all.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    if(x>twoD_minX) dist_right.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    if(x<twoD_minX) dist_left.push_back(sqrt(pow(twoD_minX-x, 2)+pow(twoD_minY-y, 2)));
    xValues.push_back(x);
    yValues.push_back(y);
  }

  // #################################################################################################
  // find index ######################################################################################
  if(debug) cout<< "\n2D Chi2 All - find index\n";

  double index_max_l_s = ReturnIndex_High(dist_left);
  double index_max_r_s = ReturnIndex_High(dist_right);
  double index_min_l_s = ReturnIndex_Low(dist_left);
  double index_min_r_s = ReturnIndex_Low(dist_right);

  double index_max_l;
  double index_max_r;
  double index_min_l;
  double index_min_r;

  // Compare ------------------------------------------------------------------------------------------
  for(unsigned int i=0; i<dist_all.size(); i++){
    if(dist_left[index_max_l_s]==dist_all[i]) index_max_l = i;
    if(dist_right[index_max_r_s]==dist_all[i]) index_max_r = i;
    if(dist_left[index_min_l_s]==dist_all[i]) index_min_l = i;
    if(dist_right[index_min_r_s]==dist_all[i]) index_min_r = i;
  }
  cout << index_max_l << " " << index_max_r << " " << index_min_l << " " << index_min_r << " " << endl;

  Ps_max.push_back(points[index_max_l]);
  Ps_max.push_back(points[index_max_r]);
  Ps_max.push_back(points[index_min_l]);
  Ps_max.push_back(points[index_min_r]);

  cout << "Number Max Points: " << Ps_max.size() << endl;
  for(unsigned int point=0; point<Ps_max.size(); point++) cout << "(" << setw(9) << Ps_max[point][0] << ", " << setw(9) << Ps_max[point][1] << ")\n";
  cout << "Number Min Points: " << Ps_min.size() << endl;
  for(unsigned int point=0; point<Ps_min.size(); point++) cout << "(" << setw(9) << Ps_min[point][0] << ", " << setw(9) << Ps_min[point][1] << ")\n";

  // #################################################################################################
  // in file #########################################################################################
  if(debug) cout<< "\n2D Chi2 All - in file\n";

  jec_txt.open(save_path_general+"/jec_factor_all.txt", ios::out | ios::app);
  jec_txt << "\nP1_up:   " << "(" << setw(6) << Round(Ps_max[0][0], 3) << ", " << setw(6) << Round(Ps_max[0][1], 3) << ")" <<  endl;
  jec_txt << "P2_up:   " << "(" << setw(6) << Round(Ps_max[1][0], 3) << ", " << setw(6) << Round(Ps_max[1][1], 3) << ")" <<  endl;
  jec_txt << "P1_down: " << "(" << setw(6) << Round(Ps_max[2][0], 3) << ", " << setw(6) << Round(Ps_max[2][1], 3) << ")" <<  endl;
  jec_txt << "P2_down: " << "(" << setw(6) << Round(Ps_max[3][0], 3) << ", " << setw(6) << Round(Ps_max[3][1], 3) << ")" <<  endl;
  jec_txt.close();

  // #################################################################################################
  // Settings ########################################################################################
  if(debug) cout<< "\n2D Chi2 All - settings\n";
  TGraph2D *extreme_points = new TGraph2D();
  for(unsigned int i=0; i<Ps_min.size(); i++) extreme_points->SetPoint(i, Ps_min[i][0], Ps_min[i][1], Ps_min[i][2]);
  for(unsigned int i=0; i<Ps_max.size(); i++) extreme_points->SetPoint(i, Ps_max[i][0], Ps_max[i][1], Ps_max[i][2]);

  extreme_points->SetMarkerColor(kBlue);
  extreme_points->SetMarkerStyle(kFullCircle);
  extreme_points->SetMarkerSize(0.3);


  // #################################################################################################
  // Plot ############################################################################################
  if(debug) cout<< "\n2D Chi2 All - plot\n";
  full_chi2_function->SetContour(50);
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
