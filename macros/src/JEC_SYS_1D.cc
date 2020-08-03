#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  /*
  Explanation:
  1) This script loops many times over bins. The bin number starts at 1 (e.g. for the function GetBinContent)
  On the other hand, the bin content or titles are stored in a vector. The vector starts at 0.
  Mostly, the loops start at bin=0. Therefor, using the function GetBinContent the bin needs to be bin+1
  */
  // #################################################################################################
  // Declare different variables used in the code ####################################################
  print_seperater();

  bool debug = false;
  bool print = false;
  bool wmass_peak = false; // changing number of bins
  TString MC_uncert = "central"; // central avg nominal
  cout.precision(6);

  if(argc != 4){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number> <selection>\n" << "rebin_number: hist->Rebin(rebin_number)\n" << "Selection: match, btagcut, btag_sel, compare, min_mass\n";
    return 0;
  }

  TString reconst = argv[3];

  // Rebin -------------------------------------------------------------------------------------------
  /*
  Right now every Data bin will be excluded, if the BinContent is zero. This leads to a problem,
  since in the first bins (at a binning of 1 GeV) the bin Content is (e.g.) 000001000.
  Therefor, the rebinning has to be taken into account.
  Rebin  4: exclude bin 1-4
  Rebin  5: exclude bin 1-2
  Rebin 10: exclude bin 1
  */
  int rebin_hist = atoi(argv[2]); // atoi() c-string representing integer into int
  int number_bins = 180/rebin_hist;
  TString str_number_bins = to_string(number_bins);

  // Range -------------------------------------------------------------------------------------------
  int peak_width = 20; int w_peak = 80; int max_mass = 200; // Only looking at the wmass peak. Boundaries are abitrary
  int lower_limit = w_peak-peak_width; int upper_limit = w_peak+peak_width;
  int first_bins = 0; int last_bins = 0;

  if(wmass_peak){
    first_bins = lower_limit/rebin_hist;
    last_bins = (max_mass-upper_limit)/rebin_hist;
    // number_bins -= last_bins+first_bins; // number of bins only in the mass peak range
    number_bins -= last_bins;
  }
  else{
    /*
    exclude bins with DataBinContent = 0 -- HARDCODED
    */
    if     (rebin_hist==4)  first_bins = 4;
    else if(rebin_hist==5)  first_bins = 3; // for btag_sel == 2 !!!!!!!!!
    else if(rebin_hist==10) first_bins = 1;
    else throw runtime_error("At this point, only the rebinning 4, 5 and 10 works");
    /*
    The reason why the first empty bins are included is, because in later histograms
    the empty bins are shown for a better comparison of the results. They could be
    excluded, but it is not done yet.
    */
  }

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false;
  if     (strcmp(year, "2016")==0) is16 = true;
  else if(strcmp(year, "2017")==0) is17 = true;
  else if(strcmp(year, "2018")==0) is18 = true;
  else throw runtime_error("Give me the correct year please (2016, 2017 or 2018)");

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

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"+year+"/"+reconst;

  if(mkdir(save_path_general+"/rebin"+str_number_bins,0777) != 0) cout << "couldn't creat rebin folder!\n";
  else cout << save_path_general+"/rebin"+str_number_bins+" created\n";
  save_path_general += "/rebin"+str_number_bins; // New save_path

  if(mkdir(save_path_general+"/JEC",0777) != 0) cout << "couldn't creat rebin/JEC folder!\n";
  else cout << save_path_general+"/JEC created\n";
  if(mkdir(save_path_general+"/JEC/single_bins",0777) != 0) cout << "couldn't creat /JEC/single_bins folder!\n";
  else cout << save_path_general+"/JEC/single_bins created\n";
  TString save_path_JEC = save_path_general+"/JEC"; // New save_path

  if(mkdir(save_path_general+"/XCone",0777) != 0) cout << "couldn't creat rebin/XCone folder!\n";
  else cout << save_path_general+"/XCone created\n";
  if(mkdir(save_path_general+"/XCone/single_bins",0777) != 0) cout << "couldn't creat /XCone/single_bins folder!\n";
  else cout << save_path_general+"/XCone/single_bins created\n";
  TString save_path_xcone = save_path_general+"/XCone"; // New save_path

  // if(wmass_peak){
  //   if(mkdir(save_path_general+"/masspeak_width"+to_string(peak_width),0777) != 0) cout << "couldn't creat masspeak folder!\n";
  //   else cout << save_path_general+"/masspeak_width"+to_string(peak_width)+" created\n";
  //   if(mkdir(save_path_general+"/masspeak_width"+to_string(peak_width)+"/single_bins",0777) != 0) cout << "couldn't creat masspeak/single_bins folder!\n";
  //   else cout << save_path_general+"/masspeak_width"+to_string(peak_width)+"/single_bins created\n";
  //   save_path += "/masspeak_width"+to_string(peak_width);
  // }

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
  if(reconst=="min_mass") w_mass = hist_class+"wmass_min";                      // minimal mass from subjet combination
  if(reconst=="btag")     w_mass = hist_class+"wmass_match";                    // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
  if(reconst=="btag_cut") w_mass = hist_class+"wmass_btagcut";                  // same mass as above mit with btag_high > 0.7
  if(reconst=="compare")  w_mass = hist_class+"wmass_compare";                  // mass closest to Wmass from subjet combination
  if(reconst=="btag_sel") w_mass = hist_class+"wmass_btagcut_one_btag_subjet";  // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
  //                                                                               Selection: btag_high>0.7; ONE subjet with dr(ak4, subjet)<0.4; ONE high btag

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

  TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file      = new TFile(dir+year+"/muon/"+data_path);

  TH1F *data            = (TH1F*)data_file->Get(w_mass);
  TH1F *data_rebin      = rebin(data, rebin_hist);
  TH1F *data_norm       = normalize(data);
  TH1F *data_rebin_norm = normalize(data_rebin);

  // Plot Setting ------------------------------------------------------------------------------------
  data_rebin_norm->SetMarkerStyle(8);  // data hist style
  data_rebin_norm->SetMarkerColor(kBlack);

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file  = new TFile(dir+year+"/muon/"+ttbar_path);

  TH1F *ttbar            = (TH1F*)ttbar_file->Get(w_mass);
  ttbar->Add(bkg, 1);
  TH1F* ttbar_norm       = normalize(ttbar);
  TH1F* ttbar_rebin      = rebin(ttbar, rebin_hist);
  TH1F* ttbar_rebin_norm = rebin(ttbar_norm, rebin_hist);

  // Plot Setting ------------------------------------------------------------------------------------
  ttbar_rebin_norm->SetTitle("");
  if(!wmass_peak) ttbar_rebin_norm->GetXaxis()->SetRangeUser(0, 180);
  ttbar_rebin_norm->GetYaxis()->SetRangeUser(0, ttbar_rebin_norm->GetMaximum()*1.2);
  ttbar_rebin_norm->GetXaxis()->SetNdivisions(505);
  ttbar_rebin_norm->GetYaxis()->SetNdivisions(505);
  ttbar_rebin_norm->GetXaxis()->SetTitleSize(0.05);
  ttbar_rebin_norm->GetYaxis()->SetTitleSize(0.04);
  ttbar_rebin_norm->GetXaxis()->SetTitleOffset(0.9);
  ttbar_rebin_norm->GetYaxis()->SetTitleOffset(1.5);
  ttbar_rebin_norm->GetXaxis()->SetTitle("m_{Wjet}");
  ttbar_rebin_norm->GetYaxis()->SetTitle("");
  ttbar_rebin_norm->SetLineWidth(2);  // ttbar hist style
  ttbar_rebin_norm->SetLineColor(kRed);

  // #################################################################################################
  // Get SYS #########################################################################################
  if(debug) cout << "JEC" << endl;

  TFile *JECup_file   = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
  TFile *JECdown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");

  TH1F *JECup   = (TH1F*)JECup_file->Get(w_mass);
  TH1F *JECdown = (TH1F*)JECdown_file->Get(w_mass);

  // Add hists from FSR and background ---------------------------------------------------------------
  if(debug) cout << "JEC+Background" << endl;

  JECup->Add(bkg, 1);
  JECdown->Add(bkg, 1);

  TH1F* JECup_norm         = normalize(JECup);
  TH1F* JECdown_norm       = normalize(JECdown);
  TH1F* JECup_rebin        = rebin(JECup, rebin_hist);
  TH1F* JECdown_rebin      = rebin(JECdown, rebin_hist);
  TH1F* JECup_rebin_norm   = rebin(JECup_norm, rebin_hist);
  TH1F* JECdown_rebin_norm = rebin(JECdown_norm, rebin_hist);

  // ------------------------------------------------------------------------------------------------
  if(debug) cout << "XCone" << endl;

  TFile *XConeup_file   = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConeup.root");
  TFile *XConedown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConedown.root");

  TH1F *XConeup   = (TH1F*)JECup_file->Get(w_mass);
  TH1F *XConedown = (TH1F*)JECdown_file->Get(w_mass);

  // Add hists from FSR and background ---------------------------------------------------------------
  if(debug) cout << "XCone+Background" << endl;

  XConeup->Add(bkg, 1);
  XConedown->Add(bkg, 1);

  TH1F* XConeup_norm         = normalize(XConeup);
  TH1F* XConedown_norm       = normalize(XConedown);
  TH1F* XConeup_rebin        = rebin(XConeup, rebin_hist);
  TH1F* XConedown_rebin      = rebin(XConedown, rebin_hist);
  TH1F* XConeup_rebin_norm   = rebin(XConeup_norm, rebin_hist);
  TH1F* XConedown_rebin_norm = rebin(XConedown_norm, rebin_hist);

  vector<TH1F*> corrections_up_v   = {JECup_rebin_norm, XConeup_rebin_norm};
  vector<TH1F*> corrections_down_v = {JECdown_rebin_norm, XConedown_rebin_norm};
  vector<TString> save_path_v      = {save_path_JEC, save_path_xcone};

  // #################################################################################################
  // Empty Data Bins #################################################################################

  vector<int> empty_bins_v;
  for(int bin=0; bin < number_bins; bin++){
    if(!(abs(data_rebin_norm->GetBinContent(bin+1))>0)){
      empty_bins_v.push_back(bin+1);
      if(debug) cout << "Empty Bins: " << empty_bins_v[bin] << "\n";
    }
  }
  cout << "\n";

  // #################################################################################################
  // Chi2 declarations ###############################################################################

  TString str_chi2_function_jec, str_chi2_function_xcone;
  vector<double> parameters_all;
  int NParams = 4;
  int number_empty_bins = empty_bins_v.size();

  for(int correction=0; correction<2; correction++){
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

    TString save_path    = save_path_v[correction];
    TH1F* CorrectionUp   = corrections_up_v[correction];
    TH1F* CorrectionDown = corrections_down_v[correction];

    TLegend *leg;
    if(!wmass_peak){
      CorrectionUp->SetLineWidth(2);
      CorrectionUp->SetLineColor(kBlue);

      CorrectionDown->SetLineStyle(2);
      CorrectionDown->SetLineColor(kBlue);

      TCanvas *A = new TCanvas("A", "A", 600, 600);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      ttbar_rebin_norm->Draw("HIST");
      CorrectionDown->Draw("SAME HIST");
      CorrectionUp->Draw("SAME HIST");
      data_rebin_norm->Draw("SAME P");
      leg = new TLegend(0.6,0.65,0.8,0.85);
      leg->SetTextSize(0.3);
      leg->AddEntry(ttbar_rebin_norm,"nominal","l");
      if(correction==0){
        leg->AddEntry(CorrectionUp,"JECup","l");
        leg->AddEntry(CorrectionDown,"JECdown","l");
      }
      if(correction==1){
        leg->AddEntry(CorrectionUp,"XConeup","l");
        leg->AddEntry(CorrectionDown,"XConedown","l");
      }
      leg->AddEntry(data_rebin_norm,"Data","p");
      leg->SetTextSize(0.05);
      leg->Draw();
      gPad->RedrawAxis();
      A->SaveAs(save_path+"/Wjet_mass_sensitivity_JEC.pdf");
      delete A;
      leg->Clear();
    }

    /*
    ███████ ██ ████████
    ██      ██    ██
    █████   ██    ██
    ██      ██    ██
    ██      ██    ██
    */
    /*
    Looking at each bin seperatly. Getting the bin content of the nominal file and
    the variations and making a fit through the values. The variation JECup and JECdown
    assigned to the values +1 and -1.
    Afterward the fit parameters will be stored in vectors for later purpose.
    */

    if(debug) print_seperater();
    if(debug) cout << "Start: Dependency in WJet bins" << endl;

    // #################################################################################################
    // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
    if(debug) cout << "Title for bins" << endl;
    vector<TString> mass_bin;
    if(debug && wmass_peak) cout << number_bins-first_bins << '\n';
    for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin*rebin_hist)+" < m_{Wjet} < "+to_string((bin+1)*rebin_hist));

    if(debug){
      cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << rebin_hist << "\n\n";
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

    vector<double> jec_factor, events;
    vector<double> jec_factor_err, events_err;
    vector<double> events_rel_err;

    int n_factors = 3; // JECdown, ttbar, JECup
    jec_factor = {-1, 0, 1};
    jec_factor_err = {0, 0, 0};

    vector<vector<double>> ordered_sys_norm_err;
    vector<TH1F*> ordered_sys_norm = {CorrectionDown, ttbar_rebin_norm, CorrectionUp};
    vector<TH1F*> ordered_sys = {JECdown_rebin, ttbar_rebin, JECup_rebin}; // for norm uncertainty
    /* ordered_sys_norm_err includes ttbar_bin_err, JECup_bin_err and JECdown_bin_err */
    for(unsigned int factor=0; factor<n_factors; factor++) ordered_sys_norm_err.push_back(normalize_error(ordered_sys[factor]));

    // for the chi2 Method -----------------------------------------------------------------------------
    if(debug) cout << "Vectors for Chi2" << endl;
    vector<double> data_bin_content;
    for(int bin=1; bin < number_bins+1; bin++) data_bin_content.push_back(data_rebin_norm->GetBinContent(bin));
    vector<double> data_bin_error     =normalize_error(data_rebin);;
    vector<double> XConeup_bin_error  =normalize_error(XConeup_rebin);
    vector<double> XConedown_bin_error=normalize_error(XConedown_rebin);

    // #################################################################################################
    // Cout Bin Error ##################################################################################
    if(debug){
      cout << "" << '\n' << fixed;
      cout << "|======================================================================================|" << endl;
      cout << "| ------------------------------------ Bin Errors -------------------------------------|" << endl;
      cout << "| -------------------------------------------------------------------------------------|" << endl;
      cout << "|   Bin  |    Data    |  JEC down  |nominal(cen)|   JEC up   |  Average   |  cen/avg   |" << endl;
      cout << "| -------|------------|------------|------------|------------|------------|------------|" << endl;
      /*
      In this loop no additional treatment is  necessary. The +first_bins includes the correct bins
      and the vectors only consists of number_bins entries.
      */
      for(unsigned int i=0; i<number_bins; i++){
        TString bin, data_str, Jd_str, tt_str, Ju_str, avg_str, cenavg_str;
        double avg = (ordered_sys_norm_err[0][i]+ordered_sys_norm_err[1][i]+ordered_sys_norm_err[2][i])/3;
        double cenavg;
        if(abs(avg) > 0) cenavg = ordered_sys_norm_err[1][i]/avg;                 // make sure bin error is not 0
        else             cenavg = avg;
        if(i+1<10)     bin = "|    "+to_string(i+1)+"   |";                       // i+1 to get bin number
        if(i+1>=10)    bin = "|   "+to_string(i+1)+"   |";                        // +first_bins to skip first bins
        data_str   = "  "+to_string(data_bin_content[i])+"  |";
        Jd_str     = "  "+to_string(ordered_sys_norm_err[0][i])+"  |";
        tt_str     = "  "+to_string(ordered_sys_norm_err[1][i])+"  |";
        Ju_str     = "  "+to_string(ordered_sys_norm_err[2][i])+"  |";
        avg_str    = "  "+to_string(avg)+"  |";
        cenavg_str = "  "+to_string(cenavg)+"  |";

        TString table_row = bin+data_str+Jd_str+tt_str+Ju_str+avg_str+cenavg_str;
        cout << table_row << '\n';
      }
      cout << "| -------------------------------------------------------------------------------------|" << endl;
      cout << "" << endl;
    }

    // #################################################################################################
    // Bin Error one sigma estimation ##################################################################
    vector<double> bin_error_one_sigma_avg, bin_error_one_sigma_central;
    for(int bin=0; bin<number_bins; bin++){
      bin_error_one_sigma_central.push_back(ordered_sys_norm_err[1][bin]);
      bin_error_one_sigma_avg.push_back((ordered_sys_norm_err[0][bin]+ordered_sys_norm_err[1][bin]+ordered_sys_norm_err[2][bin])/3);
    }

    // #################################################################################################
    // Canvas for all bins in one ######################################################################
    if(debug) cout << "Star: filling axis" << endl;

    TCanvas *E = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    if(rebin_hist==4) E->Divide(5,9);
    else if(rebin_hist==5) E->Divide(4,9);
    else if(rebin_hist==10) E->Divide(3,6);
    else cout << "Not the correct bin number - " << number_bins << '\n';
    E->SetCanvasSize(800, 1200);
    E->SetWindowSize(800, 1200);

    // #################################################################################################
    // set first fit function ##########################################################################
    TF1 *fit = new TF1("fit", "[0] + [1]*x");
    vector<TF1*> fit_functions, fit_functions_norm;
    vector<double> a, a_err, b, b_err, chi2;

    /*
    ██       ██████   ██████  ██████
    ██      ██    ██ ██    ██ ██   ██
    ██      ██    ██ ██    ██ ██████
    ██      ██    ██ ██    ██ ██
    ███████  ██████   ██████  ██
    */

    for(int bin=0; bin < number_bins; bin++){
      /*
      For the mass peak, the first bins are excluded, since there is no comparison necessary.
      This means that the bin iterator in the loop needs to be added with the first bins.
      */
      if(print) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()){
        a.push_back(0);
        a_err.push_back(0);
        b.push_back(0);
        b_err.push_back(0);
        chi2.push_back(0);
        continue;
      }

      // Getting BinContent & Error for each Factor for this bin ---------------------
      for(int factor=0; factor<n_factors; factor++){
        events.push_back(ordered_sys_norm[factor]->GetBinContent(bin+1)); // bin+1, since bin number starts at 1
        events_err.push_back(ordered_sys_norm_err[factor][bin]); // bin, since vector starts at 0
      }

      // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
      // #################################################################################################
      // Fit #############################################################################################
      if(debug) cout << "Start: TGraphError - Single - bin " << bin+1 << endl;
      TGraphErrors *single_bin;
      // TCanvas *A = new TCanvas("A","A", 600, 600);
      single_bin = new TGraphErrors(n_factors, &jec_factor[0], &events[0], &jec_factor_err[0], &events_err[0]);
      single_bin->Fit("fit", "Q"); // set first fit function (search for it)

      // #################################################################################################
      // All in One ######################################################################################
      E->cd(bin+1);
      if(debug) cout << "Start: TGraphError - all in one - bin " << bin+1 << endl;

      single_bin->SetTitle(mass_bin[bin]);
      single_bin->SetMarkerStyle(2);
      single_bin->GetXaxis()->SetTitle("m_{Wjet}");
      single_bin->Draw("AP");

      events = {}; events_err = {}; // clear vectors

      /*
      Get values of the parameter of the fit
      Values du not depend on scale of the xAxis, but only on the scale of the yAxis ((not-)norm)
      Therefore the vectors for the Value are filled once
      */
      TF1 *fit = single_bin->GetFunction("fit");
      a.push_back(fit->GetParameter(0));
      a_err.push_back(fit->GetParError(0));
      b.push_back(fit->GetParameter(1));
      b_err.push_back(fit->GetParError(1));
      chi2.push_back(fit->GetChisquare());
      fit_functions_norm.push_back(fit);
      cout << "" << endl;
    }

    E->SaveAs(save_path+"/mass_bin_all.pdf");

    if(print){
      TString bin, a_str, a_err_str, b_str, b_err_str, chi2_str;
      cout << "" << '\n' << fixed;
      cout << "|=================================================================|" << endl;
      cout << "| ------------------------- Fit Values ---------------------------|" << endl;
      cout << "| ----------------------------------------------------------------|" << endl;
      cout << "|   Bin  |    a     |  a_err   |     b     |   b_err  |   chi2    |" << endl;
      cout << "| -------|----------|----------|-----------|----------|-----------|" << endl;
      /*
      In this loop no additional treatment is  necessary. The +first_bins includes the correct bins
      and the vectors only consists of number_bins entries.
      */

      for(unsigned int i=0; i<number_bins; i++){
        /* The parameters are stored in vetors. To get the correct iterator, the first bin
        needs to be counted +1 for each bin with data bin content = 0. The iterator then needs to be subtracted
        by first_bins.*/
        if(find(empty_bins_v.begin(), empty_bins_v.end(), i+1) != empty_bins_v.end()) continue; // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)

        if(i+1<10)     bin = "|     "+to_string(i+1)+"  |"; // i+1 to get bin number
        if(i+1>=10)    bin = "|    "+to_string(i+1)+"  |";  // +first_bins to skip first bins
        if(chi2[i]<10) chi2_str = "  "+to_string(chi2[i])+" |";
        if(chi2[i]>10) chi2_str = " "+to_string(chi2[i])+" |";
        a_str = " "+to_string(a[i])+" |";
        a_err_str = " "+to_string(a_err[i])+" |";
        if(b[i]>0)     b_str = "  "+to_string(b[i])+" |";
        if(b[i]<0)     b_str = " "+to_string(b[i])+" |";
        b_err_str = " "+to_string(b_err[i])+" |";

        TString table_row = bin+a_str+a_err_str+b_str+b_err_str+chi2_str;
        cout << table_row << '\n';
      }
      cout << "| ----------------------------------------------------------------|" << endl;
      cout << "" << endl;
    }

    /*
    .██████ ██   ██ ██ ██████
    ██      ██   ██ ██      ██
    ██      ███████ ██  █████
    ██      ██   ██ ██ ██
    .██████ ██   ██ ██ ███████
    */
    if(debug) cout << "Chi2 - Parameters\n";

    int parameter = 0;
    TString str_chi2_function = ""; TString str_parameter;
    while(parameter<NParams*(number_bins-number_empty_bins)){ // first bins does not include Data
      // TString chi2_formula = "(([0]-[1]-[2]*x)/[3])^2";
      str_chi2_function += "+(([" + to_string(parameter) + "]";
      parameter++;
      str_chi2_function += "-[" + to_string(parameter) + "]";
      parameter++;
      str_chi2_function += "-[" + to_string(parameter) + "]*x)";
      parameter++;
      str_chi2_function += "/[" + to_string(parameter) + "])^2";
      parameter++;
    }
    if(debug) cout << str_chi2_function << endl;

    TF1 *chi2_function = new TF1("chi2_function", str_chi2_function, -10, 10);

    // if(correction==0) str_chi2_function_jec   = str_chi2_function;
    // else{
    //   str_chi2_function_xcone = str_chi2_function;
    //   str_chi2_function_xcone.ReplaceAll('x', 'y');
    //   parameter = 0;
    //   while(parameter<NParams*(number_bins-number_empty_bins)){
    //     str_chi2_function_xcone.ReplaceAll("[" + to_string(parameter) + "]", "[" + to_string(parameter+NParams*(number_bins-number_empty_bins)) + "]");
    //     parameter++;
    //   }
    // }

    // #################################################################################################
    // get parameters for terms ########################################################################
    vector<double> parameters;
    for(int bin=0; bin<number_bins; bin++){
      double data_bin_error_new = data_bin_error[bin];
      if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue;
      /* bin content vector contains every bin, also the first bins.
      Therefor "bin+first_bins". The parameter vectors only consists of the used bins for the fit,
      excluding the first_bins. Therefor the vetor is number_bins-first_bins long and "bin" is used */
      parameters.push_back(data_bin_content[bin]);
      parameters.push_back(a[bin]);
      parameters.push_back(b[bin]);
      parameters.push_back(data_bin_error_new+bin_error_one_sigma_central[bin]);

      // parameters_all.push_back(data_bin_content[bin]);
      // parameters_all.push_back(a[bin]);
      // parameters_all.push_back(b[bin]);
      // parameters_all.push_back(data_bin_error_new+bin_error_one_sigma_central[bin]);
    }

    // #################################################################################################
    // Set all parameters ##############################################################################
    for(int ipar=0; ipar<NParams*(number_bins-number_empty_bins); ipar++) chi2_function->SetParameter(ipar, parameters[ipar]);
    chi2_function->SetTitle("");

    // #################################################################################################
    // Get Minimum of chi2 #############################################################################
    if(debug) cout << "\nChi2 - Minimum" << endl;
    double minY        = chi2_function->GetMinimum();
    double minX        = chi2_function->GetX(minY, -10, 10);
    double X_sigmaup   = chi2_function->GetX(minY+1, minX, 10);
    double X_sigmadown = chi2_function->GetX(minY+1, -10, minX);
    double sigmaup     = X_sigmaup - minX; double sigmadown = minX - X_sigmadown;

    // #################################################################################################
    // Print and save Factor ###########################################################################
    cout << "\u03c7          = " << minY << endl;
    cout << "JEC factor =  " << minX << " \u00b1" << sigmaup << endl;
    fstream jec_txt;
    jec_txt.open(save_path+"/jec_factor.txt", ios::out);
    jec_txt << number_bins << " bins [0, 180] GeV\n";
    jec_txt << minX << /*" \u00b1"*/ "  err: " << sigmaup << endl;
    jec_txt.close();
    print_seperater();

    // #################################################################################################
    // unified x axis interval  ########################################################################
    if(debug) cout << "Chi2 - shifted Function" << endl;
    TF1 *chi2_function_shifted = new TF1("chi2_function_shifted", str_chi2_function, minX-3, minX+3);

    // Set all parameters ------------------------------------------------------------------------------
    for(int ipar=0; ipar<NParams*(number_bins-number_empty_bins); ipar++) chi2_function_shifted->SetParameter(ipar, parameters[ipar]);
    chi2_function_shifted->SetTitle("#chi^{2} of "+year);

    TCanvas *B = new TCanvas("B","B", 600, 600);
    chi2_function_shifted->Draw();
    B->SaveAs(save_path+"/chi2_shifted.pdf");
    B->Clear();

    /*
    .██████  ██████  ██████  ██████  ███████ ████████ ███████ ██████      ██   ██ ██ ███████ ████████
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██ ██         ██
    ██      ██    ██ ██████  ██████  █████      ██    █████   ██   ██     ███████ ██ ███████    ██
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██      ██    ██
    .██████  ██████  ██   ██ ██   ██ ███████    ██    ███████ ██████      ██   ██ ██ ███████    ██
    */
    if(debug) cout << "\nNew Histogram with correct Minimum from Chi2" << endl;
    cout << '\n';

    vector<double> JEC_min, JEC_min_up_err, JEC_min_down_err;
    for(unsigned int bin=0; bin<number_bins; bin++){
      /* No dummy necessary, since a and b are already zero*/
      double a_par = a[bin];
      double b_par = b[bin];
      double sigmaup_value = a_par + b_par * X_sigmaup; //
      double sigmadown_value = a_par + b_par * X_sigmadown; //
      double nominal_value = a_par + b_par * minX;
      JEC_min.push_back(nominal_value);
      if(sigmaup_value<nominal_value && nominal_value<sigmadown_value){
        if(debug) cout << "Bin " << bin << ": up(down) is down(up)" << endl;
        JEC_min_up_err.push_back(abs(JEC_min[bin]-sigmadown_value)); // Here and next line: (JEC_min_up & _down)
        JEC_min_down_err.push_back(abs(JEC_min[bin]-sigmaup_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else if(sigmadown_value<nominal_value && nominal_value<sigmaup_value){
        if(debug) cout << "Bin " << bin << ": up(down) is up(down)" << endl;
        JEC_min_up_err.push_back(abs(JEC_min[bin]-sigmaup_value)); // Here and next line: (JEC_min_up & _down)
        JEC_min_down_err.push_back(abs(JEC_min[bin]-sigmadown_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else if(sigmaup_value==0 && sigmadown_value==0){
        JEC_min_up_err.push_back(0);
        JEC_min_down_err.push_back(0);
      }
      else throw runtime_error("Something is wrong with the bin value for the FSR variation using the calculated minimum. both variations are smaller or greater than the nominal value");
    }
    if(debug) cout << "\n";

    TH1F *h_JEC_min;
    if(wmass_peak){
      h_JEC_min = new TH1F("JEC_min", "", number_bins-first_bins, lower_limit, upper_limit);
      for(unsigned int bin=0; bin<number_bins-first_bins;bin++) h_JEC_min->SetBinContent(bin+1, JEC_min[bin]);
    }
    else{
      h_JEC_min = new TH1F("JEC_min", "", number_bins, 0, 180);
      for(unsigned int bin=0; bin<number_bins;bin++) h_JEC_min->SetBinContent(bin+1, JEC_min[bin]);
    }

    // #################################################################################################
    // Settings ########################################################################################
    double JEC_min_integral = h_JEC_min->Integral();
    cout << "\nIntegral of calculated JEC_min distribution: " << JEC_min_integral << endl;

    h_JEC_min->SetTitle(year);
    if(wmass_peak) h_JEC_min->GetXaxis()->SetRangeUser(lower_limit, upper_limit);
    else           h_JEC_min->GetXaxis()->SetRangeUser(0, 180);
    h_JEC_min->GetYaxis()->SetRangeUser(0, h_JEC_min->GetMaximum()*1.2);
    h_JEC_min->GetXaxis()->SetNdivisions(505);
    h_JEC_min->GetYaxis()->SetNdivisions(505);
    h_JEC_min->GetXaxis()->SetTitleSize(0.05);
    h_JEC_min->GetYaxis()->SetTitleSize(0.04);
    h_JEC_min->GetXaxis()->SetTitleOffset(0.9);
    h_JEC_min->GetYaxis()->SetTitleOffset(1.5);
    h_JEC_min->GetXaxis()->SetTitle("w_{Wjet}");
    h_JEC_min->GetYaxis()->SetTitle("");
    h_JEC_min->SetLineWidth(1);
    h_JEC_min->SetLineColor(kBlack);

    CorrectionDown->SetLineStyle(2);

    vector<double> xbins, xbins_err; // = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
    //vector<double> xbins_err = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    int bin_steps = rebin_hist;
    int bin_position = bin_steps/2;
    if(wmass_peak) number_bins -= first_bins;
    for(unsigned int bin=0; bin<number_bins; bin++){ // center of the bins
      xbins.push_back(bin_position+bin_steps*bin);
      xbins_err.push_back(bin_position);
    }
    if(debug) cout << "Number Bins of corrected Histogram: " << xbins.size() << '\n';

    auto errors = new TGraphAsymmErrors(number_bins, &xbins[0], &JEC_min[0], &xbins_err[0], &xbins_err[0], &JEC_min_down_err[0], &JEC_min_up_err[0]);
    errors->SetFillColor(kGray);
    errors->SetFillStyle(1001); // 1001: solid - same as nothing
    errors->SetTitle(year);

    // #################################################################################################
    // Draw Hist #######################################################################################
    if(debug) cout << "Draw New Histogram with correct Minimum from Chi2" << endl;
    TCanvas *C = new TCanvas("C", "C", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    h_JEC_min->Draw("HIST");
    errors->Draw("SAME E2");
    h_JEC_min->Draw("SAME HIST"); // Draw again, because it will be overdrawn otherwise
    ttbar_rebin_norm->Draw("SAME HIST");
    data_rebin_norm->Draw("SAME P");

    leg = new TLegend(0.20,0.65,0.40,0.85);
    leg->AddEntry(ttbar_rebin_norm,"nominal","l");
    leg->AddEntry(h_JEC_min, "JEC best fit value","l");
    leg->AddEntry(data_rebin_norm, "data","ple"); // point + error
    leg->SetTextSize(0.02);
    leg->Draw();
    gPad->RedrawAxis();
    C->SaveAs(save_path+"/JEC_minimum_factor.pdf");
    CorrectionUp->Draw("SAME HIST");
    CorrectionDown->Draw("SAME HIST");
    leg->AddEntry(CorrectionUp, "JECup","l");
    leg->AddEntry(CorrectionDown, "JECdown","l");
    C->SaveAs(save_path+"/JEC_minimum_factor_compare.pdf");
    leg->Clear();
    delete C;
  }

  // /*
  // .██████  ██████       ██████ ██   ██ ██ ██████
  // .     ██ ██   ██     ██      ██   ██ ██      ██
  // . █████  ██   ██     ██      ███████ ██  █████
  // .██      ██   ██     ██      ██   ██ ██ ██
  // .███████ ██████       ██████ ██   ██ ██ ███████
  // */
  //
  // /*
  // The used method here is not correct, but kept as a template.
  // In this case the two corrections only depend on x or y. In addition, the data
  // is considered twice, which is not correct.
  // */
  //
  // if(debug) cout<< "\n2D Chi2\n";
  // // #################################################################################################
  // // JEC and XCone Chi2 in one #######################################################################
  // TF2 *chi2_function_combined = new TF2("chi2_combined", str_chi2_function_jec+str_chi2_function_xcone, -3, 3, -3, 3);
  // if(debug) cout << str_chi2_function_jec+str_chi2_function_xcone << endl;
  //
  // // #################################################################################################
  // // Set all paramters ###############################################################################
  // if(debug) cout<< "\n2D Chi2 - parameters\n";
  // for(int ipar=0; ipar<2*NParams*(number_bins-number_empty_bins); ipar++) chi2_function_combined->SetParameter(ipar, parameters_all[ipar]);
  // chi2_function_combined->SetTitle("");
  // chi2_function_combined->GetXaxis()->SetTitle("JEC");
  // chi2_function_combined->GetYaxis()->SetTitle("XCone");
  //
  // // #################################################################################################
  // // Minimum #########################################################################################
  // /*
  // Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447.
  // */
  // double twoD_minZ_alt = chi2_function_combined->GetMinimum(); // y=0 -> Found out numerically
  // if(debug) cout << "z_min (->GetMinimum() - y=0):   " << twoD_minZ_alt << endl;
  // double twoD_minX, twoD_minY;
  // double twoD_minZ = chi2_function_combined->GetMinimumXY(twoD_minX,twoD_minY);
  // cout  << "\nz_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n' << endl;
  // int number_para = chi2_function_combined->GetNpar();
  //
  // fstream jec_txt;
  // jec_txt.open(save_path_general+"/jec_factor.txt", ios::out);
  // jec_txt << number_bins << " bins [0, 180] GeV\n\n";
  // jec_txt << "z_min: " << twoD_minZ << "\nx    :  " << twoD_minX << "\ny    :  " << twoD_minY << endl;
  // jec_txt.close();
  //
  // // #################################################################################################
  // // Evaluate 1sigma #################################################################################
  //
  // // Draw Points -------------------------------------------------------------------------------------
  // vector<TPoint*> all_points_draw;
  // TGraph *zmin_point = new TGraph();
  // vector<vector<double>> points = FindXY(chi2_function_combined, twoD_minZ+1, 0, 2, 0, 2, 5000, 0.001);
  // if(print) cout << "Number Points at 1sigma: " << points.size() << endl;
  // // for(unsigned int i=0; i<points.size(); i++) cout << "x: " << points[i][0] << " | y: " << points[i][1] << " | z: " << points[i][2] << endl;
  //
  // zmin_point->SetPoint(0, twoD_minX, twoD_minY);
  // zmin_point->SetMarkerColor(kGray);
  // zmin_point->SetMarkerStyle(kFullCircle);
  // zmin_point->SetMarkerSize(0.5);
  // TGraph *sigma_points = new TGraph();
  // for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1]);
  // sigma_points->SetMarkerColor(kRed);
  // sigma_points->SetMarkerStyle(kFullCircle);
  // sigma_points->SetMarkerSize(0.1);
  // TMultiGraph *all_points = new TMultiGraph();
  // all_points->Add(zmin_point);
  // all_points->Add(sigma_points);
  //
  // // Draw Ellipse ------------------------------------------------------------------------------------
  // double r_sigma_single; double rx_sigma=0; double ry_sigma=0;
  // int count_rx=0; int count_ry=0;
  // for(unsigned int i=0; i<points.size(); i++){
  //   /* 1sigma region is a ellipse and not a circle. To get the parameters for the
  //   ellipse, only the extrempoints (closest and furthest points) are considerd.
  //   Therefor, an intervall is estimated and the average of the points within the
  //   intervall are taken as the radius.*/
  //   if(points[i][0]<twoD_minX+0.1 && points[i][0]>twoD_minX-0.1){
  //     count_rx++;
  //     rx_sigma+=sqrt(pow((points[i][0]-twoD_minX), 2)+pow((points[i][1]-twoD_minY), 2));
  //   }
  //   if(points[i][1]<twoD_minY+0.1 && points[i][1]>twoD_minY-0.1){
  //     count_ry++;
  //     ry_sigma+=sqrt(pow((points[i][0]-twoD_minX), 2)+pow((points[i][1]-twoD_minY), 2));
  //   }
  // }
  //
  // rx_sigma = rx_sigma/count_rx++;
  // ry_sigma = ry_sigma/count_ry++;
  // if(print) cout << "rx_sigma: " << rx_sigma << "(Number Points:" << count_rx << ")\n";
  // if(print) cout << "ry_sigma: " << ry_sigma << "(Number Points:" << count_ry << ")\n";
  // TEllipse *sigma_area = new TEllipse(twoD_minX, twoD_minY, ry_sigma, rx_sigma);
  // sigma_area->SetLineWidth(1);
  // sigma_area->SetLineColor(kRed);
  // sigma_area->SetFillStyle(0);
  //
  // // #################################################################################################
  // // Plot ############################################################################################
  // if(debug) cout<< "\n2D Chi2 - plot\n";
  // TCanvas *D = new TCanvas("D","D", 600, 600);
  // chi2_function_combined->Draw("COLZ");
  // D->SaveAs(save_path_general+"/chi2_combined_colz.pdf");
  // sigma_points->Draw("SAME P");
  // zmin_point->Draw("SAME P");
  // D->SaveAs(save_path_general+"/chi2_combined_points_colz.pdf");
  // D->Clear();
  //
  // chi2_function_combined->Draw("COLZ");
  // sigma_area->Draw("SAME");
  // zmin_point->Draw("SAME P");
  // D->SaveAs(save_path_general+"/chi2_combined_ellipse_colz.pdf");
}
