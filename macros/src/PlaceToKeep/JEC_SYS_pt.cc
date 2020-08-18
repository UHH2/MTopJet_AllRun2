#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
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

  bool debug = true;
  if(argc != 3){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number>\n" << "rebin_number: hist->Rebin(rebin_number)\n";
    return 0;
  }


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

  /*
  The reason why the first empty bins are included is, because in later histograms
  the empty bins are shown for a better comparison of the results. They could be
  excluded, but it is not done yet.
  */

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

  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"+year+"/btag/";

  // rmdir(save_path+"/rebin_"+str_number_bins);
  if(mkdir(save_path+"/rebin_"+str_number_bins,0777) != 0) cout << "couldn't creat rebin folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+" created\n";
  if(mkdir(save_path+"/rebin_"+str_number_bins+"/single_bins",0777) != 0) cout << "couldn't creat single_bins folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+"/single_bins created\n";

  if(mkdir(save_path+"/rebin_"+str_number_bins+"/pt",0777) != 0) cout << "couldn't creat rebin/pt folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+"/pt created\n";
  if(mkdir(save_path+"/rebin_"+str_number_bins+"/pt/single_bins",0777) != 0) cout << "couldn't creat pt/single_bins folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+"/pt/single_bins created\n";
  save_path += "/rebin_"+str_number_bins+"/pt"; // New save_path

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */
  fstream jec_txt;
  jec_txt.open(save_path+"/jec_factor_"+str_number_bins+"_all.txt", ios::out);

  vector<TF1*>   chi2_functions, chi2_functions_grounded; // To show all Chi2 functions in one plot later
  vector<double> jec_factors, jec_errors;

  TString hist_class = "comparison_topjet_xcone_pass_rec/";
  TString w_mass_low = "wmass_match_ptbin_low";
  TString w_mass_midlow = "wmass_match_ptbin_midlow";
  TString w_mass_midhigh = "wmass_match_ptbin_midhigh";
  TString w_mass_high = "wmass_match_ptbin_high";
  vector<TString> pt_range = {"    p_{T}<200", "200<p_{T}<300", "300<p_{T}<400", "400<p_{T}"};
  vector<TString> pt_bins = {"low", "midlow", "midhigh", "high"};
  vector<TString> w_mass = {hist_class+w_mass_low, hist_class+w_mass_midlow, hist_class+w_mass_midhigh, hist_class+w_mass_high};

  for(unsigned int pt_bin=0; pt_bin<pt_bins.size(); pt_bin++){
    cout << "################################################################### "+pt_bins[pt_bin]+"\n";
    cout << "###################################################################\n";
    cout << "###################################################################\n";
    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    vector<TFile*> file_bkg_v;
    vector<TH1F*> hists_bkg_v;
    vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
    for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+year+"/muon/"+path_bkg_v[i]));
    for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass[pt_bin]));
    TH1F *bkg = AddHists(hists_bkg_v, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    TFile *data_file      = new TFile(dir+year+"/muon/"+data_path);

    TH1F *data            = (TH1F*)data_file->Get(w_mass[pt_bin]);
    TH1F *data_rebin      = rebin(data, rebin_hist);
    TH1F *data_norm       = normalize(data);
    TH1F *data_rebin_norm = normalize(data_rebin);

    // Empty Bins --------------------------------------------------------------------------------------
    if(debug){
      for(int bin=0; bin<number_bins; bin++){
        if(data_rebin_norm->GetBinContent(bin+1)==0) cout << "Data bin "+to_string(bin)+" is empty\n";
      }
    }

    // #################################################################################################
    // Get TTbar #######################################################################################
    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
    TFile  *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);

    if(debug) cout << "TTbar low" << endl;
    TH1F*ttbar  = (TH1F*)ttbar_file->Get(w_mass[pt_bin]);
    ttbar->Add(bkg, 1);
    TH1F* ttbar_norm = normalize(ttbar);
    TH1F* ttbar_rebin = rebin(ttbar, rebin_hist);
    TH1F* ttbar_rebin_norm = rebin(ttbar_norm, rebin_hist);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC" << endl;

    TFile *JECup_file   = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
    TFile *JECdown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");

    TH1F *JECup   = (TH1F*)JECup_file->Get(w_mass[pt_bin]);
    TH1F *JECdown = (TH1F*)JECdown_file->Get(w_mass[pt_bin]);

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
    if(debug) cout << "JEC+Background - " << JECdown_rebin_norm->GetEntries() << endl;
    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone" << endl;

    TFile *XConeup_file   = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConeup.root");
    TFile *XConedown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConedown.root");

    TH1F *XConeup   = (TH1F*)JECup_file->Get(w_mass[pt_bin]);
    TH1F *XConedown = (TH1F*)JECdown_file->Get(w_mass[pt_bin]);

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

    TLegend *leg;
    ttbar_rebin_norm->SetTitle("");
    ttbar_rebin_norm->GetXaxis()->SetRangeUser(0, 180);
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

    data_rebin_norm->SetMarkerStyle(8);  // data hist style
    data_rebin_norm->SetMarkerColor(kBlack);

    JECup_rebin_norm->SetLineWidth(2);
    JECup_rebin_norm->SetLineColor(kBlue);

    JECdown_rebin_norm->SetLineStyle(2);
    JECdown_rebin_norm->SetLineColor(kBlue);

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    ttbar_rebin_norm->Draw("HIST");
    JECdown_rebin_norm->Draw("SAME HIST");
    JECup_rebin_norm->Draw("SAME HIST");
    data_rebin_norm->Draw("SAME P");
    leg = new TLegend(0.6,0.65,0.8,0.85);
    leg->SetTextSize(0.3);
    leg->AddEntry(ttbar_rebin_norm,"nominal","l");
    leg->AddEntry(JECup_rebin_norm,"JECup","l");
    leg->AddEntry(JECdown_rebin_norm,"JECdown","l");
    leg->AddEntry(data_rebin_norm,"Data","p");
    leg->SetTextSize(0.05);
    leg->Draw();
    gPad->RedrawAxis();
    A->SaveAs(save_path+"/wmass_JEC_rebin"+str_number_bins+"_"+pt_bins[pt_bin]+".pdf");
    delete A;
    leg->Clear();

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

    if(debug) cout << "-----------------------------------------------------------" << endl;
    if(debug) cout << "Start: Dependency in WJet bins" << endl;

    // #################################################################################################
    // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
    if(debug) cout << "Title for bins" << endl;
    vector<TString> mass_bin;
    for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin*rebin_hist)+" < m_{Wjet} < "+to_string((bin+1)*rebin_hist));

    if(debug){
      cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << rebin_hist << '\n';
      cout << "number_bins " << number_bins-first_bins << " | number titles " << mass_bin.size()-first_bins << " | bin width " << rebin_hist << '\n';
      for(unsigned int bin=first_bins; bin<mass_bin.size(); bin++) cout << mass_bin[bin] << '\n';
    }

    vector<TString> number_bin;
    for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

    // #################################################################################################
    // Bin Content and error ###########################################################################
    if(debug) cout << "Vector for Content and Error" << endl;

    vector<double> jec_factor, events;
    vector<double> jec_factor_err, events_err;
    vector<double> events_rel_err;

    int n_factors = 3; // JECdown, ttbar, JECup
    jec_factor = {-1, 0, 1};
    jec_factor_err = {0, 0, 0};

    vector<vector<double>> ordered_sys_norm_err;
    vector<TH1F*> ordered_sys_norm = {JECdown_rebin_norm, ttbar_rebin_norm, JECup_rebin_norm};
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
      cout.precision(6);
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

    // Bin Error one sigma estimation --------------------------------------------------------------------
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
      cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      if(!(abs(data_rebin_norm->GetBinContent(bin+1))>0)) continue; // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)

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
      single_bin->Fit("fit"); // set first fit function (search for it)


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

    E->SaveAs(save_path+"/mass_bin_all_rebin"+str_number_bins+"_"+pt_bins[pt_bin]+".pdf");

    cout.precision(6);
    TString bin_str, a_str, a_err_str, b_str, b_err_str, chi2_str;
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
    int empty_bins=0;
    for(unsigned int bin=0; bin<number_bins; bin++){
      /* The parameters are stored in vetors. To get the correct iterator, the first bin
      needs to be counted +1 for each bin with data bin content = 0. The iterator then needs to be subtracted
      by first_bins.*/
      if(data_rebin_norm->GetBinContent(bin+1)==0){
        empty_bins = empty_bins+1;
        continue; // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      }

      if(bin+1<10)  bin_str = "|     "+to_string(bin+1)+"  |";
      if(bin+1>=10) bin_str = "|    "+to_string(bin+1)+"  |";
      if(chi2[bin-empty_bins]<10) chi2_str = "  "+to_string(chi2[bin-empty_bins])+" |";
      if(chi2[bin-empty_bins]>10) chi2_str = " "+to_string(chi2[bin-empty_bins])+" |";
      a_str = " "+to_string(a[bin-empty_bins])+" |";
      a_err_str = " "+to_string(a_err[bin-empty_bins])+" |";
      if(b[bin-empty_bins]>0) b_str = "  "+to_string(b[bin-empty_bins])+" |";
      if(b[bin-empty_bins]<0) b_str = " "+to_string(b[bin-empty_bins])+" |";
      b_err_str = " "+to_string(b_err[bin-empty_bins])+" |";

      TString table_row = bin_str+a_str+a_err_str+b_str+b_err_str+chi2_str;
      cout << table_row << '\n';
    }
    cout << "| ----------------------------------------------------------------|" << endl;
    cout << "" << endl;

    /*
    .██████ ██   ██ ██ ██████
    ██      ██   ██ ██      ██
    ██      ███████ ██  █████
    ██      ██   ██ ██ ██
    .██████ ██   ██ ██ ███████
    */
    if(debug) cout << "Chi2 - Parameters" << endl;

    int NParams = 4; int parameter = 0;
    TString str_chi2_function = ""; TString str_parameter;
    while(parameter<NParams*(number_bins-empty_bins)){
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

    // #################################################################################################
    // get parameters for terms ########################################################################
    vector<double> parameters;
    empty_bins = 0;
    for(int bin=0; bin<number_bins; bin++){
      if(data_rebin_norm->GetBinContent(bin+1)==0){
        empty_bins = empty_bins+1;
        continue; // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      }
      /* bin content vector contains every bin, also the first bins.
      Therefor "bin+first_bins". The parameter vectors only consists of the used bins for the fit,
      excluding the first_bins. Therefor the vetor is number_bins-first_bins long and "bin" is used */
      parameters.push_back(data_bin_content[bin]);
      parameters.push_back(a[bin-empty_bins]);
      parameters.push_back(b[bin-empty_bins]);
      parameters.push_back(data_bin_error[bin]+bin_error_one_sigma_central[bin]);
    }

    // #################################################################################################
    // Set all parameters ##############################################################################
    for(int ipar=0; ipar<NParams*(number_bins-empty_bins); ipar++) chi2_function->SetParameter(ipar, parameters[ipar]);
    chi2_function->SetTitle("");
    chi2_functions.push_back(chi2_function);

    // #################################################################################################
    // Get Minimum of chi2 #############################################################################
    if(debug) cout << "Chi2 - Minimum" << endl;
    double minY = chi2_function->GetMinimum();
    double minX = chi2_function->GetX(minY, -10, 10);
    double X_sigmaup = chi2_function->GetX(minY+1, minX, 10);
    double X_sigmadown = chi2_function->GetX(minY+1, -10, minX);

    double sigmaup = X_sigmaup - minX; double sigmadown = minX - X_sigmadown;

    // #################################################################################################
    // Print and save Factor ###########################################################################
    cout << "" << endl;
    cout << "JEC factor = " << minX << " \u00b1" << sigmaup << endl;
    jec_factors.push_back(minX);
    jec_errors.push_back(sigmaup);

    jec_txt << "----------------------------------------------------- "+pt_range[pt_bin]+"\n\n";
    jec_txt << number_bins << " bins [0, 180] GeV\n";
    jec_txt << minX << /*" \u00b1"*/ "  err: " << sigmaup << "\n\n";
    print_seperater();

    // #################################################################################################
    // #################################################################################################
    // Groundline for Chi2 #############################################################################
    if(debug) cout << "Chi2 - grounded Function" << endl;
    TF1 *chi2_function_grounded = new TF1("chi2_function_grounded", str_chi2_function+"-"+to_string(minY), minX-5, minX+5);

    // #################################################################################################
    // Print and save Factor ###########################################################################
    for(int ipar=0; ipar<NParams*(number_bins-empty_bins); ipar++) chi2_function_grounded->SetParameter(ipar, parameters[ipar]);
    chi2_function_grounded->SetTitle("");
    chi2_functions_grounded.push_back(chi2_function_grounded);

    TCanvas *B = new TCanvas("B","B", 600, 600);
    chi2_function_grounded->Draw();
    B->SaveAs(save_path+"/chi2_rebin_grounded_"+pt_bins[pt_bin]+".pdf");
  }

  jec_txt.close();

  /*
  .█████  ██      ██          ██ ███    ██      ██████  ███    ██ ███████
  ██   ██ ██      ██          ██ ████   ██     ██    ██ ████   ██ ██
  ███████ ██      ██          ██ ██ ██  ██     ██    ██ ██ ██  ██ █████
  ██   ██ ██      ██          ██ ██  ██ ██     ██    ██ ██  ██ ██ ██
  ██   ██ ███████ ███████     ██ ██   ████      ██████  ██   ████ ███████
  */
  if(debug) cout << "Plot all chi2 in one\n";

  TLegend *leg;
  vector<int> color = {kRed, kBlue, kGreen, kBlack};
  TCanvas *C = new TCanvas("C", "C", 600, 600);
  for(unsigned int function=0; function<chi2_functions.size(); function++) chi2_functions[function]->SetLineColor(color[function]);
  chi2_functions[0]->Draw();
  for(unsigned int function=1; function<chi2_functions.size(); function++) chi2_functions[function]->Draw("SAME");
  leg = new TLegend(0.4,0.50,0.7,0.85);
  leg->SetTextSize(0.01);
  leg->AddEntry((TObject*)0, "p_{T} in GeV", "");
  for(unsigned int function=0; function<chi2_functions.size(); function++){
    if(function==0) leg->AddEntry(chi2_functions[function], "    "+pt_range[function],"l");
    else            leg->AddEntry(chi2_functions[function], pt_range[function],"l");
  }
  leg->Draw();
  gPad->RedrawAxis();
  C->SaveAs(save_path+"/comparison_chi2_all_ptbins.pdf");
  delete C;
  leg->Clear();

  // #################################################################################################
  // Groundline for Chi2 functions ###################################################################
  if(debug) cout << "Plot all grounded chi2 in one\n";

  C = new TCanvas("C", "C", 600, 600);
  chi2_functions_grounded[0]->GetYaxis()->SetRangeUser(0, 10);
  chi2_functions_grounded[0]->GetXaxis()->SetRangeUser(-2.5, 3);
  for(unsigned int function=0; function<chi2_functions_grounded.size(); function++) chi2_functions_grounded[function]->SetLineColor(color[function]);

  chi2_functions_grounded[0]->Draw();
  for(unsigned int function=1; function<chi2_functions_grounded.size(); function++) chi2_functions_grounded[function]->Draw("SAME");
  leg = new TLegend(0.45,0.6,0.65,0.85);
  leg->SetTextSize(0.03);
  leg->AddEntry((TObject*)0, "p_{T} in GeV", "");
  for(unsigned int function=0; function<chi2_functions_grounded.size(); function++){
    if(function==0) leg->AddEntry(chi2_functions_grounded[function], "    "+pt_range[function],"l");
    else            leg->AddEntry(chi2_functions_grounded[function], pt_range[function],"l");
  }
  leg->Draw();
  gPad->RedrawAxis();
  C->SaveAs(save_path+"/comparison_chi2_all_ptbins_grounded.pdf");
  delete C;
  leg->Clear();

  // #################################################################################################
  // JEC factors in one ##############################################################################
  if(debug) cout << "Plot all jec factors and errors\n";

  vector<TH1F*> h_jec_factors;
  TH1F *h_JEC_pt;
  for(unsigned int bin=0; bin<pt_bins.size();bin++){
    // h_JEC_pt = new TH1F("JEC_factor_"+pt_bins[bin], "", 1, bin, bin+1);
    // h_JEC_pt->SetBinContent(1, jec_factors[bin]);
    // h_JEC_pt->SetBinError(1, jec_errors[bin]);

    h_JEC_pt = new TH1F("JEC_factor_"+pt_bins[bin], "", pt_bins.size(), 0, pt_bins.size());
    h_JEC_pt->SetBinContent(bin+1, jec_factors[bin]);
    h_JEC_pt->SetBinError(bin+1, jec_errors[bin]);
    h_jec_factors.push_back(h_JEC_pt);
  }

  h_jec_factors[0]->SetTitle("");
  h_jec_factors[0]->GetXaxis()->SetRangeUser(0, 4);
  h_jec_factors[0]->GetYaxis()->SetRangeUser(-1, 2);
  h_jec_factors[0]->GetXaxis()->SetNdivisions(505);
  h_jec_factors[0]->GetYaxis()->SetNdivisions(505);
  h_jec_factors[0]->GetXaxis()->SetTitleSize(0.05);
  h_jec_factors[0]->GetYaxis()->SetTitleSize(0.04);
  h_jec_factors[0]->GetXaxis()->SetTitleOffset(0.9);
  h_jec_factors[0]->GetYaxis()->SetTitleOffset(1.5);
  h_jec_factors[0]->GetXaxis()->SetTitle("");
  h_jec_factors[0]->GetYaxis()->SetTitle("");

  for(unsigned int bin=0; bin<pt_bins.size();bin++){
    h_jec_factors[bin]->SetMarkerStyle(8);  // data hist style
    h_jec_factors[bin]->SetMarkerColor(color[bin]);
  }

  if(debug) cout << "Draw all central JEC factors in one plot" << endl;
  C = new TCanvas("C", "C", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  h_jec_factors[0]->Draw("HIST P E1");
  for(unsigned int bin=1; bin<pt_bins.size();bin++) h_jec_factors[bin]->Draw("SAME HIST P E1"); // Draw again, because it will be overdrawn otherwise

  leg = new TLegend(0.70,0.65,0.90,0.85);
  leg->AddEntry(h_jec_factors[0], "    "+pt_range[0],"ple"); // point + error
  for(unsigned int bin=1; bin<pt_bins.size();bin++) leg->AddEntry(h_jec_factors[bin], pt_range[bin],"ple"); // point + error
  leg->SetTextSize(0.02);
  leg->Draw();
  gPad->RedrawAxis();
  C->SaveAs(save_path+"/central_JEC_factor_all.pdf");
  leg->Clear();
  delete C;

}
