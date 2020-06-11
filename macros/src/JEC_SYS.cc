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
  1) This script loops many times over bins. The bin number start at 1 (e.g. for the function GetBinContent)
     On the other hand, the bin content or titles are stored within a vector. The vector starts at 0.
     Mostly, the loops start at bin=0. Therefor, using the function GetBinContent the bin needs to be bin+1
  */
  // #################################################################################################
  // Declare different variables used in the code ####################################################

  bool debug = true;
  bool rel_error = false;
  bool wmass_peak = false;

  if(argc != 3){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number>\n" << "rebin_number: hist->Rebin(rebin_number)\n";
    return 0;
  }

  // Rebin -------------------------------------------------------------------------
  /*
  Right now every Data bin will be excluded, if the BinContent is zero. This leads to a problem,
  since in the first bin the bin Content is (e.g.) 000001000. Therefor, the rebinning has to be
  taken into account.
  Rebin  4: exclude bin 1-4
  Rebin  5: exclude bin 1-2
  Rebin 10: exclude bin 1
  */
  int rebin_hist = atoi(argv[2]); // atoi() c-string representing integer into int
  int number_bins = 200/rebin_hist;
  TString str_number_bins = to_string(number_bins);

  // exclude bins with DataBinContent = 0 -- HARDCODED
  // number_bins-20/rebin_hist -> number_bins for 200 GeV
  //                              20 to get 180 GeV
  //                              20/rebin_hist to get the number of bins to cut of
  //                              taking minus to the the upper limit
  int first_bins = 0; int last_bins = 0; // last_bins for cut at m=180 GeV
  if     (rebin_hist==4)  first_bins = 4;
  else if(rebin_hist==5)  first_bins = 2;
  else if(rebin_hist==10) first_bins = 1;
  else throw runtime_error("At this point, only the rebinning 4, 5 and 10 works");
  number_bins -= 20/rebin_hist; // Only looking ab to 180 GeV

  // Year --------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false;
  if     (strcmp(year, "2016")==0) is16 = true;
  else if(strcmp(year, "2017")==0) is17 = true;
  else if(strcmp(year, "2018")==0) is18 = true;
  else throw runtime_error("Give me the correct year please (2016, 2017 or 2018)");

  // creating subdirectories. Necessary, because different binning examined
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

  TString reconst = "btag_sel"; // compare btag btag_sel min_mass
  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"+year+"/"+reconst;

  rmdir(save_path+"/rebin_"+str_number_bins);
  if(mkdir(save_path+"/rebin_"+str_number_bins,0777) != 0) cout << "couldn't creat rebin folder!\n";
  if(mkdir(save_path+"/rebin_"+str_number_bins+"/single_bins",0777) != 0) cout << "couldn't creat single_bins folder!\n";
  save_path += "/rebin_"+str_number_bins; // New save_path

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
  if(reconst=="min_mass") w_mass = "XCone_cor_subjets/min_mass_Wjet";           // minimal mass from subjet combination
  if(reconst=="btag")     w_mass = hist_class+"wmass";                          // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
  if(reconst=="compare")  w_mass = hist_class+"wmass_comparison";               // mass closest to Wmass from subjet combination
  if(reconst=="btag_sel") w_mass = hist_class+"wmass_btagcut_discrepancy_many_jets"; // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
  //                                                                                    Selection: btag_high>0.7; ONE subjet with dr(ak4, subjet)<0.4; ONE high btag

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

  TString data_path = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file = new TFile(dir+year+"/muon/"+data_path);
  TH1F *data = (TH1F*)data_file->Get(w_mass);

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar = (TH1F*)ttbar_file->Get(w_mass);
  ttbar->Add(bkg, 1);

  // #################################################################################################
  // Get JEC #########################################################################################
  if(debug) cout << "JEC" << endl;

  TFile *JECup_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
  TFile *JECdown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");

  TH1F *JECup = (TH1F*)JECup_file->Get(w_mass);
  TH1F *JECdown = (TH1F*)JECdown_file->Get(w_mass);

  // Add hists from FSR and background ---------------------------------------------------------------
  if(debug) cout << "FSR+Background" << endl;

  JECup->Add(bkg, 1);
  JECdown->Add(bkg, 1);

  // #################################################################################################
  // Normalize Hists #################################################################################
  if(debug) cout << "Normalize Hists" << endl;

  TH1F* ttbar_norm = normalize(ttbar);
  TH1F* data_norm = normalize(data);
  TH1F* JECup_norm = normalize(JECup);
  TH1F* JECdown_norm = normalize(JECdown);

  // #################################################################################################
  // Rebin Hists #####################################################################################
  if(debug) cout << "Rebin Hists" << endl;

  TH1F* ttbar_rebin_norm = rebin(ttbar_norm, rebin_hist);
  TH1F* data_rebin_norm = rebin(data_norm, rebin_hist);
  TH1F* JECup_rebin_norm = rebin(JECup_norm, rebin_hist);
  TH1F* JECdown_rebin_norm = rebin(JECdown_norm, rebin_hist);

  // #################################################################################################
  // Hists into Vector ###############################################################################
  if(debug) cout << "Hists into vector" << endl;

  vector<TH1F*> ttbar_v = {ttbar_norm, ttbar_rebin_norm};
  vector<TH1F*> data_v = {data_norm, data_rebin_norm};
  vector<TH1F*> JECup_v = {JECup_norm, JECup_rebin_norm};
  vector<TH1F*> JECdown_v = {JECdown_norm, JECdown_rebin_norm};

  /*
  ██████  ██       ██████  ████████ ███████
  ██   ██ ██      ██    ██    ██    ██
  ██████  ██      ██    ██    ██    ███████
  ██      ██      ██    ██    ██         ██
  ██      ███████  ██████     ██    ███████
  */
  if(debug) cout << "Plots" << endl;

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
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
  leg->AddEntry(data_rebin_norm,"Data","l");
  leg->SetTextSize(0.05);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path+"/Wjet_mass_sensitivity_JEC_rebin"+str_number_bins+".pdf");
  delete A;
  leg->Clear();


  /*
  ███████ ██ ████████
  ██      ██    ██
  █████   ██    ██
  ██      ██    ██
  ██      ██    ██
  */
  if(debug) cout << "-----------------------------------------------------------" << endl;
  if(debug) cout << "Start: Dependency in WJet bins" << endl;

  // set first fit function --------------------------------------------------------
  TF1 *fit = new TF1("fit", "[0] + [1]*x");
  vector<TF1*> fit_functions, fit_functions_norm;
  vector<double> a, a_err, b, b_err, chi2;

  // vectors for values (Bin Content, Errors, etc) ---------------------------------
  vector<double> jec_factor, events;
  vector<double> data_bin_content, data_bin_error;
  vector<double> ttbar_bin_content;
  vector<double> jec_factor_err, events_err;
  vector<double> events_rel_err;

  vector<TH1F*> ordered_jec = {JECdown_rebin_norm, ttbar_rebin_norm, JECup_rebin_norm};

  jec_factor = {-1, 0, 1};
  jec_factor_err = {0, 0, 0};

  // get the mass bins title for plots with all fits on one page - e.g. "10 < m_{Wjet} < 20"
  vector<TString> mass_bin;
  for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin*rebin_hist)+" < m_{Wjet} < "+to_string((bin+1)*rebin_hist));
  if(debug){
    cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << rebin_hist << '\n';
    for(unsigned int bin=0; bin<mass_bin.size(); bin++) cout << mass_bin[bin] << '\n';
  }

  vector<TString> number_bin;
  for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

  int n_factors = 3; // JECdown, ttbar, JECup
  vector<vector<double>> ordered_jec_norm_err;
  vector<TH1F*> ordered_jec_normal = {rebin(JECdown, rebin_hist), rebin(ttbar, rebin_hist), rebin(JECup, rebin_hist)}; // for norm uncertainty
  for(unsigned int factor=0; factor<n_factors; factor++) ordered_jec_norm_err.push_back(normalize_error(ordered_jec_normal[factor]));

  // -------------------------------------------------------------------------------
  if(debug) cout << "Star: filling axis" << endl;

  TCanvas *E = new TCanvas(); // All Graphs in one
  if(number_bins==50) E->Divide(5,9);
  else if(number_bins==40) E->Divide(4,9);
  else if(number_bins==20) E->Divide(3,6);
  else cout << "Not the correct bin number - " << number_bins << '\n';
  E->SetCanvasSize(800, 1200);
  E->SetWindowSize(800, 1200);

  // #################################################################################################
  // #################################################################################################
  // Start for loop for fits #########################################################################
  for(unsigned int bin=0; bin < number_bins; bin++){
    cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
    if(bin<first_bins) continue; // Exclude bins which are empty (for Data) or greater than 180 GeV

    // for the chi2 Method ---------------------------------------------------------
    data_bin_content.push_back(data_rebin_norm->GetBinContent(bin+1)); // Explanation 1) (search for it)
    data_bin_error.push_back(data_rebin_norm->GetBinError(bin+1));
    ttbar_bin_content.push_back(ttbar_rebin_norm->GetBinContent(bin+1));

    // Getting BinContent & Error for each Factor for this bin ---------------------
    for(int factor=0; factor<n_factors; factor++){
      events.push_back(ordered_jec[factor]->GetBinContent(bin+1));
      events_err.push_back(ordered_jec[factor]->GetBinError(bin+1));
    }

    // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
    if(debug) cout << "Start: TGraphError - Single - bin " << bin << endl;
    TGraphErrors *single_bin;
    // each Graph seperatly---------------------------------------------------------
    TCanvas *A = new TCanvas("A","A", 600, 600);
    single_bin = new TGraphErrors(n_factors, &jec_factor[0], &events[0], &jec_factor_err[0], &events_err[0]);

    single_bin->SetTitle(mass_bin[bin]);
    single_bin->GetXaxis()->SetTitleOffset(0.9);
    single_bin->GetXaxis()->SetTitle("m_{Wjet}");

    single_bin->SetMarkerStyle(2);
    single_bin->Fit("fit"); // set first fit function (search for it)
    single_bin->Draw("APE");
    A->SaveAs(save_path+"/single_bins/Wjet_mass_sensitivity_JEC_rebin"+str_number_bins+"_"+number_bin[bin]+".pdf");
    delete A;

    // Now plot all graphs into one file ----------------------------------------
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

  E->SaveAs(save_path+"/mass_bin_all_rebin"+str_number_bins+".pdf");

  cout.precision(6);
  TString bin, a_str, a_err_str, b_str, b_err_str, chi2_str;
  cout << "" << '\n' << fixed;
  cout << "|=================================================================|" << endl;
  cout << "| ------------------------- Fit Values ---------------------------|" << endl;
  cout << "| ----------------------------------------------------------------|" << endl;
  cout << "|   Bin  |    a     |  a_err   |     b     |   b_err  |   chi2    |" << endl;
  cout << "| -------|----------|----------|-----------|----------|-----------|" << endl;
  for(unsigned int i=0; i<number_bins-first_bins; i++){
    if(i+1+first_bins<10)     bin = "|     "+to_string(i+1+first_bins)+"  |"; // i+1 to get bin number
    if(i+1+first_bins>=10)    bin = "|    "+to_string(i+1+first_bins)+"  |";  // +first_bins to skip first bins
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
    while(parameter<NParams*(number_bins-first_bins)){ // first bins does not include Data
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

    // get parameters for terms ------------------------------------------------------
    vector<double> parameters;
    for(int bin=0; bin<number_bins-first_bins; bin++){ // first bins does not include Data
      parameters.push_back(data_bin_content[bin]);
      parameters.push_back(a[bin]);
      parameters.push_back(b[bin]);
      parameters.push_back(data_bin_error[bin]);
    }

    // Set all parameters ------------------------------------------------------------
    for(int ipar=0; ipar<NParams*(number_bins-first_bins); ipar++) chi2_function->SetParameter(ipar, parameters[ipar]);
    chi2_function->SetTitle("#chi^{2} of "+year);

    // Get Minimum of chi2 -----------------------------------------------------------
    if(debug) cout << "Chi2 - Minimum" << endl;
    double minY = chi2_function->GetMinimum();
    double minX = chi2_function->GetX(minY, -10, 10);
    double X_sigmaup = chi2_function->GetX(minY+1, minX, 10);
    double X_sigmadown = chi2_function->GetX(minY+1, -10, minX);

    double sigmaup = X_sigmaup - minX; double sigmadown = minX - X_sigmadown;

    // Print and save Factor ---------------------------------------------------------
    cout << "" << endl;
    cout << " JEC factor = " << minX << " +" << sigmaup << " -" << sigmadown << endl;

    fstream jec_txt;
    jec_txt.open(save_path+"/jec_factor_"+str_number_bins+".txt", ios::out);
    jec_txt << "JEC factor = " << minX << " +" << sigmaup << " -" << sigmadown << endl;
    jec_txt.close();
    print_seperater();

    // -------------------------------------------------------------------------------
    // unified x axis interval -------------------------------------------------------
    if(debug) cout << "Chi2 - shifted Function" << endl;
    TF1 *chi2_function_shifted = new TF1("chi2_function_shifted", str_chi2_function, minX-3, minX+3);

    // Set all parameters ------------------------------------------------------------
    for(int ipar=0; ipar<NParams*(number_bins-first_bins); ipar++) chi2_function_shifted->SetParameter(ipar, parameters[ipar]);
    chi2_function_shifted->SetTitle("#chi^{2} of "+year);

    TCanvas *B = new TCanvas("B","B", 600, 600);
    chi2_function_shifted->Draw();
    B->SaveAs(save_path+"/chi2_rebin"+str_number_bins+"_shifted.pdf");

    /*
    .██████  ██████  ██████  ██████  ███████ ████████ ███████ ██████      ██   ██ ██ ███████ ████████
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██ ██         ██
    ██      ██    ██ ██████  ██████  █████      ██    █████   ██   ██     ███████ ██ ███████    ██
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██      ██    ██
    .██████  ██████  ██   ██ ██   ██ ███████    ██    ███████ ██████      ██   ██ ██ ███████    ██
    */
    if(debug) cout << "New Histogram with correct Minimum from Chi2" << endl;
    cout << '\n';

    vector<double> JEC_min, JEC_min_up_err, JEC_min_down_err;
    for(int bin=0; bin<first_bins; bin++){
      JEC_min.push_back(0);
      JEC_min_up_err.push_back(0);
      JEC_min_down_err.push_back(0);
    }

    for(unsigned int bin=0; bin<(number_bins-first_bins); bin++){
      double a_par = a[bin];
      double b_par = b[bin];
      double sigmaup_value = a_par + b_par * X_sigmaup; //
      double sigmadown_value = a_par + b_par * X_sigmadown; //
      double nominal_value = a_par + b_par * minX;
      JEC_min.push_back(nominal_value);
      if(sigmaup_value<nominal_value && nominal_value<sigmadown_value){
        if(debug) cout << "Bin " << bin << ": up(down) is down(up)" << endl;
        JEC_min_up_err.push_back(abs(JEC_min[bin+first_bins]-sigmadown_value)); // Here and next line: (JEC_min_up & _down)
        JEC_min_down_err.push_back(abs(JEC_min[bin+first_bins]-sigmaup_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else if(sigmadown_value<nominal_value && nominal_value<sigmaup_value){
        if(debug) cout << "Bin " << bin << ": up(down) is up(down)" << endl;
        JEC_min_up_err.push_back(abs(JEC_min[bin+first_bins]-sigmaup_value)); // Here and next line: (JEC_min_up & _down)
        JEC_min_down_err.push_back(abs(JEC_min[bin+first_bins]-sigmadown_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else throw runtime_error("Something is wrong with the bin value for the FSR variation using the calculated minimum. both variations are smaller or greater than the nominal value");
    }


    TH1F *h_JEC_min = new TH1F("JEC_min", "", number_bins, 0, 180);
    for(unsigned int bin=1; bin<number_bins+1;bin++) h_JEC_min->SetBinContent(bin, JEC_min[bin-1]);

    // Integral-----------------------------------------------------------------------------------------------
    double JEC_min_integral = h_JEC_min->Integral();
    cout << "Integral of calculated JEC_min distribution: " << JEC_min_integral << endl;

    h_JEC_min->SetTitle(year);
    h_JEC_min->GetXaxis()->SetRangeUser(0, 180);
    h_JEC_min->GetYaxis()->SetRangeUser(0, h_JEC_min->GetMaximum()*1.5);
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

    JECdown_rebin_norm->SetLineStyle(2);

    vector<double> xbins, xbins_err;                                              // = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
    //vector<double> xbins_err = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    int bin_steps = rebin_hist;                                              // 180 because we are looking at m_Wjet from 0 to 180 GeV
    int bin_position = bin_steps/2;                                               // center of the bins
    for(unsigned int bin=0; bin<number_bins; bin++){
      xbins.push_back(bin_position+bin_steps*bin);
      xbins_err.push_back(bin_position);
    }

    cout << '\n' << xbins.size() << '\n';
    auto errors = new TGraphAsymmErrors(number_bins, &xbins[0], &JEC_min[0], &xbins_err[0], &xbins_err[0], &JEC_min_down_err[0], &JEC_min_up_err[0]);
    errors->SetFillColor(kGray);
    errors->SetFillStyle(1001); // 1001: solid - same as nothing
    errors->SetTitle(year);

    if(debug) cout << "Draw New Histogram with correct Minimum from Chi2" << endl;
    TCanvas *C = new TCanvas("C", "C", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    h_JEC_min->Draw("HIST");
    errors->Draw("SAME A2");
    h_JEC_min->Draw("SAME HIST"); // Draw again, because it will be overdrawn otherwise
    ttbar_rebin_norm->Draw("SAME HIST");
    data_rebin_norm->Draw("SAME P");

    leg = new TLegend(0.20,0.65,0.40,0.85);
    leg->AddEntry(ttbar_rebin_norm,"nominal","l");
    leg->AddEntry(h_JEC_min, "JEC best fit value","l");
    leg->AddEntry(data_rebin_norm, "data","l");
    leg->SetTextSize(0.02);
    leg->Draw();
    gPad->RedrawAxis();
    C->SaveAs(save_path+"/JEC_minimum_factor_rebin"+str_number_bins+".pdf");
    JECup_rebin_norm->Draw("SAME HIST");
    JECdown_rebin_norm->Draw("SAME HIST");
    leg->AddEntry(JECup_rebin_norm, "JECup","l");
    leg->AddEntry(JECdown_rebin_norm, "JECdown","l");
    C->SaveAs(save_path+"/JEC_minimum_factor_compare_rebin"+str_number_bins+".pdf");
    leg->Clear();
    delete C;
}
