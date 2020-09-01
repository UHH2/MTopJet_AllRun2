#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  /*
  This script is equal to JEC_SYS, but apply the same analysis for the mTop variations.

  Explanation:
  1) This script loops many times over bins. The bin number starts at 1 (e.g. for the function GetBinContent)
  On the other hand, the bin content or titles are stored in a vector. The vector starts at 0.
  Mostly, the loops start at bin=0. Therefor, using the function GetBinContent the bin needs to be bin+1
  */
  // #################################################################################################
  // Declare different variables used in the code ####################################################

  bool debug = false;
  TString reconst = "btag"; // btag btag_cut btag_sel compare min_mass

  if(argc != 3){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number>\n" << "rebin_number: hist->Rebin(rebin_number)\n";
    return 0;
  }

  // Rebin -------------------------------------------------------------------------
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

  // Range -------------------------------------------------------------------------
  int peak_width = 20; int w_peak = 80; int max_mass = 180; // Only looking at the wmass peak. Boundaries are abitrary
  int lower_limit = w_peak-peak_width; int upper_limit = w_peak+peak_width;
  int first_bins = 0; int last_bins = 0;

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

  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"+year+"/"+reconst;

  // rmdir(save_path+"/rebin_"+str_number_bins);
  if(mkdir(save_path+"/rebin_"+str_number_bins,0777) != 0) cout << "couldn't creat rebin folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+" created\n";
  if(mkdir(save_path+"/rebin_"+str_number_bins+"/single_bins",0777) != 0) cout << "couldn't creat rebin/single_bins folder!\n";
  else cout << save_path+"/rebin_"+str_number_bins+"/single_bins created\n";
  save_path += "/rebin_"+str_number_bins; // New save_path

  if(mkdir(save_path+"/mtop",0777) != 0) cout << "couldn't creat rebin/mtop folder!\n";
  else cout << save_path+"/mtop created\n";
  save_path += "/mtop/";

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
  if(reconst=="btag_cut") w_mass = hist_class+"wmass_btagcut";                  // same mass as above mit with btag_high > 0.7
  if(reconst=="compare")  w_mass = hist_class+"wmass_comparison";               // mass closest to Wmass from subjet combination
  if(reconst=="btag_sel") w_mass = hist_class+"wmass_btagcut_one_high_btag_one_subjet"; // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
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
  /* Necessary for comparison plots */
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar = (TH1F*)ttbar_file->Get(w_mass);
  ttbar->Add(bkg, 1);


  // #################################################################################################
  // Get mTop ########################################################################################
  if(debug) cout << "mTop" << endl;
  vector<TString> mtop_variations = {"1665", "1695", "1715", "1735", "1755", "1785"};
  /* For 1695 and 1755 the JEC variations are also included */

  vector<TFile*> file_mtop_v;
  vector<TH1F*> hists_mtop_v;
  TString mtop_path;
  for(unsigned int variation=0; variation<mtop_variations.size(); variation++){
    if(debug) cout << to_string(variation)+'\n';
    mtop_path = ttbar_path;
    mtop_path.Insert(34, "_mtop"+mtop_variations[variation]);
    file_mtop_v.push_back(new TFile(dir+year+"/muon/"+mtop_path));
    hists_mtop_v.push_back((TH1F*)file_mtop_v[variation]->Get(w_mass));
    hists_mtop_v[variation]->Add(bkg, 1);
  }

  // JEC ---------------------------------------------------------------------------------------------

  // #################################################################################################
  // Normalize Hists #################################################################################
  if(debug) cout << "Normalize Hists" << endl;

  TH1F* ttbar_norm = normalize(ttbar);
  TH1F* data_norm = normalize(data);
  vector<TH1F*> mtop_norm = normalize(hists_mtop_v);

  // #################################################################################################
  // Rebin Hists #####################################################################################
  if(debug) cout << "Rebin Hists" << endl;

  TH1F* ttbar_rebin_norm = rebin(ttbar_norm, rebin_hist);
  TH1F* data_rebin_norm = rebin(data_norm, rebin_hist);
  vector<TH1F*> mtop_rebin_norm = rebin(mtop_norm, rebin_hist);

  // #################################################################################################
  // Hists into Vector ###############################################################################
  if(debug) cout << "Hists into vector" << endl;

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
  // #################################################################################################
  // mTop & nominal Comparison #######################################################################
  if(debug) cout << "Plots mTop" << endl;
  vector<int> mtop_color = {kBlack, kBlue, kGreen, kGray, kOrange, kViolet};

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

  for(unsigned int variation=0; variation<mtop_variations.size(); variation++){
    mtop_rebin_norm[variation]->SetLineWidth(2);
    mtop_rebin_norm[variation]->SetLineColor(mtop_color[variation]);

    TCanvas *AA = new TCanvas("AA", "AA", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    ttbar_rebin_norm->Draw("HIST");
    mtop_rebin_norm[variation]->Draw("SAME HIST");
    leg = new TLegend(0.6,0.65,0.8,0.85);
    leg->SetTextSize(0.3);
    leg->AddEntry(ttbar_rebin_norm,"nominal","l");
    leg->AddEntry(mtop_rebin_norm[variation],"m_{top}="+mtop_variations[variation],"l");
    leg->SetTextSize(0.05);
    leg->Draw();
    gPad->RedrawAxis();
    AA->SaveAs(save_path+"/mtop"+mtop_variations[variation]+"_rebin"+str_number_bins+".pdf");
    delete AA;
    leg->Clear();
  }

  // #################################################################################################
  // All in one ######################################################################################
  TCanvas *AB = new TCanvas("AB", "AB", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar_rebin_norm->Draw("HIST");
  mtop_rebin_norm[1]->Draw("SAME HIST");
  mtop_rebin_norm[4]->Draw("SAME HIST");
  leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->SetTextSize(0.3);
  leg->AddEntry(ttbar_rebin_norm,"nominal","l");
  leg->AddEntry(mtop_rebin_norm[1],"m_{top}="+mtop_variations[1],"l");
  leg->AddEntry(mtop_rebin_norm[4],"m_{top}="+mtop_variations[4],"l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->RedrawAxis();
  AB->SaveAs(save_path+"/mtop_main"+"_rebin"+str_number_bins+".pdf");
  mtop_rebin_norm[0]->Draw("SAME HIST");
  mtop_rebin_norm[2]->Draw("SAME HIST");
  mtop_rebin_norm[3]->Draw("SAME HIST");
  mtop_rebin_norm[5]->Draw("SAME HIST");
  leg->AddEntry(mtop_rebin_norm[0],"m_{top}="+mtop_variations[0],"l");
  leg->AddEntry(mtop_rebin_norm[2],"m_{top}="+mtop_variations[2],"l");
  leg->AddEntry(mtop_rebin_norm[3],"m_{top}="+mtop_variations[3],"l");
  leg->AddEntry(mtop_rebin_norm[5],"m_{top}="+mtop_variations[5],"l");
  leg->Draw();
  AB->SaveAs(save_path+"mtop_all"+"_rebin"+str_number_bins+".pdf");
  delete AB;
  leg->Clear();

  /*
  Looking at each bin seperatly. Getting the bin content of the nominal file and
  the variations and making a fit through the values. The variation JECup and JECdown
  assigned to the values +1 and -1.
  Afterward the fit parameters will be stored in vectors for later purpose.
  */
  for(unsigned int var=0; var<mtop_variations.size(); var++){
    if(var!=1 && var!=4) continue; // Only 1695 and 1755 since they also have JEC files.
    cout << "\n###############################################################################\n";
    cout << "####################### "+mtop_variations[var]+" ##################################################\n";
    cout << "###############################################################################\n\n";

    if(mkdir(save_path+mtop_variations[var],0777) != 0) cout << "couldn't creat rebin/mtop/"+mtop_variations[var]+" folder!\n";
    else cout << save_path+mtop_variations[var]+" created\n";

    if(debug) cout << "-----------------------------------------------------------" << var << endl;
    if(debug) cout << "Start: Dependency in WJet bins" << endl;

    // #################################################################################################
    // Getting JEC hists ###############################################################################
    TFile *JECup_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_variations[var]+"_JECup.root");
    TFile *JECdown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_variations[var]+"_JECdown.root");

    TH1F *JECup = (TH1F*)JECup_file->Get(w_mass);
    TH1F *JECdown = (TH1F*)JECdown_file->Get(w_mass);

    if(debug) cout << "JEC+Background" << endl;
    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);

    if(debug) cout << "JEC rebin and norm" << endl;
    TH1F *JECup_rebin_norm = rebin(normalize(JECup), rebin_hist);
    TH1F *JECdown_rebin_norm = rebin(normalize(JECdown), rebin_hist);

    // #################################################################################################
    // Mass Sensitivity ################################################################################
    if(debug) cout << "Plots mTop JEC sensitivity\n";
    mtop_rebin_norm[var]->SetTitle("");
    mtop_rebin_norm[var]->GetXaxis()->SetTitle("m_{Wjet}");
    mtop_rebin_norm[var]->SetLineColor(kRed);

    data_rebin_norm->SetMarkerStyle(8);  // data hist style
    data_rebin_norm->SetMarkerColor(kBlack);

    JECup_rebin_norm->SetLineWidth(2);
    JECup_rebin_norm->SetLineColor(kBlue);
    JECdown_rebin_norm->SetLineWidth(2);
    JECdown_rebin_norm->SetLineStyle(2);
    JECdown_rebin_norm->SetLineColor(kBlue);

    TCanvas *AC = new TCanvas("AC", "AC", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    mtop_rebin_norm[var]->Draw("HIST");
    JECup_rebin_norm->Draw("SAME HIST");
    JECdown_rebin_norm->Draw("SAME HIST");
    data_rebin_norm->Draw("SAME P");
    leg = new TLegend(0.6,0.65,0.8,0.85);
    leg->SetTextSize(0.3);
    leg->AddEntry(mtop_rebin_norm[var],"m_{top}"+mtop_variations[var],"l");
    leg->AddEntry(JECup_rebin_norm,"JECup","l");
    leg->AddEntry(JECdown_rebin_norm,"JECdown","l");
    leg->AddEntry(data_rebin_norm,"Data","p");
    leg->SetTextSize(0.03);
    leg->Draw();
    gPad->RedrawAxis();
    AC->SaveAs(save_path+mtop_variations[var]+"/JEC_sensitivity"+"_rebin"+str_number_bins+"_mtop"+mtop_variations[var]+".pdf");

    /*
    ███████ ██ ████████
    ██      ██    ██
    █████   ██    ██
    ██      ██    ██
    ██      ██    ██
    */

    // #################################################################################################
    // set first fit function ##########################################################################
    TF1 *fit = new TF1("fit", "[0] + [1]*x");
    vector<TF1*> fit_functions, fit_functions_norm;
    vector<double> a, a_err, b, b_err, chi2;

    // #################################################################################################
    // vectors for values (Bin Content, Errors, etc) ###################################################
    vector<double> jec_factor, events;
    vector<double> data_bin_content, data_bin_error;
    vector<double> mtop_bin_content;
    vector<double> jec_factor_err, events_err;
    vector<double> events_rel_err;

    vector<TH1F*> ordered_jec = {JECdown_rebin_norm, mtop_rebin_norm[var], JECup_rebin_norm};

    jec_factor = {-1, 0, 1};
    jec_factor_err = {0, 0, 0};

    // #################################################################################################
    // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
    vector<TString> mass_bin;
    for(int bin=0; bin<number_bins; bin++) mass_bin.push_back(to_string(bin*rebin_hist)+" < m_{Wjet} < "+to_string((bin+1)*rebin_hist));

    if(debug){
      cout << "number_bins " << number_bins << " | number titles " << mass_bin.size() << " | bin width " << rebin_hist << '\n';
      cout << "number_bins " << number_bins-first_bins << " | number titles " << mass_bin.size()-first_bins << " | bin width " << rebin_hist << '\n';
      for(unsigned int bin=first_bins; bin<mass_bin.size(); bin++) cout << mass_bin[bin] << '\n';
    }

    // #################################################################################################
    // number bin ######################################################################################
    vector<TString> number_bin;
    for(int bin=0; bin<number_bins; bin++) number_bin.push_back(to_string(bin+1));

    // #################################################################################################
    // Hists + Errs ####################################################################################
    int n_factors = 3; // JECdown, ttbar, JECup
    vector<vector<double>> ordered_jec_norm_err;
    /* for normalized uncertainty. Function is defined in HistogramUtils.h */
    vector<TH1F*> ordered_jec_normal = {rebin(JECdown, rebin_hist), rebin(ttbar, rebin_hist), rebin(JECup, rebin_hist)};
    for(unsigned int factor=0; factor<n_factors; factor++) ordered_jec_norm_err.push_back(normalize_error(ordered_jec_normal[factor]));

    // #################################################################################################
    // Canvas for als fits in one page #################################################################
    if(debug) cout << "Star: filling axis" << endl;

    TCanvas *E = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    if(rebin_hist==4) E->Divide(5,9);
    else if(rebin_hist==5) E->Divide(4,9);
    else if(rebin_hist==10) E->Divide(3,6);
    else cout << "Not the correct bin number - " << number_bins << '\n';
    E->SetCanvasSize(800, 1200);
    E->SetWindowSize(800, 1200);

    // #################################################################################################
    // #################################################################################################
    // #################################################################################################
    // Start for loop for fits #########################################################################
    for(int bin=0; bin < number_bins; bin++){
      /*
      For the mass peak, the first bins are excluded, since there is no comparison necessary.
      This means that the bin iterator in the loop needs to be added with the first bins.
      */
      // cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      if(bin<first_bins) continue; // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)

      // for the chi2 Method ---------------------------------------------------------------------------
      data_bin_content.push_back(data_rebin_norm->GetBinContent(bin+1)); // Explanation 1) (search for it)
      data_bin_error.push_back(data_rebin_norm->GetBinError(bin+1));
      mtop_bin_content.push_back(mtop_rebin_norm[var]->GetBinContent(bin+1));

      // Getting BinContent & Error for each Factor for this bin ---------------------------------------
      for(int factor=0; factor<n_factors; factor++){
        events.push_back(ordered_jec[factor]->GetBinContent(bin+1));
        events_err.push_back(ordered_jec[factor]->GetBinError(bin+1));
      }

      // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
      if(debug) cout << "Start: TGraphError - Single - bin " << bin+1 << endl;
      TGraphErrors *single_bin;

      // each Graph seperatly---------------------------------------------------------------------------
      // TCanvas *A = new TCanvas("A","A", 600, 600);
      single_bin = new TGraphErrors(n_factors, &jec_factor[0], &events[0], &jec_factor_err[0], &events_err[0]);
      single_bin->Fit("fit", "Q"); // set first fit function (search for it) --- Option Q: No print

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

    E->SaveAs(save_path+mtop_variations[var]+"/mass_bin_all_rebin"+str_number_bins+"_mtop"+mtop_variations[var]+".pdf");

    cout.precision(6);
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
    jec_txt.open(save_path+mtop_variations[var]+"/jec_factor_"+str_number_bins+".txt", ios::out);
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
    B->SaveAs(save_path+mtop_variations[var]+"/chi2_rebin"+str_number_bins+"_shifted_mtop"+mtop_variations[var]+".pdf");
    delete B;
  }
}
