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
  bool print = true;
  bool into_latex = true;
  TString reconst   = "btag";    // match btag_cut btag_sel compare min_mass
  TString MC_uncert = "central"; // central avg nominal
  cout.precision(6);

  if(argc != 4){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number> <only>\n" << "rebin_number: hist->Rebin(rebin_number)\n" << "Selection: match, btag_cut, btag_sel, compare, min_mass\n";
    cout << "Only one fit: 0(false), 1(true)\n";
    return 0;
  }

  // #################################################################################################
  // Only one fit for all bins #######################################################################
  bool only = stob(argv[3]);
  int same_fit;
  TString addition="";
  if(only){
    cout << "0: xy | 1: xy2 | 2: xyy2 | 3: x2y2 | 4: xx2yy2\n";
    cout << "Which fit function to use?: ";
    cin  >> same_fit;
    if     (same_fit==0) addition = "_xy";
    else if(same_fit==1) addition = "_xy2";
    else if(same_fit==2) addition = "_xyy2";
    else if(same_fit==3) addition = "_x2y2";
    else if(same_fit==4) addition = "_xx2yy2";
    else throw runtime_error("Not valid fit!");
  }

  // #################################################################################################
  // Print Options ###################################################################################
  Int_t oldLevel = gErrorIgnoreLevel;
  // Set by: gErrorIgnoreLevel = ...

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
  int peak_width = 40; int w_peak = 80; int max_mass = 180; // Only looking at the wmass peak. Boundaries are abitrary
  int lower_limit = w_peak-peak_width; int upper_limit = w_peak+peak_width;
  int first_bins = 0; int last_bins = 0;

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0) is16   = true;
  else if(strcmp(year, "2017")==0) is17   = true;
  else if(strcmp(year, "2018")==0) is18   = true;
  else if(strcmp(year, "all")==0)  isAll  = true;
  else if(strcmp(year, "1718")==0) is1718 = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018, all or 1718)");

  vector<TString> all_years;
  if(is16 || is17 || is18) all_years = {year};
  else if(isAll) all_years = {"2016", "2017", "2018"};
  else if(is1718) all_years = {"2017", "2018"};
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018, all or 1718)");

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

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/";

  if(mkdir(save_path_general+year,0777) != 0) cout << "couldn't creat "+year+" folder!\n";
  else cout << save_path_general+year+" created\n";
  save_path_general += year; // New save_path

  if(mkdir(save_path_general+"/"+reconst,0777) != 0) cout << "couldn't creat "+reconst+" folder!\n";
  else cout << save_path_general+"/"+reconst+" created\n";
  save_path_general += "/"+reconst; // New save_path

  if(mkdir(save_path_general+"/rebin"+str_number_bins,0777) != 0) cout << "couldn't creat rebin folder!\n";
  else cout << save_path_general+"/rebin"+str_number_bins+" created\n";
  if(mkdir(save_path_general+"/rebin"+str_number_bins+"/single_bins",0777) != 0) cout << "couldn't creat /rebin/single_bins folder!\n";
  else cout << save_path_general+"/rebin"+str_number_bins+"/single_bins created\n";
  save_path_general += "/rebin"+str_number_bins; // New save_path

  if(mkdir(save_path_general+"/projection",0777) != 0) cout << "couldn't creat projection folder!\n";
  else cout << save_path_general+"/projection created\n";

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  vector<TH1F*> allyears_bkg,  allyears_data,  allyears_ttbar;
  vector<TH1F*> allyears_JECup, allyears_JECdown, allyears_XCup, allyears_XCdown;

  cout << '\n';
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
  TH1F *data_rebin      = rebin(data, rebin_hist);
  TH1F *data_norm       = normalize(data);
  TH1F *data_rebin_norm = normalize(data_rebin);

  vector<double> data_rebin_norm_err = normalize_error(data_rebin);

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
  TH1F* ttbar_norm       = normalize(ttbar);
  TH1F* ttbar_rebin      = rebin(ttbar, rebin_hist);
  TH1F* ttbar_rebin_norm = rebin(ttbar_norm, rebin_hist);

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

  TH1F  *JECup   = AddHists(allyears_JECup, 1);
  TH1F  *JECdown = AddHists(allyears_JECdown, 1);

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

  TH1F  *XConeup   = AddHists(allyears_XCup, 1);
  TH1F  *XConedown = AddHists(allyears_XCdown, 1);

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

  vector<double> XConeup_rebin_norm_err   = normalize_error(XConeup_rebin);
  vector<double> XConedown_rebin_norm_err = normalize_error(XConedown_rebin);

  vector<TH1F*> corrections_up_v = {JECup_rebin_norm, XConeup_rebin_norm};
  vector<TH1F*> corrections_down_v = {JECdown_rebin_norm, XConedown_rebin_norm};

  // #################################################################################################
  // Vectors from above into one Hist ################################################################
  if(debug) cout << "All Year Vector into one Hist" << endl;

  // #################################################################################################
  // Hists for Fits ##################################################################################
  TH1I used_fits = TH1I("used_fits", "used_fits", number_bins, 1, number_bins+1);
  TH1I number_fits = TH1I("number_fits", "number_fits", 5, 0, 5);

  // #################################################################################################
  // Empty Data Bins #################################################################################

  vector<int> empty_bins_v;
  for(int bin=0; bin < number_bins; bin++){
    if(!(abs(data_rebin_norm->GetBinContent(bin+1))>0)){
      empty_bins_v.push_back(bin+1);
      if(debug) cout << "Empty Bins: " << empty_bins_v[bin] << "\n";
    }
  }
  int number_empty_bins = empty_bins_v.size();
  if(debug) cout << "Number Empty Bins: " << number_empty_bins << "\n";

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

  // Data --------------------------------------------------------------------------------------------
  data_rebin_norm->SetMarkerStyle(8);  // data hist style
  data_rebin_norm->SetMarkerColor(kBlack);

  // Legend ------------------------------------------------------------------------------------------
  TLegend *leg;

  /*
  ███    ███  █████  ███████ ███████     ██████  ██       ██████  ████████ ███████
  ████  ████ ██   ██ ██      ██          ██   ██ ██      ██    ██    ██    ██
  ██ ████ ██ ███████ ███████ ███████     ██████  ██      ██    ██    ██    ███████
  ██  ██  ██ ██   ██      ██      ██     ██      ██      ██    ██    ██         ██
  ██      ██ ██   ██ ███████ ███████     ██      ███████  ██████     ██    ███████
  */

  for(int correction=0; correction<2; correction++){
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetLegendBorderSize(0);
    /*
    Plotting the histograms included above.
    */
    if(debug) cout << "Mass Plots" << endl;

    corrections_up_v[correction]->SetLineWidth(2);
    corrections_up_v[correction]->SetLineColor(kBlue);

    corrections_down_v[correction]->SetLineStyle(2);
    corrections_down_v[correction]->SetLineColor(kBlue);

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    ttbar_rebin_norm->Draw("HIST");
    corrections_down_v[correction]->Draw("SAME HIST");
    corrections_up_v[correction]->Draw("SAME HIST");
    data_rebin_norm->Draw("SAME P");
    leg = new TLegend(0.6,0.65,0.8,0.85);
    leg->SetTextSize(0.2);
    leg->AddEntry(ttbar_rebin_norm,"nominal","l");
    if(correction==0){
      leg->AddEntry(corrections_up_v[correction],"JECup","l");
      leg->AddEntry(corrections_down_v[correction],"JECdown","l");
    }
    if(correction==1){
      leg->AddEntry(corrections_up_v[correction],"XConeup","l");
      leg->AddEntry(corrections_down_v[correction],"XConedown","l");
    }
    leg->AddEntry(data_rebin_norm,"Data","p");
    leg->SetTextSize(0.05);
    leg->Draw();
    gPad->RedrawAxis();
    if(correction==0) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_JEC.pdf");
    if(correction==1) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_XCone.pdf");
    delete A;
    leg->Clear();
  }


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

  vector<TH1F*> XCone_ordered = {XConedown_rebin_norm, XConeup_rebin_norm};
  vector<TH1F*> JEC_ordered = {JECdown_rebin_norm, JECup_rebin_norm};

  // for the chi2 Method -----------------------------------------------------------------------------
  if(debug) cout << "Vectors for Chi2" << endl;
  vector<double> data_rebin_norm_bin_content;
  for(int bin=1; bin < number_bins+1; bin++) data_rebin_norm_bin_content.push_back(data_rebin_norm->GetBinContent(bin));

  // #################################################################################################
  // Cout Bin Error ##################################################################################
  if(print){
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
      Xd   = XConeup_rebin_norm->GetBinContent(i+1);
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
  }

  // #################################################################################################
  // Cout Bin Error ##################################################################################
  if(print){
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
      if(i+1>=10)    bin = to_string(i+1);                        // +first_bins to skip first bins
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
      if(i+1>=10)    bin = to_string(i+1);                        // +first_bins to skip first bins
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
  TString str_fit_xyy2 = "[0] + [1]*x + [2]*y + [3]*y*y"; // 4
  TString str_fit_quad = "[0] + [1]*x*x + [2]*y*y"; // 3
  TString str_fit_poly = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y"; // 5
  vector<TString> str_fits  = {str_fit_lin, str_fit_xy2, str_fit_xyy2, str_fit_quad, str_fit_poly};
  vector<TString> name_fits = {"linear", "mixed (xy2)", "mixed (xyy2)", "quadratic", "polynomial of order 2"};
  double npar=0;
  vector<int> ndfs  = {2, 2, 1, 2, 0}; // ndf = n_points - n_parameters = 5 - n_para
  vector<int> n_fit_para  = {3, 3, 4, 3, 5};
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

  /*
  ███████ ██ ████████     ██       ██████   ██████  ██████
  ██      ██    ██        ██      ██    ██ ██    ██ ██   ██
  █████   ██    ██        ██      ██    ██ ██    ██ ██████
  ██      ██    ██        ██      ██    ██ ██    ██ ██
  ██      ██    ██        ███████  ██████   ██████  ██
  */

  if(debug) cout << "Start: looping over bins" << endl;
  for(int bin=0; bin < number_bins; bin++){
    cout << "\n";
    if(print || debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
    if(find(empty_bins_v.begin(), empty_bins_v.end(), bin+1) != empty_bins_v.end()) continue;
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
    /* The additional points are necessary to scale the z axis properly. One can not zoom out of the z axis without
    cutting of the errors above the highest- and below the lowest point. The trick is to add two additional points
    to streth the z axis. Afterward one can zoom in (and not out).*/
    bin_fit = new TGraph2DErrors(5, &factor_x[0], &factor_y[0],&bin_content[0],&error_x[0],&error_y[0],&bin_content_err[0]);
    bin_fit->SetName(number_bin[bin]); // To avoid Warning: Replacing existing TGraph2D

    // #################################################################################################
    // set fit function ################################################################################
    if(debug) cout << "\nDefine fit function" << endl;
    // TF2 *fit_con   = new TF2("fit_con",  "[0]");
    TF2 *fit_lin  = new TF2("fit_lin",  str_fits[0]);
    TF2 *fit_xy2  = new TF2("fit_xy2",  str_fits[1]);
    TF2 *fit_xyy2 = new TF2("fit_xyy2", str_fits[2]);
    TF2 *fit_quad = new TF2("fit_quad", str_fits[3]);
    TF2 *fit_poly = new TF2("fit_poly", str_fits[4]);
    vector<TF2*> fits = {fit_lin, fit_xy2, fit_xyy2, fit_quad, fit_poly};

    // #################################################################################################
    // Fit Bin #########################################################################################
    if(debug) cout << "Start: Fitting" << endl;
    double chi2_fit;
    TF2 *used_fit;
    int order_index = 0;
    for(unsigned int order=0; order<fits.size(); order++){
      if(order!=same_fit && only) continue;
      if(debug) cout << "\nfit a "+name_fits[order]+" function ------------\n";
      if(debug) bin_fit->Fit(fits[order]);
      else      bin_fit->Fit(fits[order], "Q");
      chi2_fit = fits[order]->GetChisquare();
      double prob = TMath::Prob(chi2_fit, ndfs[order]);
      cout << "ndfs: " << ndfs[order] << " | chi2: " << chi2_fit <<" | prob: " << prob << endl;

      if(prob>0.05 || ndfs[order]==0){
        cout << "Bin " << number_bin[bin] << ": A " << GREEN << name_fits[order] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
        used_fit = new TF2(*fits[order]); // in if clause
        str_all_fits.push_back(str_fits[order]);
        all_fits.push_back(fits[order]);
        npar += n_fit_para[order];
        all_n_para.push_back(n_fit_para[order]);
        order_index=order;
        if(!only){
          cout << order << endl;
          number_fits.AddBinContent(order+1);
          used_fits.AddBinContent(bin+1, order+1);
        }
        break;
      }
      else{
        cout << "Bin " << number_bin[bin] << ": "<< RED << name_fits[order] << RESET << " fit is rejected (" << prob << ")" << endl;
        if(only){
          used_fit = new TF2(*fits[order]); // in if clause
          str_all_fits.push_back(str_fits[order]);
          all_fits.push_back(fits[order]);
          npar += n_fit_para[order];
          all_n_para.push_back(n_fit_para[order]);
          order_index=order;
        }
      }
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
    TCanvas *B = new TCanvas(number_bin[bin]+"B", "B", 600, 600);
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
    for(int ipar=0; ipar<n_fit_para[order_index]; ipar++){
      chi2_parameters.push_back(used_fit->GetParameter(ipar));
      stacker.push_back(used_fit->GetParameter(ipar));
      stacker_err.push_back(used_fit->GetParError(ipar));
    }
    stacker.push_back(data_rebin_norm_err[bin]);
    chi2_parameters.push_back(data_rebin_norm_err[bin]);

    fits_parameters.push_back(stacker);
    fits_parameters_err.push_back(stacker_err);

  }

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
  TCanvas *H = new TCanvas("H", "H", 600, 600);
  used_fits.Draw("P");
  H->SaveAs(save_path_general+"/used_fits_per_bin.pdf");


  number_fits.SetTitle("");
  number_fits.GetXaxis()->SetTitle("Fit");
  number_fits.GetYaxis()->SetTitle("Number");
  TCanvas *H1 = new TCanvas("H1", "H1", 600, 600);
  number_fits.Draw("Hist");
  leg = new TLegend(0.4,0.65,0.7,0.85);
  for(unsigned int i=0; i<str_fits.size();i++) leg->AddEntry((TObject*)0,to_string(i)+": "+str_fits[i],"");
  leg->SetTextSize(0.025);
  leg->Draw();
  gPad->RedrawAxis();
  H1->SaveAs(save_path_general+"/number_fits.pdf");
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
  if(str_all_fits.size()!= number_bins-number_empty_bins) throw runtime_error("Something is wrong with the number of fits");

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
  TF2 *chi2_function = new TF2("chi2_function", str_chi2_function, -3, 3, -3, 3);
  for(unsigned int ipar=0; ipar<chi2_parameters.size(); ipar++) chi2_function->SetParameter(ipar, chi2_parameters[ipar]);
  chi2_function->SetTitle("");
  chi2_function->GetXaxis()->SetTitle("JEC");
  chi2_function->GetYaxis()->SetTitle("XCone");
  chi2_function->SetFillStyle(1000);
  chi2_function->SetLineWidth(1);

  // #################################################################################################
  // Get Formula with parameters #####################################################################
  TString chi2_function_parameters = str_chi2_function;
  fstream formula_txt;
  formula_txt.open(save_path_general+"/formula"+addition+".txt", ios::out);
  for(unsigned int ipar=0; ipar<chi2_parameters.size(); ipar++){
    chi2_function_parameters.ReplaceAll("["+to_string(ipar)+"]", "["+to_string(chi2_function->GetParameter(ipar))+"]");
  }
  formula_txt << chi2_function_parameters << endl;
  formula_txt.close();

  // #################################################################################################
  // Minimum #########################################################################################
  /*
  Difference in the two functions: https://root.cern.ch/doc/master/TF2_8cxx_source.html#l00345 line 407 and 447.
  */
  double twoD_minX, twoD_minY;
  double twoD_minZ = chi2_function->GetMinimumXY(twoD_minX,twoD_minY);
  cout  << "z_min (->GetMinimumXY()): " << twoD_minZ << " | x value: " << twoD_minX << " | y value: " << twoD_minY << '\n' << endl;
  int number_para = chi2_function->GetNpar();

  fstream jec_txt;
  jec_txt.open(save_path_general+"/jec_factor"+addition+".txt", ios::out);
  jec_txt << number_bins << " bins [0, 180] GeV\n\n";
  jec_txt << "z_min: " << twoD_minZ << "\nx    :  " << twoD_minX << "\ny    :  " << twoD_minY << endl;
  jec_txt.close();

  // #################################################################################################
  // Evaluate 1sigma #################################################################################

  // Draw Points -------------------------------------------------------------------------------------
  vector<TPoint*> all_points_draw;
  TGraph2D *zmin_point = new TGraph2D();

  vector<vector<double>> points;
  auto start = high_resolution_clock::now(); // Calculation time - start
  points = FindXY(chi2_function, twoD_minZ+2.3, twoD_minX-2, twoD_minX+2, twoD_minY-2, twoD_minY+2, 4000, 0.001);
  auto stop = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "Numeric solution for 1sigma estimation took " << GREEN << duration.count() << "ms" << RESET << endl;
  if(print) cout << "Number Points at 1sigma: " << points.size() << endl;
  // for(unsigned int i=0; i<points.size(); i++) cout << "x: " << points[i][0] << " | y: " << points[i][1] << " | z: " << points[i][2] << endl;

  zmin_point->SetPoint(0, twoD_minX, twoD_minY, twoD_minZ);
  zmin_point->SetMarkerColor(kGray);
  zmin_point->SetMarkerStyle(kFullCircle);
  zmin_point->SetMarkerSize(0.5);
  TGraph2D *sigma_points = new TGraph2D();
  for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
  // for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1]);

  sigma_points->SetMarkerColor(kRed);
  sigma_points->SetMarkerStyle(kFullCircle);
  sigma_points->SetMarkerSize(0.1);

  // #################################################################################################
  // Plot ############################################################################################
  if(debug) cout<< "\n2D Chi2 - plot\n";
  gErrorIgnoreLevel = kWarning; // suppress TCanvas output
  TString option = "cont4z";
  TString option_add = "";
  // Contours ----------------------------------------------------------------------------------------
  chi2_function->SetContour(50);
  gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
  TColor::InvertPalette();
  TCanvas *D = new TCanvas("D","D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  chi2_function->Draw("cont4z");
  D->SaveAs(save_path_general+"/chi2_cont4z"+addition+".pdf");
  D->SetTheta(90);
  D->SetPhi(0);
  sigma_points->Draw("SAME P");
  zmin_point->Draw("SAME P");
  D->SaveAs(save_path_general+"/chi2_cont4z_points"+addition+".pdf");
  D->Clear();

  D = new TCanvas("D","D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  chi2_function->Draw("surf2z");
  D->SaveAs(save_path_general+"/chi2_surf2z"+option_add+addition+".pdf");
  sigma_points->Draw("SAME P");
  zmin_point->Draw("SAME P");
  D->SaveAs(save_path_general+"/chi2_surf2z"+option_add+"_points"+addition+".pdf");
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
  // vector<TString>  slice_constant = {to_string(0), to_string(0)};
  vector<TString>  slice_constant;
  TString constant; // zero, min

  // #################################################################################################
  // All in One ######################################################################################
  if(debug) cout << "All fits in one" << endl;
  cout << endl;

  // #################################################################################################
  // Canvas for all bins in one ######################################################################

  // double number_rows = (number_bins-number_empty_bins)/4.;
  // number_rows  = ceil(number_rows);
  // int divison  = number_rows;
  // int divison2 = divison/9*1200;
  // if(rebin_hist==4) divison2 *= 1.5;
  // if(rebin_hist==3) divison2 *= 2.5;
  // if(rebin_hist==2) divison2 *= 4.5;
  // if(rebin_hist==1) divison2 *= 5.5;

  // TCanvas *E1 = new TCanvas();
  // #################################################################################################
  // Constant x/y ####################################################################################
  for(int con=0; con<2; con++){
    if(con==0){
      slice_constant = {to_string(twoD_minX), to_string(twoD_minY)};
      constant = "min";
    } else if(con==1){
      slice_constant = {to_string(0), to_string(0)};
      constant = "zero";
    }

    // #################################################################################################
    // Which correction ################################################################################
    for(int factor=0; factor<2; factor++){
      vector<double> oneD_factor     = {-1,0,1};
      vector<double> oneD_factor_err = { 0,0,0};

      // E1->Divide(4,divison);
      // E1->SetCanvasSize(800, divison2);
      // E1->SetWindowSize(800, divison2);

      // #################################################################################################
      // Bins ############################################################################################
      int empty_bins = 0; // to skip empty bins
      for(int bin=0; bin<number_bins; bin++){

        if(debug) cout << "------------------------------------------------------- " << bin+1 << endl;
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
        if(debug) cout << setw(8) << oneD_events[0] << "   " << setw(8) << oneD_events[1] << "   " << setw(8) << oneD_events[2] << "   \n";
        if(debug) cout << setw(8) << oneD_events_err[0] << "   " << setw(8) << oneD_events_err[1] << "   " << setw(8) << oneD_events_err[2] << "   \n";

        // XCone -----------------------------------------------------------------------------------------
        TString str_fit = str_all_fits[used_bin];
        if(factor==0) str_fit.ReplaceAll("y", slice_constant[factor]);
        else if(factor==1){
          str_fit.ReplaceAll("x", slice_constant[factor]);
          str_fit.ReplaceAll("y", "x");
        }

        TF1 *fit = new TF1("fit", str_fit, -2, 2);
        TGraphErrors *graph = new TGraphErrors(3, &oneD_factor[0], &oneD_events[0], &oneD_factor_err[0], &oneD_events_err[0]);
        for(int ipar=0; ipar<all_n_para[used_bin]; ipar++) fit->SetParameter(ipar, all_fits[used_bin]->GetParameter(ipar));

        // #################################################################################################
        // Plot ############################################################################################
        graph->SetTitle(mass_bin[bin]);

        if(con==1){
          gErrorIgnoreLevel = kWarning; // suppress TCanvas output
          TCanvas *E2 = new TCanvas(number_bin[bin]+"E2_"+factor_name[factor], "E2", 600, 600);
          if(factor==0) graph->GetXaxis()->SetTitle("JEC");
          else          graph->GetXaxis()->SetTitle("XCone");
          graph->Draw("AP");
          fit->Draw("SAME");
          E2->SaveAs(save_path_general+"/projection/Bin"+number_bin[bin]+"_"+factor_name[factor]+addition+".pdf");
          gErrorIgnoreLevel = oldLevel; // Reset TCanvas output
        }

        // E1->cd(used_bin+1); // No blank spaces in PDF
        // if(factor==0) graph->GetXaxis()->SetTitle("JEC");
        // else          graph->GetXaxis()->SetTitle("XCone");
        // graph->Draw("AP");
        // fit->Draw("SAME");

      }
      // E1->Print(save_path_general+"/fits_all_bins_"+constant+"_"+factor_name[factor]+addition+".pdf", "pdf");
      // E1->Clear();// E2->Clear(); E3->Clear(); E4->Clear(); E5->Clear(); E6->Clear();
    }
    if(into_latex) system("python single_bins_latex.py "+year+" "+to_string(number_bins));
  }
}
