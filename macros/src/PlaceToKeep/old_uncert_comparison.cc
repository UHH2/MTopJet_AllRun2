#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){

  TString chi2_str_hh, chi2_str_hl, chi2_str_lh, chi2_str_ll;
  vector<double> chi2_parameters_hh, chi2_parameters_hl, chi2_parameters_lh, chi2_parameters_ll;
  print_seperater();

  bool debug        = false;
  cout.precision(6);

  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_comparison_uncert <year>";
    return 0;
  }

  // #################################################################################################
  // Print Options ###################################################################################
  Int_t oldLevel = gErrorIgnoreLevel;
  // Set by: gErrorIgnoreLevel = ...

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
  save_path_general = creat_folder_and_path(save_path_general, "pt_bins"); // CHANGE_PT

  // #################################################################################################
  // creat subdirectories ############################################################################
  save_path_general = creat_folder_and_path(save_path_general, year);

  // #################################################################################################
  // get the mass bin title for the plots - e.g. "10 < m_{Wjet} < 20" ################################
  int number_bins = 45;
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
  for(int ptbin=0; ptbin<5; ptbin++){ // CHANGE_PT
    TString ptbin_str = to_string(ptbin);
    if(ptbin==0) addition="_hh";
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="_full";
    creat_folder(save_path_general, addition);

    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */

    cout << '\n'; // CHANGE_PT
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
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile *ttbar_file  = new TFile(dir+year+"/muon/"+ttbar_path);
    TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);
    ttbar->Add(bkg, 1);
    TH1F* ttbar_norm              = normalize(ttbar);
    vector<double> ttbar_norm_err = normalize_error(ttbar);

    // #################################################################################################
    // Get TTbar 1695 ##################################################################################
    if(debug) cout << "1695" << endl;

    TString ttbar_path_1695 = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root"; // CHANGE_MTOP
    TFile *ttbar_file_1695  = new TFile(dir+year+"/muon/"+ttbar_path_1695);
    TH1F  *ttbar_1695       = (TH1F*)ttbar_file_1695->Get(w_mass);
    ttbar_1695->Add(bkg, 1);
    TH1F  *ttbar_norm_1695  = normalize(ttbar_1695);
    vector<double> ttbar_norm_1695_err = normalize_error(ttbar_1695);

    // #################################################################################################
    // Get TTbar 1755 ##################################################################################
    if(debug) cout << "1755" << endl;

    TString ttbar_path_1755 = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root"; // CHANGE_MTOP
    TFile *ttbar_file_1755  = new TFile(dir+year+"/muon/"+ttbar_path_1755);
    TH1F  *ttbar_1755       = (TH1F*)ttbar_file_1755->Get(w_mass);
    ttbar_1755->Add(bkg, 1);
    TH1F  *ttbar_norm_1755  = normalize(ttbar_1755);
    vector<double> ttbar_norm_1755_err = normalize_error(ttbar_1755);

    /*
    ███████ ██ ████████     ██       ██████   ██████  ██████
    ██      ██    ██        ██      ██    ██ ██    ██ ██   ██
    █████   ██    ██        ██      ██    ██ ██    ██ ██████
    ██      ██    ██        ██      ██    ██ ██    ██ ██
    ██      ██    ██        ███████  ██████   ██████  ██
    */
    if(debug) cout << "Start: looping over bins" << endl;

    // TCanvas *E0 = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    // E0->SetCanvasSize(800, 1200);
    // E0->Divide(4,8);
    //
    // TCanvas *E1 = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    // E1->SetCanvasSize(800, 1200);
    // E1->Divide(4,8);
    //
    // TCanvas *E2 = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    // E2->SetCanvasSize(800, 1200);
    // E2->Divide(4,8);
    //
    // TCanvas *E3 = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    // E3->SetCanvasSize(800, 1200);
    // E3->Divide(4,8);
    //
    // TCanvas *E4 = new TCanvas(); // All fits on one page. Fit is done for a mass up to 180 GeV
    // E4->SetCanvasSize(800, 1200);
    // E4->Divide(4,8);

    TLegend *leg;

    for(int bin=0; bin < number_bins; bin++){
      if(debug) cout << "---------------------------------------------------------------------------------------- " << bin+1 << endl;
      if(debug) cout << mass_bin[bin] << endl;
      // Exclude bins which are empty (for Data) or out of the massrange (wmass_peak)
      // #################################################################################################
      // Fill Vectors ####################################################################################
      if(debug) cout << "Get Bin Content" << endl;

      vector<double> bin_content, bin_content_err;
      bin_content.push_back(ttbar_norm->GetBinContent(bin+1));
      bin_content.push_back(ttbar_norm_1695->GetBinContent(bin+1));
      bin_content.push_back(ttbar_norm_1755->GetBinContent(bin+1));

      bin_content_err.push_back(ttbar_norm_err[bin]);
      bin_content_err.push_back(ttbar_norm_1695_err[bin]);
      bin_content_err.push_back(ttbar_norm_1755_err[bin]);

      // #################################################################################################
      // Set TGraph ######################################################################################
      if(debug) cout << "TGraphError - Single bin " << bin+1 << endl;
      TGraphErrors *bin_fit_nom = new TGraphErrors();
      bin_fit_nom->SetTitle(mass_bin[bin]);
      bin_fit_nom->SetName(number_bin[bin]);
      bin_fit_nom->SetPoint(0, 2, bin_content[0]);
      bin_fit_nom->SetPointError(0, 0, bin_content_err[0]);
      bin_fit_nom->SetMarkerStyle(kFullCircle);
      bin_fit_nom->SetMarkerColor(kBlack);
      bin_fit_nom->SetMarkerSize(0.5);

      TGraphErrors *bin_fit_69 = new TGraphErrors();
      bin_fit_69->SetName(number_bin[bin]);
      bin_fit_69->SetPoint(0, 1, bin_content[1]);
      bin_fit_69->SetPointError(0, 0, bin_content_err[1]);
      bin_fit_69->SetMarkerStyle(kFullCircle);
      bin_fit_69->SetMarkerColor(kBlue);
      bin_fit_69->SetMarkerSize(0.5);

      TGraphErrors *bin_fit_75 = new TGraphErrors();
      bin_fit_75->SetName(number_bin[bin]);
      bin_fit_75->SetPoint(0, 3, bin_content[2]);
      bin_fit_75->SetPointError(0, 0, bin_content_err[2]);
      bin_fit_75->SetMarkerStyle(kFullCircle);
      bin_fit_75->SetMarkerColor(kRed);
      bin_fit_75->SetMarkerSize(0.5);

      /*
      ██████  ██       ██████  ████████ ███████
      ██   ██ ██      ██    ██    ██    ██
      ██████  ██      ██    ██    ██    ███████
      ██      ██      ██    ██    ██         ██
      ██      ███████  ██████     ██    ███████
      */
      if(debug) cout << "Plots" << endl;
      TCanvas *A = new TCanvas("A", "A", 600, 600);
      // if     (0  <=bin&&bin<36)  E0->cd(bin+1);
      // else if(36 <=bin&&bin<72)  E1->cd(bin+1-36);
      // else if(72 <=bin&&bin<108) E2->cd(bin+1-72);
      // else if(108<=bin&&bin<144) E3->cd(bin+1-108);
      // else if(144<=bin&&bin<180) E3->cd(bin+1-144);
      bin_fit_nom->SetTitle(mass_bin[bin]);
      bin_fit_nom->GetXaxis()->SetLimits(0, 4);
      bin_fit_nom->GetYaxis()->SetRangeUser(0, bin_content[0]*2.5);

      if(debug) cout << "Changed Canvas" << endl;
      gStyle->SetPadTickY(1);
      gStyle->SetPadTickX(1);
      gStyle->SetOptStat(kFALSE);
      gStyle->SetLegendBorderSize(0);

      if(debug) cout << "Set Draw" << endl;
      bin_fit_nom->Draw("AP");
      bin_fit_69->Draw("same P");
      bin_fit_75->Draw("same P");

      if(debug) cout << "Set Leg" << endl;
      leg = new TLegend(0.20,0.65,0.40,0.85);
      leg->AddEntry(bin_fit_nom,"nominal","p");
      leg->AddEntry(bin_fit_69,"1695","p");
      leg->AddEntry(bin_fit_75,"1755","p");
      leg->Draw();
      leg->Clear();

      A->SaveAs(save_path_general+"/"+addition+"/compare_bin_"+number_bin[bin]+".pdf");
    }

    /*
    ██       █████  ████████ ███████ ██   ██
    ██      ██   ██    ██    ██       ██ ██
    ██      ███████    ██    █████     ███
    ██      ██   ██    ██    ██       ██ ██
    ███████ ██   ██    ██    ███████ ██   ██
    */

    cout << "\n####################################################################\n";
    cout << "############################ Latex #################################\n";
    cout << "####################################################################\n";
    // suppress output with &> /dev/null

    // system("python python/single_bins_latex_uncert.py "+save_path_general+" "+addition);
    system("python python/single_bins_latex_uncert.py "+save_path_general+" "+addition+" single_bins &> /dev/null");

  }
  if(debug) cout << "Save Canvas" << endl;
  // E0->SaveAs(save_path_general+"/"+addition+"/compare_uncert_0.pdf");
  // E1->SaveAs(save_path_general+"/"+addition+"/compare_uncert_1.pdf");
  // E2->SaveAs(save_path_general+"/"+addition+"/compare_uncert_2.pdf");
  // E3->SaveAs(save_path_general+"/"+addition+"/compare_uncert_3.pdf");
  // E4->SaveAs(save_path_general+"/"+addition+"/compare_uncert_4.pdf");

}
