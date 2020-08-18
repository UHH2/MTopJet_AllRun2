#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"
using namespace std;

int main(int argc, char* argv[]){
  bool debug = true;
  bool rel_error = false;

  if(argc != 2){
    cout << "Usage: ./Radiation_SYS <year>" << endl;
    return 0;
  }

  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];

  bool is16=false; bool is17=false; bool is18=false; bool is1718=false;
  if(strcmp(year, "2016")==0) is16 = true;
  else if(strcmp(year, "2017")==0) is17 = true;
  else if(strcmp(year, "2018")==0) is18 = true;
  else if(strcmp(year, "combined")==0) is1718 = true;
  else{
    cout << "" << endl;
    cout << "WHAT IS THE ************* YEAR?" << endl;
    cout << "... Sorry, had a long day, but can you please give me the correct year? (2016, 2017, 2018, combined -(17&18))" << endl;
    cout << "" << endl;
    throw runtime_error("By the way ... it is line 17");
  }

  if(is1718) year = "2017"; // this will be changed to "combined" again

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString tau32 = "comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau32";
  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/";


  // Difference between 2016 and 17&18
  std::vector<int> color;
  vector<TString> FSR_strings, FSRup_strings, FSRdown_strings;
  if(is16){
    color = {kBlue};
    FSR_strings = {"FSR_2"};
    FSRup_strings = {"FSRup_2"};
    FSRdown_strings = {"FSRdown_2"};
  }
  if(is17 || is18 || is1718){
    color = {kOrange, kBlue, kGreen};
    FSR_strings = {"FSR_sqrt2", "FSR_2", "FSR_4"};
    FSRup_strings = {"FSRup_sqrt2", "FSRup_2", "FSRup_4"};
    FSRdown_strings = {"FSRdown_sqrt2", "FSRdown_2", "FSRdown_4"};
  }

  // #################################################################################################
  // Get Background ##################################################################################
  if(debug) cout << "Background" << endl;

  vector<TFile*>  file_v;
  vector<TH1F*>   hist_v;
  vector<TString> path_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
  for(unsigned int i=0; i<path_v.size(); i++){
    file_v.push_back(new TFile(dir+year+"/muon/"+path_v[i]));
    cout << file_v.size() << endl;
  }
  for(unsigned int i=0; i<file_v.size(); i++){
    hist_v.push_back((TH1F*)file_v[i]->Get(tau32));
    cout << hist_v.size() << endl;
  }
  if(debug) cout << "Background6" << endl;
  TH1F *bkg = AddHists(hist_v, 1);
  if(debug) cout << "Background" << endl;
  if(is1718) bkg = add_second_year("2018", bkg, dir, path_v, tau32);

  // #################################################################################################
  // Get Data ########################################################################################
  if(debug) cout << "Data" << endl;

  TString data_path = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file = new TFile(dir+year+"/muon/"+data_path);
  TH1F *data = (TH1F*)data_file->Get(tau32);
  if(is1718) data = add_second_year("2018", data, dir, data_path, tau32);

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar = (TH1F*)ttbar_file->Get(tau32);
  ttbar->Add(bkg, 1);
  if(is1718) ttbar = add_second_year("2018", ttbar, dir, ttbar_path, tau32);

  // #################################################################################################
  // Get FSR #########################################################################################
  if(debug) cout << "FSR" << endl;

  tau32 = "AK8_tau_pass_rec_masscut_140/tau32_hadjet";
  vector<TFile*> FSRup_file, FSRdown_file;
  vector<TH1F*> FSRup, FSRdown;

  if(is17 || is18 || is1718){
    for(unsigned int i=0; i<FSRup_strings.size();i++){
      FSRup_file.push_back(new TFile(dir+year+"/muon/"+FSRup_strings[i]+"/"+ttbar_path));
      FSRdown_file.push_back(new TFile(dir+year+"/muon/"+FSRdown_strings[i]+"/"+ttbar_path));

      FSRup.push_back((TH1F*)FSRup_file[i]->Get(tau32));
      FSRdown.push_back((TH1F*)FSRdown_file[i]->Get(tau32));

      if(is1718){
        FSRup[i] = add_second_year("2018", FSRup[i], dir, FSRup_strings[i]+"/"+ttbar_path, tau32);
        FSRdown[i] = add_second_year("2018", FSRdown[i], dir, FSRdown_strings[i]+"/"+ttbar_path, tau32);
      }
    }
  }

  if(is16){
    TFile *FSRup_2016 = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root");
    TFile *FSRdown_2016 = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root");

    FSRup.push_back((TH1F*)FSRup_2016->Get(tau32));
    FSRdown.push_back((TH1F*)FSRdown_2016->Get(tau32));
  }

  // Add hists from FSR and background ---------------------------------------------------------------
  if(debug) cout << "FSR+Background" << endl;
  for(unsigned int i=0; i<FSRup.size();i++){
    FSRup[i]->Add(bkg, 1);
    FSRdown[i]->Add(bkg, 1);
  }

  // #################################################################################################
  // Rebin ###########################################################################################

  // Rebin 10 ----------------------------------------------------------------------------------------
  if(debug) cout << "Rebin 10" << endl;
  TH1F *ttbar_rebin10 = rebin(ttbar, 5);
  TH1F *data_rebin10 = rebin(data, 5);

  vector<TH1F*> FSRup_rebin10 = rebin(FSRup, 5);
  vector<TH1F*> FSRdown_rebin10 = rebin(FSRdown, 5);

  // For tau bin dependency
  vector<TH1F*> FSRup_rebin10_norm = normalize(FSRup_rebin10);
  vector<TH1F*> FSRdown_rebin10_norm = normalize(FSRdown_rebin10);

  // #################################################################################################
  // Vector for compact Plots ########################################################################
  if(debug) cout << "Compact Plots" << endl;

  vector<TH1F*> ttbar_all = {ttbar, ttbar_rebin10};
  vector<TH1F*> data_all = {data, data_rebin10};

  vector<vector<TH1F*>> FSRup_all = {FSRup, FSRup_rebin10};
  vector<vector<TH1F*>> FSRdown_all = {FSRdown, FSRdown_rebin10};

  // #################################################################################################
  // Normalize Vectors ###############################################################################
  if(debug) cout << "Normalized Vector" << endl;

  // Vectors have the same order like the ones above
  vector<TH1F*> ttbar_all_norm = normalize(ttbar_all);
  vector<TH1F*> data_all_norm = normalize(data_all);

  vector<vector<TH1F*>> FSRup_all_norm = normalize(FSRup_all);
  vector<vector<TH1F*>> FSRdown_all_norm = normalize(FSRdown_all);

  // #################################################################################################
  // Combine Non normalized and normalized ###########################################################
  if(debug) cout << "Norm and not norm Vectors" << endl;

  vector<vector<TH1F*>> ttbar_general = {ttbar_all, ttbar_all_norm};
  vector<vector<TH1F*>> data_general = {data_all, data_all_norm};

  vector<vector<vector<TH1F*>>> FSRup_general = {FSRup_all, FSRup_all_norm};
  vector<vector<vector<TH1F*>>> FSRdown_general = {FSRdown_all, FSRdown_all_norm};

  /*
  ██████  ██       ██████  ████████ ███████
  ██   ██ ██      ██    ██    ██    ██
  ██████  ██      ██    ██    ██    ███████
  ██      ██      ██    ██    ██         ██
  ██      ███████  ██████     ██    ███████
  */
  if(debug) cout << "Plots" << endl;
  if(is1718) year = "combined";


  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  TLegend *leg;

  for(unsigned int a=0; a<ttbar_general.size();a++){ //norm and not norm
    for(unsigned int b=0; b<ttbar_general[a].size();b++){ // not rebin and rebin
      // ttbar hist style
      ttbar_general[a].at(b)->SetLineWidth(2);
      ttbar_general[a].at(b)->SetLineColor(kRed);
      // data hist style
      data_general[a].at(b)->SetMarkerStyle(8);
      data_general[a].at(b)->SetMarkerColor(kBlack);
    }
  }

  // #############################################################################################
  // Single Plots ################################################################################
  if(debug) cout << "Sinlge Plots" << endl;

  for(unsigned int a=0; a<FSRup_general.size();a++){ // not norm and norm (ttbar, data and variations)
    for(unsigned int b=0; b<FSRup_general.at(a).size();b++){ // not rebin and rebin (ttbar, data and variations)
      for(unsigned int c=0; c<FSRup_general.at(a).at(b).size();c++){ // variation factors (sqrt_2, 2, 4)
        FSRup_general[a].at(b).at(c)->SetTitle(FSR_strings[c]);
        FSRup_general[a].at(b).at(c)->GetXaxis()->SetRangeUser(0, 400);
        FSRup_general[a].at(b).at(c)->GetYaxis()->SetRangeUser(0, FSRup_general[a].at(b).at(c)->GetMaximum()*1.2);
        FSRup_general[a].at(b).at(c)->GetXaxis()->SetNdivisions(505);
        FSRup_general[a].at(b).at(c)->GetYaxis()->SetNdivisions(505);
        FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitleSize(0.05);
        FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitleSize(0.04);
        FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitleOffset(0.9);
        FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitleOffset(1.5);
        FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitle("#tau_{32}");
        FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitle("");
        FSRup_general[a].at(b).at(c)->SetLineWidth(2);
        FSRup_general[a].at(b).at(c)->SetLineColor(color[c]);

        FSRdown_general[a].at(b).at(c)->SetLineWidth(2);
        FSRdown_general[a].at(b).at(c)->SetLineStyle(2);
        FSRdown_general[a].at(b).at(c)->SetLineColor(color[c]);

        TCanvas *A = new TCanvas("A", "A", 600, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        FSRup_general[a].at(b).at(c)->Draw("HIST");
        FSRdown_general[a].at(b).at(c)->Draw("SAME HIST");
        ttbar_general[a].at(b)->Draw("SAME HIST");
        data_general[a].at(b)->Draw("SAME P");
        if(a==0){ // not norm
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR_strings[c]+".pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR_strings[c]+"_Rebin10.pdf");
        }
        if(a==1){ // norm
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR_strings[c]+"_norm.pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR_strings[c]+"_Rebin10_norm.pdf");
        }
        delete A;
      }
    }
  }

  // #############################################################################################
  // All in one Plot #############################################################################
  cout << "" << endl;

  if(is17 || is18 || is1718){
    if(debug) cout << "All Plots in one" << endl;

    for(unsigned int a=0; a<FSRup_general.size();a++){ // not norm and norm
      for(unsigned int b=0; b<FSRup_general[a].size();b++){ // not rebin and rebin
        ttbar_general[a].at(b)->SetTitle("");
        ttbar_general[a].at(b)->GetXaxis()->SetRangeUser(0, 400);
        ttbar_general[a].at(b)->GetYaxis()->SetRangeUser(0, get_highest_peak(FSRup_general[a].at(b))*1.2);
        ttbar_general[a].at(b)->GetXaxis()->SetNdivisions(505);
        ttbar_general[a].at(b)->GetYaxis()->SetNdivisions(505);
        ttbar_general[a].at(b)->GetXaxis()->SetTitleSize(0.05);
        ttbar_general[a].at(b)->GetYaxis()->SetTitleSize(0.04);
        ttbar_general[a].at(b)->GetXaxis()->SetTitleOffset(0.9);
        ttbar_general[a].at(b)->GetYaxis()->SetTitleOffset(1.5);
        ttbar_general[a].at(b)->GetXaxis()->SetTitle("#tau_{32}");
        ttbar_general[a].at(b)->GetYaxis()->SetTitle("");

        for(unsigned int c=0; c<FSRup_general[a].at(b).size(); c++){
          FSRup_general[a].at(b).at(c)->SetLineWidth(2);
          FSRup_general[a].at(b).at(c)->SetLineColor(color[c]);

          FSRdown_general[a].at(b).at(c)->SetLineWidth(2);
          FSRdown_general[a].at(b).at(c)->SetLineColor(color[c]);
        }

        TCanvas *A = new TCanvas("A", "A", 600, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        ttbar_general[a].at(b)->Draw("HIST");
        data_general[a].at(b)->Draw("SAME P");
        for(unsigned int c=0; c<FSRup_general[a].at(b).size();c++){
          FSRup_general[a].at(b).at(c)->Draw("SAME HIST");
          FSRdown_general[a].at(b).at(c)->Draw("SAME HIST");
        }
        leg = new TLegend(0.25,0.65,0.45,0.85);
        leg->SetTextSize(0.5);
        leg->AddEntry(ttbar_general[a].at(b),"nominal","l");
        for(unsigned int c=0; c<FSRup_general[a].at(b).size();c++) leg->AddEntry(FSRup_general[a].at(b).at(c),FSR_strings[c],"l");
        leg->SetTextSize(0.05);
        leg->Draw();
        gPad->RedrawAxis();
        if(a==0){
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_all.pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin10.pdf");
        }
        if(a==1){
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_all_norm.pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin10_norm.pdf");
        }
        leg->Clear();
        delete A;
      }
    }
  }


  /*
  ███████ ██ ████████
  ██      ██    ██
  █████   ██    ██
  ██      ██    ██
  ██      ██    ██
  */


  /*
  **********************************************************************************************
  **************************** FIT  ************************************************************
  **********************************************************************************************
  */

  // set first fit function
  TF1 *fit = new TF1("fit", "[0] + [1]*log(x)");

  /*
  OVERVIEW ---------------------------------------------------------------------
  xxx - ttbar; data; FSRxx;
  FSRxx - FSRup; FSRdown;
  xxx_general = {xxx_all, xxx_all_norm}
  xxx_all     = {xxx, xxx_rebin10}
  FSRxx       = {FSRxx_sqrt2, FSRxx_2, FSRxx_4}
  */

  if(debug) cout << "-----------------------------------------------------------" << endl;
  if(debug) cout << "Start: Dependency in tau32 bins" << endl;

  // At first only for Masscut_140 and Rebin10
  int n_factors;
  vector<TH1F*> ordered_fsr, ordered_fsr_norm;

  vector<vector<double>> ordered_fsr_norm_err;
  vector<vector<double>> fit_bin_values, fit_bin_values_norm; // chi2, parameters and errors --- {a, a_err, b, b_err, chi2}
  vector<TF1*> fit_functions, fit_functions_norm;
  vector<double> a, a_err, b, b_err, chi2;
  int fit_value_size = 5;

  // first index: factor (order: 1/4, 1/2, 1/sqrt2, 1, sqrt2, 2, 4) --- second index: bin
  vector<double> mu, mu2, events, events_norm;
  vector<double> data_norm_bin_content, data_norm_bin_error, ttbar_norm_bin_content;
  vector<double> mu_err, events_err, events_norm_err;
  vector<double> events_rel_err, events_norm_rel_err;

  if(is16){
    n_factors = 3;
    ordered_fsr = {FSRdown_rebin10[0], ttbar_all[1], FSRup_rebin10[0]};
    ordered_fsr_norm = {FSRdown_rebin10_norm[0], ttbar_all_norm[1], FSRup_rebin10_norm[0]};
    mu = {0.5, 1, 2};
    mu_err = {0, 0, 0};
  }

  if(is17 || is18 || is1718){
    n_factors = 7;
    ordered_fsr = {FSRdown_rebin10[2], FSRdown_rebin10[1], FSRdown_rebin10[0], ttbar_all[1], FSRup_rebin10[0], FSRup_rebin10[1], FSRup_rebin10[2]};
    ordered_fsr_norm = {FSRdown_rebin10_norm[2], FSRdown_rebin10_norm[1], FSRdown_rebin10_norm[0], ttbar_all_norm[1], FSRup_rebin10_norm[0], FSRup_rebin10_norm[1], FSRup_rebin10_norm[2]};
    mu = {0.25, 0.5, 1/sqrt(2), 1, sqrt(2), 2, 4};
    mu_err = {0, 0, 0, 0, 0, 0, 0};
  }

  mu2 = square_vector(mu);
  vector<TString> tau_bin = {"0 < #tau_{32} < 0.1", "0.1 < #tau_{32} < 0.2", "0.2 < #tau_{32} < 0.3", "0.3 < #tau_{32} < 0.4", "0.4 < #tau_{32} < 0.5", "0.5 < #tau_{32} < 0.6", "0.6 < #tau_{32} < 0.7", "0.7 < #tau_{32} < 0.8", "0.8 < #tau_{32} < 0.9", "0.9 < #tau_{32} < 1"};
  vector<TString> number_bin = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  for(unsigned int factor=0; factor<n_factors; factor++) ordered_fsr_norm_err.push_back(normalize_error(ordered_fsr[factor]));

  for(unsigned int isnorm=0; isnorm<2; isnorm++){ // isnorm=0: normal --- isnorm=1: norm
    if(debug) cout << "Star: filling axis" << endl;

    cout << "" << endl;
    cout << "******************************************************************************************************************************" << endl;
    if(isnorm==0) cout << "*********** yAxis: normal --- xAxis: log *************************************************************************************" << endl;
    if(isnorm==1) cout << "*********** yAxis: norm   --- xAxis: log *************************************************************************************" << endl;
    cout << "******************************************************************************************************************************" << endl;
    cout << "" << endl;

    TCanvas *E = new TCanvas(); // All Graphs in one
    E->Divide(3,3);
    E->SetCanvasSize(800, 1200);
    E->SetWindowSize(800, 1200);

    // START FOR LOOP FOR FITS -------------------------------------------------
    for(unsigned int bin=3; bin < 11; bin++){ // The first two Bins are empty (for Data), therefor they are excluded
      cout << "------------------------------------------------------------------------------------------------" << endl;
      if(isnorm == 1){ // data_norm_rebin10 - for the chi2 Method
        data_norm_bin_content.push_back(data_all_norm[1]->GetBinContent(bin));
        data_norm_bin_error.push_back(data_all_norm[1]->GetBinError(bin));
        ttbar_norm_bin_content.push_back(ttbar_all_norm[1]->GetBinContent(bin));
      }
      for(int factor=0; factor<n_factors; factor++){
        events.push_back(ordered_fsr[factor]->GetBinContent(bin));
        events_err.push_back(ordered_fsr[factor]->GetBinError(bin));
        events_norm.push_back(ordered_fsr_norm[factor]->GetBinContent(bin));
        events_norm_err.push_back(ordered_fsr_norm_err[factor].at(bin-1));
      }

      // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
      if(debug) cout << "Start: TGraphError - Single - bin " << bin << endl;
      TGraphErrors *single_bin;

      TCanvas *A = new TCanvas("A","A", 600, 600); // each Graph seperatly
      if(isnorm==0) single_bin = new TGraphErrors(n_factors, &mu2[0], &events[0], &mu_err[0], &events_err[0]);
      if(isnorm==1) single_bin = new TGraphErrors(n_factors, &mu2[0], &events_norm[0], &mu_err[0], &events_norm_err[0]);

      single_bin->SetTitle(tau_bin[bin-1]);
      single_bin ->GetXaxis()->SetTitleOffset(0.9);
      gPad->SetLogx();
      single_bin->GetXaxis()->SetTitle("log(#mu)");

      single_bin->SetMarkerStyle(2);
      single_bin->Fit("fit");
      single_bin->Draw("APE");
      if(isnorm==0) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+"/single_bins/tau32_bin_"+number_bin[bin-1]+".pdf");
      if(isnorm==1) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+"/single_bins/tau32_bin_"+number_bin[bin-1]+"_norm.pdf");
      delete A;

      // Now plot all graphs into one file
      E->cd(bin-1);
      if(debug) cout << "Start: TGraphError - all in one - bin " << bin << endl;

      single_bin->SetTitle(tau_bin[bin-1]);
      single_bin->SetMarkerStyle(2);
      gPad->SetLogx();
      single_bin->GetXaxis()->SetTitle("log(#mu)");
      single_bin->Draw("AP");

      events = {}; events_err = {}; events_norm = {}; events_norm_err = {}; // clear vectors

      if(isnorm==1){
        /*
        Get values of the parameter of the fit
        Values du not depend on scale of ((not-)log) xAxis, but only on the scale of the yAxis ((not-)norm)
        Therefore the vectors for the Value are filled once
        */
        TF1 *fit = single_bin->GetFunction("fit");
        a.push_back(fit->GetParameter(0));
        a_err.push_back(fit->GetParError(0));
        b.push_back(fit->GetParameter(1));
        b_err.push_back(fit->GetParError(1));
        chi2.push_back(fit->GetChisquare());
        fit_functions_norm.push_back(fit);
      }

      cout << "" << endl;
    }

    if(isnorm==0) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+"/tau_bin_all.pdf");
    if(isnorm==1) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+"/tau_bin_all_norm.pdf");

  }

  cout.precision(6);
  cout << "" << '\n' << fixed;
  cout << "|================================================================|" << endl;
  cout << "| ------------------------- Fit Values --------------------------|" << endl;
  cout << "| ---------------------------------------------------------------|" << endl;
  cout << "|   Bin  |    a     |  a_err   |     b     |   b_err  |   chi2   |" << endl;
  cout << "| -------|----------|----------|-----------|----------|----------|" << endl;
  for(unsigned int i=0; i<8; i++){
    if(is17){
      if(i==0) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==1) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==2) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==3) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==4) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==5) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==6) cout<<"|     "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==7) cout<<"|    "<<i+3<<"  | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<< " | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
    }
    if(is18){
      if(i==0) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" |"<<chi2.at(i)<<" |"<< endl;
      if(i==1) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" |"<<chi2.at(i)<<" |"<< endl;
      if(i==2) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==3) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==4) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==5) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==6) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==7) cout<<"|   "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<< " | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
    }
    if(is16){
      if(i==0) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==1) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==2) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==3) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==4) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==5) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==6) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==7) cout<<"|   "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<< " | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
    }
    if(is1718){
      if(i==0) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" |"<<chi2.at(i)<<" |"<< endl;
      if(i==1) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" |" <<chi2.at(i)<<" |"<< endl;
      if(i==2) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" | "<<b.at(i)<<" | "<<b_err.at(i)<<" |" <<chi2.at(i)<<" |"<< endl;
      if(i==3) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" |" <<chi2.at(i)<<" |"<< endl;
      if(i==4) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==5) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==6) cout<<"|    "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<<" | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
      if(i==7) cout<<"|   "<<i+3<<"   | "<<a.at(i)<<" | "<<a_err.at(i)<<" |  "<<b.at(i)<< " | "<<b_err.at(i)<<" | "<<chi2.at(i)<<" |"<< endl;
    }
  }
  cout << "| ---------------------------------------------------------------|" << endl;
  cout << "" << endl;

  /*
  .██████ ██   ██ ██ ██████
  ██      ██   ██ ██      ██
  ██      ███████ ██  █████
  ██      ██   ██ ██ ██
  .██████ ██   ██ ██ ███████
  */
  if(debug) cout << "Chi2 - Parameters" << endl;
  // At first the chi2 fit is done for the normilized case.
  // [0]=data_bin_content | [1]=a | [2]=b | [3]=error (at first only data_bin_error)

  int NParams = 4;
  vector<vector<double>> StartParameters;
  // TString chi2_formula = "(([0]-[1]-[2]*log(x))/[3])^2";

  // for 2017, 2018
  // Chi2 term for all bins
  TF1 *b3 = new TF1("b3", "(([0]-[1]-[2]*log(x))/[3])^2");
  TF1 *b4 = new TF1("b4", "(([4]-[5]-[6]*log(x))/[7])^2");
  TF1 *b5 = new TF1("b5", "(([8]-[9]-[10]*log(x))/[11])^2");
  TF1 *b6 = new TF1("b6", "(([12]-[13]-[14]*log(x))/[15])^2");
  TF1 *b7 = new TF1("b7", "(([16]-[17]-[18]*log(x))/[19])^2");
  TF1 *b8 = new TF1("b8", "(([20]-[21]-[22]*log(x))/[23])^2");
  TF1 *b9 = new TF1("b9", "(([24]-[25]-[26]*log(x))/[27])^2");
  TF1 *b10 = new TF1("b10", "(([28]-[29]-[30]*log(x))/[31])^2");

  TF1 *chi2_function;
  if(is16) chi2_function = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 0.1, 10);
  if(is17) chi2_function = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 10, 200);
  if(is18) chi2_function = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 1, 50);
  if(is1718) chi2_function = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 1, 50);

  // get parameters for terms
  vector<double> parameters;
  for(int bin=0; bin<8; bin++){
    parameters.push_back(data_norm_bin_content[bin]);
    parameters.push_back(a[bin]);
    parameters.push_back(b[bin]);
    parameters.push_back(data_norm_bin_error[bin]);
  }

  // Set all parameters
  for(int ipar=0; ipar<36; ipar++) chi2_function->SetParameter(ipar, parameters[ipar]);
  chi2_function->SetTitle("#chi^{2} of "+year);

  if(debug) cout << "Chi2 - Function" << endl;
  TCanvas *B = new TCanvas("B","B", 600, 600);
  chi2_function->Draw();
  gPad->SetLogx();
  B->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+"/chi2.pdf");

  // Get Minimum of chi2
  if(debug) cout << "Chi2 - Minimum" << endl;
  double minY = chi2_function->GetMinimum();
  double minX2 = chi2_function->GetX(minY, 0, 100);
  double minX = sqrt(minX2);
  double X2_sigmaup = chi2_function->GetX(minY+1, minX2, 200);
  double X2_sigmadown = chi2_function->GetX(minY+1, 0.01, minX2);

  double X_sigmaup = sqrt(X2_sigmaup); double X_sigmadown = sqrt(X2_sigmadown);
  double sigmaup = X_sigmaup - minX;   double sigmadown = minX - X_sigmadown;

  // \u03BC = Mu
  cout << "" << endl;
  cout << " \u03BC = " << minX << " +" << sigmaup << " -" << sigmadown << endl;
  print_seperater();

  fstream fsr_txt;
  fsr_txt.open(save_path+year+"/fsr_factor.txt", ios::out);
  fsr_txt << "ymin     : " << minY << "\nxmin     : " << minX << "\nxmin_up  :  " << sigmaup << "\nxmin_down:  " << sigmadown  << endl;
  fsr_txt.close();

  /*
  .██████  ██████  ██████  ██████  ███████ ████████ ███████ ██████      ██   ██ ██ ███████ ████████
  ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██ ██         ██
  ██      ██    ██ ██████  ██████  █████      ██    █████   ██   ██     ███████ ██ ███████    ██
  ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██      ██    ██
  .██████  ██████  ██   ██ ██   ██ ███████    ██    ███████ ██████      ██   ██ ██ ███████    ██
  */
  if(debug) cout << "New Histogram with correct Minimum from Chi2" << endl;
  cout << '\n';
  /*
  Data is described mainly through the FSRup(_4) variation. Therefore, the calculated variation
  is only used for the up variation. The error for the calculated factor is therefore used for the up variation.
  */

  vector<double> FSRup_min = {0, 0}; vector<double> FSRup_min_up_err = {0, 0}; vector<double> FSRup_min_down_err = {0, 0};
  for(unsigned int bin=0; bin<8;bin++){ // Starting at Bin 3. For Completeness the first two bins are set to 0.
    double a_par = a[bin];
    double b_par = b[bin];
    double sigmaup_value = a_par + b_par * log(X2_sigmaup); // pow(minX*sqrt(2),2)
    double sigmadown_value = a_par + b_par * log(X2_sigmadown); // pow(minX*(1/sqrt(2)),2)
    double nominal_value = a_par + b_par * log(minX2);
    FSRup_min.push_back(nominal_value);
    if(sigmaup_value<nominal_value && nominal_value<sigmadown_value){
      if(debug) cout << "Bin " << bin << ": up(down) is down(up)" << endl;
      FSRup_min_up_err.push_back(abs(FSRup_min[bin+2]-sigmadown_value)); // Here and next line: (FSRup_min_up & _down)
      FSRup_min_down_err.push_back(abs(FSRup_min[bin+2]-sigmaup_value)); // take minus to get error for TGraphAsymmErrors (below)
    }
    else if(sigmadown_value<nominal_value && nominal_value<sigmaup_value){
      if(debug) cout << "Bin " << bin << ": up(down) is up(down)" << endl;
      FSRup_min_up_err.push_back(abs(FSRup_min[bin+2]-sigmaup_value)); // Here and next line: (FSRup_min_up & _down)
      FSRup_min_down_err.push_back(abs(FSRup_min[bin+2]-sigmadown_value)); // take minus to get error for TGraphAsymmErrors (below)
    }
    else throw runtime_error("Something is wrong with the bin value for the FSR variation using the calculated minimum. both variations are smaller or greater than the nominal value");
  }


  TH1F *h_FSRup_min = new TH1F("FSRup_min", "", 10, 0, 1);
  for(unsigned int bin=1; bin<11;bin++) h_FSRup_min->SetBinContent(bin, FSRup_min[bin-1]);

  // Integral-----------------------------------------------------------------------------------------------
  double fsr_min_integral = h_FSRup_min->Integral();
  cout << "Integral of calculated \u03C4_32 distribution: " << fsr_min_integral << endl;

  h_FSRup_min->SetTitle(year);
  h_FSRup_min->GetXaxis()->SetRangeUser(0, 400);
  h_FSRup_min->GetYaxis()->SetRangeUser(0, h_FSRup_min->GetMaximum()*1.5);
  h_FSRup_min->GetXaxis()->SetNdivisions(505);
  h_FSRup_min->GetYaxis()->SetNdivisions(505);
  h_FSRup_min->GetXaxis()->SetTitleSize(0.05);
  h_FSRup_min->GetYaxis()->SetTitleSize(0.04);
  h_FSRup_min->GetXaxis()->SetTitleOffset(0.9);
  h_FSRup_min->GetYaxis()->SetTitleOffset(1.5);
  h_FSRup_min->GetXaxis()->SetTitle("#tau_{32}");
  h_FSRup_min->GetYaxis()->SetTitle("");
  h_FSRup_min->SetLineWidth(1);
  h_FSRup_min->SetLineColor(kBlack);

  if(is17 || is18 || is1718) FSRdown_general[1].at(1).at(2)->SetLineStyle(2);
  if(is16) FSRdown_general[1].at(1).at(0)->SetLineStyle(2);

  vector<double> xbins = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  vector<double> xbins_err = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  auto errors = new TGraphAsymmErrors(10, &xbins[0], &FSRup_min[0], &xbins_err[0], &xbins_err[0], &FSRup_min_down_err[0], &FSRup_min_up_err[0]);
  errors->SetFillColor(kGray);
  errors->SetFillStyle(1001); // 1001: solid - same as nothing
  errors->SetTitle(year);

  if(debug) cout << "Draw New Histogram with correct Minimum from Chi2" << endl;
  TCanvas *C = new TCanvas("C", "C", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  h_FSRup_min->Draw("HIST");
  errors->Draw("SAME A2");
  h_FSRup_min->Draw("SAME HIST"); // Draw again, because it will be over over otherwise
  ttbar_all_norm[1]->Draw("SAME HIST");
  data_all_norm[1]->Draw("SAME P");

  leg = new TLegend(0.20,0.65,0.40,0.85);
  leg->AddEntry(ttbar_all_norm[1],"nominal","l");
  leg->AddEntry(h_FSRup_min, "FSR best fit value","l");
  leg->AddEntry(data_all_norm[1], "data","l");
  leg->SetTextSize(0.02);
  leg->Draw();
  gPad->RedrawAxis();
  C->SaveAs(save_path+year+"/tau32_minimum_factor.pdf");
  // C->SaveAs(save_path+year+"/Addition_err/tau32_minimum_factor_err_factorsqrt2.pdf");
  if(is17 || is18 || is1718) {
    FSRup_general[1].at(1).at(2)->Draw("SAME HIST");
    FSRdown_general[1].at(1).at(2)->Draw("SAME HIST");
    leg->AddEntry(FSRup_general[1].at(1).at(2), "FSRup 4","l");
    leg->AddEntry(FSRdown_general[1].at(1).at(2), "FSRdown 4","l");
  }
  if(is16) {
    FSRup_general[1].at(1).at(0)->Draw("SAME HIST");
    FSRdown_general[1].at(1).at(0)->Draw("SAME HIST");
    leg->AddEntry(FSRup_general[1].at(1).at(0), "FSRup 2","l");
    leg->AddEntry(FSRdown_general[1].at(1).at(0), "FSRdown 2","l");
  }
  C->SaveAs(save_path+year+"/tau32_minimum_factor_compare.pdf");
  // C->SaveAs(save_path+year+"/Addition_err/tau32_minimum_factor_compare_err_factorsqrt2.pdf");
  leg->Clear();
  delete C;
}
