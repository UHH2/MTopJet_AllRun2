#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(int argc, char* argv[]){
  bool debug = false;
  bool with_bkg = true;
  bool rel_error = true;

  if(argc != 2){
    /* <Background True/False> : True  - (ttbar + bkg) & False - (data - bkg) */
    cout << "Usage: ./Radiation_SYS <year>" << endl;
    return 0;
  }

  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];

  cout << "" << endl;
  cout << "###################################################################" << endl;
  if(with_bkg)  cout << "#################### WITH BKG #####################################" << endl;
  else          cout << "################### WITHOUT BKG ###################################" << endl;
  cout << "###################################################################" << endl;
  cout << "" << endl;

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString tau32 = "AK8_tau_pass_rec_masscut_140/tau32_hadjet";
  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/";

  // #################################################################################################
  // Get Background ##################################################################################
  if(debug) cout << "Background" << endl;

  vector<TFile*> file_v;
  vector<TH1F*> hist_v;
  vector<TString> path_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.ST.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
  for(unsigned int i=0; i<path_v.size(); i++) file_v.push_back(new TFile(dir+year+"/muon/"+path_v[i]));
  for(unsigned int i=0; i<file_v.size(); i++) hist_v.push_back((TH1F*)file_v[i]->Get(tau32));
  TH1F *bkg = AddHists(hist_v, 1);

  // #################################################################################################
  // Get Data ########################################################################################
  if(debug) cout << "Data" << endl;

  TString data_path = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file = new TFile(dir+year+"/muon/"+data_path);
  TH1F *data = (TH1F*)data_file->Get(tau32);
  if(!with_bkg) data->Add(bkg, -1);

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar = (TH1F*)ttbar_file->Get(tau32);
  if(with_bkg) ttbar->Add(bkg, 1);

  // #################################################################################################
  // Get FSR #########################################################################################
  if(debug) cout << "FSR" << endl;

  vector<TString> FSR_strings = {"FSR_sqrt2", "FSR_2", "FSR_4"};
  vector<TString> FSRup_strings = {"FSRup_sqrt2", "FSRup_2", "FSRup_4"};
  vector<TString> FSRdown_strings = {"FSRdown_sqrt2", "FSRdown_2", "FSRdown_4"};

  vector<TFile*> FSRup_file, FSRdown_file;
  vector<TH1F*> FSRup, FSRdown;

  for(unsigned int i=0; i<FSRup_strings.size();i++){
    FSRup_file.push_back(new TFile(dir+year+"/muon/"+FSRup_strings[i]+"/"+ttbar_path));
    FSRdown_file.push_back(new TFile(dir+year+"/muon/"+FSRdown_strings[i]+"/"+ttbar_path));

    FSRup.push_back((TH1F*)FSRup_file[i]->Get(tau32));
    FSRdown.push_back((TH1F*)FSRdown_file[i]->Get(tau32));
  }

  // Add hists from FSR and background ---------------------------------------------------------------
  if(with_bkg){
    if(debug) cout << "FSR+Background" << endl;
    for(unsigned int i=0; i<FSRup.size();i++){
      FSRup[i]->Add(bkg, 1);
      FSRdown[i]->Add(bkg, 1);
    }
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
  ██████  ██ ███    ██
  ██   ██ ██ ████   ██
  ██████  ██ ██ ██  ██
  ██   ██ ██ ██  ██ ██
  ██████  ██ ██   ████
  */

  /*
  OVERVIEW ---------------------------------------------------------------------
  xxx - ttbar; data; FSRxx;
  FSRxx - FSRup; FSRdown;
  xxx_general = {xxx_all, xxx_all_norm}
  xxx_all     = {xxx, xxx_rebin10}
  FSRxx       = {FSRxx_sqrt2, FSRxx_2, FSRxx_4}
  */

  TString folder_rel_err;
  if(rel_error) folder_rel_err = "/rel_err";
  else          folder_rel_err = "/normal_err";

  if(debug) cout << "-----------------------------------------------------------" << endl;
  if(debug) cout << "Start: Dependency in tau32 bins" << endl;

  int n_factors = 7;
  // At first only for Masscut_140and Rebin10
  vector<TH1F*> ordered_fsr = {FSRdown_rebin10[2], FSRdown_rebin10[1], FSRdown_rebin10[0], ttbar_all[1], FSRup_rebin10[0], FSRup_rebin10[1], FSRup_rebin10[2]};
  vector<TH1F*> ordered_fsr_norm = {FSRdown_rebin10_norm[2], FSRdown_rebin10_norm[1], FSRdown_rebin10_norm[0], ttbar_all_norm[1], FSRup_rebin10_norm[0], FSRup_rebin10_norm[1], FSRup_rebin10_norm[2]};
  vector<TString> tau_bin = {"0 < #tau_{32} < 0.1", "0.1 < #tau_{32} < 0.2", "0.2 < #tau_{32} < 0.3", "0.3 < #tau_{32} < 0.4", "0.4 < #tau_{32} < 0.5", "0.5 < #tau_{32} < 0.6", "0.6 < #tau_{32} < 0.7", "0.7 < #tau_{32} < 0.8", "0.8 < #tau_{32} < 0.9", "0.9 < #tau_{32} < 1"};
  vector<TString> number_bin = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};

  vector<vector<float>> ordered_fsr_norm_err;
  // first index: factor (order: 1/4, 1/2, 1/sqrt2, 1, sqrt2, 2, 4) --- second index: bin
  for(unsigned int factor=0; factor<n_factors; factor++) ordered_fsr_norm_err.push_back(normalize_error(ordered_fsr[factor]));
  vector<float> mu, events, events_norm;
  vector<float> mu_err, events_err, events_norm_err;
  vector<float> events_rel_err, events_norm_rel_err;

  mu = {0.25, 0.5, 1/sqrt(2), 1, sqrt(2), 2, 4};
  mu_err = {0, 0, 0, 0, 0, 0, 0};

  for(unsigned int isnorm=0; isnorm<2; isnorm++){ // isnorm=0: normal --- isnorm=1: norm
    for(unsigned int islog=0; islog<2; islog++){ // islog=0: normal --- islog=1: log
      if(debug) cout << "Star: filling axis" << endl;

      for(unsigned int bin=1; bin < 11; bin++){
        for(Int_t factor=0; factor<n_factors; factor++){
          events.push_back(ordered_fsr[factor]->GetBinContent(bin));
          events_err.push_back(ordered_fsr[factor]->GetBinError(bin));
          events_norm.push_back(ordered_fsr_norm[factor]->GetBinContent(bin));
          events_norm_err.push_back(ordered_fsr_norm_err[factor].at(bin-1));
        }

        if(rel_error){
          events_rel_err = relative_error(events, events_err);
          events_norm_rel_err = relative_error(events_norm, events_norm_err);
        }
        // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
        if(debug) cout << "Start: TGraphError - Single - bin " << bin << endl;
        TGraphErrors *single_bin;

        TCanvas *A = new TCanvas("A","A", 600, 600); // each Graph seperatly
        if(isnorm==0) single_bin = new TGraphErrors(n_factors, &mu[0], &events[0], &mu_err[0], &events_err[0]);
        if(isnorm==1) single_bin = new TGraphErrors(n_factors, &mu[0], &events_norm[0], &mu_err[0], &events_norm_err[0]);

        if(rel_error){
          if(isnorm==0) single_bin = new TGraphErrors(n_factors, &mu[0], &events[0], &mu_err[0], &events_rel_err[0]);
          if(isnorm==1) single_bin = new TGraphErrors(n_factors, &mu[0], &events_norm[0], &mu_err[0], &events_norm_rel_err[0]);
        }

        single_bin->SetTitle(tau_bin[bin-1]);
        if(islog==0) single_bin->GetXaxis()->SetTitle("#mu");
        if(islog==1){
          single_bin->GetXaxis()->SetTitle("log(#mu)");
          gPad->SetLogx();
        }
        single_bin->SetMarkerStyle(2);
        single_bin->Draw("APE");
        if(isnorm==0 && islog==0) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/single_bins/tau32_bin_"+number_bin[bin-1]+".pdf");
        if(isnorm==1 && islog==0) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/single_bins/tau32_bin_"+number_bin[bin-1]+"_norm.pdf");
        if(isnorm==0 && islog==1) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/single_bins/tau32_bin_"+number_bin[bin-1]+"_log.pdf");
        if(isnorm==1 && islog==1) A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/single_bins/tau32_bin_"+number_bin[bin-1]+"_norm_log.pdf");
        delete A;
        events = {}; events_err = {}; events_norm = {}; events_norm_err = {}; // clear vectors
      }

      // At this point only for-loop with norm
      if(debug) cout << "Start: TGraphError - all in one" << endl;

      TCanvas *E = new TCanvas(); // All Graphs in one
      E->Divide(3,4);
      E->SetCanvasSize(800, 1200);
      E->SetWindowSize(800, 1200);

      for(unsigned int bin=1; bin < 11; bin++){
        E->cd(bin);

        if(debug) cout << "Start: TGraphError - all in one - bin " << bin << endl;
        for(int factor=0; factor<n_factors; factor++){
          events.push_back(ordered_fsr[factor]->GetBinContent(bin));
          events_err.push_back(ordered_fsr[factor]->GetBinError(bin));
          events_norm.push_back(ordered_fsr_norm[factor]->GetBinContent(bin));
          events_norm_err.push_back(ordered_fsr_norm_err[factor].at(bin-1));
        }

        if(rel_error){
          events_rel_err = relative_error(events, events_err);
          events_norm_rel_err = relative_error(events_norm, events_norm_err);
        }

        TGraphErrors *all_bins;
        if(isnorm==0) all_bins = new TGraphErrors(n_factors, &mu[0], &events[0], &mu_err[0], &events_err[0]);
        if(isnorm==1) all_bins = new TGraphErrors(n_factors, &mu[0], &events_norm[0], &mu_err[0], &events_norm_err[0]);

        if(rel_error){
          if(isnorm==0) all_bins = new TGraphErrors(n_factors, &mu[0], &events[0], &mu_err[0], &events_rel_err[0]);
          if(isnorm==1) all_bins = new TGraphErrors(n_factors, &mu[0], &events_norm[0], &mu_err[0], &events_norm_rel_err[0]);
        }

        all_bins->SetTitle(tau_bin[bin-1]);
        all_bins->SetMarkerStyle(2);
        all_bins->GetXaxis()->SetTitle("#mu");
        if(islog==1){
          all_bins->GetXaxis()->SetTitle("log(#mu)");
          gPad->SetLogx();
        }
        all_bins->Draw("AP");

        events = {}; events_err = {}; events_norm = {}; events_norm_err = {}; // clear vectors
      }

      if(isnorm==0 && islog==0) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/tau_bin_all.pdf");
      if(isnorm==1 && islog==0) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/tau_bin_all_norm.pdf");
      if(isnorm==0 && islog==1) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/tau_bin_all_log.pdf");
      if(isnorm==1 && islog==1) E->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/"+year+folder_rel_err+"/tau_bin_all_norm_log.pdf");
    }
  }
}
