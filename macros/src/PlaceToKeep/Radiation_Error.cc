#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(int argc, char* argv[]){
  bool debug = false;
  bool with_bkg = true;

  if(argc != 2){
    /* <Background True/False> : True  - (ttbar + bkg) & False - (data - bkg) */
    cout << "Usage: ./Radiation_SYS <year>" << endl;
    return 0;
  }

  if(debug) cout << "Getting Input Variables" << endl;
  TString year = argv[1];

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
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar = (TH1F*)ttbar_file->Get(tau32);
  ttbar->Add(bkg, 1);

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

  vector<TH1F*> FSRup_rebin10 = rebin(FSRup, 5);
  vector<TH1F*> FSRdown_rebin10 = rebin(FSRdown, 5);

  // For tau bin dependency
  vector<TH1F*> FSRup_rebin10_norm = normalize(FSRup_rebin10);
  vector<TH1F*> FSRdown_rebin10_norm = normalize(FSRdown_rebin10);

  // #################################################################################################
  // Normalize Vectors ###############################################################################
  if(debug) cout << "Normalized Vector" << endl;

  vector<TH1F*> ttbar_all = {ttbar, ttbar_rebin10};
  // Vectors have the same order like the ones above
  vector<TH1F*> ttbar_all_norm = normalize(ttbar_all);



  /*
  ███████ ██████  ██████   ██████  ██████
  ██      ██   ██ ██   ██ ██    ██ ██   ██
  █████   ██████  ██████  ██    ██ ██████
  ██      ██   ██ ██   ██ ██    ██ ██   ██
  ███████ ██   ██ ██   ██  ██████  ██   ██
  */

  /*
  OVERVIEW ---------------------------------------------------------------------
  xxx - ttbar; data; FSRxx;
  FSRxx - FSRup; FSRdown;
  xxx_general = {xxx_all, xxx_all_norm}
  xxx_all     = {xxx, xxx_rebin10}
  FSRxx       = {FSRxx_sqrt2, FSRxx_2, FSRxx_4}
  */


  // REBIN in 10 Bins and MASSCUT 140
  vector<TH1F*> ordered_fsr = {FSRdown_rebin10[2], FSRdown_rebin10[1], FSRdown_rebin10[0], ttbar_all[1], FSRup_rebin10[0], FSRup_rebin10[1], FSRup_rebin10[2]};
  vector<TH1F*> ordered_fsr_norm = {FSRdown_rebin10_norm[2], FSRdown_rebin10_norm[1], FSRdown_rebin10_norm[0], ttbar_all_norm[1], FSRup_rebin10_norm[0], FSRup_rebin10_norm[1], FSRup_rebin10_norm[2]};

  cout << "-------------- Not normalized ------------------------" << endl;
  auto [content, error] = get_bin_content_and_error(ordered_fsr[2]);
  vector<double> rel_error = relative_error(content, error);
  for(unsigned int i=0; i<rel_error.size(); i++) cout << "Bin " << i+1 << " --- content: " << content[i] << " | error: " << error[i] << " | rel error: " << error[i]/content[i] << " | rel error function: " << rel_error[i] << endl;

  cout << "-------------- normalized ----------------------------" << endl;
  vector<double> content_norm = get_bin_content(ordered_fsr_norm[2]);
  vector<double> error_norm = normalize_error(ordered_fsr[2]);
  vector<double> rel_error_norm = relative_error(content_norm, error_norm);
  for(unsigned int i=0; i<rel_error_norm.size(); i++) cout << "Bin " << i+1 << " --- content: " << content_norm[i] << " | error: " << error_norm[i] << " | rel error: " << error_norm[i]/content_norm[i] << " | rel error function: " << rel_error_norm[i] << endl;

  cout << "-------------- Comparison ----------------------------" << endl;
  for(unsigned int i=0; i<content_norm.size(); i++) cout << "Bin " << i+1 << " --- normal: " << rel_error[i] << " | normalized: " << rel_error_norm[i] << endl;
}
