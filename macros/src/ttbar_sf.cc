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
  bool debug = false;

  print_seperater();
  cout.precision(6);

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  vector<TString> years={"2016", "2017", "2018"};
  cout << '\n'; // CHANGE_PT
  TString w_mass; TString hist_class = "XCone_cor/";
  w_mass = hist_class+"M_jet1_";

  for(auto year: years){
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
    TH1F  *data           = (TH1F*)data_file->Get(w_mass);


    /*
    .████████ ████████ ██████   █████  ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ███████ ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ██   ██ ██   ██
    */

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_nom_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile *ttbar_nom_file  = new TFile(dir+year+"/muon/"+ttbar_nom_path);
    TH1F  *ttbar_nom       = (TH1F*)ttbar_nom_file->Get(w_mass);


    /*
    .███████ ███████
    .██      ██
    .███████ █████
    .     ██ ██
    .███████ ██
    */

    double number_ttbar       = ttbar_nom->Integral();
    double number_mc          = bkg->Integral();
    double number_simulations = number_ttbar+number_mc;
    double number_data        = data->Integral();
    double ratio              = number_data/number_simulations;
    double ratio_tt           = number_ttbar/number_simulations;
    double ratio_other        = (number_data-number_mc)/number_ttbar;
    double contribution       = number_ttbar/(number_ttbar+number_mc);


    cout << "Year - " << year << " | SF-Minus- " << ratio_other << " | tt contribution - " << contribution << endl;
  }

  // #################################################################################################
  // #################################################################################################
  // #################################################################################################
  // #################################################################################################

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  hist_class = "PreSel04_Muon/";
  w_mass = hist_class+"charge";

  for(auto year: years){
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
    TH1F  *data           = (TH1F*)data_file->Get(w_mass);


    /*
    .████████ ████████ ██████   █████  ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ███████ ██████
    .   ██       ██    ██   ██ ██   ██ ██   ██
    .   ██       ██    ██████  ██   ██ ██   ██
    */

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_nom_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile *ttbar_nom_file  = new TFile(dir+year+"/muon/"+ttbar_nom_path);
    TH1F  *ttbar_nom       = (TH1F*)ttbar_nom_file->Get(w_mass);

    // #################################################################################################
    // Get Numbers #####################################################################################
    double number_ttbar       = ttbar_nom->Integral();
    double number_mc          = bkg->Integral();
    double number_simulations = number_ttbar+number_mc;
    double number_data        = data->Integral();
    double ratio              = number_data/number_simulations;
    double ratio_tt           = number_ttbar/number_simulations;
    double ratio_other        = (number_data-number_mc)/number_ttbar;
    double contribution       = number_ttbar/(number_ttbar+number_mc);

    cout << "Year - " << year << " | SF-Minus- " << ratio_tt << endl;

  }
}
