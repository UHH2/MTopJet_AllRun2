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

  vector<TString> years={"2016", "2017", "2018", "combined"};
  cout << '\n'; // CHANGE_PT
  TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
  w_mass = hist_class+"ak4_distance_xcone_btag";

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

    TString data_path = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    TFile  *data_file = new TFile(dir+year+"/muon/"+data_path);
    TH1F   *data      = (TH1F*)data_file->Get(w_mass);


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

    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile  *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
    TH1F   *ttbar      = (TH1F*)ttbar_file->Get(w_mass);


    /*
    ███    ██ ██    ██ ███    ███ ██████  ███████ ██████  ███████
    ████   ██ ██    ██ ████  ████ ██   ██ ██      ██   ██ ██
    ██ ██  ██ ██    ██ ██ ████ ██ ██████  █████   ██████  ███████
    ██  ██ ██ ██    ██ ██  ██  ██ ██   ██ ██      ██   ██      ██
    ██   ████  ██████  ██      ██ ██████  ███████ ██   ██ ███████
    */


    double inside_subjet_tt   = ttbar->Integral(1, 20);
    double outside_subjet_tt  = ttbar->Integral(21, 51);
    double inside_subjet_bkg  = bkg->Integral(1, 20);
    double outside_subjet_bkg = bkg->Integral(21, 51);
    double inside_subjet_da   = data->Integral(1, 20);
    double outside_subjet_da  = data->Integral(21, 51);

    double ratio_tt  = inside_subjet_tt/(inside_subjet_tt+outside_subjet_tt);
    double ratio_da  = inside_subjet_da/(inside_subjet_da+outside_subjet_da);
    double ratio_all = (inside_subjet_tt+inside_subjet_bkg)/(inside_subjet_tt+outside_subjet_tt+inside_subjet_bkg+outside_subjet_bkg);

    cout << setw(9) << year << ": Data - " << ratio_da << " | ttbar - " << ratio_tt << " | tt+bkg - " << ratio_all << endl;
  }

}
