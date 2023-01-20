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
  bool debug   = false;
  TString save_path = get_save_path();
  TString year = "2017";

  print_seperater();
  cout.precision(6);

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString hist_class = "XCone_cor/";
  TString hist = hist_class+"M_jet1_";

  // // #################################################################################################
  // // Get Background ##################################################################################
  // if(debug) cout << "Background" << endl;
  //
  // vector<TFile*> file_bkg_v;
  // vector<TH1F*> hists_bkg_v;
  // vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
  // for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+year+"/muon/"+path_bkg_v[i]));
  // for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(hist));
  // TH1F *bkg = AddHists(hists_bkg_v, 1);

  // #################################################################################################
  // Get Data ########################################################################################
  if(debug) cout << "Data" << endl;

  TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file      = new TFile(dir+year+"/muon/"+data_path);
  TH1F  *data           = (TH1F*)data_file->Get(hist);


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
  TH1F  *ttbar_nom       = (TH1F*)ttbar_nom_file->Get(hist);

  TString ttbar_up_path = "FSRup_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
  TFile *ttbar_up_file  = new TFile(dir+year+"/muon/"+ttbar_up_path);
  TH1F  *ttbar_up       = (TH1F*)ttbar_up_file->Get(hist);

  TString ttbar_down_path = "FSRdown_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
  TFile *ttbar_down_file  = new TFile(dir+year+"/muon/"+ttbar_down_path);
  TH1F  *ttbar_down       = (TH1F*)ttbar_down_file->Get(hist);

  cout << "Nom: " << ttbar_nom->Integral() << endl;
  cout << "Up4: " << ttbar_up->Integral() << endl;
  cout << "Do4: " << ttbar_down->Integral() << endl;

  ttbar_nom->SetTitle("");
  ttbar_nom->GetXaxis()->SetRangeUser(0, 400);
  ttbar_nom->GetYaxis()->SetRangeUser(0, ttbar_nom->GetMaximum()*1.2);
  ttbar_nom->GetXaxis()->SetNdivisions(505);
  ttbar_nom->GetYaxis()->SetNdivisions(505);
  ttbar_nom->GetXaxis()->SetTitleSize(0.05);
  ttbar_nom->GetYaxis()->SetTitleSize(0.04);
  ttbar_nom->GetXaxis()->SetTitleOffset(0.9);
  ttbar_nom->GetYaxis()->SetTitleOffset(1.5);
  ttbar_nom->GetXaxis()->SetTitle("m_{jet} [GeV]");
  ttbar_nom->GetYaxis()->SetTitle("Events");

  ttbar_nom->SetLineColor(kRed);
  ttbar_up->SetLineColor(kGreen+2);
  ttbar_down->SetLineColor(kBlue);

  ttbar_nom->SetLineWidth(2);
  ttbar_up->SetLineWidth(2);
  ttbar_down->SetLineWidth(2);

  TCanvas *A= new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLegend *leg = new TLegend();
  ttbar_nom->Draw("HIST");
  ttbar_up->Draw("same hist");
  ttbar_down->Draw("same hist");
  A->SaveAs(save_path+"/Plots/ComparisonMass/FSRup4_2017.pdf");

}
