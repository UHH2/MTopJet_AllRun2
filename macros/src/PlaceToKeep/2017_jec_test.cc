#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  // #################################################################################################
  // Declare different variables used in the code ####################################################
  bool debug = true;
  if(debug) cout << "start" << endl;


  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */
  TString dir = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2017/muon";
  if(debug) cout << "Get Classes" << endl;
  cout << '\n';
  // TString hist_class = "03_Lepton_Muon/";
  TString hist_class = "comparison_topjet_xcone_pass_rec/";
  TString w_mass = hist_class+"wmass_match";                    // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet

  TString year = "2017";
  if(debug) cout << "TTbar" << endl;

  // TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file  = new TFile(dir+"/"+ttbar_path);
  TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);



  // #################################################################################################
  // Get SYS #########################################################################################
  if(debug) cout << "JEC" << endl;

  TFile *JECup_file   = new TFile(dir+"/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");
  TFile *JECdown_file = new TFile(dir+"/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
  TH1F  *JECup        = (TH1F*)JECup_file->Get(w_mass);
  TH1F  *JECdown      = (TH1F*)JECdown_file->Get(w_mass);

  // ------------------------------------------------------------------------------------------------
  // if(debug) cout << "XCone" << endl;
  // TFile *XConeup_file   = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConeup.root");
  // TFile *XConedown_file = new TFile(dir+year+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConedown.root");
  // TH1F  *XConeup        = (TH1F*)XConeup_file->Get(w_mass);
  // TH1F  *XConedown      = (TH1F*)XConedown_file->Get(w_mass);

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
  /* The vectors are only used here */

  ttbar->SetTitle("");
  if(debug) cout << "Plots Sensitivity" << endl;
  ttbar->GetXaxis()->SetRangeUser(0, 400);
  ttbar->GetYaxis()->SetRangeUser(0, JECup->GetMaximum()*1.2);
  ttbar->GetXaxis()->SetNdivisions(505);
  ttbar->GetYaxis()->SetNdivisions(505);
  ttbar->GetXaxis()->SetTitleSize(0.05);
  ttbar->GetYaxis()->SetTitleSize(0.04);
  ttbar->GetXaxis()->SetTitleOffset(0.9);
  ttbar->GetYaxis()->SetTitleOffset(1.5);
  ttbar->GetXaxis()->SetTitle("m_{Wjet}");
  ttbar->GetYaxis()->SetTitle("");
  ttbar->SetLineWidth(2);  // ttbar hist style
  ttbar->SetLineColor(kRed);
  if(debug) cout << "Plots Sensitivity" << endl;

  JECup->SetLineWidth(2);  // ttbar hist style
  JECup->SetLineColor(kBlue);
  if(debug) cout << "Plots Sensitivity" << endl;

  JECdown->SetLineStyle(2);  // ttbar hist style
  JECdown->SetLineWidth(2);  // ttbar hist style
  JECdown->SetLineColor(kBlue);

  // Legend ------------------------------------------------------------------------------------------
  TLegend *leg;
  if(debug) cout << "Plots Sensitivity" << endl;

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar->Draw("Hist");
  JECup->Draw("same hist");
  JECdown->Draw("same hist");
  A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/test.pdf");
}
