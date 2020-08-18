#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  bool debug = false;
  TString w_mass = "comparison_topjet_xcone_pass_rec/wmass_min";

  // // #################################################################################################
  // // Get Background ##################################################################################
  // if(debug) cout << "Background" << endl;
  // vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
  // for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+path_bkg_v[i]));
  // for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
  // TH1F *bkg = AddHists(hists_bkg_v, 1);

  // TString dir  = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2017/muon/";
  TString dir  = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2017/muon/";
  TString save = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS";

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar_Mtt0000to0700_SemiLep_2017v2.root ";
  TFile *ttbar_file  = new TFile(dir+ttbar_path);
  TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);

  // ------------------------------------------------------------------------------------------------
  if(debug) cout << "XCone" << endl;
  TFile *XCup_file   = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_Mtt0000to0700_SemiLep_XConeup_2017v2.root ");
  TFile *XCdown_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_Mtt0000to0700_SemiLep_XConedown_2017v2.root ");
  TH1F  *XCup        = (TH1F*)XCup_file->Get(w_mass);
  TH1F  *XCdown      = (TH1F*)XCdown_file->Get(w_mass);


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


  ttbar->SetTitle("");
  ttbar->GetXaxis()->SetRangeUser(0, 180);
  ttbar->GetYaxis()->SetRangeUser(0, ttbar->GetMaximum()*1.2);
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



  // Legend ------------------------------------------------------------------------------------------
  TLegend *leg;

  /*
  ███    ███  █████  ███████ ███████     ██████  ██       ██████  ████████ ███████
  ████  ████ ██   ██ ██      ██          ██   ██ ██      ██    ██    ██    ██
  ██ ████ ██ ███████ ███████ ███████     ██████  ██      ██    ██    ██    ███████
  ██  ██  ██ ██   ██      ██      ██     ██      ██      ██    ██    ██         ██
  ██      ██ ██   ██ ███████ ███████     ██      ███████  ██████     ██    ███████
  */

  TLine *line = new TLine(0, 200, 180, 200);
  line->SetLineColor(kGray);
  line->SetLineWidth(1);

  TString norm_str="";
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  /*
  Plotting the histograms included above.
  */
  if(debug) cout << "Mass Plots" << endl;

  XCup->SetLineWidth(2);
  XCup->SetLineColor(kBlue);

  XCdown->SetLineStyle(2);
  XCdown->SetLineColor(kBlue);

  TCanvas *A = new TCanvas("A", "A", 1000, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar->Draw("HIST");
  XCdown->Draw("SAME HIST");
  XCup->Draw("SAME HIST");
  leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->SetTextSize(0.2);
  leg->AddEntry(ttbar,"nominal","l");
  leg->AddEntry(XCup,"XConeup","l");
  leg->AddEntry(XCdown,"XConedown","l");
  leg->SetTextSize(0.05);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save+"/test_had_x.pdf");
  delete A;
  leg->Clear();
}
