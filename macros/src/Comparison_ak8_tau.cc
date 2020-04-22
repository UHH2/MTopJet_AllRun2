#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){
  vector<TString> year = {"2017", "2018"};

  TString ttbar_string = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString data_string = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/Comparison_ak8/";

  vector<TString> masscut = {"no_masscut", "masscut_120", "masscut_130", "masscut_140", "masscut_150"};
  vector<TString> classes = {"AK8_tau_pass_rec", "AK8_tau_pass_rec_masscut_120", "AK8_tau_pass_rec_masscut_130", "AK8_tau_pass_rec_masscut_140", "AK8_tau_pass_rec_masscut_150"};
  vector<TH1F*> hists;

  TFile *file_ttbar;

  // Main-Option for Plots
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  TLegend *leg;
  // TCanvas *A = new TCanvas("A", "A", 600, 600);

  // At first only ttbar
  for(unsigned int i=0; i<year.size(); i++){

    file_ttbar = new TFile(dir+year[i]+"/muon/"+ttbar_string); // file for each year
    for(unsigned int j=0; j<classes.size(); j++) hists.push_back((TH1F*)file_ttbar->Get(classes[j]+"/tau32_hadjet")); // push_back Histogramm for N-subjettiness for each masscut

    for(unsigned int j=0; j<hists.size(); j++){
      if(j==0) hists[j]->SetTitle(masscut[j]);
      if(j>0)  hists[j]->SetTitle(masscut[j]+" GeV"); // The vector masscut should have the same length like hists/classes
      hists[j]->GetXaxis()->SetRangeUser(0, 400);
      hists[j]->GetYaxis()->SetRangeUser(0, hists[j]->GetMaximum()*1.2);
      hists[j]->GetXaxis()->SetNdivisions(505);
      hists[j]->GetYaxis()->SetNdivisions(505);
      hists[j]->GetXaxis()->SetTitleSize(0.05);
      hists[j]->GetYaxis()->SetTitleSize(0.04);
      hists[j]->GetXaxis()->SetTitleOffset(0.9);
      hists[j]->GetYaxis()->SetTitleOffset(1.5);
      hists[j]->GetXaxis()->SetTitle("#tau_{32}");
      hists[j]->GetYaxis()->SetTitle("");
      hists[j]->SetLineWidth(2);
      hists[j]->SetLineColor(kRed);

      TCanvas *A = new TCanvas("A", "A", 600, 600);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      hists[j]->Draw("HIST");
      A->SaveAs(save_path+year[i]+"/tau32_"+masscut[j]+"_"+year[i]+".pdf");
      delete A;



      /////////////////////////////////////////////////////////////////////////

      if(j==0){ // Comparison Hists
        hists[j]->SetTitle("");
        for(unsigned int k=1;k<hists.size();k++){
          hists[k]->SetLineWidth(2);
          hists[k]->SetLineColor(kBlue);

          TCanvas *A = new TCanvas("A", "A", 600, 600);
          hists[0]->Draw("HIST");
          hists[k]->Draw("SAME HIST");
          leg = new TLegend(0.25,0.65,0.45,0.85);
          leg->AddEntry(hists[0],"no masscut","l");
          leg->AddEntry(hists[k],masscut[k],"l");
          leg->SetTextSize(0.05);
          leg->Draw();
          gPad->RedrawAxis();
          A->SaveAs(save_path+year[i]+"/Comparison_tau32_"+masscut[k]+"_"+year[i]+".pdf");
          leg->Clear();
          delete A;
        }
      }
    }
    hists = {}; // Empty Hist at the end of the loop to avoid Segmentation Violation
  }
}



// later for ratio
// "tau32_hadjet_fullymerged", "tau32_hadjetsemimerged", "tau32_hadjet_notmerged"
