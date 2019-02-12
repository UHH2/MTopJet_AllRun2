#include "../include/CentralInclude.h"


using namespace std;


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TString> histnames = {"pt_subjets_had1", "pt_subjets_had2", "pt_subjets_had3", "pt_subjets_had"};
  vector<TString> filenames = {"pt1_gen", "pt2_gen", "pt3_gen", "pt_gen"};
  vector<TH1F*> gen;

  for(auto name: histnames){
    gen.push_back( (TH1F*)TT_file->Get( "XCone_GEN_GenOnly/"+name) );
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  for(unsigned int i=0; i<gen.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    gen[i]->SetTitle(" ");
    gen[i]->GetXaxis()->SetTitle("subjet p_{T}^{gen}");
    gen[i]->GetYaxis()->SetTitle("events");
    gen[i]->GetYaxis()->SetTitleSize(0.06);
    gen[i]->GetXaxis()->SetTitleSize(0.05);
    gen[i]->GetZaxis()->SetTitleSize(0.05);
    gen[i]->GetXaxis()->SetTitleOffset(0.9);
    gen[i]->GetYaxis()->SetTitleOffset(1.1);
    gen[i]->GetZaxis()->SetTitleOffset(0.9);
    gen[i]->GetXaxis()->SetNdivisions(505);
    gen[i]->GetYaxis()->SetNdivisions(505);
    gen[i]->SetFillColor(kRed);
    gen[i]->SetLineColor(1);
    gen[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/SubjetPT/"+filenames[i]+".pdf");
    delete a;
  }

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
