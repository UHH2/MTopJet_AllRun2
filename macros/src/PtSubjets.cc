#include "../include/CentralInclude.h"


using namespace std;


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TString> histnames = {"pt_had_subjet1", "pt_had_subjet2", "pt_had_subjet3", "pt_had_subjets"};
  vector<TString> filenames = {"pt1_rec", "pt2_rec", "pt3_rec", "pt_rec"};
  vector<TH1F*> rec;

  for(auto name: histnames){
    rec.push_back( (TH1F*)TT_file->Get( "XCone_jec_subjets/"+name) );
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  for(unsigned int i=0; i<rec.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    rec[i]->SetTitle(" ");
    rec[i]->GetXaxis()->SetTitle("subjet p_{T}^{rec}");
    rec[i]->GetYaxis()->SetTitle("events");
    rec[i]->GetYaxis()->SetTitleSize(0.06);
    rec[i]->GetXaxis()->SetTitleSize(0.05);
    rec[i]->GetZaxis()->SetTitleSize(0.05);
    rec[i]->GetXaxis()->SetTitleOffset(0.9);
    rec[i]->GetYaxis()->SetTitleOffset(1.1);
    rec[i]->GetZaxis()->SetTitleOffset(0.9);
    rec[i]->GetXaxis()->SetNdivisions(505);
    rec[i]->GetYaxis()->SetNdivisions(505);
    rec[i]->SetFillColor(kRed);
    rec[i]->SetLineColor(1);
    rec[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/SubjetPT/"+filenames[i]+".pdf");
    delete a;
  }

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
