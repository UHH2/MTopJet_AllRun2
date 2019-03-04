#include "../include/CentralInclude.h"


using namespace std;


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TFile*> files;
  files.push_back(new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1715.root"));
  files.push_back(new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  files.push_back(new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1735.root"));
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<vector<TH1F*>> hist;

  vector<TString> histnames;
  histnames.push_back("XCone_cor_subjets/pt_had_subjet1");
  histnames.push_back("XCone_cor_subjets/pt_had_subjet2");
  histnames.push_back("XCone_cor_subjets/pt_had_subjet3");
  histnames.push_back("MuonHists/pt_1");
  histnames.push_back("XCone_cor/M_jet1_");
  histnames.push_back("XCone_cor/M_jet1");
  histnames.push_back("XCone_cor/M_jet2");
  histnames.push_back("XCone_cor/pt_jet1");
  histnames.push_back("XCone_cor/pt_jet2");
  histnames.push_back("EventHists/MET");
  histnames.push_back("GenParticles_GenOnly/hadtop_pt");
  histnames.push_back("GenParticles_GenOnly/leptop_pt");
  histnames.push_back("GenParticles_GenOnly/hadtop_mass");
  histnames.push_back("GenParticles_GenOnly/leptop_mass");

  for(auto name: histnames){
    vector<TH1F*> h;
    for(auto file: files){
      h.push_back((TH1F*)file->Get( name) );
    }
    hist.push_back(h);
  }


  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<hist.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    hist[i][0]->SetLineColor(kRed+1);
    hist[i][1]->SetLineColor(13);
    hist[i][2]->SetLineColor(kAzure+7);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    for(unsigned int j=0; j<hist[i].size(); j++){
      hist[i][j]->SetLineWidth(4);
    }
    hist[i][0]->GetYaxis()->SetRangeUser(0, 1.3 * hist[i][0]->GetMaximum());
    TString title = hist[i][0]->GetTitle();
    hist[i][0]->GetXaxis()->SetTitle(title);
    hist[i][0]->SetTitle(" ");
    hist[i][0]->Draw("HIST");
    hist[i][1]->Draw("HIST SAME");
    hist[i][2]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.63,0.6,0.88,0.88);
    leg->AddEntry(hist[i][0], "m_{top} = 171.5", "l");
    leg->AddEntry(hist[i][1], "m_{top} = 172.5", "l");
    leg->AddEntry(hist[i][2], "m_{top} = 173.5", "l");
    leg->Draw();
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/DifferencesMassSamples/"+histnames[i]+".pdf");
    delete a;
  }

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
