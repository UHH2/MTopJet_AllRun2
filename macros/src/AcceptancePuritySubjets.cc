#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  // declare files
  TFile *f_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_fractions.root");

  // set binning
  int n_bins = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};

  TH1F * acceptance = new TH1F("acceptance"," ",n_bins, bins);
  TH1F * purity = new TH1F("purity"," ",n_bins, bins);

  // now get plots from file
  vector<TH1F*> MatchGen, MatchRec, TotalGen, TotalRec;

  std::string dir = "RecGenHists_subjets_noJEC/";

  std::string name_matchGen;
  std::string name_matchRec;
  std::string name_TotalGen;
  std::string name_TotalRec;

  for(int ptbin = 1; ptbin <= n_bins; ptbin++){
    name_matchGen = dir + "MatchedGen_" + std::to_string(ptbin);
    name_matchRec = dir + "MatchedRec_" + std::to_string(ptbin);
    name_TotalGen = dir + "TotalGen_" + std::to_string(ptbin);
    name_TotalRec = dir + "TotalRec_" + std::to_string(ptbin);

    MatchGen.push_back( (TH1F*)f_file->Get(name_matchGen.c_str()) );
    MatchRec.push_back( (TH1F*)f_file->Get(name_matchRec.c_str()) );
    TotalGen.push_back( (TH1F*)f_file->Get(name_TotalGen.c_str()) );
    TotalRec.push_back( (TH1F*)f_file->Get(name_TotalRec.c_str()) );
  }

  double number_recjets = 0;
  for(unsigned int i=0; i<TotalRec.size(); i++){
    int bin=i+1;
    number_recjets += TotalRec[i]->GetBinContent(1);
    double ratio = MatchRec[i]->GetBinContent(1) / TotalRec[i]->GetBinContent(1);
    // double error_1 = sqrt(MatchRec[i]->GetBinContent(1)) / TotalRec[i]->GetBinContent(1);
    // double error_2 = sqrt(TotalRec[i]->GetBinContent(1)) * ratio / TotalRec[i]->GetBinContent(1);
    // double error = error_1 + error_2;
    double error = (1-ratio)/ratio;
    purity->SetBinContent(bin,ratio);
    purity->SetBinError(bin, error);
  }
  double number_genjets = 0;
  for(unsigned int i=0; i<TotalGen.size(); i++){
    int bin=i+1;
    number_genjets += TotalGen[i]->GetBinContent(1);
    double ratio = MatchGen[i]->GetBinContent(1) / TotalGen[i]->GetBinContent(1);
    // double error_1 = sqrt(MatchGen[i]->GetBinContent(1)) / TotalGen[i]->GetBinContent(1);
    // double error_2 = sqrt(TotalGen[i]->GetBinContent(1)) * ratio / TotalGen[i]->GetBinContent(1);
    // double error = error_1 + error_2;
    double error = (1-ratio)/ratio;
    acceptance->SetBinContent(bin,ratio);
    acceptance->SetBinError(bin, error);
  }

  cout << "number of rec jets: " << number_recjets << endl;
  cout << "number of gen jets: " << number_genjets << endl;
  cout << "max acceptance = " << acceptance->GetMaximum() << endl;
  cout << "max purity     = " <<     purity->GetMaximum() << endl;

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // start plotting
  TLine *line = new TLine(0, 1, 600, 1);
  line->SetLineColor(15);
  line->SetLineWidth(3);
  line->SetLineStyle(7);

  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  purity->GetYaxis()->SetRangeUser(0.9, 1.05);
  purity->GetXaxis()->SetNdivisions(505);
  purity->GetYaxis()->SetNdivisions(505);
  purity->GetXaxis()->SetTitleSize(0.05);
  purity->GetYaxis()->SetTitleSize(0.04);
  purity->GetXaxis()->SetTitleOffset(0.9);
  purity->GetYaxis()->SetTitleOffset(1.5);
  purity->GetXaxis()->SetTitle("p_{T}^{rec} or p_{T}^{gen}");
  purity->GetYaxis()->SetTitle("#frac{matched}{total}");
  purity->SetLineColor(kBlack);
  purity->SetMarkerColor(kBlack);
  purity->SetMarkerStyle(8);
  purity->SetMarkerSize(1);

  acceptance->GetYaxis()->SetRangeUser(0.9, 1.05);
  acceptance->GetXaxis()->SetNdivisions(505);
  acceptance->GetYaxis()->SetNdivisions(505);
  acceptance->GetXaxis()->SetTitleSize(0.05);
  acceptance->GetYaxis()->SetTitleSize(0.04);
  acceptance->GetXaxis()->SetTitleOffset(0.9);
  acceptance->GetYaxis()->SetTitleOffset(1.5);
  acceptance->GetXaxis()->SetTitle("p_{T}^{rec} or p_{T}^{gen}");
  acceptance->GetYaxis()->SetTitle("#frac{matched}{total}");
  acceptance->SetLineColor(kRed);
  acceptance->SetMarkerColor(kRed);
  acceptance->SetMarkerStyle(21);
  acceptance->SetMarkerSize(1);

  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------

  TCanvas *c0 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  purity->Draw("E1");
  acceptance->Draw("SAME E1");
  line->Draw("SAME");
  TLegend *leg0 = new TLegend(0.45,0.85,0.80,0.65);
  leg0->AddEntry(acceptance,"acceptance (vs p_{T}^{gen})","l");
  leg0->AddEntry(purity,"purity (vs p_{T}^{rec})","l");
  leg0->SetTextSize(0.05);
  leg0->Draw();
  TLatex text0;
  text0.SetNDC(kTRUE);
  text0.SetTextFont(43);
  text0.SetTextSize(18);
  text0.DrawLatex(.3,.2, "#DeltaR(gen,reco)<0.2");
  gPad->RedrawAxis();
  c0->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/AcceptancePurity/AcceptancePurity.pdf");


  return 0;
}
