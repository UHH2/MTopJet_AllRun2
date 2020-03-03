#include "../include/CentralInclude.h"

using namespace std;

int main(){

  TString uhh = "2016/muon/uhh2.AnalysisModuleRunner.MC.";
  vector<TString> name = {"nominal", "1695", "1715", "1735", "1755"};
  vector<int> colors = {kRed, kBlue, kGreen, kYellow, kBlack};

  TFile *f_tt = new TFile(dir+uhh+"TTbar.root");
  TFile *f_1695 = new TFile(dir+uhh+"TTbar_mtop1695_2016v3.root");
  TFile *f_1715 = new TFile(dir+uhh+"TTbar_mtop1715_2016v3.root");
  TFile *f_1735 = new TFile(dir+uhh+"TTbar_mtop1735_2016v3.root");
  TFile *f_1755 = new TFile(dir+uhh+"TTbar_mtop1755_2016v3.root");
  vector<TFile*> files = {f_tt, f_1695, f_1715, f_1735, f_1755};

  TString hist_perp = "LeptonicTop/MuonPtPerpFull";
  TString hist_boost = "LeptonicTop/MuonPtBoostFull";

  vector<TH1F*> pt_perp, pt_boost, pt_perp_rebin_5, pt_boost_rebin_5;
  vector<TH1F*> pt_perp_rebin_10, pt_boost_rebin_10;
  TH1F *hist;

  for(unsigned int i=0; i<files.size(); i++){
    // --- Boost ----------------------------
    hist = (TH1F*)files[i]->Get(hist_boost);
    pt_boost.push_back(hist);
    // --- Rebin 5 --------------------------
    hist->Rebin(5);
    pt_boost_rebin_5.push_back(hist);
    // --- Rebin 10 -------------------------
    hist->Rebin();
    pt_boost_rebin_10.push_back(hist);

    // --- Perp -----------------------------
    hist = (TH1F*)files[i]->Get(hist_perp);
    pt_perp.push_back(hist);
    // --- Rebin ----------------------------
    hist->Rebin(5);
    pt_perp_rebin_5.push_back(hist);
    // --- Rebin 10 -------------------------
    hist->Rebin();
    pt_perp_rebin_10.push_back(hist);
  };


  /*
  ██████  ██       ██████  ████████
  ██   ██ ██      ██    ██    ██
  ██████  ██      ██    ██    ██
  ██      ██      ██    ██    ██
  ██      ███████  ██████     ██
  */

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  pt_perp[0]->GetXaxis()->SetRangeUser(0, 100);
  pt_perp[0]->GetYaxis()->SetRangeUser(0, (pt_perp[0]->GetMaximum())*1.2);
  pt_perp[0]->GetXaxis()->SetNdivisions(505);
  pt_perp[0]->GetYaxis()->SetNdivisions(505);
  pt_perp[0]->GetXaxis()->SetTitleSize(0.05);
  pt_perp[0]->GetYaxis()->SetTitleSize(0.04);
  pt_perp[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_perp[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_perp[0]->GetXaxis()->SetTitle("p_{T,#perp}");
  pt_perp[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_perp[i]->SetLineWidth(2);
    pt_perp[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *c = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_perp[0]->Draw("HIST");
  pt_perp[1]->Draw("SAME HIST");
  pt_perp[4]->Draw("SAME HIST");
  TLegend *leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->AddEntry(pt_perp[0],"nominal","l");
  leg->AddEntry(pt_perp[1],"1695","l");
  leg->AddEntry(pt_perp[4],"1755","l");
  leg->SetTextSize(0.05);
  leg->Draw();
  gPad->RedrawAxis();
  c->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_main_2016.pdf");

  pt_perp[2]->Draw("SAME HIST");
  pt_perp[3]->Draw("SAME HIST");
  leg->AddEntry(pt_perp[2],"1715","l");
  leg->AddEntry(pt_perp[3],"1735","l");
  c->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_all_2016.pdf");
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  pt_boost[0]->GetXaxis()->SetRangeUser(0, 1000);
  pt_boost[0]->GetYaxis()->SetRangeUser(0, (pt_boost[4]->GetMaximum())*1.2);
  pt_boost[0]->GetXaxis()->SetNdivisions(505);
  pt_boost[0]->GetYaxis()->SetNdivisions(505);
  pt_boost[0]->GetXaxis()->SetTitleSize(0.05);
  pt_boost[0]->GetYaxis()->SetTitleSize(0.04);
  pt_boost[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_boost[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_boost[0]->GetXaxis()->SetTitle("p_{T}");
  pt_boost[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_boost[i]->SetLineWidth(2);
    pt_boost[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *d = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_boost[0]->Draw("HIST");
  pt_boost[1]->Draw("SAME HIST");
  pt_boost[4]->Draw("SAME HIST");
  TLegend *leg1 = new TLegend(0.65,0.65,0.85,0.85);
  leg1->AddEntry(pt_boost[0],"nominal","l");
  leg1->AddEntry(pt_boost[1],"1695","l");
  leg1->AddEntry(pt_boost[4],"1755","l");
  leg1->SetTextSize(0.05);
  leg1->Draw();
  gPad->RedrawAxis();
  d->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_main_2016.pdf");

  pt_boost[2]->Draw("SAME HIST");
  pt_boost[3]->Draw("SAME HIST");
  leg1->AddEntry(pt_boost[2],"1715","l");
  leg1->AddEntry(pt_boost[3],"1735","l");
  d->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_all_2016.pdf");

  /*
  ██████  ███████ ██████  ██ ███    ██
  ██   ██ ██      ██   ██ ██ ████   ██
  ██████  █████   ██████  ██ ██ ██  ██
  ██   ██ ██      ██   ██ ██ ██  ██ ██
  ██   ██ ███████ ██████  ██ ██   ████
  */
  /*
  ███████
  ██
  ███████
  .....██
  ███████
  */

  pt_perp_rebin_5[0]->GetXaxis()->SetRangeUser(0, 100);
  pt_perp_rebin_5[0]->GetYaxis()->SetRangeUser(0, (pt_perp_rebin_5[0]->GetMaximum())*1.2);
  pt_perp_rebin_5[0]->GetXaxis()->SetNdivisions(505);
  pt_perp_rebin_5[0]->GetYaxis()->SetNdivisions(505);
  pt_perp_rebin_5[0]->GetXaxis()->SetTitleSize(0.05);
  pt_perp_rebin_5[0]->GetYaxis()->SetTitleSize(0.04);
  pt_perp_rebin_5[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_perp_rebin_5[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_perp_rebin_5[0]->GetXaxis()->SetTitle("p_{T,#perp}");
  pt_perp_rebin_5[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_perp_rebin_5[i]->SetLineWidth(2);
    pt_perp_rebin_5[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *e = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_perp_rebin_5[0]->Draw("HIST");
  pt_perp_rebin_5[1]->Draw("SAME HIST");
  pt_perp_rebin_5[4]->Draw("SAME HIST");
  TLegend *leg2 = new TLegend(0.65,0.65,0.85,0.85);
  leg2->AddEntry(pt_perp_rebin_5[0],"nominal","l");
  leg2->AddEntry(pt_perp_rebin_5[1],"1695","l");
  leg2->AddEntry(pt_perp_rebin_5[4],"1755","l");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  gPad->RedrawAxis();
  e->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_main_rebin5_2016.pdf");

  pt_perp_rebin_5[2]->Draw("SAME HIST");
  pt_perp_rebin_5[3]->Draw("SAME HIST");
  leg2->AddEntry(pt_perp_rebin_5[2],"1715","l");
  leg2->AddEntry(pt_perp_rebin_5[3],"1735","l");
  e->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_all_rebin5_2016.pdf");

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  pt_boost_rebin_5[0]->GetXaxis()->SetRangeUser(0, 1000);
  pt_boost_rebin_5[0]->GetYaxis()->SetRangeUser(0, (pt_boost_rebin_5[4]->GetMaximum())*1.2);
  pt_boost_rebin_5[0]->GetXaxis()->SetNdivisions(505);
  pt_boost_rebin_5[0]->GetYaxis()->SetNdivisions(505);
  pt_boost_rebin_5[0]->GetXaxis()->SetTitleSize(0.05);
  pt_boost_rebin_5[0]->GetYaxis()->SetTitleSize(0.04);
  pt_boost_rebin_5[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_boost_rebin_5[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_boost_rebin_5[0]->GetXaxis()->SetTitle("p_{T}");
  pt_boost_rebin_5[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_boost_rebin_5[i]->SetLineWidth(2);
    pt_boost_rebin_5[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *f = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_boost_rebin_5[0]->Draw("HIST");
  pt_boost_rebin_5[1]->Draw("SAME HIST");
  pt_boost_rebin_5[4]->Draw("SAME HIST");
  TLegend *leg3 = new TLegend(0.65,0.65,0.85,0.85);
  leg3->AddEntry(pt_boost_rebin_5[0],"nominal","l");
  leg3->AddEntry(pt_boost_rebin_5[1],"1695","l");
  leg3->AddEntry(pt_boost_rebin_5[4],"1755","l");
  leg3->SetTextSize(0.05);
  leg3->Draw();
  gPad->RedrawAxis();
  f->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_main_rebin5_2016.pdf");

  pt_boost_rebin_5[2]->Draw("SAME HIST");
  pt_boost_rebin_5[3]->Draw("SAME HIST");
  leg3->AddEntry(pt_boost_rebin_5[2],"1715","l");
  leg3->AddEntry(pt_boost_rebin_5[3],"1735","l");
  f->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_all_rebin5_2016.pdf");

  /*
  .██  ██████
  ███ ██  ████
  .██ ██ ██ ██
  .██ ████  ██
  .██  ██████
  */

  pt_perp_rebin_10[0]->GetXaxis()->SetRangeUser(0, 100);
  pt_perp_rebin_10[0]->GetYaxis()->SetRangeUser(0, (pt_perp_rebin_10[0]->GetMaximum())*1.2);
  pt_perp_rebin_10[0]->GetXaxis()->SetNdivisions(505);
  pt_perp_rebin_10[0]->GetYaxis()->SetNdivisions(505);
  pt_perp_rebin_10[0]->GetXaxis()->SetTitleSize(0.05);
  pt_perp_rebin_10[0]->GetYaxis()->SetTitleSize(0.04);
  pt_perp_rebin_10[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_perp_rebin_10[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_perp_rebin_10[0]->GetXaxis()->SetTitle("p_{T,#perp}");
  pt_perp_rebin_10[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_perp_rebin_10[i]->SetLineWidth(2);
    pt_perp_rebin_10[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *g = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_perp_rebin_10[0]->Draw("HIST");
  pt_perp_rebin_10[1]->Draw("SAME HIST");
  pt_perp_rebin_10[4]->Draw("SAME HIST");
  TLegend *leg4 = new TLegend(0.65,0.65,0.85,0.85);
  leg4->AddEntry(pt_perp_rebin_10[0],"nominal","l");
  leg4->AddEntry(pt_perp_rebin_10[1],"1695","l");
  leg4->AddEntry(pt_perp_rebin_10[4],"1755","l");
  leg4->SetTextSize(0.05);
  leg4->Draw();
  gPad->RedrawAxis();
  g->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_main_rebin10_2016.pdf");

  pt_perp_rebin_10[2]->Draw("SAME HIST");
  pt_perp_rebin_10[3]->Draw("SAME HIST");
  leg4->AddEntry(pt_perp_rebin_10[2],"1715","l");
  leg4->AddEntry(pt_perp_rebin_10[3],"1735","l");
  g->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/PerpendicularPt_Comparison_all_rebin10_2016.pdf");

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  pt_boost_rebin_10[0]->GetXaxis()->SetRangeUser(0, 1000);
  pt_boost_rebin_10[0]->GetYaxis()->SetRangeUser(0, (pt_boost_rebin_10[4]->GetMaximum())*1.2);
  pt_boost_rebin_10[0]->GetXaxis()->SetNdivisions(505);
  pt_boost_rebin_10[0]->GetYaxis()->SetNdivisions(505);
  pt_boost_rebin_10[0]->GetXaxis()->SetTitleSize(0.05);
  pt_boost_rebin_10[0]->GetYaxis()->SetTitleSize(0.04);
  pt_boost_rebin_10[0]->GetXaxis()->SetTitleOffset(0.9);
  pt_boost_rebin_10[0]->GetYaxis()->SetTitleOffset(1.5);
  pt_boost_rebin_10[0]->GetXaxis()->SetTitle("p_{T}");
  pt_boost_rebin_10[0]->GetYaxis()->SetTitle("");
  for(unsigned int i=0; i<files.size();i++){
    pt_boost_rebin_10[i]->SetLineWidth(2);
    pt_boost_rebin_10[i]->SetLineColor(colors[i]);
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  TCanvas *h = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  pt_boost_rebin_10[0]->Draw("HIST");
  pt_boost_rebin_10[1]->Draw("SAME HIST");
  pt_boost_rebin_10[4]->Draw("SAME HIST");
  TLegend *leg5 = new TLegend(0.65,0.65,0.85,0.85);
  leg5->AddEntry(pt_boost_rebin_10[0],"nominal","l");
  leg5->AddEntry(pt_boost_rebin_10[1],"1695","l");
  leg5->AddEntry(pt_boost_rebin_10[4],"1755","l");
  leg5->SetTextSize(0.05);
  leg5->Draw();
  gPad->RedrawAxis();
  h->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_main_rebin10_2016.pdf");

  pt_boost_rebin_10[2]->Draw("SAME HIST");
  pt_boost_rebin_10[3]->Draw("SAME HIST");
  leg5->AddEntry(pt_boost_rebin_10[2],"1715","l");
  leg5->AddEntry(pt_boost_rebin_10[3],"1735","l");
  h->SaveAs("/afs/desy.de/user/p/paaschal/Plots/LeptonPlots/BoostPt_Comparison_all_rebin10_2016.pdf");


}
