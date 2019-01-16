#include "../include/CentralInclude.h"

using namespace std;

// compile with:
// g++ -o resolution_subjets_flavorJEC resolution_subjets_flavorJEC.cc `root-config --cflags --glibs`

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown);
vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error);

int main(int argc, char* argv[]){

  bool use_median = false;
  if(argc >1 && strcmp(argv[1], "median") == 0) use_median = true;

  // declare files
  TFile *f_std = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_std_jec.root");
  TFile *f_match = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_match.root");
  TFile *f_leadingpt = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_leadingpt.root");
  TFile *f_lightonly = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_lightonly.root");
  TFile *f_btag = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_btag.root");
  TFile *f_fractions = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_flavor_fractions.root");

  // set binning
  int n_ptbin = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TString bin_strings[] = {"0", "50", "80", "120", "170", "220", "270", "320", "370", "420", "600"};

  // this is for naming of pdf files
  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  // now get plots from file
  vector<TH1F*> reso_match, reso_match_noJEC;
  vector<TH1F*> reso_leadingpt, reso_leadingpt_noJEC;
  vector<TH1F*> reso_lightonly, reso_lightonly_noJEC;
  vector<TH1F*> reso_btag, reso_btag_noJEC;
  vector<TH1F*> reso_fractions, reso_fractions_noJEC;
  vector<TH1F*> reso_std, reso_std_noJEC, reso_std_match;

  vector<TH1F*> ptrec_lightonly, ptrec_btag, ptrec_leadingpt, ptrec_match, ptrec_fractions;

  std::string dir_jec = "RecGenHists_subjets/";
  std::string dir_raw = "RecGenHists_subjets_noJEC/";
  std::string dir_match = "RecGenHists_subjets_matched/";

  std::string name_jec;
  std::string name_raw;
  std::string name_match;
  std::string name_ptrec;

  for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){

    name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
    name_raw = dir_raw + "PtResolution_" + std::to_string(ptbin);
    name_match = dir_match + "PtResolution_" + std::to_string(ptbin);
    name_ptrec = dir_jec + "PtRec_" + std::to_string(ptbin);

    reso_match.push_back( (TH1F*)f_match->Get(name_jec.c_str()) );
    reso_match_noJEC.push_back( (TH1F*)f_match->Get(name_raw.c_str()) );
    ptrec_match.push_back( (TH1F*)f_match->Get(name_ptrec.c_str()) );
    reso_leadingpt.push_back( (TH1F*)f_leadingpt->Get(name_jec.c_str()) );
    reso_leadingpt_noJEC.push_back( (TH1F*)f_leadingpt->Get(name_raw.c_str()) );
    ptrec_leadingpt.push_back( (TH1F*)f_leadingpt->Get(name_ptrec.c_str()) );
    reso_lightonly.push_back( (TH1F*)f_lightonly->Get(name_jec.c_str()) );
    reso_lightonly_noJEC.push_back( (TH1F*)f_lightonly->Get(name_raw.c_str()) );
    ptrec_lightonly.push_back( (TH1F*)f_lightonly->Get(name_ptrec.c_str()) );
    reso_btag.push_back( (TH1F*)f_btag->Get(name_jec.c_str()) );
    reso_btag_noJEC.push_back( (TH1F*)f_btag->Get(name_raw.c_str()) );
    ptrec_btag.push_back( (TH1F*)f_btag->Get(name_ptrec.c_str()) );
    reso_fractions.push_back( (TH1F*)f_fractions->Get(name_jec.c_str()) );
    reso_fractions_noJEC.push_back( (TH1F*)f_fractions->Get(name_raw.c_str()) );
    ptrec_fractions.push_back( (TH1F*)f_fractions->Get(name_ptrec.c_str()) );
    reso_std.push_back( (TH1F*)f_std->Get(name_jec.c_str()) );
    reso_std_noJEC.push_back( (TH1F*)f_std->Get(name_raw.c_str()) );
    reso_std_match.push_back( (TH1F*)f_std->Get(name_match.c_str()) );
  }

  // fit histograms and extract mean/median value
  TH1F * resolution_match  = GetResoPlot(reso_match, use_median);
  TH1F * resolution_match_noJEC = GetResoPlot(reso_match_noJEC, use_median);
  TH1F * resolution_leadingpt  = GetResoPlot(reso_leadingpt, use_median);
  TH1F * resolution_leadingpt_noJEC = GetResoPlot(reso_leadingpt_noJEC, use_median);
  TH1F * resolution_lightonly  = GetResoPlot(reso_lightonly, use_median);
  TH1F * resolution_lightonly_noJEC = GetResoPlot(reso_lightonly_noJEC, use_median);
  TH1F * resolution_btag  = GetResoPlot(reso_btag, use_median);
  TH1F * resolution_btag_noJEC = GetResoPlot(reso_btag_noJEC, use_median);
  TH1F * resolution_fractions  = GetResoPlot(reso_fractions, use_median);
  TH1F * resolution_fractions_noJEC = GetResoPlot(reso_fractions_noJEC, use_median);
  TH1F * resolution_std  = GetResoPlot(reso_std, use_median);
  TH1F * resolution_std_noJEC = GetResoPlot(reso_std_noJEC, use_median);
  TH1F * resolution_std_match = GetResoPlot(reso_std_match, use_median);
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // calculate non-closure as function of ptrec
  // arguments: (hists, use median?, do ptrec?, mean or error)
  vector<double> ptrec = GetMeans(ptrec_fractions, use_median, true, "mean");
  vector<double> ptrec_err = GetMeans(ptrec_fractions, use_median, true, "error");
  vector<double> MeanForUncert = GetMeans(reso_fractions, use_median, false, "mean");
  vector<double> MeanForUncert_err = GetMeans(reso_fractions, use_median, false, "error");
  ////
  TGraphErrors* non_closure = new TGraphErrors(n_ptbin, &ptrec[0], &MeanForUncert[0], &ptrec_err[0], &MeanForUncert_err[0]);
  // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
  // do same to reach ptrec = 600
  non_closure->SetPoint(non_closure->GetN(), 0, MeanForUncert[0]);
  int n_last = MeanForUncert.size()-1;
  non_closure->SetPoint(non_closure->GetN(), 600, MeanForUncert[n_last]);

  TGraph* area = AreaFromEnvelope(non_closure, "area");
  TGraph* upvar = AreaFromEnvelope(non_closure, "up");
  TGraph* downvar = AreaFromEnvelope(non_closure, "down");

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // non-closure in root file
  TString rootname = "files/FlavorJEC_nonClosure.root";
  TFile * outfile = new TFile(rootname,"RECREATE");;
  upvar->Write("Up");
  downvar->Write("Down");
  outfile->Close();
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // start plotting

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  resolution_std->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_std->GetXaxis()->SetNdivisions(505);
  resolution_std->GetYaxis()->SetNdivisions(505);
  resolution_std->GetXaxis()->SetTitleSize(0.05);
  resolution_std->GetYaxis()->SetTitleSize(0.04);
  resolution_std->GetXaxis()->SetTitleOffset(0.9);
  resolution_std->GetYaxis()->SetTitleOffset(1.5);
  resolution_std->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_std->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_std->SetLineWidth(4);
  resolution_std->SetLineColor(kAzure+7);
  resolution_std_noJEC->SetLineWidth(4);
  resolution_std_noJEC->SetLineColor(kOrange+1);
  resolution_std_match->SetLineWidth(4);
  resolution_std_match->SetLineColor(kGreen);

  resolution_match->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_match->GetXaxis()->SetNdivisions(505);
  resolution_match->GetYaxis()->SetNdivisions(505);
  resolution_match->GetXaxis()->SetTitleSize(0.05);
  resolution_match->GetYaxis()->SetTitleSize(0.04);
  resolution_match->GetXaxis()->SetTitleOffset(0.9);
  resolution_match->GetYaxis()->SetTitleOffset(1.5);
  resolution_match->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_match->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_match->SetLineWidth(4);
  resolution_match->SetLineColor(kRed+1);
  resolution_match_noJEC->SetLineWidth(4);
  resolution_match_noJEC->SetLineColor(kOrange+1);

  resolution_leadingpt->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_leadingpt->GetXaxis()->SetNdivisions(505);
  resolution_leadingpt->GetYaxis()->SetNdivisions(505);
  resolution_leadingpt->GetXaxis()->SetTitleSize(0.05);
  resolution_leadingpt->GetYaxis()->SetTitleSize(0.04);
  resolution_leadingpt->GetXaxis()->SetTitleOffset(0.9);
  resolution_leadingpt->GetYaxis()->SetTitleOffset(1.5);
  resolution_leadingpt->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_leadingpt->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_leadingpt->SetLineWidth(4);
  resolution_leadingpt->SetLineColor(kRed+1);
  resolution_leadingpt_noJEC->SetLineWidth(4);
  resolution_leadingpt_noJEC->SetLineColor(kOrange+1);

  resolution_lightonly->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_lightonly->GetXaxis()->SetNdivisions(505);
  resolution_lightonly->GetYaxis()->SetNdivisions(505);
  resolution_lightonly->GetXaxis()->SetTitleSize(0.05);
  resolution_lightonly->GetYaxis()->SetTitleSize(0.04);
  resolution_lightonly->GetXaxis()->SetTitleOffset(0.9);
  resolution_lightonly->GetYaxis()->SetTitleOffset(1.5);
  resolution_lightonly->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_lightonly->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_lightonly->SetLineWidth(4);
  resolution_lightonly->SetLineColor(kRed+1);
  resolution_lightonly_noJEC->SetLineWidth(4);
  resolution_lightonly_noJEC->SetLineColor(kOrange+1);

  resolution_btag->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_btag->GetXaxis()->SetNdivisions(505);
  resolution_btag->GetYaxis()->SetNdivisions(505);
  resolution_btag->GetXaxis()->SetTitleSize(0.05);
  resolution_btag->GetYaxis()->SetTitleSize(0.04);
  resolution_btag->GetXaxis()->SetTitleOffset(0.9);
  resolution_btag->GetYaxis()->SetTitleOffset(1.5);
  resolution_btag->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_btag->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_btag->SetLineWidth(4);
  resolution_btag->SetLineColor(kRed+1);
  resolution_btag_noJEC->SetLineWidth(4);
  resolution_btag_noJEC->SetLineColor(kOrange+1);

  resolution_fractions->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_fractions->GetXaxis()->SetNdivisions(505);
  resolution_fractions->GetYaxis()->SetNdivisions(505);
  resolution_fractions->GetXaxis()->SetTitleSize(0.05);
  resolution_fractions->GetYaxis()->SetTitleSize(0.04);
  resolution_fractions->GetXaxis()->SetTitleOffset(0.9);
  resolution_fractions->GetYaxis()->SetTitleOffset(1.5);
  resolution_fractions->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_fractions->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_fractions->SetLineWidth(4);
  resolution_fractions->SetLineColor(kRed+1);
  resolution_fractions_noJEC->SetLineWidth(4);
  resolution_fractions_noJEC->SetLineColor(kOrange+1);

  non_closure->GetXaxis()->SetRangeUser(0, 600);
  non_closure->GetYaxis()->SetRangeUser(-10, 10);
  non_closure->GetXaxis()->SetNdivisions(505);
  non_closure->GetYaxis()->SetNdivisions(505);
  non_closure->GetXaxis()->SetTitleSize(0.05);
  non_closure->GetYaxis()->SetTitleSize(0.04);
  non_closure->GetXaxis()->SetTitleOffset(0.9);
  non_closure->GetYaxis()->SetTitleOffset(1.5);
  non_closure->GetXaxis()->SetTitle("p_{T}^{rec}");
  non_closure->GetYaxis()->SetTitle("non-closure uncertainty [%]");
  non_closure->SetTitle(" ");
  non_closure->SetMarkerStyle(20);
  non_closure->SetMarkerSize(0.8);
  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------

  TCanvas *c0 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_std->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_std_noJEC->Draw("SAME E1");
  // resolution_std_match->Draw("SAME E1");
  TLegend *leg0 = new TLegend(0.45,0.85,0.80,0.65);
  leg0->AddEntry(resolution_std,"Standard JEC applied","l");
  leg0->AddEntry(resolution_std_noJEC,"no JEC applied","l");
  leg0->SetTextSize(0.05);
  leg0->Draw();
  TLatex text0;
  text0.SetNDC(kTRUE);
  text0.SetTextFont(43);
  text0.SetTextSize(18);
  text0.DrawLatex(.2,.2, "standard JEC");
  gPad->RedrawAxis();
  c0->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_"+mean_median+"_before.pdf");

  TCanvas *c1 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_match->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_match->Draw("SAME E1");
  resolution_match_noJEC->Draw("SAME E1");
  TLegend *leg1 = new TLegend(0.45,0.85,0.80,0.65);
  leg1->AddEntry(resolution_match,"Flavor JEC applied","l");
  leg1->AddEntry(resolution_std,"Standard JEC applied","l");
  leg1->AddEntry(resolution_match_noJEC,"no JEC applied","l");
  leg1->SetTextSize(0.05);
  leg1->Draw();
  TLatex text1;
  text1.SetNDC(kTRUE);
  text1.SetTextFont(43);
  text1.SetTextSize(18);
  text1.DrawLatex(.2,.2, "flavor JEC, matched to gen particle");
  gPad->RedrawAxis();
  c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorMatch/pt_"+mean_median+"_before.pdf");

  TCanvas *c2 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_leadingpt->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_leadingpt->Draw("SAME E1");
  resolution_leadingpt_noJEC->Draw("SAME E1");
  TLegend *leg2 = new TLegend(0.45,0.85,0.80,0.65);
  leg2->AddEntry(resolution_leadingpt,"Flavor JEC applied","l");
  leg2->AddEntry(resolution_std,"Standard JEC applied","l");
  leg2->AddEntry(resolution_leadingpt_noJEC,"no JEC applied","l");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  TLatex text2;
  text2.SetNDC(kTRUE);
  text2.SetTextFont(43);
  text2.SetTextSize(18);
  text2.DrawLatex(.2,.2, "flavor JEC, leading jet assumed to be b-jet");
  gPad->RedrawAxis();
  c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorLeadingPT/pt_"+mean_median+"_before.pdf");

  TCanvas *c3 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_lightonly->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_lightonly->Draw("SAME E1");
  resolution_lightonly_noJEC->Draw("SAME E1");
  TLegend *leg3 = new TLegend(0.45,0.85,0.80,0.65);
  leg3->AddEntry(resolution_lightonly,"Flavor JEC applied","l");
  leg3->AddEntry(resolution_std,"Standard JEC applied","l");
  leg3->AddEntry(resolution_lightonly_noJEC,"no JEC applied","l");
  leg3->SetTextSize(0.05);
  leg3->Draw();
  TLatex text3;
  text3.SetNDC(kTRUE);
  text3.SetTextFont(43);
  text3.SetTextSize(18);
  text3.DrawLatex(.2,.2, "flavor JEC, only light");
  gPad->RedrawAxis();
  c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorLightOnly/pt_"+mean_median+"_before.pdf");

  TCanvas *c4 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_btag->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_btag->Draw("SAME E1");
  resolution_btag_noJEC->Draw("SAME E1");
  TLegend *leg4 = new TLegend(0.45,0.85,0.80,0.65);
  leg4->AddEntry(resolution_btag,"Flavor JEC applied","l");
  leg4->AddEntry(resolution_std,"Standard JEC applied","l");
  leg4->AddEntry(resolution_btag_noJEC,"no JEC applied","l");
  leg4->SetTextSize(0.05);
  leg4->Draw();
  TLatex text4;
  text4.SetNDC(kTRUE);
  text4.SetTextFont(43);
  text4.SetTextSize(18);
  text4.DrawLatex(.2,.2, "flavor JEC, matched to b-tag");
  gPad->RedrawAxis();
  c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorBTag/pt_"+mean_median+"_before.pdf");

  TCanvas *c5 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_fractions->Draw("E1");
  zero_line->Draw("SAME");
  resolution_std->Draw("SAME E1");
  resolution_fractions->Draw("SAME E1");
  resolution_fractions_noJEC->Draw("SAME E1");
  TLegend *leg5 = new TLegend(0.45,0.85,0.80,0.65);
  leg5->AddEntry(resolution_fractions,"Flavor JEC applied","l");
  leg5->AddEntry(resolution_std,"Standard JEC applied","l");
  leg5->AddEntry(resolution_fractions_noJEC,"no JEC applied","l");
  leg5->SetTextSize(0.05);
  leg5->Draw();
  TLatex text5;
  text5.SetNDC(kTRUE);
  text5.SetTextFont(43);
  text5.SetTextSize(18);
  text5.DrawLatex(.2,.2, "flavor JEC, flavor fractions");
  gPad->RedrawAxis();
  c5->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorFractions/pt_"+mean_median+"_before.pdf");

  for(unsigned int i=0; i<reso_fractions.size(); i++){
    TCanvas *f5 = new TCanvas(" ", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    double mean = resolution_fractions->GetBinContent(i+1);
    TLine* line = new TLine(mean, 0, mean, reso_fractions[i]->GetMaximum() );
    line->SetLineWidth(4);
    line->SetLineColor(kRed+1);
    f5->cd(i+1);
    reso_fractions[i]->SetTitle(bin_strings[i] + " GeV < p_{T} < " + bin_strings[i+1] + " GeV");
    reso_fractions[i]->GetXaxis()->SetTitle("#frac{ p_{T}^{rec} - p_{T}^{gen} }{ p_{T}^{gen} }");
    reso_fractions[i]->GetYaxis()->SetTitle("Events");
    reso_fractions[i]->GetXaxis()->SetNdivisions(505);
    reso_fractions[i]->GetYaxis()->SetNdivisions(505);
    reso_fractions[i]->GetXaxis()->SetTitleSize(0.05);
    reso_fractions[i]->GetYaxis()->SetTitleSize(0.04);
    reso_fractions[i]->GetXaxis()->SetTitleOffset(0.9);
    reso_fractions[i]->GetYaxis()->SetTitleOffset(1.5);
    reso_fractions[i]->SetLineColor(kBlack);
    reso_fractions[i]->SetLineWidth(3);
    reso_fractions[i]->SetFillColor(15);
    reso_fractions[i]->Draw("HIST");
    line->Draw("SAME");
    f5->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/FlavorFractions/Fits_"+mean_median+"_ptbin-" + std::to_string(i)+ ".pdf");
    delete f5;
  }

  TCanvas *c5a = new TCanvas(" ", " ", 1600, 900);
  TLatex t;
  t.SetNDC(kTRUE);
  t.SetTextFont(43);
  t.SetTextSize(18);
  c5a->Divide(2,2);
  c5a->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_match->Draw("E1");
  zero_line->Draw("SAME");
  resolution_match->Draw("E1 SAME");
  resolution_match_noJEC->Draw("E1 SAME");
  leg1->Draw();
  t.DrawLatex(.2,.2, "match to gen particle");
  gPad->RedrawAxis();
  c5a->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_leadingpt->Draw("E1");
  zero_line->Draw("SAME");
  resolution_leadingpt->Draw("E1 SAME");
  resolution_leadingpt_noJEC->Draw("E1 SAME");
  leg2->Draw();
  t.DrawLatex(.2,.2, "leading p_{T}");
  gPad->RedrawAxis();
  c5a->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_lightonly->Draw("E1");
  zero_line->Draw("SAME");
  resolution_lightonly->Draw("E1 SAME");
  resolution_lightonly_noJEC->Draw("E1 SAME");
  leg3->Draw();
  t.DrawLatex(.2,.2, "only light flavor");
  gPad->RedrawAxis();
  c5a->cd(4);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_btag->Draw("E1");
  zero_line->Draw("SAME");
  resolution_btag->Draw("E1 SAME");
  resolution_btag_noJEC->Draw("E1 SAME");
  leg4->Draw();
  t.DrawLatex(.2,.2, "match to AK4 b-tag");
  gPad->RedrawAxis();
  c5a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/Flavor_comparison_"+mean_median+".pdf");

  TCanvas *c5b = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_btag->SetLineColor(kAzure+7);
  resolution_lightonly->SetLineColor(kRed+1);
  resolution_leadingpt->SetLineColor(13);
  resolution_match->SetLineColor(kGreen);
  resolution_btag->Draw("E1");
  zero_line->Draw("SAME");
  resolution_btag->Draw("SAME E1");
  resolution_lightonly->Draw("SAME E1");
  resolution_leadingpt->Draw("SAME E1");
  resolution_match->Draw("SAME E1");
  resolution_btag_noJEC->Draw("SAME E1");
  TLegend *leg5b = new TLegend(0.33,0.85,0.85,0.60);
  leg5b->AddEntry(resolution_btag,"matched with AK4 b-tag","l");
  leg5b->AddEntry(resolution_lightonly,"only light flavor","l");
  leg5b->AddEntry(resolution_leadingpt,"p_{T} leading jet as b-jet","l");
  leg5b->AddEntry(resolution_match,"matched to generator particles","l");
  leg5b->AddEntry(resolution_btag_noJEC,"no JEC applied","l");
  leg5b->SetTextSize(0.05);
  leg5b->Draw();
  TLatex text5b;
  text5b.SetNDC(kTRUE);
  text5b.SetTextFont(43);
  text5b.SetTextSize(18);
  text5b.DrawLatex(.2,.2, "lepton + jets");
  gPad->RedrawAxis();
  c5b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/Flavor_comparison2_"+mean_median+".pdf");



  TCanvas *c6 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  TGraphErrors* non_closure_percent = (TGraphErrors*) non_closure->Clone();
  for (int i=0;i<non_closure_percent->GetN();i++) non_closure_percent->GetY()[i] *= 100;
  TGraphErrors* area_percent = (TGraphErrors*) area->Clone();
  for (int i=0;i<area_percent->GetN();i++) area_percent->GetY()[i] *= 100;
  non_closure_percent->RemovePoint(n_last+2); // remove points at x=600 for drawing
  non_closure_percent->RemovePoint(n_last+1); // remove points at x=0 for drawing
  non_closure_percent->Draw("AP");
  area_percent->SetFillColor(16);
  area_percent->Draw("f SAME");
  zero_line->SetLineColor(kRed);
  zero_line->Draw("SAME");
  non_closure_percent->Draw("P SAME");
  TLegend *leg6 = new TLegend(0.45,0.85,0.80,0.65);
  leg6->AddEntry(non_closure_percent, mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]","pl");
  leg6->AddEntry(area_percent, "non-closure", "f");
  leg6->Draw();
  gPad->RedrawAxis();
  c6->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/Flavor_nonClosure_"+mean_median+".pdf");



  return 0;
}
//------------------------------------------------------------------------------
/*
██████  ███████ ████████     ██████  ███████ ███████  ██████
██       ██         ██        ██   ██ ██      ██      ██    ██
██   ███ █████      ██        ██████  █████   ███████ ██    ██
██    ██ ██         ██        ██   ██ ██           ██ ██    ██
██████  ███████    ██        ██   ██ ███████ ███████  ██████
*/

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median){
  int n_bins = hists.size();
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TH1F * reso = new TH1F("reso"," ",n_bins, bins);
  vector<double> mean, error;

  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1 = new TF1("f1","gaus",-0.2,0.2);
    hists[i]->Fit("f1", "QR");
    double m = f1->GetParameter(1);
    double e = f1->GetParError(1);
    if(!use_median) mean.push_back(m);
    error.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  TAxis *xaxis = reso->GetXaxis();
  for(unsigned int i=0; i<n_bins; i++){
    reso->Fill(xaxis->GetBinCenter(i+1), mean[i]);
    reso->SetBinError(i+1, error[i]);
  }
  return reso;
}
//------------------------------------------------------------------------------
/*
██████  ███████ ████████     ███    ███ ███████  █████  ███    ██ ███████
██       ██         ██        ████  ████ ██      ██   ██ ████   ██ ██
██   ███ █████      ██        ██ ████ ██ █████   ███████ ██ ██  ██ ███████
██    ██ ██         ██        ██  ██  ██ ██      ██   ██ ██  ██ ██      ██
██████  ███████    ██        ██      ██ ███████ ██   ██ ██   ████ ███████
*/

vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error){
  int n_bins = hists.size();
  vector<double> mean, err;

  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1;
    if(do_ptrec) f1 = new TF1("f1","gaus");
    else         f1 = new TF1("f1","gaus", -0.2, 0.2);
    TString options = "QR";
    if(do_ptrec) options = "Q";
    hists[i]->Fit("f1", options);
    double m = f1->GetParameter(1);
    if(!use_median) mean.push_back(m);
    double e = f1->GetParError(1);
    err.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  if(error == "error") return err;
  else return mean;
}

//------------------------------------------------------------------------------
/*
███    ██  ██████  ███    ██        ██████ ██       ██████  ███████ ██    ██ ██████  ███████
████   ██ ██    ██ ████   ██       ██      ██      ██    ██ ██      ██    ██ ██   ██ ██
██ ██  ██ ██    ██ ██ ██  ██ █████ ██      ██      ██    ██ ███████ ██    ██ ██████  █████
██  ██ ██ ██    ██ ██  ██ ██       ██      ██      ██    ██      ██ ██    ██ ██   ██ ██
██   ████  ██████  ██   ████        ██████ ███████  ██████  ███████  ██████  ██   ██ ███████
*/
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown){
  int steps = 150;
  double stepsize = 700/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  // first generate new TGraph with all positive values points
  const int Npoints = uncert->GetN();
  double xPoint[Npoints];
  double yPoint[Npoints];
  for(unsigned int i=0; i<Npoints; i++){
    double y;
    uncert->GetPoint(i, xPoint[i], y);
    yPoint[i] = sqrt(y*y); // take abs for y-value
  }
  TGraph* uncert_abs = new TGraph( Npoints, xPoint, yPoint );

  //

  for(int i=0; i<steps; i++){
    x = stepsize * i;
    up[i] = sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    down[i] = - sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    xvalues[i] = x;
  }
  TGraph *uncertfit;
  if(updown == "area"){
    uncertfit = new TGraph(2*steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
      uncertfit->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
    }
  }
  else if(updown == "up"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
    }
  }
  else if(updown == "down"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],down[i]);
    }
  }
  return uncertfit;
}
