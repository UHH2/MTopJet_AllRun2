#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  // TString year = "2016";
  //
  // TFile* files_FSRup_2 = new TFile(dir+"2016/muon/FSRup_2/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  // TFile* files_FSRdown_2 = new TFile(dir+"2016/muon/FSRdown_2/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  // TFile* files_FSRnominal = new TFile(dir+"2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  //
  // TH1F* mass_FSRup_2 = (TH1F*)files_FSRup_2->Get("XCone_cor/M_jet1");
  // TH1F* mass_FSRdown_2 = (TH1F*)files_FSRdown_2->Get("XCone_cor/M_jet1");
  // TH1F* mass_FSRnominal = (TH1F*)files_FSRnominal->Get("XCone_cor/M_jet1");
  //
  // TFile* files_FSRup_4 = new TFile(dir+"2017/muon/test/uhh2.AnalysisModuleRunner.MC.TTbar_Mtt0000to0700_SemiLep_2017v2.root ");
  // TH1F* weight_1 = (TH1F*)files_FSRup_4->Get("XCone_GEN_GenOnly/Weight_Massbin_145");
  // TH1F* weight_1_max_y100 = (TH1F*)files_FSRup_4->Get("XCone_GEN_GenOnly/Weight_Massbin_145");
  // TH1F* weight_weight = (TH1F*)files_FSRup_4->Get("XCone_GEN_GenOnly/Weight_Massbin_145_weight");
  //
  //
  //
  // /*
  // ██████  ██       ██████  ████████
  // ██   ██ ██      ██    ██    ██
  // ██████  ██      ██    ██    ██
  // ██      ██      ██    ██    ██
  // ██      ███████  ██████     ██
  // */
  // gStyle->SetPadTickY(1);
  // gStyle->SetPadTickX(1);
  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetLegendBorderSize(0);
  // TLegend *leg;
  //
  // mass_FSRup_2->GetXaxis()->SetRangeUser(50, 300);
  // mass_FSRup_2->GetYaxis()->SetRangeUser(0, (mass_FSRup_2->GetMaximum())*1.2);
  // mass_FSRup_2->GetXaxis()->SetNdivisions(505);
  // mass_FSRup_2->GetYaxis()->SetNdivisions(505);
  // mass_FSRup_2->GetXaxis()->SetTitleSize(0.05);
  // mass_FSRup_2->GetYaxis()->SetTitleSize(0.04);
  // mass_FSRup_2->GetXaxis()->SetTitleOffset(0.9);
  // mass_FSRup_2->GetYaxis()->SetTitleOffset(1.5);
  // mass_FSRup_2->GetXaxis()->SetTitle("M_{jet}");
  // mass_FSRup_2->GetYaxis()->SetTitle("Mass Distribution");
  // mass_FSRup_2->SetLineWidth(2);
  // mass_FSRup_2->SetLineColor(kRed);
  //
  // mass_FSRdown_2->SetLineWidth(2);
  // mass_FSRdown_2->SetLineColor(kBlue);
  //
  // mass_FSRnominal->SetLineWidth(2);
  // mass_FSRnominal->SetLineColor(kYellow);
  //
  // TCanvas *JEC = new TCanvas();
  // gPad->SetLeftMargin(0.15);
  // gPad->SetBottomMargin(0.12);
  // mass_FSRup_2->Draw("HIST");
  // mass_FSRdown_2->Draw("HIST same");
  // mass_FSRnominal->Draw("HIST same");
  // leg = new TLegend(0.20,0.55,0.40,0.75);
  // leg->AddEntry(mass_FSRup_2,"FSRup 2","l");
  // leg->AddEntry(mass_FSRdown_2,"FSRdown 2","l");
  // leg->AddEntry(mass_FSRnominal,"FSRnominal","l");
  // leg->SetTextSize(0.05);
  // leg->Draw();
  // gPad->RedrawAxis();
  // JEC->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/2016_FSRup2.pdf");
  // delete JEC;
  // leg->Clear();
  //
  // weight_1->GetYaxis()->SetRangeUser(0, (weight_1->GetMaximum())*1.2);
  // weight_1->GetXaxis()->SetNdivisions(505);
  // weight_1->GetYaxis()->SetNdivisions(505);
  // weight_1->GetXaxis()->SetTitleSize(0.05);
  // weight_1->GetYaxis()->SetTitleSize(0.04);
  // weight_1->GetXaxis()->SetTitleOffset(0.9);
  // weight_1->GetYaxis()->SetTitleOffset(1.5);
  // weight_1->GetXaxis()->SetTitle("weights");
  // weight_1->SetLineWidth(2);
  // weight_1->SetLineColor(kRed);
  //
  // TCanvas *NW = new TCanvas();
  // gPad->SetLeftMargin(0.15);
  // gPad->SetBottomMargin(0.12);
  // weight_1->Draw("HIST");
  // leg = new TLegend(0.20,0.55,0.40,0.75);
  // leg->AddEntry(weight_1,"no weight","l");
  // leg->SetTextSize(0.05);
  // leg->Draw();
  // gPad->RedrawAxis();
  // NW->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/2017_noweight_gen_fsrup4.pdf");
  // delete NW;
  // leg->Clear();
  //
  // weight_1_max_y100->GetYaxis()->SetRangeUser(0, 100);
  // weight_1_max_y100->GetXaxis()->SetNdivisions(505);
  // weight_1_max_y100->GetYaxis()->SetNdivisions(505);
  // weight_1_max_y100->GetXaxis()->SetTitleSize(0.05);
  // weight_1_max_y100->GetYaxis()->SetTitleSize(0.04);
  // weight_1_max_y100->GetXaxis()->SetTitleOffset(0.9);
  // weight_1_max_y100->GetYaxis()->SetTitleOffset(1.5);
  // weight_1_max_y100->GetXaxis()->SetTitle("weights");
  // weight_1_max_y100->SetLineWidth(2);
  // weight_1_max_y100->SetLineColor(kRed);
  //
  // TCanvas *NW100 = new TCanvas();
  // gPad->SetLeftMargin(0.15);
  // gPad->SetBottomMargin(0.12);
  // weight_1_max_y100->Draw("HIST");
  // leg = new TLegend(0.20,0.55,0.40,0.75);
  // leg->AddEntry(weight_1_max_y100,"no weight","l");
  // leg->SetTextSize(0.05);
  // leg->Draw();
  // gPad->RedrawAxis();
  // NW100->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/2017_noweight_gen_max_y100_fsrup4.pdf");
  // delete NW100;
  // leg->Clear();
  //
  // weight_weight->GetYaxis()->SetRangeUser(0, (weight_weight->GetMaximum())*1.2);
  // weight_weight->GetXaxis()->SetNdivisions(505);
  // weight_weight->GetYaxis()->SetNdivisions(505);
  // weight_weight->GetXaxis()->SetTitleSize(0.05);
  // weight_weight->GetYaxis()->SetTitleSize(0.04);
  // weight_weight->GetXaxis()->SetTitleOffset(0.9);
  // weight_weight->GetYaxis()->SetTitleOffset(1.5);
  // weight_weight->GetXaxis()->SetTitle("weights");
  // weight_weight->SetLineWidth(2);
  // weight_weight->SetLineColor(kRed);
  //
  // TCanvas *W = new TCanvas();
  // gPad->SetLeftMargin(0.15);
  // gPad->SetBottomMargin(0.12);
  // weight_weight->Draw("HIST");
  // leg = new TLegend(0.20,0.55,0.40,0.75);
  // leg->AddEntry(weight_weight,"weight","l");
  // leg->SetTextSize(0.05);
  // leg->Draw();
  // gPad->RedrawAxis();
  // W->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/2017_weight_gen_fsrup4.pdf");
  // delete W;
  // leg->Clear();


  TString fsrup_4_2017 = dir+"2017/muon/FSRup_4/uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString fsrup_4_2018 = dir+"2018/muon/FSRup_4/uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString fsrdown_4_2017 = dir+"2017/muon/FSRdown_4/uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString fsrdown_4_2018 = dir+"2018/muon/FSRdown_4/uhh2.AnalysisModuleRunner.MC.TTbar.root";
  cout << "debug" << endl;
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  vector<TString> names = {"fsrup_4_2017", "fsrup_4_2018", "fsrdown_4_2017", "fsrdown_4_2018"};
  vector<TString> files = {fsrup_4_2017, fsrup_4_2018, fsrdown_4_2017, fsrdown_4_2018};
  TFile* file;
  TH1F* hist1, *hist2;
  cout << "debug" << endl;
  for(unsigned int i=0; i<files.size(); i++){
    cout << "debug"<< i << endl;
    file = new TFile(files[i]);
    hist1 = (TH1F*)file->Get("Gen_Weight_300/Events_noweight");
    hist2 = (TH1F*)file->Get("Gen_Weight_gensel_300/Events_noweight");
    cout << "debug"<< i << endl;
    hist1->GetXaxis()->SetRangeUser(50, 300);
    hist1->GetYaxis()->SetRangeUser(0, (hist1->GetMaximum())*1.2);
    hist1->GetXaxis()->SetNdivisions(505);
    hist1->GetYaxis()->SetNdivisions(505);
    hist1->GetXaxis()->SetTitleSize(0.05);
    hist1->GetYaxis()->SetTitleSize(0.04);
    hist1->GetXaxis()->SetTitleOffset(0.9);
    hist1->GetYaxis()->SetTitleOffset(1.5);
    hist1->GetXaxis()->SetTitle("M_{jet}");
    hist1->GetYaxis()->SetTitle("Mass Distribution");
    hist1->SetLineWidth(2);
    hist1->SetLineColor(kRed);
    cout << "debug"<< i << endl;
    hist2->GetXaxis()->SetRangeUser(50, 300);
    hist2->GetYaxis()->SetRangeUser(0, (hist2->GetMaximum())*1.2);
    hist2->GetXaxis()->SetNdivisions(505);
    hist2->GetYaxis()->SetNdivisions(505);
    hist2->GetXaxis()->SetTitleSize(0.05);
    hist2->GetYaxis()->SetTitleSize(0.04);
    hist2->GetXaxis()->SetTitleOffset(0.9);
    hist2->GetYaxis()->SetTitleOffset(1.5);
    hist2->GetXaxis()->SetTitle("M_{jet}");
    hist2->GetYaxis()->SetTitle("Mass Distribution");
    hist2->SetLineWidth(2);
    hist2->SetLineColor(kRed);
    cout << "debug"<< i << endl;

    /*
    ██████  ██       ██████  ████████
    ██   ██ ██      ██    ██    ██
    ██████  ██      ██    ██    ██
    ██      ██      ██    ██    ██
    ██      ███████  ██████     ██
    */

    TCanvas *JEC = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    hist1->Draw("HIST");
    gPad->RedrawAxis();
    JEC->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/number_event_"+names[i]+".pdf");
    delete JEC;
    cout << "debug"<< i << endl;
    TCanvas *JEC1 = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    hist1->Draw("HIST");
    gPad->RedrawAxis();
    JEC1->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/number_event_"+names[i]+"_genweight.pdf");
    delete JEC1;

  }

}
