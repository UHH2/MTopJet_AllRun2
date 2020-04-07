#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc != 3){
    cout << "Usage: ./MassComparison_SYS rec/gen cut:true/false" << endl;
    return 0;
  }

  bool isCut = false;
  if(strcmp(argv[2], "true") == 0) isCut = true;
  string cut;
  if(isCut){
    cout << "At which value did you cut the gen weight? ";
    cin >> cut;
  }

  TString level;
  bool isRec = true;
  if(argc > 1 && strcmp(argv[1], "rec") == 0) level = "rec";
  else if(argc > 1 && strcmp(argv[1], "gen") == 0){
    level = "gen";
    isRec = false;
  }
  else{
    cout << "Usage: ./MassComparison_SYS rec/gen" << endl;
    return 0;
  }

  vector<TString> year = {"2016", "2017", "2018"};

  vector<TFile*> files_nominal, files_JECup, files_JECdown;
  vector<TFile*> files_FSRup_sqrt2, files_FSRdown_sqrt2, files_FSRup_2, files_FSRdown_2, files_FSRup_4, files_FSRdown_4;
  vector<TFile*> files_ISRup_sqrt2, files_ISRdown_sqrt2, files_ISRup_2, files_ISRdown_2, files_ISRup_4, files_ISRdown_4;

  // Empty Histogram to keep the code clean. Since 2016 does not have FSR- & ISR-4,
  // the code gets dirty due to different vector sizes.
  TH1F *empty_hist;
  TFile *empty_file;

  for(unsigned int i=0; i<year.size(); i++){
    files_nominal.push_back(new TFile(dir+year[i]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root"));

    files_JECup.push_back(new TFile(dir+year[i]+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
    files_JECdown.push_back(new TFile(dir+year[i]+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));

    files_FSRup_2.push_back(new TFile(dir+year[i]+"/muon/FSRup_2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
    files_FSRdown_2.push_back(new TFile(dir+year[i]+"/muon/FSRdown_2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));

    files_ISRup_2.push_back(new TFile(dir+year[i]+"/muon/ISRup_2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
    files_ISRdown_2.push_back(new TFile(dir+year[i]+"/muon/ISRdown_2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
    if(i==0){
      files_FSRup_sqrt2.push_back(empty_file);
      files_FSRdown_sqrt2.push_back(empty_file);

      files_ISRup_sqrt2.push_back(empty_file);
      files_ISRdown_sqrt2.push_back(empty_file);

      files_FSRup_4.push_back(empty_file);
      files_FSRdown_4.push_back(empty_file);

      files_ISRup_4.push_back(empty_file);
      files_ISRdown_4.push_back(empty_file);
    }
    if(i>0){
      files_FSRup_sqrt2.push_back(new TFile(dir+year[i]+"/muon/FSRup_sqrt2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
      files_FSRdown_sqrt2.push_back(new TFile(dir+year[i]+"/muon/FSRdown_sqrt2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));

      files_ISRup_sqrt2.push_back(new TFile(dir+year[i]+"/muon/ISRup_sqrt2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
      files_ISRdown_sqrt2.push_back(new TFile(dir+year[i]+"/muon/ISRdown_sqrt2/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
      files_FSRup_4.push_back(new TFile(dir+year[i]+"/muon/FSRup_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
      files_FSRdown_4.push_back(new TFile(dir+year[i]+"/muon/FSRdown_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"));

      files_ISRup_4.push_back(new TFile(dir+year[i]+"/muon/ISRup_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
      files_ISRdown_4.push_back(new TFile(dir+year[i]+"/muon/ISRdown_4/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
    }
  }

  TString hist_mass;
  if(isRec) hist_mass = "XCone_cor/M_jet1";
  else      hist_mass = "XCone_GEN_GenOnly/Mass_HadJet33";
  vector<TH1F*> mass_nominal, mass_JECup, mass_JECdown;
  vector<TH1F*> mass_FSRup_sqrt2, mass_FSRdown_sqrt2, mass_FSRup_2, mass_FSRdown_2, mass_FSRup_4, mass_FSRdown_4;
  vector<TH1F*> mass_ISRup_sqrt2, mass_ISRdown_sqrt2, mass_ISRup_2, mass_ISRdown_2, mass_ISRup_4, mass_ISRdown_4;

  for(unsigned int i=0; i<year.size(); i++){
    mass_nominal.push_back((TH1F*)files_nominal[i]->Get(hist_mass));

    // --- JEC -----------------------------------------------------------------
    mass_JECup.push_back((TH1F*)files_JECup[i]->Get(hist_mass));
    mass_JECdown.push_back((TH1F*)files_JECdown[i]->Get(hist_mass));
    // --- FSR2 -----------------------------------------------------------------
    mass_FSRup_2.push_back((TH1F*)files_FSRup_2[i]->Get(hist_mass));
    mass_FSRdown_2.push_back((TH1F*)files_FSRdown_2[i]->Get(hist_mass));
    // --- ISR2 -----------------------------------------------------------------
    mass_ISRup_2.push_back((TH1F*)files_ISRup_2[i]->Get(hist_mass));
    mass_ISRdown_2.push_back((TH1F*)files_ISRdown_2[i]->Get(hist_mass));
    if(i==0){
      // --- FSRsqrt2 -----------------------------------------------------------------
      mass_FSRup_sqrt2.push_back(empty_hist);
      mass_FSRdown_sqrt2.push_back(empty_hist);
      // --- ISRsqrt2 -----------------------------------------------------------------
      mass_ISRup_sqrt2.push_back(empty_hist);
      mass_ISRdown_sqrt2.push_back(empty_hist);
      // --- FSR4 -----------------------------------------------------------------
      mass_FSRup_4.push_back(empty_hist);
      mass_FSRdown_4.push_back(empty_hist);
      // --- ISR4 -----------------------------------------------------------------
      mass_ISRup_4.push_back(empty_hist);
      mass_ISRdown_4.push_back(empty_hist);
    }
    if(i>0){
      // --- FSRsqrt2 -----------------------------------------------------------------
      mass_FSRup_sqrt2.push_back((TH1F*)files_FSRup_sqrt2[i]->Get(hist_mass));
      mass_FSRdown_sqrt2.push_back((TH1F*)files_FSRdown_sqrt2[i]->Get(hist_mass));
      // --- ISRsqrt2 -----------------------------------------------------------------
      mass_ISRup_sqrt2.push_back((TH1F*)files_ISRup_sqrt2[i]->Get(hist_mass));
      mass_ISRdown_sqrt2.push_back((TH1F*)files_ISRdown_sqrt2[i]->Get(hist_mass));
      // --- FSR4 -----------------------------------------------------------------
      mass_FSRup_4.push_back((TH1F*)files_FSRup_4[i]->Get(hist_mass));
      mass_FSRdown_4.push_back((TH1F*)files_FSRdown_4[i]->Get(hist_mass));
      // --- ISR4 -----------------------------------------------------------------
      mass_ISRup_4.push_back((TH1F*)files_ISRup_4[i]->Get(hist_mass));
      mass_ISRdown_4.push_back((TH1F*)files_ISRdown_4[i]->Get(hist_mass));
    }
  }
  vector<vector<TH1F*>> UPs = {mass_FSRup_sqrt2, mass_FSRup_2, mass_FSRup_4, mass_ISRup_sqrt2, mass_ISRup_2, mass_ISRup_4};
  vector<vector<TH1F*>> DOWNs = {mass_FSRdown_sqrt2, mass_FSRdown_2, mass_FSRdown_4, mass_ISRdown_sqrt2, mass_ISRdown_2, mass_ISRdown_4};

  std::vector<TString> file_name = {"FSR_sqrt2_", "FSR_2_", "FSR_4_","ISR_sqrt2_", "ISR_2_", "ISR_4_"};
  std::vector<TString> mass_name_up = {"FSRup_sqrt2", "FSRup_2", "FSRup_4","ISRup_sqrt2", "ISRup_2", "ISRup_4"};
  std::vector<TString> mass_name_down = {"FSRdown_sqrt2", "FSRdown_2", "FSRdown_4","ISRdown_sqrt2",  "ISRdown_2", "ISRdown_4"};

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
  TLegend *leg;

  // i -> iterates through years
  // j -> iterates through variations
  for(unsigned int i=0; i<year.size(); i++){
    mass_nominal[i]->SetTitle("");
    mass_nominal[i]->GetXaxis()->SetRangeUser(50, 300);
    mass_nominal[i]->GetYaxis()->SetRangeUser(0, (mass_nominal[i]->GetMaximum())*1.2);
    mass_nominal[i]->GetXaxis()->SetNdivisions(505);
    mass_nominal[i]->GetYaxis()->SetNdivisions(505);
    mass_nominal[i]->GetXaxis()->SetTitleSize(0.05);
    mass_nominal[i]->GetYaxis()->SetTitleSize(0.04);
    mass_nominal[i]->GetXaxis()->SetTitleOffset(0.9);
    mass_nominal[i]->GetYaxis()->SetTitleOffset(1.5);
    mass_nominal[i]->GetXaxis()->SetTitle("M_{jet}");
    mass_nominal[i]->GetYaxis()->SetTitle("Events");
    mass_nominal[i]->SetLineWidth(2);
    mass_nominal[i]->SetLineColor(kRed);
    mass_JECup[i]->SetLineWidth(2);
    mass_JECup[i]->SetLineColor(kBlue);
    mass_JECdown[i]->SetLineWidth(2);
    mass_JECdown[i]->SetLineColor(kGreen);
    for(unsigned int j=0; j<DOWNs.size(); j++){ // DOWNs and UPs (should) have the same size
      // skip 2016: FSR- & ISR-4
      if(i==0 && (j==0 || j==2 || j==3 || j==5)){
        continue; // "continue" skips to next iteration.
      }
      UPs[j].at(i)->SetLineWidth(2);
      UPs[j].at(i)->SetLineColor(kBlue);
      // ---------------------------------------------------------------------------
      DOWNs[j].at(i)->SetLineWidth(2);
      DOWNs[j].at(i)->SetLineColor(kGreen);
    }


    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    TCanvas *JEC = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    mass_nominal[i]->Draw("HIST");
    mass_JECup[i]->Draw("SAME HIST");
    mass_JECdown[i]->Draw("SAME HIST");
    leg = new TLegend(0.20,0.55,0.40,0.75);
    leg->AddEntry(mass_nominal[i],"nominal","l");
    leg->AddEntry(mass_JECup[i],"JEC up","l");
    leg->AddEntry(mass_JECdown[i],"JEC down","l");
    leg->SetTextSize(0.05);
    leg->Draw();
    gPad->RedrawAxis();
    if(isCut) JEC->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/"+level+"/gencut/"+cut+"/MassComparison_JEC_"+level+"_" +year[i]+".pdf");
    else JEC->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/"+level+"/MassComparison_JEC_"+level+"_" +year[i]+".pdf");
    delete JEC;
    leg->Clear();
    for(unsigned int j=0; j<DOWNs.size(); j++){
      // skip 2016: FSR- & ISR-4
      if(i==0 && (j==0 || j==2 || j==3 || j==5)) continue; // "continue" skips to next iteration. "break" leaves the loop
      TCanvas *SYS = new TCanvas();
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      mass_nominal[i]->Draw("HIST");
      UPs[j].at(i)->Draw("SAME HIST");
      DOWNs[j].at(i)->Draw("SAME HIST");
      leg = new TLegend(0.20,0.55,0.40,0.75);
      leg->AddEntry(mass_nominal[i],"nominal","l");
      leg->AddEntry(UPs[j].at(i), mass_name_up[j],"l");
      leg->AddEntry(DOWNs[j].at(i), mass_name_down[j],"l");
      leg->SetTextSize(0.05);
      leg->Draw();
      gPad->RedrawAxis();
      if(isCut) SYS->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/"+level+"/gencut/"+cut+"/MassComparison_"+file_name[j]+level+"_"+year[i]+".pdf");
      else SYS->SaveAs("/afs/desy.de/user/p/paaschal/Plots/ComparisonMass/"+level+"/MassComparison_"+file_name[j]+level+"_"+year[i]+".pdf");
      delete SYS;
      leg->Clear();
    }
  }
}
