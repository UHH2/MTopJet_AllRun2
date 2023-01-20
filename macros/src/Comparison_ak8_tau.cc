#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){
  vector<TString> year = {"2017", "2018"};

  TString ttbar_string = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString data_string = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TString save_path = get_save_path()+"/Plots/Comparison_ak8/";

  vector<TString> masscut = {"no_masscut", "masscut_120", "masscut_130", "masscut_140", "masscut_150"};
  vector<TString> classes = {"AK8_tau_pass_rec", "AK8_tau_pass_rec_masscut_120", "AK8_tau_pass_rec_masscut_130", "AK8_tau_pass_rec_masscut_140", "AK8_tau_pass_rec_masscut_150"};
  vector<TH1F*> hists_all, hists_fm, hists_sm, hists_nm;

  TFile *file_ttbar;

  TString tau32_had_all = "/tau32_hadjet";
  TString tau32_had_fm = "/tau32_hadjet_fullymerged";
  TString tau32_had_sm = "/tau32_hadjet_semimerged";
  TString tau32_had_nm = "/tau32_hadjet_notmerged";

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
    for(unsigned int j=0; j<classes.size(); j++){
      hists_all.push_back((TH1F*)file_ttbar->Get(classes[j]+tau32_had_all)); // push_back Histogramm for N-subjettiness for each masscut
      hists_fm.push_back((TH1F*)file_ttbar->Get(classes[j]+tau32_had_fm));
      hists_sm.push_back((TH1F*)file_ttbar->Get(classes[j]+tau32_had_sm));
      hists_nm.push_back((TH1F*)file_ttbar->Get(classes[j]+tau32_had_nm));
    }

    for(unsigned int j=0; j<hists_all.size(); j++){
      if(j==0) hists_all[j]->SetTitle(masscut[j]);
      if(j>0)  hists_all[j]->SetTitle(masscut[j]+" GeV"); // The vector masscut should have the same length like hists_all/classes
      hists_all[j]->GetXaxis()->SetRangeUser(0, 400);
      hists_all[j]->GetYaxis()->SetRangeUser(0, hists_all[j]->GetMaximum()*1.2);
      hists_all[j]->GetXaxis()->SetNdivisions(505);
      hists_all[j]->GetYaxis()->SetNdivisions(505);
      hists_all[j]->GetXaxis()->SetTitleSize(0.05);
      hists_all[j]->GetYaxis()->SetTitleSize(0.04);
      hists_all[j]->GetXaxis()->SetTitleOffset(0.9);
      hists_all[j]->GetYaxis()->SetTitleOffset(1.5);
      hists_all[j]->GetXaxis()->SetTitle("#tau_{32}");
      hists_all[j]->GetYaxis()->SetTitle("");
      hists_all[j]->SetLineWidth(2);
      hists_all[j]->SetLineColor(kRed);

      TCanvas *A = new TCanvas("A", "A", 600, 600);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      hists_all[j]->Draw("HIST");
      A->SaveAs(save_path+year[i]+"/tau32_"+masscut[j]+"_"+year[i]+".pdf");
      delete A;



      /////////////////////////////////////////////////////////////////////////

      if(j==0){ // Comparison Hists
        hists_all[j]->SetTitle("");
        for(unsigned int k=1;k<hists_all.size();k++){
          hists_all[k]->SetLineWidth(2);
          hists_all[k]->SetLineColor(kBlue);

          hists_fm[k]->SetLineStyle(1);
          hists_fm[k]->SetLineWidth(1);
          hists_fm[k]->SetLineColor(kTeal);

          hists_sm[k]->SetLineStyle(10);
          hists_sm[k]->SetLineWidth(1);
          hists_sm[k]->SetLineColor(kTeal);

          hists_nm[k]->SetLineStyle(2);
          hists_nm[k]->SetLineWidth(1);
          hists_nm[k]->SetLineColor(kTeal);

          TCanvas *A = new TCanvas("A", "A", 600, 600);
          hists_all[0]->Draw("HIST");
          hists_all[k]->Draw("SAME HIST");
          leg = new TLegend(0.2,0.65,0.4,0.85);
          leg->AddEntry(hists_all[0],"no masscut","l");
          leg->AddEntry(hists_all[k],masscut[k],"l");
          leg->SetTextSize(0.05);
          leg->Draw();
          gPad->RedrawAxis();
          A->SaveAs(save_path+year[i]+"/Comparison_tau32_"+masscut[k]+"_"+year[i]+".pdf");
          hists_fm[k]->Draw("SAME HIST");
          hists_sm[k]->Draw("SAME HIST");
          hists_nm[k]->Draw("SAME HIST");
          leg->AddEntry(hists_fm[k],masscut[k]+" fm","l");
          leg->AddEntry(hists_sm[k],masscut[k]+" sm","l");
          leg->AddEntry(hists_nm[k],masscut[k]+" nm","l");
          leg->Draw();
          A->SaveAs(save_path+year[i]+"/Comparison_tau32_"+masscut[k]+"_merged_"+year[i]+".pdf");
          leg->Clear();
          delete A;
        }
      }
    }
    hists_all = {}; // Empty Hist at the end of the loop to avoid Segmentation Violation
  }
}



// later for ratio
// "tau32_hadjet_fullymerged", "tau32_hadjetsemimerged", "tau32_hadjet_notmerged"
