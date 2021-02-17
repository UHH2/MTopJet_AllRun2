#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc != 2){
    cout << "Usage: ./MassComparison_mTop <year>" << endl;
    return 0;
  }

  TString save_path = get_save_path();
  TString year = argv[1];

  TString uhh = year+"/muon/uhh2.AnalysisModuleRunner.MC.";
  vector<TString> name = {"nominal", "1665", "1695", "1715", "1735", "1755", "1785"};
  vector<int> colors = {kRed, kGray, kBlue, kGreen, kYellow, kBlack, kOrange};

  TFile *f_tt = new TFile(dir+uhh+"TTbar.root");
  TFile *f_1665 = new TFile(dir+uhh+"TTbar_mtop1665.root");
  TFile *f_1695 = new TFile(dir+uhh+"TTbar_mtop1695.root");
  TFile *f_1715 = new TFile(dir+uhh+"TTbar_mtop1715.root");
  TFile *f_1735 = new TFile(dir+uhh+"TTbar_mtop1735.root");
  TFile *f_1755 = new TFile(dir+uhh+"TTbar_mtop1755.root");
  TFile *f_1785 = new TFile(dir+uhh+"TTbar_mtop1785.root");
  vector<TFile*> files = {f_tt, f_1665, f_1695, f_1715, f_1735, f_1755, f_1785};
  vector<TString> hists_name = {"nominal 172.5", "mTop 166.5", "mTop 169.5", "mTop 171.5", "mTop 173.5", "mTop 175.5", "mTop 178.5"};

  TString hist_jetmass = "XCone_cor/M_jet1";

  vector<TH1F*> hists, hists_norm;
  for(unsigned int i=0; i<files.size(); i++) hists.push_back((TH1F*)files[i]->Get(hist_jetmass));
  hists_norm = normalize(hists);

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

  hists[0]->SetTitle("");
  hists[0]->GetXaxis()->SetRangeUser(0, 400);
  hists[0]->GetYaxis()->SetRangeUser(0, get_highest_peak(hists)*1.2);
  hists[0]->GetXaxis()->SetNdivisions(505);
  hists[0]->GetYaxis()->SetNdivisions(505);
  hists[0]->GetXaxis()->SetTitleSize(0.05);
  hists[0]->GetYaxis()->SetTitleSize(0.04);
  hists[0]->GetXaxis()->SetTitleOffset(0.9);
  hists[0]->GetYaxis()->SetTitleOffset(1.5);
  hists[0]->GetXaxis()->SetTitle("M_{jet}^{had}");
  hists[0]->GetYaxis()->SetTitle("");

  for(unsigned int i=0; i<hists.size(); i++){
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colors[i]);
  }

  // 0- nominal; 1- 66; 2- 69; 3- 71; 4- 73; 5- 75; 6- 78;
  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hists[0]->Draw("HIST");
  hists[2]->Draw("HIST SAME");
  hists[5]->Draw("HIST SAME");
  leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->AddEntry(hists[0],"nominal 172.5","l");
  leg->AddEntry(hists[2],"mTop 169.5","l");
  leg->AddEntry(hists[5],"mTop 175.5","l");
  leg->SetTextSize(0.03);
  leg->Draw();
  leg->Clear();
  gPad->RedrawAxis();
  A->SaveAs(save_path+"/Plots/MassSensitivity/HadjetMass_main_"+year+".pdf");
  hists[1]->Draw("SAME HIST");
  hists[3]->Draw("SAME HIST");
  hists[4]->Draw("SAME HIST");
  hists[6]->Draw("SAME HIST");
  leg->AddEntry(hists[0],"nominal 172.5","l");
  leg->AddEntry(hists[1],"mTop 166.5","l");
  leg->AddEntry(hists[2],"mTop 169.5","l");
  leg->AddEntry(hists[3],"mTop 171.5","l");
  leg->AddEntry(hists[4],"mTop 173.5","l");
  leg->AddEntry(hists[5],"mTop 175.5","l");
  leg->AddEntry(hists[6],"mTop 178.5","l");
  A->SaveAs(save_path+"/Plots/MassSensitivity/HadjetMass_full_"+year+".pdf");
  delete A;
  leg->Clear();

  /*
  ███    ██  ██████  ██████  ███    ███
  ████   ██ ██    ██ ██   ██ ████  ████
  ██ ██  ██ ██    ██ ██████  ██ ████ ██
  ██  ██ ██ ██    ██ ██   ██ ██  ██  ██
  ██   ████  ██████  ██   ██ ██      ██
  */

  hists_norm[0]->SetTitle("");
  hists_norm[0]->GetXaxis()->SetRangeUser(0, 400);
  hists_norm[0]->GetYaxis()->SetRangeUser(0, get_highest_peak(hists_norm)*1.2);
  hists_norm[0]->GetXaxis()->SetNdivisions(505);
  hists_norm[0]->GetYaxis()->SetNdivisions(505);
  hists_norm[0]->GetXaxis()->SetTitleSize(0.05);
  hists_norm[0]->GetYaxis()->SetTitleSize(0.04);
  hists_norm[0]->GetXaxis()->SetTitleOffset(0.9);
  hists_norm[0]->GetYaxis()->SetTitleOffset(1.5);
  hists_norm[0]->GetXaxis()->SetTitle("M_{jet}^{had}");
  hists_norm[0]->GetYaxis()->SetTitle("");

  for(unsigned int i=0; i<hists_norm.size(); i++){
    hists_norm[i]->SetLineWidth(2);
    hists_norm[i]->SetLineColor(colors[i]);
  }

  // 0- nominal; 1- 66; 2- 69; 3- 71; 4- 73; 5- 75; 6- 78;
  A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hists_norm[0]->Draw("HIST");
  hists_norm[2]->Draw("HIST SAME");
  hists_norm[5]->Draw("HIST SAME");
  leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->AddEntry(hists_norm[0],"nominal 172.5","l");
  leg->AddEntry(hists_norm[2],"mTop 169.5","l");
  leg->AddEntry(hists_norm[5],"mTop 175.5","l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path+"/Plots/MassSensitivity/HadjetMass_main_"+year+"_norm.pdf");
  hists_norm[1]->Draw("SAME HIST");
  hists_norm[3]->Draw("SAME HIST");
  hists_norm[4]->Draw("SAME HIST");
  hists_norm[6]->Draw("SAME HIST");
  leg->AddEntry(hists_norm[1],"mTop 166.5","l");
  leg->AddEntry(hists_norm[3],"mTop 171.5","l");
  leg->AddEntry(hists_norm[4],"mTop 173.5","l");
  leg->AddEntry(hists_norm[6],"mTop 178.5","l");
  A->SaveAs(save_path+"/Plots/MassSensitivity/HadjetMass_full_"+year+"_norm.pdf");
  delete A;
  leg->Clear();

}
