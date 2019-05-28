#include "../include/CentralInclude.h"


using namespace std;

void PlotEfficiency(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TString histname);


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *DATA_file = new TFile(dir+"uhh2.AnalysisModuleRunner.DATA.DATA.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TString> histnames = {"number_smalldR_all", "number_smalldR_pass", "number_largedR_all", "number_largedR_pass"};
  vector<TH1F*> rec, rec_data;

  for(auto name: histnames){
    rec.push_back( (TH1F*)TT_file->Get( "XCone_cor/"+name) );
    rec_data.push_back( (TH1F*)DATA_file->Get( "XCone_cor/"+name) );
  }

  double effi_tt_iso = rec[1]->Integral()/rec[0]->Integral();
  double effi_tt_noiso = rec[3]->Integral()/rec[2]->Integral();
  double effi_data_iso = rec_data[1]->Integral()/rec_data[0]->Integral();
  double effi_data_noiso = rec_data[3]->Integral()/rec_data[2]->Integral();


  TGraphAsymmErrors* h_effi_tt_iso = new TGraphAsymmErrors(rec[1], rec[0],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_data_iso = new TGraphAsymmErrors(rec_data[1], rec_data[0],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_tt_noiso = new TGraphAsymmErrors(rec[3], rec[2],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_data_noiso = new TGraphAsymmErrors(rec_data[3], rec_data[2],"cl=0.683 b(1,1) mode");

  PlotEfficiency(h_effi_data_iso, h_effi_data_noiso, "_data");
  PlotEfficiency(h_effi_tt_iso, h_effi_tt_noiso, "_tt");


  return 0;
}

void PlotEfficiency(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TString histname){
  // if(histname == "effi_pt") cout << h_data->Eval(25) << endl;
  // if(histname == "effi_pt") cout << h_mc->Eval(25) << endl;

  h_data->SetTitle(" ");
  h_data->GetXaxis()->SetTitle(" ");
  h_data->GetYaxis()->SetTitle("efficiency");
  h_data->GetYaxis()->SetTitleOffset(1.6);
  h_data->GetXaxis()->SetTitleOffset(1.3);
  h_data->GetXaxis()->SetNdivisions(505);
  h_data->GetYaxis()->SetNdivisions(505);
  h_data->GetYaxis()->SetRangeUser(0.0, .5);

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);

  h_mc->SetLineColor(kRed);
  h_mc->SetMarkerColor(kRed);
  h_mc->SetMarkerStyle(8);
  h_mc->SetMarkerSize(1);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h_data->Draw("AP");
  h_mc->Draw("P SAME");
  TLegend *leg = new TLegend(0.33,0.20,0.66,0.33);
  leg->AddEntry(h_data,"isolated","pl");
  leg->AddEntry(h_mc,"not isolated","pl");
  leg->Draw("");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/BTagEffi"+histname+".pdf");
  delete A;
  return;
}
