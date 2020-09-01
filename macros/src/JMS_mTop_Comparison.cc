#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <TMath.h>
#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  TString hist = "XCone_cor/M_jet1";

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file  = new TFile(dir+"2017/muon/"+ttbar_path);
  TH1F  *ttbar       = (TH1F*)ttbar_file->Get(hist);

  ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root";
  TFile *ttbar_file_1695  = new TFile(dir+"2017/muon/"+ttbar_path);
  TH1F  *ttbar_1695       = (TH1F*)ttbar_file->Get(hist);

  ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root";
  TFile *ttbar_file_1755  = new TFile(dir+"2017/muon/"+ttbar_path);
  TH1F  *ttbar_1755       = (TH1F*)ttbar_file->Get(hist);


  ttbar_1695->SetTitle("");
  ttbar_1695->GetXaxis()->SetRangeUser(0, 500);
  ttbar_1695->GetYaxis()->SetRangeUser(0, ttbar_1695->GetMaximum()*1.2);
  ttbar_1695->GetXaxis()->SetNdivisions(505);
  ttbar_1695->GetYaxis()->SetNdivisions(505);
  ttbar_1695->GetXaxis()->SetTitleSize(0.05);
  ttbar_1695->GetYaxis()->SetTitleSize(0.04);
  ttbar_1695->GetXaxis()->SetTitleOffset(0.9);
  ttbar_1695->GetYaxis()->SetTitleOffset(1.5);
  ttbar_1695->GetXaxis()->SetTitle("m_{Wjet}");
  ttbar_1695->GetYaxis()->SetTitle("");
  ttbar_1695->SetLineWidth(2);  // ttbar hist style
  ttbar_1695->SetLineColor(kRed);

  ttbar_1755->SetLineWidth(2);  // ttbar hist style
  ttbar_1755->SetLineColor(kBlack);

  ttbar->SetLineWidth(2);  // ttbar hist style
  ttbar->SetLineColor(kBlue);

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar_1695->Draw("Hist");
  ttbar->Draw("Hist same");
  ttbar_1755->Draw("Hist same");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/ellipse_compare.pdf");


  // vector<double> nominal_min = {0.493, 0.208};
  // vector<double> nominal_uu  = {0.667, 0.276};
  // vector<double> nominal_dd  = {0.314, 0.149};
  // vector<double> nominal_du  = {0.241, 0.88};
  // vector<double> nominal_ud  = {0.739, -0.463};
  //
  // vector<double> min_min = {0.480, 0.298};
  // vector<double> min_uu  = {0.608, 0.348};
  // vector<double> min_dd  = {0.345, 0.240};
  // vector<double> min_du  = {0.339, 0.623};
  // vector<double> min_ud  = {0.591, -0.021};
  //
  // vector<double> max_min = {0.568, -0.065};
  // vector<double> max_uu  = {0.716, 0.0};
  // vector<double> max_dd  = {0.425, -0.147};
  // vector<double> max_du  = {0.421, -0.298};
  // vector<double> max_ud  = {0.718, -0.429};
  //
  // vector<vector<double>> points_min = {nominal_min, min_min, max_min};
  // vector<vector<double>> points_uu  = {nominal_uu, min_uu, max_uu};
  // vector<vector<double>> points_dd  = {nominal_dd, min_dd, max_dd};
  // vector<vector<double>> points_du  = {nominal_du, min_du, max_du};
  // vector<vector<double>> points_ud  = {nominal_ud, min_ud, max_ud};
  // vector<vector<double>> distance   = {{0.1875, 0.7165}, {0.142, 0.345}, {0.1635, 0.393}};
  //
  //
  // vector<TGraph*> graphs;
  // vector<TEllipse*> ellipses;
  //
  // gStyle->SetPadTickY(1);
  // gStyle->SetPadTickX(1);
  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetLegendBorderSize(0);
  //
  // TCanvas *A = new TCanvas("A", "A", 600, 600);
  // gPad->SetLeftMargin(0.15);
  // gPad->SetBottomMargin(0.12);
  // for(int i =0; i<3; i++){
  //   TGraph *graph = new TGraph();
  //   graphs.push_back(graph);
  //   graphs[i]->SetPoint(0, points_min[i][0], points_min[i][1]);
  //   if(i==0) graphs[i]->SetMarkerColor(kBlack);
  //   if(i==1) graphs[i]->SetMarkerColor(kBlue);
  //   if(i==2) graphs[i]->SetMarkerColor(kRed);
  //   graphs[i]->SetMarkerStyle(kFullCircle);
  //   graphs[i]->SetMarkerSize(0.4);
  //   graphs[i]->GetYaxis()->SetRangeUser(-3, 3);
  //   graphs[i]->GetXaxis()->SetRangeUser(-3, 3);
  //
  //   double scalar   = abs(points_uu[i][0]*1+points_uu[i][1]*0); // look at vector (1,0)
  //   double length_v = sqrt(points_uu[i][0]*points_uu[i][0]+points_uu[i][1]*points_uu[i][1]);
  //   double theta= acos(scalar/length_v);
  //   cout << theta << endl;
  //   TEllipse *ellipse = new TEllipse(points_min[i][0], points_min[i][1], distance[i][0], distance[i][1], 0, 360, theta);
  //   ellipses.push_back(ellipse);
  //   if(i==0) graphs[i]->SetLineColor(kBlack);
  //   if(i==1) graphs[i]->SetLineColor(kBlue);
  //   if(i==2) graphs[i]->SetLineColor(kRed);
  //   if(i==0) graphs[i]->Draw("P");
  //   else     graphs[i]->Draw("SAME P");
  //   ellipses[i]->Draw("SAME");
  // }
  // gPad->RedrawAxis();
  // A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/ellipse_compare.pdf");

}
