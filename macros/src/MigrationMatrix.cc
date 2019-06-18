#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  // declare files
  TFile *file = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Histograms_combine.root");

  // get matrix
  TH2F* matrix = (TH2F*)file->Get("mc_matrix");

  // set up lines to seperate regions
  vector<TString> labels_rec, labels_gen;
  labels_rec.push_back("#splitline{measurement}{(3 p_{T} bins)}");
  labels_rec.push_back("#splitline{p_{T} < 400}{(3 p_{T} bins)}");
  labels_rec.push_back("p_{T}^{subjet} < 30");
  labels_rec.push_back("mass cut");
  labels_rec.push_back("p_{T}^{#mu} < 60");
  labels_rec.push_back("b-tag");
  labels_gen.push_back("#splitline{measurement}{(2 p_{T} bins)}");
  labels_gen.push_back("#splitline{p_{T} < 400}{(2 p_{T} bins)}");
  labels_gen.push_back("p_{T}^{subjet} < 30");
  labels_gen.push_back("mass cut");
  labels_gen.push_back("p_{T}^{#mu} < 60");



  vector<double> bins_rec = {0.5, 78.5, 120.5, 144.5, 168.5, 174.5, 200.5};
  vector<double> subbins_rec = {26.5, 52.5, 92.5, 106.5};
  vector<double> bins_gen = {0.5, 14.5, 28.5, 35.5, 48.5, 55.5};
  vector<double> subbins_gen = {7.5, 21.5};

  vector<TLine*> lines_rec, sublines_rec;
  vector<TLine*> lines_gen, sublines_gen;

  for(auto bin: bins_rec){
    double xmin = -2.5;
    double xmax = bins_gen[bins_gen.size()-1];
    double ymin = bin;
    double ymax = bin;
    lines_rec.push_back(new TLine(xmin, ymin, xmax, ymax));
  }
  for(auto bin: subbins_rec){
    double xmin = -1.25;
    double xmax = bins_gen[bins_gen.size()-1];
    double ymin = bin;
    double ymax = bin;
    sublines_rec.push_back(new TLine(xmin, ymin, xmax, ymax));
  }

  for(auto bin: bins_gen){
    double xmin = bin;
    double xmax = bin;
    double ymin = -12;
    double ymax = bins_rec[bins_rec.size()-1];
    lines_gen.push_back(new TLine(xmin, ymin, xmax, ymax));
  }
  for(auto bin: subbins_gen){
    double xmin = bin;
    double xmax = bin;
    double ymin = -6;
    double ymax = bins_rec[bins_rec.size()-1];
    sublines_gen.push_back(new TLine(xmin, ymin, xmax, ymax));
  }


  // draw
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *c= new TCanvas("c","",2000,2000);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);
  // gPad->SetLogz();
  // matrix->Reset();
  matrix->SetTitle(" ");
  matrix->GetXaxis()->SetTitleOffset(1.6);
  matrix->GetYaxis()->SetTitleOffset(2.3);
  matrix->GetXaxis()->SetTitleSize(0.03);
  matrix->GetYaxis()->SetTitleSize(0.03);
  matrix->GetXaxis()->SetTitle("generator binning");
  matrix->GetYaxis()->SetTitle("detector binning");
  matrix->Draw("COLZ");
  matrix->Draw("BOX SAME");

  for(auto line: sublines_rec){
    line->SetLineColor(kRed-2);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw("SAME");
  }
  for(auto line: sublines_gen){
    line->SetLineColor(kRed-2);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw("SAME");
  }
  lines_rec.erase(lines_rec.begin());                          // remove first line
  lines_rec.erase(lines_rec.begin() + lines_rec.size() - 1);   // remove last line
  for(auto line: lines_rec){
    line->SetLineColor(kRed);
    line->SetLineWidth(4);
    line->Draw("SAME");
  }
  lines_gen.erase(lines_gen.begin());                          // remove first line
  lines_gen.erase(lines_gen.begin() + lines_gen.size() - 1);   // remove last line
  for(auto line: lines_gen){
    line->SetLineColor(kRed);
    line->SetLineWidth(4);
    line->Draw("SAME");
  }

  // first delete bin labels (numbers)
  matrix->GetXaxis()->SetLabelSize(0.0);
  matrix->GetYaxis()->SetLabelSize(0.0);

  // custom labels
  for(unsigned int i=1; i<bins_rec.size(); i++){
    TLatex text;
    text.SetTextFont(43);
    text.SetTextSize(40);
    text.SetTextAngle(90);
    text.SetTextAlign(21);
    double position = bins_rec[i-1] + (bins_rec[i]-bins_rec[i-1]) / 2;
    TString label = labels_rec[i-1];
    text.DrawLatex(-2.5, position, label);
  }
  for(unsigned int i=1; i<bins_gen.size(); i++){
    TLatex text;
    text.SetTextFont(43);
    text.SetTextSize(40);
    text.SetTextAngle(0);
    text.SetTextAlign(21);
    double position = bins_gen[i-1] + (bins_gen[i]-bins_gen[i-1]) / 2;
    TString label = labels_gen[i-1];
    text.DrawLatex(position, -12, label);
  }
  gPad->RedrawAxis();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/MigrationMatrix_Scheme.pdf");
  delete c;

  return 0;
}
