#include "../include/CentralInclude.h"


using namespace std;
bool passed;
vector<double> *pt = new vector<double>(2);
double weight;


void fill_pt(TTree* tree, TH1F* h_pt);
void fill_pt_reweight(TTree* tree, TH1F* h_pt, TString pm);
void PlotHist(TH1F* h_hist, TString xaxis, TString histname);
void PlotHistCompare(TH1F* h1, TH1F* h2, TH1F* h3, TString xaxis, TString histname);
void PlotRatio(TH1F* h1, TH1F* h2, TString xaxis, TString histname);
void PlotRatioCompare(TH1F* h1, TH1F* h2, TH1F* h3, TString xaxis, TString histname);
double ReweightFunction(double pt, double maxdeviation, double crossing, double slope, TString pm);


int main(int argc, char* argv[]){
  TFile* file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TTree* tree = (TTree *) file->Get("AnalysisTree");

  TH1F* pt = new TH1F("ptb", "p_{T}^{b}", 100, 0, 500);
  TH1F* pt_plus = new TH1F("ptb_plus", "p_{T}^{b}", 100, 0, 500);
  TH1F* pt_minus = new TH1F("ptb_minus", "p_{T}^{b}", 100, 0, 500);
  fill_pt(tree, pt);
  fill_pt_reweight(tree, pt_plus, "plus");
  fill_pt_reweight(tree, pt_minus, "minus");
  PlotHist(pt, "p_{T}^{b}", "pt_bquarks");
  PlotHist(pt_plus, "p_{T}^{b}", "pt_bquarks_plus");
  PlotHist(pt_minus, "p_{T}^{b}", "pt_bquarks_minus");
  PlotHistCompare(pt_plus, pt_minus, pt, "p_{T}^{b}", "pt_bquarks_all");
  PlotRatio(pt_plus, pt, "p_{T}^{b}", "ratio_plus");
  PlotRatio(pt_minus, pt, "p_{T}^{b}", "ratio_minus");
  PlotRatioCompare(pt_plus, pt_minus, pt, "p_{T}^{b}", "ratio_all");

  cout << "central: " << pt->Integral() << endl;
  cout << "up:      " << pt_plus->Integral() << endl;
  cout << "down:    " << pt_minus->Integral() << endl;
}


void fill_pt(TTree* tree, TH1F* h_pt){
  if(!tree) cout << "could not read tree\n";
  else      cout << "Filling Histograms...\n";
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("bquark_pt",&pt);
  tree->SetBranchAddress("passed_measurement_gen",&passed);
  tree->SetBranchAddress("gen_weight",&weight);
  tree->SetBranchStatus("*",1);
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(passed){
      h_pt->Fill((*pt)[0], weight);
      h_pt->Fill((*pt)[1], weight);
    }
  }
  return;
}

void fill_pt_reweight(TTree* tree, TH1F* h_pt, TString pm){
  if(!tree) cout << "could not read tree\n";
  else      cout << "Filling Histograms...\n";
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("bquark_pt",&pt);
  tree->SetBranchAddress("passed_measurement_gen",&passed);
  tree->SetBranchAddress("gen_weight",&weight);
  tree->SetBranchStatus("*",1);
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    double maxdeviation = 0.1;
    double crossing = 165;
    double slope = 0.01;
    double sf1 = ReweightFunction((*pt)[0], maxdeviation, crossing, slope, pm);
    double sf2 = ReweightFunction((*pt)[1], maxdeviation, crossing, slope, pm);
    weight *= (sf1 * sf2);
    if(passed){
      h_pt->Fill((*pt)[0], weight);
      h_pt->Fill((*pt)[1], weight);
    }
  }
  return;
}

void PlotHist(TH1F* hist, TString xaxis, TString histname){
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle(xaxis);
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
  hist->SetLineColor(kBlack);
  hist->SetFillColor(13);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/BFragmentation/"+histname+".pdf");
  delete A;
  return;
}

void PlotHistCompare(TH1F* h1_, TH1F* h2_, TH1F* h3_, TString xaxis, TString histname){
  TH1F* h1 = (TH1F*) h1_->Clone();
  TH1F* h2 = (TH1F*) h2_->Clone();
  TH1F* h3 = (TH1F*) h3_->Clone();
  h3->SetTitle(" ");
  h3->GetXaxis()->SetTitle(xaxis);
  h3->GetYaxis()->SetTitle("events");
  h3->GetYaxis()->SetTitleOffset(1.6);
  h3->GetXaxis()->SetTitleOffset(1.3);
  h3->GetXaxis()->SetNdivisions(505);
  h3->GetYaxis()->SetNdivisions(505);
  h3->GetYaxis()->SetRangeUser(0, h1->GetMaximum()*1.5);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
  h3->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  h3->SetLineColor(13);
  h1->SetFillColor(0);
  h2->SetFillColor(0);
  h3->SetFillColor(0);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h3->Draw("HIST");
  h1->Draw("HIST SAME");
  h2->Draw("HIST SAME");
  TLegend *l = new TLegend(0.6, 0.65, 0.9, 0.9);
  l->AddEntry(h3, "central", "l");
  l->AddEntry(h1, "up", "l");
  l->AddEntry(h2, "down", "l");
  l->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/BFragmentation/"+histname+".pdf");
  delete A;
  return;
}


void PlotRatio(TH1F* h1, TH1F* h2, TString xaxis, TString histname){
  TH1F* hist = (TH1F*) h1->Clone();
  hist->Divide(h2);
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle(xaxis);
  hist->GetYaxis()->SetTitle("reweight/central");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetRangeUser(0.5, 1.5);
  hist->SetLineColor(kRed);
  hist->SetFillColor(0);
  hist->SetLineWidth(3);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/BFragmentation/"+histname+".pdf");
  delete A;
  return;
}

void PlotRatioCompare(TH1F* h1, TH1F* h2, TH1F* h3, TString xaxis, TString histname){
  TH1F* r1 = (TH1F*) h1->Clone();
  TH1F* r2 = (TH1F*) h2->Clone();
  TH1F* r3 = (TH1F*) h3->Clone();
  r1->Divide(h3);
  r2->Divide(h3);
  r3->Divide(h3);
  r1->SetTitle(" ");
  r1->GetXaxis()->SetTitle(xaxis);
  r1->GetYaxis()->SetTitle("reweight/central");
  r1->GetYaxis()->SetTitleOffset(1.6);
  r1->GetXaxis()->SetTitleOffset(1.3);
  r1->GetXaxis()->SetNdivisions(505);
  r1->GetYaxis()->SetNdivisions(505);
  r1->GetYaxis()->SetRangeUser(0.5, 1.5);
  r1->SetLineColor(kBlue);
  r2->SetLineColor(kRed);
  r3->SetLineColor(13);
  r1->SetFillColor(0);
  r2->SetFillColor(0);
  r3->SetFillColor(0);
  r1->SetLineWidth(3);
  r2->SetLineWidth(3);
  r3->SetLineWidth(3);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  r1->Draw("HIST");
  r2->Draw("HIST SAME");
  r3->Draw("HIST SAME");
  TLegend *l = new TLegend(0.6, 0.65, 0.9, 0.9);
  l->AddEntry(r3, "central", "l");
  l->AddEntry(r1, "up", "l");
  l->AddEntry(r2, "down", "l");
  l->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/BFragmentation/"+histname+".pdf");
  delete A;
  return;
}

double ReweightFunction(double pt, double maxdeviation, double crossing, double slope, TString pm){
  double factor = 1;
  if(pm == "plus") factor = 1 + maxdeviation * tanh( (pt-crossing)*slope );
  else             factor = 1 - maxdeviation * tanh( (pt-crossing)*slope );
  return factor;
}
