#include "../include/CentralInclude.h"


using namespace std;

TH1F* GetRatio(TH1F* h1, TH1F* h2);

int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile* file_c = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* file_u = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_fsrup.root");
  TFile* file_d = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown.root");

  TFile* file2_c = new TFile(dir_elec+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* file2_u = new TFile(dir_elec+"uhh2.AnalysisModuleRunner.MC.TTbar_fsrup.root");
  TFile* file2_d = new TFile(dir_elec+"uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TH1F*> hist_c, hist_u, hist_d;

  vector<TString> histnames;
  histnames.push_back("XCone_GEN_GenOnly_matched/Mass_HadJet33");
  histnames.push_back("XCone_GEN_GenOnly_unmatched/Mass_HadJet33");

  vector<TString> filenames = {"matched", "unmatched"};


  for(auto name: histnames){
    TH1F* c = (TH1F*)file_c->Get(name);
    TH1F* u = (TH1F*)file_u->Get(name);
    TH1F* d = (TH1F*)file_d->Get(name);
    TH1F* ec = (TH1F*)file2_c->Get(name);
    TH1F* eu = (TH1F*)file2_u->Get(name);
    TH1F* ed = (TH1F*)file2_d->Get(name);
    c->Add(ec);
    u->Add(eu);
    d->Add(ed);
    hist_c.push_back(c);
    hist_u.push_back(u);
    hist_d.push_back(d);
  }

  TH1F* r_matched_u = GetRatio(hist_u[0], hist_c[0]);
  TH1F* r_matched_d = GetRatio(hist_d[0], hist_c[0]);
  TH1F* r_unmatched_u = GetRatio(hist_u[1], hist_c[1]);
  TH1F* r_unmatched_d = GetRatio(hist_d[1], hist_c[1]);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<hist_c.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    hist_c[i]->SetLineColor(kAzure+7);
    hist_u[i]->SetLineColor(kRed+1);
    hist_d[i]->SetLineColor(kRed+1);
    hist_c[i]->SetLineWidth(3);
    hist_u[i]->SetLineWidth(3);
    hist_d[i]->SetLineWidth(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    hist_c[i]->GetYaxis()->SetRangeUser(0, 1.3 * hist_c[i]->GetMaximum());
    if(i==0)hist_c[i]->GetXaxis()->SetRangeUser(100, 250);
    TString title = hist_c[i]->GetTitle();
    hist_c[i]->GetXaxis()->SetTitle(title);
    hist_c[i]->SetTitle(" ");
    hist_c[i]->Draw("HIST");
    hist_d[i]->Draw("HIST SAME");
    hist_u[i]->Draw("HIST SAME");
    hist_c[i]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.63,0.6,0.88,0.88);
    leg->AddEntry(hist_c[i], "central", "l");
    leg->AddEntry(hist_u[i], "FSR Variations", "l");
    leg->Draw();
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/FSR_"+filenames[i]+".pdf");
    delete a;
  }

  TCanvas *b = new TCanvas("b", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  r_matched_u->SetTitle(" ");
  r_matched_u->GetYaxis()->SetTitle("#frac{variation}{nominal}");
  r_matched_u->GetXaxis()->SetNdivisions(505);
  r_matched_u->GetYaxis()->SetNdivisions(505);
  r_matched_u->GetYaxis()->SetTitleOffset(1.5);
  r_matched_u->GetXaxis()->SetTitleOffset(1.3);
  r_matched_u->GetYaxis()->SetRangeUser(0.2, 1.8);
  r_matched_u->GetXaxis()->SetRangeUser(100, 250);
  r_matched_u->SetLineColor(kAzure+7);
  r_matched_d->SetLineColor(kRed);
  r_matched_u->Draw("HIST");
  r_matched_d->Draw("HIST SAME");
  TLegend *legb = new TLegend(0.3,0.7,0.7,0.88);
  legb->AddEntry(r_matched_u, "up variation", "l");
  legb->AddEntry(r_matched_d, "down variation", "l");
  legb->Draw();
  b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/FSR_ratio_matched.pdf");
  delete b;

  TCanvas *c = new TCanvas("c", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  r_unmatched_u->SetTitle(" ");
  r_unmatched_u->GetYaxis()->SetTitle("#frac{variation}{nominal}");
  r_unmatched_u->GetXaxis()->SetNdivisions(505);
  r_unmatched_u->GetYaxis()->SetNdivisions(505);
  r_unmatched_u->GetYaxis()->SetTitleOffset(1.5);
  r_unmatched_u->GetXaxis()->SetTitleOffset(1.3);
  r_unmatched_u->GetYaxis()->SetRangeUser(0.2, 1.8);
  r_unmatched_u->SetLineColor(kAzure+7);
  r_unmatched_d->SetLineColor(kRed);
  r_unmatched_u->Draw("HIST");
  r_unmatched_d->Draw("HIST SAME");
  TLegend *legc = new TLegend(0.3,0.7,0.7,0.88);
  legc->AddEntry(r_unmatched_u, "up variation", "l");
  legc->AddEntry(r_unmatched_d, "down variation", "l");
  legc->Draw();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/FSR_ratio_unmatched.pdf");
  delete c;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}



TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double r;
    if(N2 == 0) r=1;
    else r = N1/N2;
    double error = sqrt(fabs((1-r)*r/N2));
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
      ratio->SetBinError(i, 0);
    }
    else{
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
