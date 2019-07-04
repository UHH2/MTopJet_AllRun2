#include "../include/CentralInclude.h"


using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);


int main(int argc, char* argv[]){
  TFile* file = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Histograms_combine.root");

  vector<TString> names = {"mc_gen", "SHOWER_fsrup_gen", "SHOWER_fsrdown_gen", "SCALE_upup_gen", "SHOWER_isrup_gen", "SHOWER_isrdown_gen"};
  vector<TH1F*> hists;
  for(auto n: names){
    hists.push_back((TH1F*)file->Get(n));
  }
  for(auto h: hists){
    h->SetLineWidth(4);
    h->SetLineColor(kRed+1);

  }
  hists[0]->SetLineColor(kAzure+7);

  vector<TH1F*> ratios;
  for(unsigned int i=0; i<hists.size(); i++){
    ratios.push_back(GetRatio(hists[i], hists[0]));
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=1; i<hists.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    hists[i]->GetYaxis()->SetRangeUser(0, 1.3 * hists[0]->GetMaximum());
    hists[i]->GetXaxis()->SetNdivisions(505);
    hists[i]->GetYaxis()->SetNdivisions(505);
    hists[i]->GetYaxis()->SetTitleOffset(1.6);
    hists[i]->GetXaxis()->SetTitleOffset(1.3);
    hists[i]->SetTitle(" ");
    hists[i]->GetYaxis()->SetTitle("events");
    hists[i]->Draw("HIST");
    hists[0]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.53,0.63,0.88,0.88);
    leg->AddEntry(hists[0], "nominal", "l");
    leg->AddEntry(hists[i], names[i], "l");
    leg->Draw();
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/GenComparison/"+names[i]+".pdf");
    delete a;
  }

  for(unsigned int i=1; i<ratios.size(); i++){
    TCanvas *b = new TCanvas("b", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    ratios[i]->SetTitle(" ");
    ratios[i]->GetYaxis()->SetTitle("#frac{variation}{nominal}");
    ratios[i]->GetXaxis()->SetNdivisions(505);
    ratios[i]->GetYaxis()->SetNdivisions(505);
    ratios[i]->GetYaxis()->SetTitleOffset(1.5);
    ratios[i]->GetXaxis()->SetTitleOffset(1.3);
    ratios[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratios[i]->SetLineColor(1);
    ratios[i]->Draw("E1");
    b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/GenComparison/Ratio_"+names[i]+".pdf");
    delete b;
  }

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
