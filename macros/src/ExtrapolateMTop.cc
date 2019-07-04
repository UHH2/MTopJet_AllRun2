#include "../include/CentralInclude.h"

using namespace std;

TH1F* CreateNewSample(vector<TF1*> fits, double mass, TH1F* dummy);
void PlotNewSample(TH1F* hist, TString name, TString channel);
void PlotBinHists(vector<TGraph*> binhists, vector<TF1*> fits, TString channel);

int main(int argc, char* argv[]){

  vector<TString> channels = {"muon", "elec"};

  for(auto channel: channels){

    // declare files
    TFile *file = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Histograms_"+channel+".root");

    // get truth
    vector<TString> mass_names = {"mc_mtop1695", "mc_mtop1715", "mc_mtop1735","mc_mtop1755"};
    vector<TH1F*> h_truth;
    for(auto name: mass_names) h_truth.push_back( (TH1F*)file->Get(name+"_truth"));

    // now fill hist for every bin
    vector<double> masses = {169.5, 171.5, 173.5, 175.5};
    vector<TGraph*> binhists;

    int nbins = h_truth[0]->GetXaxis()->GetNbins();
    for(int bin=1; bin<=nbins; bin++){
      vector<double> contents;
      for(unsigned int i=0; i<h_truth.size(); i++){
        contents.push_back(h_truth[i]->GetBinContent(bin));
      }
      TGraph *h = new TGraph(h_truth.size(), &masses[0], &contents[0]);
      binhists.push_back(h);
    }

    vector<TF1*> fits;
    for(auto hist: binhists){
      TF1*f1 = new TF1("f1","pol1",0,500);
      hist->Fit("f1","R");
      TF1* fit = hist->GetFunction("f1");
      fits.push_back(fit);
    }
    PlotBinHists(binhists, fits, channel);

    TH1F* mtop1705 = CreateNewSample(fits, 170.5, h_truth[0]);
    TH1F* mtop1720 = CreateNewSample(fits, 172.0, h_truth[0]);
    TH1F* mtop1730 = CreateNewSample(fits, 173.0, h_truth[0]);
    TH1F* mtop1745 = CreateNewSample(fits, 174.5, h_truth[0]);
    PlotNewSample(mtop1705, "mtop1705", channel);
    PlotNewSample(mtop1720, "mtop1720", channel);
    PlotNewSample(mtop1730, "mtop1730", channel);
    PlotNewSample(mtop1745, "mtop1745", channel);


    TFile * f_out = new TFile("files/NewMasses_"+channel+".root","RECREATE");;
    mtop1705->Write("mc_mtop1705_truth");
    mtop1720->Write("mc_mtop1720_truth");
    mtop1730->Write("mc_mtop1730_truth");
    mtop1745->Write("mc_mtop1745_truth");
    f_out->Close();

    // draw


  }
  return 0;
}

TH1F* CreateNewSample(vector<TF1*> fits, double mass, TH1F* dummy){
  TH1F* sample = (TH1F*) dummy->Clone();
  sample->Reset();
  for(int i=0; i<fits.size(); i++){
    int bin = i+1;
    double content = fits[i]->Eval(mass);
    sample->SetBinContent(bin, content);
  }
  return sample;
}

void PlotNewSample(TH1F* hist, TString name, TString channel){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *c= new TCanvas(" ","",600,600);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle("m_{t} [GeV]");
  hist->GetYaxis()->SetTitle("bin content");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->Draw("HIST");

  gPad->RedrawAxis();
  TString savename = "/afs/desy.de/user/s/schwarzd/Plots/ExtrapolateMTop/"+channel+"_";
  savename += name;
  savename += ".pdf";
  c->SaveAs(savename);
  delete c;
}

void PlotBinHists(vector<TGraph*> binhists, vector<TF1*> fits, TString channel){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(int i=0; i<binhists.size(); i++){
    int binnr = i+1;
    TCanvas *c= new TCanvas(" ","",600,600);
    gPad->SetLeftMargin(0.15);
    binhists[i]->SetTitle(" ");
    binhists[i]->GetXaxis()->SetTitle("m_{t} [GeV]");
    binhists[i]->GetYaxis()->SetTitle("bin content");
    binhists[i]->GetYaxis()->SetTitleOffset(1.5);
    binhists[i]->GetYaxis()->SetNdivisions(505);
    binhists[i]->SetLineColor(kBlack);
    binhists[i]->SetMarkerColor(kBlack);
    binhists[i]->SetMarkerStyle(20);
    binhists[i]->Draw("AP");
    fits[i]->SetLineColor(kRed);
    fits[i]->Draw("SAME");
    gPad->RedrawAxis();
    TString savename = "/afs/desy.de/user/s/schwarzd/Plots/ExtrapolateMTop/"+channel+"_bin_";
    savename += binnr;
    savename += ".pdf";
    c->SaveAs(savename);
    delete c;
  }
}
