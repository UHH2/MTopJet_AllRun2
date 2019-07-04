#include "../include/CentralInclude.h"


using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);


int main(int argc, char* argv[]){
  TFile* file_mc = new TFile(dir_combine+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* file_data = new TFile(dir_combine+"uhh2.AnalysisModuleRunner.DATA.DATA.root");

  TH2F* h_2d =(TH2F*)file_mc->Get("XCone_cor/Mass_Vertices");
  TH2F* h_2d_data =(TH2F*)file_data->Get("XCone_cor/Mass_Vertices");

  TH2F* h_2dW =(TH2F*)file_mc->Get("XCone_cor_subjets/WMass_Vertices");
  TH2F* h_2dW_data =(TH2F*)file_data->Get("XCone_cor_subjets/WMass_Vertices");

  vector<TH1F*> h_1d, h_1d_data;
  vector<TH1F*> h_1dW, h_1dW_data;
  TH1F* hh = new TH1F("hh", "hh", 50, 0, 500);

  vector<double> NPV_bins = {0, 10, 15, 20, 30, 50};
  vector<double> mjet, mjet_e, mjet_data, mjet_data_e;
  vector<double> mW, mW_e, mW_data, mW_data_e;
  vector<double> NPV, NPV_e, NPV_data_e;

  // mjet
  for(unsigned int i=1; i<NPV_bins.size(); i++){
    TH1F* h = (TH1F*) hh->Clone();
    TH1F* h_data = (TH1F*) hh->Clone();
    TH1F* hW = (TH1F*) hh->Clone();
    TH1F* hW_data = (TH1F*) hh->Clone();
    double npv_min = NPV_bins[i-1];
    double npv_max = NPV_bins[i];
    double bincenter = npv_min + (npv_max - npv_min)/2;
    double binwidth = (npv_max - npv_min);
    NPV.push_back(bincenter);
    NPV_e.push_back(binwidth/2);
    NPV_data_e.push_back(0.0);
    int xbin_min = h_2d->GetXaxis()->FindBin(npv_min);
    int xbin_max = h_2d->GetXaxis()->FindBin(npv_max) - 1; // avoid double counting with '-1'
    int ybin_min = h_2d->GetYaxis()->FindBin(120);
    int ybin_max = h_2d->GetYaxis()->FindBin(240);
    int ybinW_min = h_2d->GetYaxis()->FindBin(65);
    int ybinW_max = h_2d->GetYaxis()->FindBin(95);
    // fill mjet
    for(int ybin=ybin_min; ybin<ybin_max; ybin++){
      double entry = 0;
      double entry_data = 0;
      for(int xbin=xbin_min; xbin<xbin_max; xbin++){
        entry += h_2d->GetBinContent(xbin, ybin);
        entry_data += h_2d_data->GetBinContent(xbin, ybin);
      }
      h->SetBinContent(ybin, entry);
      h_data->SetBinContent(ybin, entry_data);
    }
    mjet.push_back(h->GetMean());
    mjet_e.push_back(h->GetMeanError());
    h_1d.push_back(h);
    mjet_data.push_back(h_data->GetMean());
    mjet_data_e.push_back(h_data->GetMeanError());
    h_1d_data.push_back(h_data);
    // fill mW
    for(int ybinW=ybinW_min; ybinW<ybinW_max; ybinW++){
      double entry = 0;
      double entry_data = 0;
      for(int xbin=xbin_min; xbin<xbin_max; xbin++){
        entry += h_2dW->GetBinContent(xbin, ybinW);
        entry_data += h_2dW_data->GetBinContent(xbin, ybinW);
      }
      hW->SetBinContent(ybinW, entry);
      hW_data->SetBinContent(ybinW, entry_data);
    }
    mW.push_back(hW->GetMean());
    mW_e.push_back(hW->GetMeanError());
    h_1dW.push_back(hW);
    mW_data.push_back(hW_data->GetMean());
    mW_data_e.push_back(hW_data->GetMeanError());
    h_1dW_data.push_back(hW_data);
  }
  TGraphErrors * mjetPU = new TGraphErrors(NPV.size(), &NPV[0], &mjet[0], &NPV_e[0], &mjet_e[0]);
  TGraphErrors * mjetPU_data = new TGraphErrors(NPV.size(), &NPV[0], &mjet_data[0], &NPV_data_e[0], &mjet_data_e[0]);
  TGraphErrors * mWPU = new TGraphErrors(NPV.size(), &NPV[0], &mW[0], &NPV_e[0], &mW_e[0]);
  TGraphErrors * mWPU_data = new TGraphErrors(NPV.size(), &NPV[0], &mW_data[0], &NPV_data_e[0], &mW_data_e[0]);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<h_1d.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_1d[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_1d[i]->GetMaximum());
    h_1d[i]->GetXaxis()->SetNdivisions(505);
    h_1d[i]->GetYaxis()->SetNdivisions(505);
    h_1d[i]->GetYaxis()->SetTitleOffset(1.6);
    h_1d[i]->GetXaxis()->SetTitleOffset(1.3);
    h_1d[i]->SetTitle(" ");
    h_1d[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass/PUbin_";
    filename += i;
    filename += "_MC.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_1d_data.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_1d_data[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_1d_data[i]->GetMaximum());
    h_1d_data[i]->GetXaxis()->SetNdivisions(505);
    h_1d_data[i]->GetYaxis()->SetNdivisions(505);
    h_1d_data[i]->GetYaxis()->SetTitleOffset(1.6);
    h_1d_data[i]->GetXaxis()->SetTitleOffset(1.3);
    h_1d_data[i]->SetTitle(" ");
    h_1d_data[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass/PUbin_";
    filename += i;
    filename += "_DATA.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_1dW.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_1dW[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_1dW[i]->GetMaximum());
    h_1dW[i]->GetXaxis()->SetNdivisions(505);
    h_1dW[i]->GetYaxis()->SetNdivisions(505);
    h_1dW[i]->GetYaxis()->SetTitleOffset(1.6);
    h_1dW[i]->GetXaxis()->SetTitleOffset(1.3);
    h_1dW[i]->SetTitle(" ");
    h_1dW[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass/PUbin_";
    filename += i;
    filename += "_MC.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_1dW_data.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_1dW_data[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_1dW_data[i]->GetMaximum());
    h_1dW_data[i]->GetXaxis()->SetNdivisions(505);
    h_1dW_data[i]->GetYaxis()->SetNdivisions(505);
    h_1dW_data[i]->GetYaxis()->SetTitleOffset(1.6);
    h_1dW_data[i]->GetXaxis()->SetTitleOffset(1.3);
    h_1dW_data[i]->SetTitle(" ");
    h_1dW_data[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass/PUbin_";
    filename += i;
    filename += "_DATA.pdf";
    a->SaveAs(filename);
    delete a;
  }

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.13);
  mjetPU_data->Draw("AP");
  mjetPU_data->GetXaxis()->SetRangeUser(0, 50);
  mjetPU_data->GetYaxis()->SetRangeUser(0, 300);
  mjetPU_data->GetXaxis()->SetTitle("number of primary vertices");
  mjetPU_data->GetYaxis()->SetTitle("mean m_{jet} [GeV]");
  mjetPU_data->SetTitle(" ");
  mjetPU_data->GetXaxis()->SetNdivisions(505);
  mjetPU_data->GetYaxis()->SetNdivisions(505);
  mjetPU_data->GetYaxis()->SetTitleOffset(1.8);
  mjetPU_data->GetXaxis()->SetTitleOffset(1.3);
  mjetPU_data->SetLineColor(kBlack);
  mjetPU_data->SetMarkerColor(kBlack);
  mjetPU_data->SetMarkerStyle(8);
  mjetPU_data->SetMarkerSize(1);
  mWPU_data->SetLineColor(14);
  mWPU_data->SetMarkerColor(14);
  mWPU_data->SetMarkerStyle(8);
  mWPU_data->SetMarkerSize(1);
  mjetPU->SetLineColor(kRed);
  mjetPU->SetFillColor(kRed);
  mjetPU->SetMarkerColor(kRed);
  mjetPU->SetMarkerStyle(8);
  mjetPU->SetMarkerSize(0);
  mWPU->SetLineColor(kOrange+7);
  mWPU->SetFillColor(kOrange+7);
  mWPU->SetMarkerColor(kOrange+7);
  mWPU->SetMarkerStyle(8);
  mWPU->SetMarkerSize(0);

  mjetPU->Draw("E2 SAME");
  mjetPU_data->Draw("P SAME");
  mWPU->Draw("E2 SAME");
  mWPU_data->Draw("P SAME");
  CMSLabel(true, 0.2, 0.85);

  TLegend* leg = new TLegend(0.50, 0.65, 0.85, 0.85);
  leg->AddEntry(mjetPU_data,"t decay data","pe");
  leg->AddEntry(mjetPU,"t decay t#bar{t}","pf");
  leg->AddEntry(mWPU_data,"W decay data","pe");
  leg->AddEntry(mWPU,"W decay t#bar{t}","pf");
  leg->Draw();
  TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass/MjetVsPU.pdf";
  a->SaveAs(filename);
  delete a;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
