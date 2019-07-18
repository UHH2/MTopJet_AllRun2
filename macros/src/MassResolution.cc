#include "../include/CentralInclude.h"


using namespace std;
void PlotHist(TH1F* hist, TF1* fit, TString histname);

int main(int argc, char* argv[]){

  bool merged = false;
  if(argc == 2){
    if(strcmp(argv[1], "merged") == 0){
      merged = true;
    }
  }

  TString histdir = "/RecGenHists_GenSelRecInfo";

  if(merged){
    cout << "Only merged jets are considered" << endl;
    histdir = "/RecGenHists_GenSelRecInfo_matched";
  }
  else cout << "All jets are considered, you can also choose 'merged'." << endl;

  TFile* file_mc = new TFile(dir_combine+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  vector<TString> PU;
  PU.push_back("_lowPU");
  PU.push_back("_midPU");
  PU.push_back("_highPU");
  PU.push_back("");

  vector<TGraphErrors*> h_resos;
  vector<double> ptbins = {400, 450, 500, 600, 700, 800, 1500};
  for(auto pu: PU){
    vector<double> pt, pt_e, reso, reso_e;
    for(int i=0; i<ptbins.size()-1; i++){
      TString hdir = histdir+pu+"/";
      double binwidth = ptbins[i+1]-ptbins[i];
      pt.push_back(ptbins[i]+binwidth/2);
      pt_e.push_back(binwidth/2);
      TString histname = "MassResolution_pt";
      histname += i+1;
      TH1F* hist = (TH1F*)file_mc->Get(hdir+histname);
      TF1* f1 = new TF1("f1","gaus",-0.2,0.2);
      hist->Fit("f1","R");
      TF1* fit = hist->GetFunction("f1");
      reso.push_back(fit->GetParameter(2));
      reso_e.push_back(fit->GetParError(2));
      if(!merged) PlotHist(hist, fit, histname+pu);
      else        PlotHist(hist, fit, histname+"_merged"+pu);
    }
    int nbins = reso.size();
    TGraphErrors* h_reso = new TGraphErrors(nbins, &pt[0], &reso[0], &pt_e[0], &reso_e[0]);
    h_resos.push_back(h_reso);
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  h_resos[0]->Draw("AP");
  h_resos[0]->GetXaxis()->SetRangeUser(405, 1495);
  h_resos[0]->GetYaxis()->SetRangeUser(0.0, 0.15);
  h_resos[0]->GetXaxis()->SetTitle("p_{T, jet}^{gen} [GeV]");
  // h_resos[0]->GetYaxis()->SetTitle("RMS #left[ #frac{m_{jet}^{rec} - m_{jet}^{gen}}{m_{jet}^{gen}} #right]");
  h_resos[0]->GetYaxis()->SetTitle("jet mass resolution");
  h_resos[0]->SetTitle(" ");
  h_resos[0]->GetXaxis()->SetNdivisions(505);
  h_resos[0]->GetYaxis()->SetNdivisions(505);
  h_resos[0]->GetYaxis()->SetTitleOffset(1.4);
  h_resos[0]->GetXaxis()->SetTitleOffset(1.1);
  h_resos[0]->GetYaxis()->SetTitleSize(0.05);
  h_resos[0]->GetXaxis()->SetTitleSize(0.05);
  h_resos[0]->GetXaxis()->SetLabelSize(0.05);
  h_resos[0]->GetYaxis()->SetLabelSize(0.05);
  Color_t col[] = {kAzure+10, 798, kRed, kBlack};
  for(int i=0; i<h_resos.size(); i++){
    h_resos[i]->SetLineColor(col[i]);
    h_resos[i]->SetMarkerColor(col[i]);
    h_resos[i]->SetMarkerStyle(8);
    h_resos[i]->SetMarkerSize(1);
  }
  for(auto h: h_resos) h->Draw("P SAME");

  CMSSimLabel(true, 0.2, 0.85);

  TLegend *leg = new TLegend(0.55, 0.55, 0.85, 0.85);
  leg->AddEntry(h_resos[0], " 0 < NPV < 10", "ple");
  leg->AddEntry(h_resos[1], "10 < NPV < 20", "ple");
  leg->AddEntry(h_resos[2], "NPV > 20", "ple");
  leg->AddEntry(h_resos[3], "all", "ple");
  leg->Draw();
  TString filename;
  if(!merged) filename = "/afs/desy.de/user/s/schwarzd/Plots/MassResolution/reso.pdf";
  else        filename = "/afs/desy.de/user/s/schwarzd/Plots/MassResolution/reso_merged.pdf";
  a->SaveAs(filename);
  delete a;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}


void PlotHist(TH1F* hist, TF1* fit, TString histname){
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle("#left[ #frac{m_{jet}^{rec} - m_{jet}^{gen}}{m_{jet}^{gen}} #right]");
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetXaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetRangeUser(-0.5, 0.5);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
  hist->SetLineColor(kBlack);
  hist->SetFillColor(15);
  fit->SetLineWidth(3);
  fit->SetLineColor(kRed);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.17);
  hist->Draw("HIST");
  fit->Draw("SAME");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassResolution/"+histname+".pdf");
  delete A;
  return;

}
