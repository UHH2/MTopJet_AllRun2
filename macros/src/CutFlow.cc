#include "../include/CentralInclude.h"
#include "../include/Utils.h"

using namespace std;

void Plot(TGraphAsymmErrors* elec, TGraphAsymmErrors* elec_stat, TGraphAsymmErrors* muon, TGraphAsymmErrors* muon_stat);
TGraphAsymmErrors* ConvertToGraph(TH1F*, TString);
TString year;

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  if(argc < 2){
    cout << "[ERROR] Specify the year!" << endl;
    cout << "./ChannelCompare <year>" << endl;
    return 1;
  }

  if(strcmp(argv[1], "2016") == 0) year = "2016";
  else if(strcmp(argv[1], "2017") == 0) year = "2017";
  else if(strcmp(argv[1], "2018") == 0) year = "2018";


  TFile* f_elec = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_data_"+year+"_elec.root");
  TFile* f_muon = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_data_"+year+"_muon.root");

  TH1F * h_elec = (TH1F*)f_elec->Get("Unfold_XS_totuncert");
  TH1F * h_elec_stat = (TH1F*)f_elec->Get("Unfold_XS_statuncert");
  TH1F * h_muon = (TH1F*)f_muon->Get("Unfold_XS_totuncert");
  TH1F * h_muon_stat = (TH1F*)f_muon->Get("Unfold_XS_statuncert");

  TGraphAsymmErrors* elec = ConvertToGraph(h_elec, "left");
  TGraphAsymmErrors* elec_stat = ConvertToGraph(h_elec_stat, "left");
  TGraphAsymmErrors* muon = ConvertToGraph(h_muon, "right");
  TGraphAsymmErrors* muon_stat = ConvertToGraph(h_muon_stat, "right");
  Plot(elec, elec_stat, muon, muon_stat);
}

void Plot(TGraphAsymmErrors* elec, TGraphAsymmErrors* elec_stat, TGraphAsymmErrors* muon, TGraphAsymmErrors* muon_stat){
  TString save_path = get_save_path();
  creat_folder(save_path+"/Plots/Unfolding_Run2/");

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  elec->Draw("AP");
  elec->SetTitle(" ");
  elec->GetXaxis()->SetRangeUser(112, 232);
  elec->GetYaxis()->SetRangeUser(0, 11);
  elec->GetXaxis()->SetTitle("#it{m}_{jet} [GeV]");
  elec->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{m}_{jet}} [#frac{fb}{GeV}]");
  elec->GetYaxis()->SetTitleOffset(1.1);
  elec->GetXaxis()->SetTitleOffset(0.9);
  elec->GetYaxis()->SetTitleSize(0.05);
  elec->GetXaxis()->SetTitleSize(0.05);
  elec->GetYaxis()->SetNdivisions(505);
  elec->SetLineColor(17);
  elec->SetFillColor(17);
  elec->SetMarkerColor(kBlack);
  elec->SetMarkerStyle(8);
  elec->SetMarkerSize(1);
  elec_stat->SetFillColor(15);
  elec_stat->SetMarkerColor(kBlack);
  elec_stat->SetMarkerStyle(8);
  elec_stat->SetMarkerSize(1);
  elec->Draw("E2 SAME");
  elec_stat->Draw("E2 SAME");
  elec_stat->Draw("PX SAME");

  muon->SetLineColor(kRed+1);
  muon->SetMarkerColor(kRed+1);
  muon->SetMarkerStyle(8);
  muon->SetMarkerSize(1);
  muon->Draw("P SAME");
  muon_stat->SetLineColor(kRed+1);
  muon_stat->SetMarkerColor(kRed+1);
  muon_stat->SetMarkerStyle(8);
  muon_stat->SetMarkerSize(1);
  muon_stat->Draw("P SAME");
  muon->Draw("P SAME");

  gStyle->SetEndErrorSize(5);
  gPad->RedrawAxis();

  CMSLabel(true, 0.2, 0.85);

  TString infotext = "35.9 fb^{-1} (13 TeV)";
  if(year == "2017")      infotext = "41.5 fb^{-1} (13 TeV)";
  else if(year == "2018") infotext = "59.7 fb^{-1} (13 TeV)";
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(0.9);
  text1->SetY(0.941);
  text1->SetTextSize(0.04);
  text1->Draw();

  TLegend *l=new TLegend(0.55,0.68,0.85,0.87);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(elec,"Electron channel","pf");
  l->AddEntry(muon,"Muon channel","pe");
  l->SetTextSize(0.03);
  l->Draw();
  c->SaveAs(save_path+"/Plots/Unfolding_Run2/ChannelComparison_"+year+".pdf");
  delete c;
}


TGraphAsymmErrors* ConvertToGraph(TH1F* hist, TString leftright){
  vector<double> x,y,xel,xer,ye;
  int nbins = hist->GetXaxis()->GetNbins();
  for(int bin = 1; bin<=nbins; bin++){
    double x_ = hist->GetBinCenter(bin);
    if(leftright == "left") x_ -= 2.0;
    if(leftright == "right") x_ += 2.0;
    // cout << x_ << ", " << hist->GetBinContent(bin) << endl;
    x.push_back(x_);
    double xel_ = hist->GetBinWidth(bin)/2;
    double xer_ = hist->GetBinWidth(bin)/2;
    if(leftright == "left"){
      xel_ -= 2.0;
      xer_ += 2.0;
    }
    if(leftright == "right"){
      // xel_ += 2.0;
      // xer_ -= 2.0;
      xel_ = 0.0;
      xer_ = 0.0;
    }
    xel.push_back(xel_);
    xer.push_back(xer_);
    y.push_back(hist->GetBinContent(bin));
    ye.push_back(hist->GetBinError(bin));
  }
  TGraphAsymmErrors *g = new TGraphAsymmErrors(nbins, &x[0], &y[0], &xel[0], &xer[0], &ye[0], &ye[0]);
  return g;
}
