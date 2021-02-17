#include "../include/CentralInclude.h"

using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);
TH1F* GetRatioUncert(TH1F* mc);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/


int main(int argc, char* argv[]){
  TString save_path = get_save_path();
  TString year = "FullRun2";
  if(argc > 1){
    if(strcmp(argv[1], "2016") == 0)      year = "2016";
    else if(strcmp(argv[1], "2017") == 0) year = "2017";
    else if(strcmp(argv[1], "2018") == 0) year = "2018";
  }

  TFile* file = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_data_"+year+"_combine.root");

  TH1F * h_data = (TH1F*)file->Get("Unfold_norm_totuncert");
  TH1F * h_data_stat = (TH1F*)file->Get("Unfold_norm_statuncert");
  TH1F * h_truth = (TH1F*)file->Get("mc_truth_norm");
  TH1F * h_truth2 = (TH1F*)file->Get("mc_mtop1695_truth_norm");
  TH1F * h_truth3 = (TH1F*)file->Get("mc_mtop1755_truth_norm");

  TH1F * h_ratio_data = GetRatioUncert(h_data);
  TH1F * h_ratio_data_stat = GetRatioUncert(h_data_stat);
  TH1F * h_ratio_powheg = GetRatio(h_truth, h_data);
  TH1F * h_ratio_1695 = GetRatio(h_truth2, h_data);
  TH1F * h_ratio_1755 = GetRatio(h_truth3, h_data);

  std::vector<TH1D*> truth, truth2;
  for(auto t: {h_truth2, h_truth, h_truth3}){
    truth.push_back( (TH1D*) t->Clone() );
    truth2.push_back( (TH1D*) t->Clone() );
  }
  vector<TString> legnames = {"#it{m}_{t} = 169.5 GeV", "#it{m}_{t} = 172.5 GeV", "#it{m}_{t} = 175.5 GeV"};




  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  bool doTicks = true;

  TString o_tot = "E1";
  TString o_stat = "E1";

  if(!doTicks){
    o_tot = "E";
    o_stat = "E1 X0";
  }


  int titlefont = 43;
  double ymax = 0.04;

  TCanvas *c = new TCanvas("c","",600,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();


  h_data->SetTitle(" ");
  h_data->GetYaxis()->SetRangeUser(0., ymax);
  h_data->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#it{m}_{jet}} #left[#frac{1}{GeV}#right]");
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetNdivisions(505);
  h_data->SetLineColor(kBlack);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);
  h_data->Draw(o_tot);
  h_data_stat->SetLineColor(kBlack);
  h_data_stat->SetMarkerColor(kBlack);
  h_data_stat->SetMarkerStyle(8);
  h_data_stat->SetMarkerSize(1);
  gStyle->SetEndErrorSize(5);
  Int_t color[] = {kRed-4, TColor::GetColor("#0059b3"), 14};
  Int_t fillcolor[] = {kRed-10, TColor::GetColor("#99ccff"), 16};
  Int_t style[] = {3, 1, 2};

  for(unsigned int i=0; i<truth.size(); i++){
    truth[i]->SetLineWidth(3);
    truth[i]->SetLineColor(color[i]);
    truth[i]->SetLineStyle(style[i]);
    truth[i]->SetMarkerStyle(0);
    truth[i]->SetFillColor(fillcolor[i]);
    truth[i]->Draw("E2 SAME");
    truth2[i]->SetLineWidth(3);
    truth2[i]->SetLineColor(color[i]);
    truth2[i]->SetLineStyle(style[i]);
    truth2[i]->Draw("HIST SAME");
  }
  truth[1]->Draw("E2 SAME");
  truth2[1]->Draw("HIST SAME");

  h_data_stat->Draw(o_stat+" SAME");
  h_data->Draw(o_tot+" SAME");

  TLegend *l=new TLegend(0.56,0.6,0.82,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(h_data,"Data","ple");
  for(unsigned int i=0; i<truth.size(); i++){
    l->AddEntry(truth[i],legnames[i],"fl");
  }
  l->SetTextSize(0.052);
  l->Draw();
  gPad->RedrawAxis();

  TLatex *cmstext = new TLatex(3.5, 24, "CMS");
  cmstext->SetNDC();
  cmstext->SetTextAlign(13);
  cmstext->SetTextFont(62);
  cmstext->SetTextSize(0.08);
  cmstext->SetX(0.24);
  cmstext->SetY(0.84);
  cmstext->Draw();

  TString simtext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, simtext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.24);
  text3->SetTextFont(52);
  text3->SetTextSize(0.06);
  text3->SetY(0.77);
  text3->Draw();


  TString infotext = "137.1 fb^{-1} (13 TeV)";
  if(year == "2016") infotext = "35.9 fb^{-1} (13 TeV)";
  else if(year == "2017") infotext = "41.5 fb^{-1} (13 TeV)";
  else if(year == "2018") infotext = "59.7 fb^{-1} (13 TeV)";

  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(0.9);
  text1->SetY(0.961);
  text1->SetTextSize(0.055);
  text1->Draw();

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  h_data->GetXaxis()->SetLabelSize(0.);
  h_data->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( 112, 0, 112, ymax, 0, ymax, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(titlefont);
  axis->SetLabelSize(21);
  axis->Draw();

  // Ratio Plot kommt in unteres pad
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.38);
  pad2->Draw();
  pad2->cd();

  h_ratio_data->GetXaxis()->SetTickLength(0.07);
  h_ratio_data->GetXaxis()->SetTitleSize(25);
  h_ratio_data->GetXaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetXaxis()->SetTitleOffset(4.0);
  h_ratio_data->GetXaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetXaxis()->SetLabelSize(21);
  h_ratio_data->GetXaxis()->SetLabelOffset(0.035);
  h_ratio_data->GetYaxis()->SetTitle("#frac{Theory}{Data}");
  h_ratio_data->GetYaxis()->CenterTitle();
  h_ratio_data->GetYaxis()->SetTitleSize(22);
  h_ratio_data->GetYaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetYaxis()->SetTitleOffset(2.2);
  h_ratio_data->GetYaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetYaxis()->SetLabelSize(19);
  h_ratio_data->GetYaxis()->SetLabelOffset(0.009);
  h_ratio_data->GetYaxis()->SetNdivisions(505);
  h_ratio_data->SetTitle(" ");
  h_ratio_data->GetYaxis()->SetRangeUser(0.2, 1.8);
  h_ratio_data->GetXaxis()->SetTitle("#it{m}_{jet} #left[GeV#right]");
  h_ratio_data->SetLineColor(kBlack);
  h_ratio_data->SetMarkerColor(kBlack);
  h_ratio_data->SetMarkerStyle(8);
  h_ratio_data->SetMarkerSize(1);
  h_ratio_data->Draw(o_tot);
  h_ratio_data_stat->SetLineColor(kBlack);
  h_ratio_data_stat->SetMarkerColor(kBlack);
  h_ratio_data_stat->SetMarkerStyle(8);
  h_ratio_data_stat->SetMarkerSize(1);

  TH1F* h_ratio_powheg_line = (TH1F*) h_ratio_powheg->Clone();
  h_ratio_powheg->SetMarkerSize(0);
  h_ratio_powheg->SetFillColor(fillcolor[1]);
  h_ratio_powheg->Draw("E2 SAME");
  h_ratio_1695->SetLineWidth(3);
  h_ratio_1695->SetLineColor(color[0]);
  h_ratio_1695->SetLineStyle(style[0]);
  h_ratio_1695->Draw("HIST SAME");
  h_ratio_1755->SetLineWidth(3);
  h_ratio_1755->SetLineColor(color[2]);
  h_ratio_1755->SetLineStyle(style[2]);
  h_ratio_1755->Draw("HIST SAME");
  h_ratio_data->Draw(o_tot+" SAME");

  gStyle->SetEndErrorSize(5);

  h_ratio_powheg_line->SetLineColor(color[1]);
  h_ratio_powheg_line->SetLineWidth(3);
  h_ratio_powheg_line->SetLineStyle(style[1]);


  // noch einmal alles in der richtigen Reihenfolge zeichnen
  h_ratio_powheg->Draw("E2 SAME");
  h_ratio_powheg_line->Draw("HIST SAME");
  h_ratio_1695->Draw("HIST SAME");
  h_ratio_data_stat->Draw(o_stat+" SAME");
  h_ratio_data->Draw(o_tot+" SAME");
  gPad->RedrawAxis();

  c->SaveAs(save_path+"/Plots/Unfolding_Run2/Unfold_norm_"+year+".pdf");
  delete c;
}


TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      // double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      double error = E1/N2; // only consider uncert from MC
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}

TH1F* GetRatioUncert(TH1F* mc){
  TH1F* mc_uncert = (TH1F*) mc->Clone();
  int Nbins = mc->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    mc_uncert->SetBinContent(i, 1.0);
    double central = mc->GetBinContent(i);
    double error = mc->GetBinError(i);
    double error_ratio = error/central;
    mc_uncert->SetBinError(i, error_ratio);
  }
  return mc_uncert;
}
