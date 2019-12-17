#include "../include/CentralInclude.h"


using namespace std;
void PlotHists(vector<TH1F*>, TString, bool, TString);
TH1F* GetRatio(TH1F* h1, TH1F* h2);
vector<TString> years = {"2016v3", "2017v2", "2018"};


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  TString prefix_ttbar = "uhh2.AnalysisModuleRunner.MC.TTbar_";
  TString prefix_data = "uhh2.AnalysisModuleRunner.DATA.SingleMu_";


  vector<TString> histnames;
  histnames.push_back("Muon_sel/pt");
  histnames.push_back("Event_sel/MET");
  histnames.push_back("Jet_sel/number");
  histnames.push_back("ProbeJet_All_Pt400/pt");
  histnames.push_back("ProbeJet_All_Pt400/mass_sub");
  histnames.push_back("ProbeJet_All_Pt400/tau32");
  histnames.push_back("ProbeJet_pt400_wp3_all_pass/mass_sub");
  histnames.push_back("ProbeJet_pt400_wp3_all_fail/mass_sub");

  vector<TString> jets = {"PUPPI", "HOTVR"};

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  for(auto jet: jets){
    TString dir;
    if(jet == "HOTVR") dir = "/nfs/dust/cms/user/schwarzd/TopTagging/PostSel_HOTVR/";
    else               dir = "/nfs/dust/cms/user/schwarzd/TopTagging/PostSel/";
    vector<TFile*> f_data, f_ttbar;
    for(auto year: years){
      f_data.push_back(new TFile(dir+prefix_data+year+".root"));
      f_ttbar.push_back(new TFile(dir+prefix_ttbar+year+".root"));
    }

    for(auto hname: histnames){
      if(jet=="HOTVR"){
        if(hname == "ProbeJet_All_Pt400/tau32") hname += "_groomed";
        if(hname.Contains("ProbeJet_pt400_wp3_all_pass")) hname = "ProbeJet_pt400_all_pass/mass_sub";
        if(hname.Contains("ProbeJet_pt400_wp3_all_fail")) hname = "ProbeJet_pt400_all_fail/mass_sub";
      }
      vector<TH1F*> h_data, h_ttbar;
      for(auto file: f_data)   h_data.push_back((TH1F*) file->Get(hname));
      for(auto file: f_ttbar)  h_ttbar.push_back((TH1F*) file->Get(hname));
      if(hname.Contains("pt") || hname.Contains("mass") || hname.Contains("MET")){
        for(auto h: h_data) h->Rebin(5);
        for(auto h: h_ttbar) h->Rebin(5);
      }
      PlotHists(h_data, "data", true, jet+"__"+hname+"__NORM");
      PlotHists(h_ttbar, "tt", true, jet+"__"+hname+"__NORM");
      PlotHists(h_data, "data", false, jet+"__"+hname);
      PlotHists(h_ttbar, "tt", false, jet+"__"+hname);
    }
  }

  return 0;
}

void PlotHists(vector<TH1F*> hists_, TString datamc, bool norm, TString filename){

  vector<TH1F*> hists;
  for(auto h: hists_) hists.push_back((TH1F*) h->Clone());

  bool isDATA = false;
  if(datamc == "data") isDATA = true;

  if(norm){
    for(auto h: hists){
      double integral = h->Integral();
      h->Scale(1/integral);
    }
  }

  double ymax = 0.0;
  for(auto h: hists){
    if(h->GetMaximum() > ymax) ymax = h->GetMaximum();
  }

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();
  // TGaxis::SetMaxDigits(3);

  hists[0]->SetTitle(" ");
  if(norm) hists[0]->GetYaxis()->SetTitle("a.u.");
  else     hists[0]->GetYaxis()->SetTitle("events");
  hists[0]->GetYaxis()->SetTitleSize(0.06);
  hists[0]->GetXaxis()->SetTitleSize(0.05);
  hists[0]->GetZaxis()->SetTitleSize(0.05);
  hists[0]->GetXaxis()->SetTitleOffset(0.9);
  hists[0]->GetYaxis()->SetTitleOffset(1.1);
  hists[0]->GetZaxis()->SetTitleOffset(0.9);
  hists[0]->GetXaxis()->SetNdivisions(505);
  hists[0]->GetYaxis()->SetNdivisions(505);
  hists[0]->GetYaxis()->SetRangeUser(0, ymax*1.2);
  if(filename.Contains("tau32")) hists[0]->GetXaxis()->SetTitle("#tau_{32}");
  if(filename.Contains("pt")){
    if(filename.Contains("Muon")) hists[0]->GetXaxis()->SetTitle("muon p_{T}");
    else                          hists[0]->GetXaxis()->SetTitle("probe jet p_{T}");
  }
  if(filename.Contains("mass")) hists[0]->GetXaxis()->SetTitle("probe jet m_{jet}");
  if(filename.Contains("MET")) hists[0]->GetXaxis()->SetTitle("p_{T}^{miss}");
  if(filename.Contains("number")) hists[0]->GetXaxis()->SetTitle("number of AK4 jets");


  Color_t color[] = {kAzure+2, kRed-4, kBlack};

  for(unsigned int i=0; i<hists.size(); i++){
    if(isDATA){
      hists[i]->SetLineColor(color[i]);
      hists[i]->SetMarkerColor(color[i]);
      hists[i]->SetMarkerStyle(8);
      hists[i]->SetMarkerSize(1);
    }
    else{
      hists[i]->SetLineColor(color[i]);
      hists[i]->SetLineWidth(3);
    }
  }

  TString drawoption = "HIST";
  if(isDATA) drawoption = "P";
  hists[0]->Draw(drawoption);
  for(auto h: hists) h->Draw(drawoption+" SAME");

  TString legoption = "l";
  if(isDATA) legoption = "pel";

  TString type = "t#bar{t}";
  if(isDATA) type = "data";

  double x1 = 0.55;
  double y1 = 0.65;

  if(filename.Contains("tau32")) x1 = 0.2;

  TLegend* leg = new TLegend(x1, y1, x1+0.3, y1+0.2);
  for(unsigned int i=0; i<hists.size(); i++){
    leg->AddEntry(hists[i], type + " " + years[i], legoption);
  }
  leg->SetBorderSize(0);
  leg->Draw();

  double ymin = 0;
  if(filename.Contains("number")) ymin = -0.5;

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  hists[0]->GetXaxis()->SetLabelSize(0.);
  hists[0]->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( ymin, 0, ymin, ymax*1.2, 0, ymax*1.2, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(43);
  axis->SetLabelSize(21);
  axis->Draw();

  // Ratio Plot kommt in unteres pad
  a->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  vector<TH1F*> ratios;
  ratios.push_back(GetRatio(hists[0], hists[2]));
  ratios.push_back(GetRatio(hists[1], hists[2]));
  ratios.push_back(GetRatio(hists[2], hists[2]));

  ratios[0]->GetXaxis()->SetTickLength(0.07);
  ratios[0]->GetXaxis()->SetTitleSize(25);
  ratios[0]->GetXaxis()->SetTitleFont(43);
  ratios[0]->GetXaxis()->SetTitleOffset(4.0);
  ratios[0]->GetXaxis()->SetLabelFont(43);
  ratios[0]->GetXaxis()->SetLabelSize(21);
  ratios[0]->GetXaxis()->SetLabelOffset(0.035);
  if(isDATA) ratios[0]->GetYaxis()->SetTitle("#frac{Data}{2018 Data}");
  else       ratios[0]->GetYaxis()->SetTitle("#frac{MC}{2018 MC}");
  ratios[0]->GetYaxis()->CenterTitle();
  ratios[0]->GetYaxis()->SetTitleSize(22);
  ratios[0]->GetYaxis()->SetTitleFont(43);
  ratios[0]->GetYaxis()->SetTitleOffset(2.2);
  ratios[0]->GetYaxis()->SetLabelFont(43);
  ratios[0]->GetYaxis()->SetLabelSize(19);
  ratios[0]->GetYaxis()->SetLabelOffset(0.009);
  ratios[0]->GetYaxis()->SetNdivisions(505);
  ratios[0]->SetTitle(" ");
  ratios[0]->GetYaxis()->SetRangeUser(0.2, 1.8);
  if(filename.Contains("tau32")) ratios[0]->GetXaxis()->SetTitle("#tau_{32}");
  if(filename.Contains("pt")){
    if(filename.Contains("Muon")) ratios[0]->GetXaxis()->SetTitle("muon p_{T}");
    else                          ratios[0]->GetXaxis()->SetTitle("probe jet p_{T}");
  }
  if(filename.Contains("mass")) ratios[0]->GetXaxis()->SetTitle("probe jet m_{jet}");
  if(filename.Contains("MET")) ratios[0]->GetXaxis()->SetTitle("p_{T}^{miss}");
  if(filename.Contains("number")) ratios[0]->GetXaxis()->SetTitle("number of AK4 jets");



  ratios[0]->Draw(drawoption);

  for(unsigned int i=0; i<ratios.size(); i++){
    if(isDATA){
      ratios[i]->SetLineColor(color[i]);
      ratios[i]->SetMarkerColor(color[i]);
      ratios[i]->SetMarkerStyle(8);
      ratios[i]->SetMarkerSize(1);
    }
    else{
      ratios[i]->SetLineColor(color[i]);
      ratios[i]->SetLineWidth(3);
    }
  }

  for(auto h: ratios) h->Draw(drawoption+" SAME");

  filename.ReplaceAll("/", "__");
  if(isDATA) filename += "__DATA";
  else       filename += "__MC";
  a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/TopTagging/YearComparison/"+filename+".pdf");
  delete a;
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
      double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
