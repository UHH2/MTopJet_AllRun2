#include "../include/CentralInclude.h"

using namespace std;

/* This File is only to compare the Distributions from 2016 and 2017 Data.
In 2017 Data and MC's show a discrepancy on the mass-axis.
If this file is used for later purposes, then it needs to be modified!
*/

void Plot(std::vector<TH1F*>, TString, bool, TString);
TH1F* GetRatio(TH1F* h1, TH1F* h2);
vector<TString> years = {"2016", "2017", "2018"};

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

int main(){
  cout << "Compare the following Distributions: \n";
  cout << "   -- Data           2016 & 17\n";
  cout << "   -- RecoLevel      2016 & 17\n";
  cout << "   -- GenLevel       2016 & 17" << endl;

  //============================================================================
  TString dir = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/";
  TString data = "/muon/uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TString ttbar = "/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root";

  vector<TFile*> f_data, f_ttbar;
  for(auto year: years){
    f_data.push_back(new TFile(dir+year+data));
    f_ttbar.push_back(new TFile(dir+year+ttbar));
  }

  //============================================================================
  vector<TString> histnames;
  histnames.push_back("MuonHists/pt");
  histnames.push_back("EventHists/MET");
  histnames.push_back("JetHits/pt_jet");
  histnames.push_back("XCone_cor_SF/M_jet1_");

  cout << histnames.size() << endl;
  cout << f_data.size() << endl;
  cout << f_ttbar.size() << endl;
  //============================================================================
  vector<TH1F*> h_data, h_ttbar;
  for(auto hname: histnames){
    cout << hname << endl;
    for(auto file: f_data) h_data.push_back((TH1F*) file->Get(hname));
    for(auto file: f_ttbar) h_ttbar.push_back((TH1F*) file->Get(hname));
    cout << h_data.size() << endl;
    cout << h_ttbar.size() << endl;
    Plot(h_data, "data", false, hname);
  }

  return 0;
}

//==============================================================================
/*
██████  ██       ██████  ████████
██   ██ ██      ██    ██    ██
██████  ██      ██    ██    ██
██      ██      ██    ██    ██
██      ███████  ██████     ██
*/


void Plot(vector<TH1F*> hists_, TString isDATA_, bool norm, TString filename){

  //--------------------- Set varibles -----------------------------------------
  vector<TH1F*> hists;
  for(auto h: hists_) hists.push_back((TH1F*) h->Clone());

  bool isDATA = false;
  if(isDATA_ == "data") isDATA = true;

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
  for(auto h: hists) cout << '1' << endl;
  //--------------------- Create plots -----------------------------------------
  cout << hists.size() << endl;
  TCanvas *a = new TCanvas("a", " ", 600, 600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  //--------------------- Set frame --------------------------------------------
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
  if(filename.Contains("pt")){
    if(filename.Contains("Muon")) hists[0]->GetXaxis()->SetTitle("muon p_{T}");
    else                          hists[0]->GetXaxis()->SetTitle("jet p_{T}");
  }
  if(filename.Contains("mass")) hists[0]->GetXaxis()->SetTitle("hadjet m_{jet}");
  if(filename.Contains("MET")) hists[0]->GetXaxis()->SetTitle("p_{T}^{miss}");

  //--------------------- Set diagramms ----------------------------------------
  Color_t color[] = {kAzure+2, kRed-4, kBlack};
  cout << "test2" << endl;
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
  cout << "test2" << endl;
  //--------------------- Set legend -- ----------------------------------------
  TLegend* leg = new TLegend(x1, y1, x1+0.3, y1+0.2);
  cout << "test3" << endl;
  cout << hists.size() << endl;
  for(unsigned int i=0; i<hists.size(); i++){
    cout << "test4" << endl;
    leg->AddEntry(hists[i], type + " " + years[i], legoption);
    cout << "test4.1" << endl;
  }
  leg->SetBorderSize(0);
  leg->Draw();
  cout << hists.size() << endl;
  double ymin = 0;

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  hists[0]->GetXaxis()->SetLabelSize(0.);
  hists[0]->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( ymin, 0, ymin, ymax*1.2, 0, ymax*1.2, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(43);
  axis->SetLabelSize(21);
  axis->Draw();
  cout << "test2" << endl;
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
  if(filename.Contains("pt")){
    if(filename.Contains("Muon")) ratios[0]->GetXaxis()->SetTitle("muon p_{T}");
    else                          ratios[0]->GetXaxis()->SetTitle("hadjet p_{T}");
  }
  if(filename.Contains("mass")) ratios[0]->GetXaxis()->SetTitle("hadjet m_{jet}");
  if(filename.Contains("MET")) ratios[0]->GetXaxis()->SetTitle("p_{T}^{miss}");

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

  filename.ReplaceAll("/", "_");
  if(isDATA) filename += "_DATA";
  else       filename += "__MC";
  a->SaveAs("/afs/desy.de/user/p/paaschal/Plots/YearComparison/"+filename+".pdf");
  delete a;
}

//==============================================================================
/*
.██████  ███████ ████████     ██████   █████  ████████ ██  ██████
██       ██         ██        ██   ██ ██   ██    ██    ██ ██    ██
██   ███ █████      ██        ██████  ███████    ██    ██ ██    ██
██    ██ ██         ██        ██   ██ ██   ██    ██    ██ ██    ██
.██████  ███████    ██        ██   ██ ██   ██    ██    ██  ██████
*/

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
