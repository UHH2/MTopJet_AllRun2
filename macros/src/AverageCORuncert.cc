#include "../include/CentralInclude.h"


using namespace std;

int main(int argc, char* argv[]){

  TString d = dir;
  TString channel = "muon";
  if(argc > 1){
    if(strcmp(argv[1], "elec") == 0){
      d = dir_elec;
      channel = "elec";
    }
  }

  TFile* file = new TFile(d+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* file_up = new TFile(d+"COR_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile* file_down = new TFile(d+"COR_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  vector<TFile*> files = {file, file_up, file_down};

  TString dirname = "XCone_cor_subjets/";
  vector<TString> histname = {"COR_factor_had", "COR_factor_lep"};

  vector<TH1F*> h1, h2;
  for(auto f: files){
    h1.push_back((TH1F*)f->Get(dirname+histname[0]));
    h2.push_back((TH1F*)f->Get(dirname+histname[1]));
  }

  vector<double> mean1, mean2;
  for(auto h: h1) mean1.push_back(h->GetMean());
  for(auto h: h2) mean2.push_back(h->GetMean());

  vector<TLine*> line1, line2;
  double ymax1 = h1[1]->GetMaximum() * 1.1;
  for(auto m: mean1) line1.push_back(new TLine(m, 0, m, ymax1));
  double ymax2 = h2[1]->GetMaximum() * 1.1;
  for(auto m: mean2) line2.push_back(new TLine(m, 0, m, ymax2));

  for(auto h: h1){
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->SetFillColor(13);
  }
  for(auto h: h2){
    h->SetLineWidth(4);
    h->SetLineColor(1);
    h->SetFillColor(13);
  }
  for(auto l: line1){
    l->SetLineWidth(4);
    l->SetLineColor(kRed+2);
    l->SetLineStyle(2);
  }
  for(auto l: line2){
    l->SetLineWidth(4);
    l->SetLineColor(kRed+2);
    l->SetLineStyle(2);
  }



gStyle->SetOptStat(kFALSE);
gStyle->SetPadTickY(1);
gStyle->SetPadTickX(1);
gStyle->SetLegendBorderSize(0);

for(unsigned int i=0; i<h1.size(); i++){
  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h1[i]->GetYaxis()->SetRangeUser(0, ymax1);
  h1[i]->GetXaxis()->SetNdivisions(505);
  h1[i]->GetYaxis()->SetNdivisions(505);
  h1[i]->GetYaxis()->SetTitleOffset(1.6);
  h1[i]->GetXaxis()->SetTitleOffset(1.3);
  h1[i]->SetTitle(" ");
  h1[i]->GetXaxis()->SetTitle("cor factor");
  h1[i]->GetYaxis()->SetTitle("events");
  h1[i]->Draw("HIST");
  line1[i]->Draw("SAME");
  TLegend *leg = new TLegend(0.53,0.63,0.88,0.88);
  TString meantitle = "mean = " + to_string(mean1[i]);
  leg->AddEntry(h1[i], "cor factor", "f");
  leg->AddEntry(line1[i], meantitle, "l");
  leg->Draw();
  TString name = "had";
  if(i==0) name += "central";
  else if(i==1) name += "up";
  else if(i==2) name += "down";
  a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction_allHad/AverageFactor/"+channel+"_"+name+".pdf");
  delete a;
}

for(unsigned int i=0; i<h2.size(); i++){
  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h2[i]->GetYaxis()->SetRangeUser(0, ymax1);
  h2[i]->GetXaxis()->SetNdivisions(505);
  h2[i]->GetYaxis()->SetNdivisions(505);
  h2[i]->GetYaxis()->SetTitleOffset(1.6);
  h2[i]->GetXaxis()->SetTitleOffset(1.3);
  h2[i]->SetTitle(" ");
  h2[i]->GetXaxis()->SetTitle("cor factor");
  h2[i]->GetYaxis()->SetTitle("events");
  h2[i]->Draw("HIST");
  line2[i]->Draw("SAME");
  TLegend *leg = new TLegend(0.53,0.63,0.88,0.88);
  TString meantitle = "mean = " + to_string(mean2[i]);
  leg->AddEntry(h2[i], "cor factor", "f");
  leg->AddEntry(line2[i], meantitle, "l");
  leg->Draw();
  TString name = "lep_";
  if(i==0) name += "central";
  else if(i==1) name += "up";
  else if(i==2) name += "down";
  a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction_allHad/AverageFactor/"+channel+"_"+name+".pdf");
  delete a;
}

//// ---------------------------------------------------------------------------------------------------------------------
//// ---------------------------------------------------------------------------------------------------------------------
return 0;
}
