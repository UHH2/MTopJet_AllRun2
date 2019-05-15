#include "../include/CentralInclude.h"

using namespace std;

// void SetHist(TH1F *hist);
void SetupGlobalStyle();

int main(int argc, char* argv[]){
  bool dosub = true;

  if(argc != 2){
    cout << "Chose 'sub' or 'fat'!" << endl;
    return 0;
  }

  if(strcmp(argv[1], "sub") == 0){
    dosub = true;
  }
  else if(strcmp(argv[1], "fat") == 0){
    dosub = false;
  }

  TFile* file = new TFile(dir_combine+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  TString histname = "Mass_HadJet33";
  TString dirall = "XCone_GEN_GenOnly/";
  TString dirmatch = "XCone_GEN_GenOnly_matched/";
  TString dirunmatch = "XCone_GEN_GenOnly_unmatched/";

  if(!dosub){
    dirmatch = "XCone_GEN_GenOnly_matched_fat/";
    dirunmatch = "XCone_GEN_GenOnly_unmatched_fat/";
  }

  TH1F * h1 = (TH1F*)file->Get(dirall+histname);
  TH1F * h1_m = (TH1F*)file->Get(dirmatch+histname);
  TH1F * h1_u = (TH1F*)file->Get(dirunmatch+histname);


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  SetupGlobalStyle();
  SetHist(h1);
  SetHist(h1_m);
  SetHist(h1_u);

  int xlow, xhigh, ylow;
  xlow = 0;
  xhigh = 400;
  ylow = 0;

  h1->SetTitle(" ");
  h1->GetXaxis()->SetTitle("m_{jet} [GeV]");
  h1->GetYaxis()->SetTitle("Events");
  h1->GetYaxis()->SetTitleOffset(1.8);
  h1->GetXaxis()->SetRangeUser(xlow, xhigh);
  h1->GetYaxis()->SetRangeUser(ylow, (h1->GetMaximum())*1.15);

  h1->SetLineColor(TColor::GetColor("#000000"));
  h1_m->SetLineColor(TColor::GetColor("#008ae6"));
  h1_u->SetLineColor(TColor::GetColor("#ff8000"));
  h1->SetLineWidth(3);
  h1_m->SetLineWidth(3);
  h1_u->SetLineWidth(3);


  TCanvas *A = new TCanvas("A", "A", 440, 400);

  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  Float_t y1, y2, y3;
  y3 = 0.99;
  y2 = y3-yplot;
  y1 = y2-yratio;
  Float_t x1, x2;
  x1 = 0.01;
  x2 = 0.99;

  TPad *pad = new TPad("pad", "Control Plots 2", x1, y1, x2, y3);
  pad->SetTopMargin(0.055);
  pad->SetBottomMargin(0.16);
  // pad->SetLeftMargin(0.19);//default
  pad->SetLeftMargin(0.19);
  //  pad->SetRightMargin(0.05);//default
  pad->SetRightMargin(0.045);

  pad->Draw();


  TLegend *leg = new TLegend(0.53,0.70,0.86,0.83, NULL,"brNDC");
  leg->AddEntry(h1,"t#bar{t} total","f");
  leg->AddEntry(h1_m,"t#bar{t} fully merged","l");
  leg->AddEntry(h1_u,"t#bar{t} not merged","l");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  h1->Draw("HIST");
  h1_u->Draw("HIST SAME");
  h1_m->Draw("HIST SAME");
  leg->Draw("");
  CMSSimLabel(false);
  LumiInfo();
  gPad->RedrawAxis();
  if(dosub) A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/GEN_match/XCone_GEN_matching_sub_combine.pdf");
  else      A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/GEN_match/XCone_GEN_matching_fat_combine.pdf");
  delete A;
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;

}
