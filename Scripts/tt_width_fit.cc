#include "TMath.h"
#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include <TTree.h>
#include <TH1.h>
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include <TStyle.h>

using namespace std;

Double_t Breit_Wigner(Double_t *x, Double_t *par){

  Double_t gamma = sqrt(par[0]*par[0]*(par[0]*par[0]+par[1]*par[1]));
  Double_t k = (2*sqrt(2)*par[0]*par[1]*gamma) / (TMath::Pi()*sqrt(par[0]*par[0]+gamma));
  Double_t value = par[2] * k/((x[0]*x[0]-par[0]*par[0])*(x[0]*x[0]-par[0]*par[0])+par[0]*par[0]*par[1]*par[1]);


  return value;
}

void tt_width_fit()
{
  TFile *TTbar_SM = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *TTbar_2G = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_2gamma.root");
  TFile *TTbar_4G = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_4gamma.root");
  TFile *TTbar_8G = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_8gamma.root");



  TH1F * SM_width = (TH1F*)TTbar_SM->Get("GenParticles_SM/hadtop_mass");
  TH1F * G2_width = (TH1F*)TTbar_2G->Get("GenParticles_newWidth/hadtop_mass");
  TH1F * G4_width = (TH1F*)TTbar_4G->Get("GenParticles_newWidth/hadtop_mass");
  TH1F * G8_width = (TH1F*)TTbar_8G->Get("GenParticles_newWidth/hadtop_mass");

  TH1F * JetMass_SM = (TH1F*)TTbar_SM->Get("XCone_cor/M_jet1");
  TH1F * JetMass_2G = (TH1F*)TTbar_2G->Get("XCone_cor/M_jet1");
  TH1F * JetMass_4G = (TH1F*)TTbar_4G->Get("XCone_cor/M_jet1");
  TH1F * JetMass_8G = (TH1F*)TTbar_8G->Get("XCone_cor/M_jet1");


  Double_t events_sm = SM_width->Integral();
  Double_t events_2G = G2_width->Integral();
  Double_t events_4G = G4_width->Integral();
  Double_t events_8G = G8_width->Integral();


  SM_width->Scale(1/events_sm);
  G2_width->Scale(1/events_2G);
  G4_width->Scale(1/events_4G);
  G8_width->Scale(1/events_8G);


  TF1 *func = new TF1("Breit_Wigner",Breit_Wigner,000,400,3);
  func->SetParameters(172.5,1.324,1);
  func->SetLineWidth(3);
  func->SetLineColor(2);
  // SM_width->Fit("Breit_Wigner", "r");
  TF1 *func_2G = new TF1("Breit_Wigner",Breit_Wigner,000,400,3);
  func_2G->SetParameters(172.5,2.648,1);
  func_2G->SetLineWidth(3);
  func_2G->SetLineColor(2);
  // G2_width->Fit("Breit_Wigner", "r");
  TF1 *func_4G = new TF1("Breit_Wigner",Breit_Wigner,000,400,3);
  func_4G->SetParameters(172.5,5.296,1);
  func_4G->SetLineWidth(3);
  func_4G->SetLineColor(2);
  // G4_width->Fit("Breit_Wigner", "r");
  TF1 *func_8G = new TF1("Breit_Wigner",Breit_Wigner,000,400,3);
  func_8G->SetParameters(172.5,10.592,1);
  func_8G->SetLineWidth(3);
  func_8G->SetLineColor(2);
  // G8_width->Fit("Breit_Wigner", "r");

  double xlow, xhigh;
  double ylow, yhigh;
  xlow = 150;
  xhigh = 200;
  ylow = 0;
  for(auto i: {SM_width, G2_width, G4_width, G8_width}){
    yhigh = (i->GetMaximum())*1.1;
    i->SetTitle(" ");
    i->GetXaxis()->SetTitle("Top Mass [GeV]");
    i->GetYaxis()->SetTitle("#Delta N / N");
    i->GetYaxis()->SetRangeUser(ylow, yhigh);
    i->GetXaxis()->SetRangeUser(xlow, xhigh);
    i->GetXaxis()->SetNdivisions(505);
    i->GetYaxis()->SetNdivisions(505);
    i->GetYaxis()->SetTitleOffset(1.7);

    i->SetFillColor(kGray);
    i->SetLineColor(kBlack);
  }
  
  for(auto i: {JetMass_SM, JetMass_2G, JetMass_4G, JetMass_8G}){
    i->SetLineWidth(3);
    i->SetTitle(" ");
    i->GetXaxis()->SetTitle("Jet Mass [GeV]");
    i->GetYaxis()->SetTitle("events");
    i->GetYaxis()->SetRangeUser(0, 1300);
    i->GetXaxis()->SetRangeUser(100, 300);
    i->GetXaxis()->SetNdivisions(505);
    i->GetYaxis()->SetNdivisions(505);
    i->GetYaxis()->SetTitleOffset(1.7);
    i->SetLineWidth(4);

  }

  JetMass_SM->SetLineColor(kGray);
  JetMass_2G->SetLineColor(kRed);
  JetMass_4G->SetLineColor(1);
  JetMass_8G->SetLineColor(kAzure+7);
  // JetMass_SM->SetLineStyle(1);
  // JetMass_2G->SetLineStyle(2);
  // JetMass_4G->SetLineStyle(1);
  // JetMass_8G->SetLineStyle(2);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLegend *leg = new TLegend(0.55,0.65,0.86,0.8);
  leg->AddEntry(SM_width,"t#bar{t} (#Gamma_{SM})","f");
  leg->AddEntry(func,"Breit Wigner","l");
  SM_width->Draw("HIST");
  func->Draw("SAME");
  leg->Draw("");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Width/BreitWignerFit_SM.pdf");

  TCanvas *B= new TCanvas("B", "B", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLegend *leg2 = new TLegend(0.55,0.65,0.86,0.8);
  leg2->AddEntry(G2_width,"t#bar{t} (2 #Gamma_{SM})","f");
  leg2->AddEntry(func_2G,"Breit Wigner","l");
  G2_width->Draw("HIST");
  func_2G->Draw("SAME");
  leg2->Draw("");
  gPad->RedrawAxis();
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Width/BreitWignerFit_2gamma.pdf");

  TCanvas *C = new TCanvas("C", "C", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLegend *leg3 = new TLegend(0.55,0.65,0.86,0.8);
  leg3->AddEntry(G4_width,"t#bar{t} (4 #Gamma_{SM})","f");
  leg3->AddEntry(func_4G,"Breit Wigner","l");
  G4_width->Draw("HIST");
  func_4G->Draw("SAME");
  leg3->Draw("");
  gPad->RedrawAxis();
  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Width/BreitWignerFit_4gamma.pdf");

  TCanvas *D = new TCanvas("D", "D", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLegend *leg4 = new TLegend(0.55,0.65,0.86,0.8);
  leg4->AddEntry(G8_width,"t#bar{t} (8 #Gamma_{SM})","f");
  leg4->AddEntry(func_8G,"Breit Wigner","l");
  G8_width->Draw("HIST");
  func_8G->Draw("SAME");
  leg4->Draw("");
  gPad->RedrawAxis();
  D->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Width/BreitWignerFit_8gamma.pdf");

  TCanvas *E = new TCanvas("E", "E", 600, 600);
  gPad->SetLeftMargin(0.15);
  JetMass_8G->Draw("HIST SAME");
  JetMass_4G->Draw("HIST SAME");
  JetMass_2G->Draw("HIST SAME");
  JetMass_SM->Draw("HIST SAME");
  TLegend *leg5 = new TLegend(0.55,0.65,0.86,0.8);
  leg5->AddEntry(JetMass_SM,"t#bar{t} (#Gamma_{SM})","l");
  leg5->AddEntry(JetMass_2G,"t#bar{t} (2 #Gamma_{SM})","l");
  leg5->AddEntry(JetMass_4G,"t#bar{t} (4 #Gamma_{SM})","l");
  leg5->AddEntry(JetMass_8G,"t#bar{t} (8 #Gamma_{SM})","l");
  leg5->Draw("");
  E->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Width/JetMass.pdf");


}





