#include "../include/CentralInclude.h"


using namespace std;


int main(int argc, char* argv[]){

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TString> histnames = {"hadtop_pt", "leptop_pt", "deltaR_hadtop_leptop", "deltaPhi_hadtop_leptop"};
  vector<TH1F*> rec, gen;

  for(auto name: histnames){
    rec.push_back( (TH1F*)TT_file->Get( "GenParticles_RecOnly/"+name) );
    gen.push_back( (TH1F*)TT_file->Get( "GenParticles_GenOnly/"+name) );
  }
  TH2F* rec_pthad_ptlep = (TH2F*)TT_file->Get( "GenParticles_RecOnly/pthadtop_ptleptop");
  TH2F* gen_pthad_ptlep = (TH2F*)TT_file->Get( "GenParticles_GenOnly/pthadtop_ptleptop");

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  rec[0]->GetXaxis()->SetTitle("had. top p_{T}");
  rec[1]->GetXaxis()->SetTitle("lep. top p_{T}");
  rec[2]->GetXaxis()->SetTitle("#Delta R(had. top, lep. top)");
  rec[3]->GetXaxis()->SetTitle("#Delta #Phi(had. top, lep. top)");

  for(unsigned int i=0; i<rec.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    rec[i]->SetTitle(" ");
    rec[i]->GetYaxis()->SetTitle("events");
    rec[i]->GetYaxis()->SetTitleSize(0.06);
    rec[i]->GetXaxis()->SetTitleSize(0.05);
    rec[i]->GetZaxis()->SetTitleSize(0.05);
    rec[i]->GetXaxis()->SetTitleOffset(0.9);
    rec[i]->GetYaxis()->SetTitleOffset(1.1);
    rec[i]->GetZaxis()->SetTitleOffset(0.9);
    rec[i]->GetXaxis()->SetNdivisions(505);
    rec[i]->GetYaxis()->SetNdivisions(505);
    rec[i]->SetFillColor(kRed);
    rec[i]->SetLineColor(1);
    rec[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/TopProperties/RecSel_"+histnames[i]+".pdf");
    delete a;
  }

  gen[0]->GetXaxis()->SetTitle("had. top p_{T}");
  gen[1]->GetXaxis()->SetTitle("lep. top p_{T}");
  gen[2]->GetXaxis()->SetTitle("#Delta R(had. top, lep. top)");
  gen[3]->GetXaxis()->SetTitle("#Delta #Phi(had. top, lep. top)");

  for(unsigned int i=0; i<gen.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    gen[i]->SetTitle(" ");
    gen[i]->GetYaxis()->SetTitle("events");
    gen[i]->GetYaxis()->SetTitleSize(0.06);
    gen[i]->GetXaxis()->SetTitleSize(0.05);
    gen[i]->GetZaxis()->SetTitleSize(0.05);
    gen[i]->GetXaxis()->SetTitleOffset(0.9);
    gen[i]->GetYaxis()->SetTitleOffset(1.1);
    gen[i]->GetZaxis()->SetTitleOffset(0.9);
    gen[i]->GetXaxis()->SetNdivisions(505);
    gen[i]->GetYaxis()->SetNdivisions(505);
    gen[i]->SetFillColor(kRed);
    gen[i]->SetLineColor(1);
    gen[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/TopProperties/GenSel_"+histnames[i]+".pdf");
    delete a;
  }


  TCanvas *b = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  // gPad->SetLogz();
  rec_pthad_ptlep->SetTitle(" ");
  rec_pthad_ptlep->GetXaxis()->SetTitle("had. top p_{T}");
  rec_pthad_ptlep->GetYaxis()->SetTitle("lep. top p_{T}");
  rec_pthad_ptlep->GetZaxis()->SetTitle("events");
  rec_pthad_ptlep->GetYaxis()->SetTitleSize(0.06);
  rec_pthad_ptlep->GetXaxis()->SetTitleSize(0.05);
  rec_pthad_ptlep->GetZaxis()->SetTitleSize(0.05);
  rec_pthad_ptlep->GetXaxis()->SetTitleOffset(0.9);
  rec_pthad_ptlep->GetYaxis()->SetTitleOffset(1.1);
  rec_pthad_ptlep->GetZaxis()->SetTitleOffset(0.9);
  rec_pthad_ptlep->GetXaxis()->SetNdivisions(505);
  rec_pthad_ptlep->GetYaxis()->SetNdivisions(505);
  rec_pthad_ptlep->Draw("COLZ");
  rec_pthad_ptlep->Draw("BOX SAME");
  b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/TopProperties/RecSel_pthad_ptlep.pdf");


  TCanvas *c = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  // gPad->SetLogz();
  gen_pthad_ptlep->SetTitle(" ");
  gen_pthad_ptlep->GetXaxis()->SetTitle("had. top p_{T}");
  gen_pthad_ptlep->GetYaxis()->SetTitle("lep. top p_{T}");
  gen_pthad_ptlep->GetZaxis()->SetTitle("events");
  gen_pthad_ptlep->GetYaxis()->SetTitleSize(0.06);
  gen_pthad_ptlep->GetXaxis()->SetTitleSize(0.05);
  gen_pthad_ptlep->GetZaxis()->SetTitleSize(0.05);
  gen_pthad_ptlep->GetXaxis()->SetTitleOffset(0.9);
  gen_pthad_ptlep->GetYaxis()->SetTitleOffset(1.1);
  gen_pthad_ptlep->GetZaxis()->SetTitleOffset(0.9);
  gen_pthad_ptlep->GetXaxis()->SetNdivisions(505);
  gen_pthad_ptlep->GetYaxis()->SetNdivisions(505);
  gen_pthad_ptlep->Draw("COLZ");
  gen_pthad_ptlep->Draw("BOX SAME");
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/TopProperties/GenSel_pthad_ptlep.pdf");
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
