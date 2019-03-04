#include "../include/CentralInclude.h"


using namespace std;


int main(int argc, char* argv[]){

  TString pdfname = "jec";
  TString mode;
  if(argc >1 && strcmp(argv[1], "raw") == 0){
    mode = "_noJEC";
    pdfname = "raw";
  }
  else if(argc >1 && strcmp(argv[1], "cor") == 0){
    mode = "_corrected";
    pdfname = "cor";
  }
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *f_jecup = new TFile(dir+"JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *f_jecdown = new TFile(dir+"JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  vector<TString> histnames = {"min_mass_Wjet", "min_mass_Wjet_zoom", "min_mass_Wjet_ptlow", "min_mass_Wjet_pthigh"};
  vector<TH1F*> rec, gen;

  for(auto name: histnames){
    rec.push_back( (TH1F*)TT_file->Get( "XCone_"+pdfname+"_subjets/"+name) );
    gen.push_back( (TH1F*)TT_file->Get( "XCone_GEN_Sel_measurement/"+name) );

  }

  TH1F* reso = (TH1F*)TT_file->Get( "RecGenHists_subjets"+mode+"/WMassResolution" );
  TH1F* reso_up = (TH1F*)f_jecup->Get( "RecGenHists_subjets"+mode+"/WMassResolution" );
  TH1F* reso_down = (TH1F*)f_jecdown->Get( "RecGenHists_subjets"+mode+"/WMassResolution" );

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  for(unsigned int i=0; i<rec.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    rec[i]->SetTitle(" ");
    rec[i]->GetXaxis()->SetTitle("min m_{ij}^{rec}");
    rec[i]->GetYaxis()->SetTitle("events");
    rec[i]->GetYaxis()->SetTitleSize(0.06);
    rec[i]->GetXaxis()->SetTitleSize(0.05);
    rec[i]->GetZaxis()->SetTitleSize(0.05);
    rec[i]->GetXaxis()->SetTitleOffset(0.9);
    rec[i]->GetYaxis()->SetTitleOffset(1.1);
    rec[i]->GetZaxis()->SetTitleOffset(0.9);
    rec[i]->GetXaxis()->SetNdivisions(505);
    rec[i]->GetYaxis()->SetNdivisions(505);
    rec[i]->GetXaxis()->SetRangeUser(0, 4);
    rec[i]->SetFillColor(kRed);
    rec[i]->SetLineColor(1);
    rec[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Wmass/REC_"+histnames[i]+"_"+pdfname+".pdf");
    delete a;
  }

  for(unsigned int i=0; i<rec.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    gen[i]->SetTitle(" ");
    gen[i]->GetXaxis()->SetTitle("min m_{ij}^{gen}");
    gen[i]->GetYaxis()->SetTitle("events");
    gen[i]->GetYaxis()->SetTitleSize(0.06);
    gen[i]->GetXaxis()->SetTitleSize(0.05);
    gen[i]->GetZaxis()->SetTitleSize(0.05);
    gen[i]->GetXaxis()->SetTitleOffset(0.9);
    gen[i]->GetYaxis()->SetTitleOffset(1.1);
    gen[i]->GetZaxis()->SetTitleOffset(0.9);
    gen[i]->GetXaxis()->SetNdivisions(505);
    gen[i]->GetYaxis()->SetNdivisions(505);
    gen[i]->GetXaxis()->SetRangeUser(0, 4);
    gen[i]->SetFillColor(13);
    gen[i]->SetLineColor(1);
    gen[i]->Draw("HIST");
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Wmass/GEN_"+histnames[i]+".pdf");
    delete a;
  }

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  reso->SetTitle(" ");
  reso->GetXaxis()->SetTitle("#frac{m^{rec}_{W} - m^{gen}_{W}}{m^{gen}_{W}}");
  reso->GetYaxis()->SetTitle("events");
  reso->GetYaxis()->SetTitleSize(0.06);
  reso->GetXaxis()->SetTitleSize(0.05);
  reso->GetZaxis()->SetTitleSize(0.05);
  reso->GetXaxis()->SetTitleOffset(1.3);
  reso->GetYaxis()->SetTitleOffset(1.1);
  reso->GetZaxis()->SetTitleOffset(0.9);
  reso->GetXaxis()->SetNdivisions(505);
  reso->GetYaxis()->SetNdivisions(505);
  reso->GetXaxis()->SetRangeUser(-1.2, 1.2);
  reso->SetFillColor(13);
  reso->SetLineColor(1);
  reso->Draw("HIST");
  a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Wmass/resolution_"+pdfname+".pdf");
  delete a;

  TCanvas *b = new TCanvas("b", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.2);
  reso->SetTitle(" ");
  reso->GetXaxis()->SetTitle("#frac{m^{rec}_{W} - m^{gen}_{W}}{m^{gen}_{W}}");
  reso->GetYaxis()->SetTitle("events");
  reso->GetYaxis()->SetTitleSize(0.06);
  reso->GetXaxis()->SetTitleSize(0.05);
  reso->GetZaxis()->SetTitleSize(0.05);
  reso->GetXaxis()->SetTitleOffset(1.3);
  reso->GetYaxis()->SetTitleOffset(1.1);
  reso->GetZaxis()->SetTitleOffset(0.9);
  reso->GetXaxis()->SetNdivisions(505);
  reso->GetYaxis()->SetNdivisions(505);
  reso->GetXaxis()->SetRangeUser(-0.5, 0.5);
  reso->SetFillColor(13);
  reso->SetLineColor(1);
  reso_down->SetLineColor(kRed+1);
  reso_up->SetLineColor(kAzure+7);
  reso->Draw("HIST");
  reso_up->Draw("HIST SAME");
  reso_down->Draw("HIST SAME");
  b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Wmass/resolution_"+pdfname+"_variations.pdf");
  delete b;


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
