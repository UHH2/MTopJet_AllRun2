#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[])
{

  TFile *f_central = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *f_nonclosure = new TFile(dir+"NONCLOSURE/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  vector<TString>  histdir = {"XCone_jec/", "XCone_jec_subjets/"};
  vector<TString>  histname = {"M_jet1_B", "pt_had_subjets_fine"};
  vector<TH1F*> h_central, h_nonclosure;
  for(unsigned int i=0; i<histname.size(); i++){
    h_central.push_back( (TH1F*)f_central->Get(histdir[i] + histname[i]) );
    h_nonclosure.push_back( (TH1F*)f_nonclosure->Get(histdir[i] + histname[i]) );
  }

  for(unsigned int i=0; i<h_central.size(); i++){

    double ylow, yhigh;
    ylow = 0;
    yhigh = h_central[i]->GetMaximum();
    if(h_nonclosure[i]->GetMaximum() > yhigh) yhigh = h_nonclosure[i]->GetMaximum();
    yhigh *= 1.1;

    h_central[i]->SetTitle(" ");
    if(i==0)       h_central[i]->GetXaxis()->SetTitle("Leading-jet m_{jet} [GeV]");
    else if (i==1) h_central[i]->GetXaxis()->SetTitle("subjet p_{T} [GeV]");
    h_central[i]->GetYaxis()->SetTitle("events");
    h_central[i]->GetYaxis()->SetRangeUser(ylow, yhigh);
    h_central[i]->GetXaxis()->SetNdivisions(505);
    h_central[i]->GetYaxis()->SetNdivisions(505);
    h_central[i]->GetYaxis()->SetTitleOffset(1.5);
    h_central[i]->GetXaxis()->SetTitleOffset(1.3);

    h_central[i]->SetLineWidth(4);
    h_central[i]->SetLineColor(kRed-4);
    h_nonclosure[i]->SetLineWidth(4);
    h_nonclosure[i]->SetLineColor(kAzure+7);

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);

    TCanvas *A = new TCanvas("A", "A", 600, 600);
    gPad->SetLeftMargin(0.1);
    TGaxis::SetMaxDigits(3);
    h_central[i]->Draw("HIST");
    h_nonclosure[i]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.55,0.85,0.80,0.65);
    leg->SetBorderSize(0);
    leg->AddEntry(h_central[i],"central","l");
    leg->AddEntry(h_nonclosure[i],"non-closure variation","l");
    leg->Draw();
    A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/NonClosure/"+histname[i]+".pdf");
    delete A;
  }
}
