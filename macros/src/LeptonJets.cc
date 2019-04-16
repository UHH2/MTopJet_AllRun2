#include "../include/CentralInclude.h"


using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);


int main(int argc, char* argv[]){
  TFile* file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_LeptonJets.root");

  vector<TString> names = {"pt_lepton_ptCut", "JetMass_ptCut", "pt_lepton_ptCut_massCut", "JetMass_ptCut_massCut"};
  vector<TString> xtitles = {"p_{T}", "m_{jet}", "p_{T}", "m_{jet}"};
  vector<TH1F*> h_mu, h_mujet;
  h_mu.push_back((TH1F*)file->Get("01a_muon/pt_lepton") );
  h_mu.push_back((TH1F*)file->Get("01a_muon_xcone/Mass_HadJet33_rebin"));
  h_mu.push_back((TH1F*)file->Get("02a_muon/pt_lepton") );
  h_mu.push_back((TH1F*)file->Get("02a_muon_xcone/Mass_HadJet33_rebin"));
  h_mujet.push_back((TH1F*)file->Get("01b_muon/pt_matchjet") );
  h_mujet.push_back((TH1F*)file->Get("01b_muon_xcone/Mass_HadJet33_rebin"));
  h_mujet.push_back((TH1F*)file->Get("02b_muon/pt_matchjet") );
  h_mujet.push_back((TH1F*)file->Get("02b_muon_xcone/Mass_HadJet33_rebin"));

  for(auto h: h_mu){
    h->SetLineWidth(4);
    h->SetLineColor(kAzure+7);
  }
  for(auto h: h_mujet){
    h->SetLineWidth(4);
    h->SetLineColor(kRed+1);
  }

  vector<TH1F*> ratio_mu;
  for(unsigned int i=0; i<h_mu.size(); i++){
    ratio_mu.push_back(GetRatio(h_mu[i], h_mujet[i]));
  }

  cout << "-----------------------------------------------------" << endl;
  cout << " Overall normalization difference: " << endl << endl;
  cout << " " << 100*fabs(h_mu[1]->Integral() - h_mujet[1]->Integral()) / h_mujet[1]->Integral() << "% after lepton pT cut" << endl;
  cout << " " << 100*fabs(h_mu[3]->Integral() - h_mujet[3]->Integral()) / h_mujet[3]->Integral() << "% after lepton pT cut and mass cut" << endl;
  cout << "-----------------------------------------------------" << endl;

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<h_mu.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_mu[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_mu[i]->GetMaximum());
    h_mu[i]->GetXaxis()->SetNdivisions(505);
    h_mu[i]->GetYaxis()->SetNdivisions(505);
    h_mu[i]->GetYaxis()->SetTitleOffset(1.6);
    h_mu[i]->GetXaxis()->SetTitleOffset(1.3);
    h_mu[i]->SetTitle(" ");
    h_mu[i]->GetXaxis()->SetTitle(xtitles[i]);
    h_mu[i]->GetYaxis()->SetTitle("events");
    h_mu[i]->Draw("HIST");
    h_mujet[i]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.53,0.63,0.88,0.88);
    leg->AddEntry(h_mu[i], "lepton from W decay", "l");
    leg->AddEntry(h_mujet[i], "lepton GenJet", "l");
    leg->Draw();
    a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/LeptonJets/"+names[i]+".pdf");
    delete a;
  }

  for(unsigned int i=0; i<ratio_mu.size(); i++){
    TCanvas *b = new TCanvas("b", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    ratio_mu[i]->GetYaxis()->SetTitle("#frac{lepton from W}{leton GenJet}");
    ratio_mu[i]->GetXaxis()->SetTitle(xtitles[i]);
    ratio_mu[i]->GetXaxis()->SetNdivisions(505);
    ratio_mu[i]->GetYaxis()->SetNdivisions(505);
    ratio_mu[i]->GetYaxis()->SetTitleOffset(1.5);
    ratio_mu[i]->GetXaxis()->SetTitleOffset(1.3);
    ratio_mu[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
    ratio_mu[i]->SetLineColor(1);
    ratio_mu[i]->Draw("E1");
    b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/LeptonJets/Ratio_"+names[i]+".pdf");
    delete b;
  }

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}

TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double r = N1/N2;
    double error = sqrt(fabs((1-r)*r/N2));
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
      ratio->SetBinError(i, 0);
    }
    else{
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
