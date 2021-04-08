#include "../include/CentralInclude.h"


using namespace std;
void PlotHist(TH1F* hist, TF1* fit, TString histname);
vector<vector<TH1F*>> get_hists(TFile* file, vector<int> bins_PV, vector<double> bins_pt);

int main(int argc, char* argv[]){

  TString directory = "/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/";
  TString prefix_mc = "uhh2.AnalysisModuleRunner.MC.";
  TString histdir = "/RecGenHists_GenSelRecInfo";

  vector<TFile*> files;
  files.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2016v3.root"));
  files.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2017v2.root"));
  files.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2018.root"));
  files.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2016v3.root"));
  files.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2017v2.root"));
  files.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2018.root"));

  vector<double> ptbins = {400, 450, 500, 600, 700, 800, 1500};
  vector<int> pvbins = {0, 10, 20, 100};

  vector<vector<TH1F*>> resohists = get_hists(files[0], pvbins, ptbins);
  for(unsigned int i=1; i<files.size(); i++){
    vector<vector<TH1F*>> hists = get_hists(files[i], pvbins, ptbins);
    for(unsigned int pvbin=0; pvbin<pvbins.size(); pvbin++){ // one more NPV bin with no cut on NPV
      for(unsigned int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
        resohists[pvbin][ptbin]->Add(hists[pvbin][ptbin]);
      }
    }
  }

  vector<TGraphErrors*> h_resos;
  for(unsigned int pvbin=0; pvbin<pvbins.size(); pvbin++){  // one more NPV bin with no cut on NPV
    vector<double> pt, pt_e, reso, reso_e;
    for(unsigned int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      double binwidth = ptbins[ptbin+1]-ptbins[ptbin];
      pt.push_back(ptbins[ptbin]+binwidth/2);
      pt_e.push_back(binwidth/2);
      TH1F* hist = (TH1F*) resohists[pvbin][ptbin]->Clone();
      TF1* f1 = new TF1("f1","gaus",-0.2,0.2);
      hist->Fit("f1","R");
      TF1* fit = hist->GetFunction("f1");
      reso.push_back(fit->GetParameter(2));
      reso_e.push_back(fit->GetParError(2));
    }
    int nbins = reso.size();
    TGraphErrors* h_reso = new TGraphErrors(nbins, &pt[0], &reso[0], &pt_e[0], &reso_e[0]);
    h_resos.push_back(h_reso);
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  h_resos[0]->Draw("AP");
  h_resos[0]->GetXaxis()->SetRangeUser(405, 1495);
  h_resos[0]->GetYaxis()->SetRangeUser(0.0, 0.15);
  h_resos[0]->GetXaxis()->SetTitle("#it{p}_{T, jet}^{gen} [GeV]");
  // h_resos[0]->GetYaxis()->SetTitle("RMS #left[ #frac{m_{jet}^{rec} - m_{jet}^{gen}}{m_{jet}^{gen}} #right]");
  h_resos[0]->GetYaxis()->SetTitle("Jet Mass Resolution");
  h_resos[0]->SetTitle(" ");
  h_resos[0]->GetXaxis()->SetNdivisions(505);
  h_resos[0]->GetYaxis()->SetNdivisions(505);
  h_resos[0]->GetYaxis()->SetTitleOffset(1.6);
  h_resos[0]->GetXaxis()->SetTitleOffset(1.1);
  h_resos[0]->GetYaxis()->SetTitleSize(0.05);
  h_resos[0]->GetXaxis()->SetTitleSize(0.05);
  h_resos[0]->GetXaxis()->SetLabelSize(0.05);
  h_resos[0]->GetYaxis()->SetLabelSize(0.05);
  Color_t col[] = {kAzure+10, 798, kRed, kBlack};
  for(int i=0; i<h_resos.size(); i++){
    h_resos[i]->SetLineColor(col[i]);
    h_resos[i]->SetMarkerColor(col[i]);
    h_resos[i]->SetMarkerStyle(8);
    h_resos[i]->SetMarkerSize(1);
  }
  for(auto h: h_resos) h->Draw("P SAME");

  CMSSimLabel(true, 0.2, 0.85);

  TLegend *leg = new TLegend(0.55, 0.55, 0.85, 0.85);
  leg->AddEntry(h_resos[0], " 0 < NPV < 10", "ple");
  leg->AddEntry(h_resos[1], "10 < NPV < 20", "ple");
  leg->AddEntry(h_resos[2], "NPV > 20", "ple");
  leg->AddEntry(h_resos[3], "all", "ple");
  leg->Draw();
  TString filename = "/afs/desy.de/user/s/schwarzd/Plots/MassResolution_Run2/reso.pdf";
  a->SaveAs(filename);
  delete a;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}


void PlotHist(TH1F* hist, TF1* fit, TString histname){
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle("#left[ #frac{m_{jet}^{rec} - m_{jet}^{gen}}{m_{jet}^{gen}} #right]");
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetXaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetRangeUser(-0.5, 0.5);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
  hist->SetLineColor(kBlack);
  hist->SetFillColor(15);
  fit->SetLineWidth(3);
  fit->SetLineColor(kRed);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.17);
  hist->Draw("HIST");
  fit->Draw("SAME");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassResolution_Run2/"+histname+".pdf");
  delete A;
  return;

}





vector<vector<TH1F*>> get_hists(TFile* file, vector<int> bins_PV, vector<double> bins_pt){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  vector<vector<TH1F*>> reso;
  int Nbins_PV = bins_PV.size()-1;
  int Nbins_pt = bins_pt.size()-1;

  TH1F* reso_dummy = new TH1F("reso", "reso", 60, -0.5, 0.5);

  for(unsigned int i=1; i<=Nbins_PV+1; i++){ // one more bin for no cut on NPV
    vector<TH1F*> dummyvec;
    for(unsigned int j=1; j<=Nbins_pt; j++){
      TString name = "bin"+to_string(i)+"-"+to_string(j);
      dummyvec.push_back( (TH1F*) reso_dummy->Clone(name));
    }
    reso.push_back(dummyvec);
  }


  Double_t mjet_rec, mjet_gen, pt_gen;
  Int_t npv;
  Bool_t passed_measurement_rec, passed_measurement_gen;
  Bool_t pass_pt350migration_rec, pass_massmigration_rec, pass_btagmigration_rec, pass_subptmigration_rec, pass_leptonptmigration_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("passed_measurement_rec", &passed_measurement_rec);
  tree->SetBranchAddress("passed_ptmigration_rec", &pass_pt350migration_rec);
  tree->SetBranchAddress("passed_massmigration_rec", &pass_massmigration_rec);
  tree->SetBranchAddress("passed_btagmigration_rec", &pass_btagmigration_rec);
  tree->SetBranchAddress("passed_subptmigration_rec", &pass_subptmigration_rec);
  tree->SetBranchAddress("passed_leptonptmigration_rec", &pass_leptonptmigration_rec);
  tree->SetBranchAddress("passed_measurement_gen", &passed_measurement_gen);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("Mass_Rec", &mjet_rec);
  tree->SetBranchAddress("Mass_Gen33", &mjet_gen);
  tree->SetBranchAddress("Pt_Gen33", &pt_gen);
  tree->SetBranchAddress("NPV", &npv);
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    weight = rec_weight * gen_weight;

    // has to pass gen sel and at least one reco sideband
    bool passedsel = passed_measurement_gen &&
    (passed_measurement_rec ||
      pass_pt350migration_rec ||
      pass_massmigration_rec ||
      pass_btagmigration_rec ||
      pass_subptmigration_rec ||
      pass_leptonptmigration_rec);



    if(passedsel){
      for(int pvbin=0; pvbin<bins_PV.size()-1; pvbin++){
        for(int ptbin=0; ptbin<bins_pt.size()-1; ptbin++){
          if(npv > bins_PV[pvbin] && npv < bins_PV[pvbin+1]){
            if(pt_gen > bins_pt[ptbin] && pt_gen < bins_pt[ptbin+1]){
              double r = (mjet_rec - mjet_gen)/mjet_gen;
              reso[pvbin][ptbin]->Fill(r,weight);
            }
          }
        }
      }
      for(int ptbin=0; ptbin<bins_pt.size()-1; ptbin++){
        if(pt_gen > bins_pt[ptbin] && pt_gen < bins_pt[ptbin+1]){
          double r = (mjet_rec - mjet_gen)/mjet_gen;
          reso[reso.size()-1][ptbin]->Fill(r,weight);
        }
      }
    }
  }
  return reso;
}
