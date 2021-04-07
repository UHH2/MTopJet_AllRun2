#include "../include/CentralInclude.h"


using namespace std;
vector<vector<TH1F*>> get_hists(TFile* file, vector<int> binning);

int main(int argc, char* argv[]){

  TString directory = "/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/";
  TString prefix_mc = "uhh2.AnalysisModuleRunner.MC.";
  TString prefix_data = "uhh2.AnalysisModuleRunner.DATA.";

  vector<TFile*> files_mc, files_data;
  files_mc.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2016v3.root"));
  files_mc.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2017v2.root"));
  files_mc.push_back(new TFile(directory+"/elec/"+prefix_mc+"TTbar_2018.root"));
  files_mc.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2016v3.root"));
  files_mc.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2017v2.root"));
  files_mc.push_back(new TFile(directory+"/muon/"+prefix_mc+"TTbar_2018.root"));

  files_data.push_back(new TFile(directory+"/elec/"+prefix_data+"DATA_2016v3.root"));
  files_data.push_back(new TFile(directory+"/elec/"+prefix_data+"DATA_2017v2.root"));
  files_data.push_back(new TFile(directory+"/elec/"+prefix_data+"DATA_2018.root"));
  files_data.push_back(new TFile(directory+"/muon/"+prefix_data+"DATA_2016v3.root"));
  files_data.push_back(new TFile(directory+"/muon/"+prefix_data+"DATA_2017v2.root"));
  files_data.push_back(new TFile(directory+"/muon/"+prefix_data+"DATA_2018.root"));

  vector<int> NPV_bins = {0, 10, 15, 20, 30, 60};

  // Fill MC hists
  vector<TH1F*> h_mjet_mc, h_mW_mc;
  for(unsigned int i=0; i<files_mc.size(); i++){
    vector<vector<TH1F*>> hists = get_hists(files_mc[i], NPV_bins);
    if(i==0){
      for(unsigned int bin=0; bin<NPV_bins.size()-1; bin++){
        h_mjet_mc.push_back((TH1F*) hists[0][bin]->Clone());
        h_mW_mc.push_back((TH1F*) hists[1][bin]->Clone());
      }
    }
    else{
      for(unsigned int bin=0; bin<NPV_bins.size()-1; bin++){
        h_mjet_mc[bin]->Add(hists[0][bin]);
        h_mW_mc[bin]->Add(hists[1][bin]);
      }
    }
  }
  // Fill Data Hists
  vector<TH1F*> h_mjet_data, h_mW_data;
  for(unsigned int i=0; i<files_data.size(); i++){
    vector<vector<TH1F*>> hists = get_hists(files_data[i], NPV_bins);
    if(i==0){
      for(unsigned int bin=0; bin<NPV_bins.size()-1; bin++){
        h_mjet_data.push_back((TH1F*) hists[0][bin]->Clone());
        h_mW_data.push_back((TH1F*) hists[1][bin]->Clone());
      }
    }
    else{
      for(unsigned int bin=0; bin<NPV_bins.size()-1; bin++){
        h_mjet_data[bin]->Add(hists[0][bin]);
        h_mW_data[bin]->Add(hists[1][bin]);
      }
    }
  }

  // Get Means
  vector<double> NPV, NPV_e, NPV_data_e, mjet, mjet_e, mjet_data, mjet_data_e, mW, mW_e, mW_data, mW_data_e;
  for(unsigned int i=0; i<NPV_bins.size()-1; i++){
    TH1F* h = (TH1F*) h_mjet_mc[i]->Clone();
    TH1F* h_data = (TH1F*) h_mjet_data[i]->Clone();
    TH1F* hW = (TH1F*) h_mW_mc[i]->Clone();
    TH1F* hW_data = (TH1F*) h_mW_data[i]->Clone();
    double npv_min = NPV_bins[i];
    double npv_max = NPV_bins[i+1];
    double bincenter = npv_min + (npv_max - npv_min)/2;
    double binwidth = (npv_max - npv_min);
    NPV.push_back(bincenter);
    NPV_e.push_back(binwidth/2);
    NPV_data_e.push_back(0.0);

    mjet.push_back(h->GetMean());
    mjet_e.push_back(h->GetMeanError());
    mjet_data.push_back(h_data->GetMean());
    mjet_data_e.push_back(h_data->GetMeanError());

    mW.push_back(hW->GetMean());
    mW_e.push_back(hW->GetMeanError());
    mW_data.push_back(hW_data->GetMean());
    mW_data_e.push_back(hW_data->GetMeanError());
  }
  TGraphErrors * mjetPU = new TGraphErrors(NPV.size(), &NPV[0], &mjet[0], &NPV_e[0], &mjet_e[0]);
  TGraphErrors * mjetPU_data = new TGraphErrors(NPV.size(), &NPV[0], &mjet_data[0], &NPV_data_e[0], &mjet_data_e[0]);
  TGraphErrors * mWPU = new TGraphErrors(NPV.size(), &NPV[0], &mW[0], &NPV_e[0], &mW_e[0]);
  TGraphErrors * mWPU_data = new TGraphErrors(NPV.size(), &NPV[0], &mW_data[0], &NPV_data_e[0], &mW_data_e[0]);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  for(unsigned int i=0; i<h_mjet_mc.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_mjet_mc[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_mjet_mc[i]->GetMaximum());
    h_mjet_mc[i]->GetXaxis()->SetNdivisions(505);
    h_mjet_mc[i]->GetYaxis()->SetNdivisions(505);
    h_mjet_mc[i]->GetYaxis()->SetTitleOffset(1.6);
    h_mjet_mc[i]->GetXaxis()->SetTitleOffset(1.3);
    h_mjet_mc[i]->SetTitle(" ");
    h_mjet_mc[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass_Run2/PUbin_";
    filename += i;
    filename += "_MC.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_mjet_data.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_mjet_data[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_mjet_data[i]->GetMaximum());
    h_mjet_data[i]->GetXaxis()->SetNdivisions(505);
    h_mjet_data[i]->GetYaxis()->SetNdivisions(505);
    h_mjet_data[i]->GetYaxis()->SetTitleOffset(1.6);
    h_mjet_data[i]->GetXaxis()->SetTitleOffset(1.3);
    h_mjet_data[i]->SetTitle(" ");
    h_mjet_data[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass_Run2/PUbin_";
    filename += i;
    filename += "_DATA.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_mW_mc.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_mW_mc[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_mW_mc[i]->GetMaximum());
    h_mW_mc[i]->GetXaxis()->SetNdivisions(505);
    h_mW_mc[i]->GetYaxis()->SetNdivisions(505);
    h_mW_mc[i]->GetYaxis()->SetTitleOffset(1.6);
    h_mW_mc[i]->GetXaxis()->SetTitleOffset(1.3);
    h_mW_mc[i]->SetTitle(" ");
    h_mW_mc[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass_Run2/PUbin_";
    filename += i;
    filename += "_MC.pdf";
    a->SaveAs(filename);
    delete a;
  }

  for(unsigned int i=0; i<h_mW_data.size(); i++){
    TCanvas *a = new TCanvas("a", " ", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.1);
    h_mW_data[i]->GetYaxis()->SetRangeUser(0, 1.1 * h_mW_data[i]->GetMaximum());
    h_mW_data[i]->GetXaxis()->SetNdivisions(505);
    h_mW_data[i]->GetYaxis()->SetNdivisions(505);
    h_mW_data[i]->GetYaxis()->SetTitleOffset(1.6);
    h_mW_data[i]->GetXaxis()->SetTitleOffset(1.3);
    h_mW_data[i]->SetTitle(" ");
    h_mW_data[i]->Draw("HIST");
    TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass_Run2/PUbin_";
    filename += i;
    filename += "_DATA.pdf";
    a->SaveAs(filename);
    delete a;
  }

  TCanvas *a = new TCanvas("a", " ", 600, 600);
  TPad *p1 = new TPad("p1","p1",0.0,0.5,1.0,1.0);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.1);
  p1->SetBottomMargin(0.0);
  p1->SetBorderMode(0);
  p1->Draw();
  TPad *p2 = new TPad("p2","p2",0.0,0.0,1.0,0.501);
  p2->SetTopMargin(0.);
  p2->SetLeftMargin(0.15);
  p1->SetRightMargin(0.1);
  p2->SetBottomMargin(0.18);
  p2->SetBorderMode(0);
  p2->Draw();

  p1->cd();
  // mjetPU_data->GetXaxis()->SetLabelSize(0);
  // mjetPU_data->GetXaxis()->SetTickLength(0);
  mjetPU_data->Draw("AP");
  mjetPU_data->GetXaxis()->SetRangeUser(0, 50);
  mjetPU_data->GetYaxis()->SetRangeUser(159, 196);
  mjetPU_data->GetXaxis()->SetTitle(" ");
  mjetPU_data->GetYaxis()->SetTitle("<#it{m}_{jet}> [GeV]");
  mjetPU_data->SetTitle(" ");
  mjetPU_data->GetXaxis()->SetNdivisions(505);
  mjetPU_data->GetYaxis()->SetNdivisions(505);
  mjetPU_data->GetYaxis()->SetTitleOffset(0.7);
  // mjetPU_data->GetXaxis()->SetTitleOffset(0.9);
  mjetPU_data->GetYaxis()->SetLabelSize(0.08);
  mjetPU_data->GetYaxis()->SetTitleSize(0.1);
  // mjetPU_data->GetXaxis()->SetTitleSize(0.05);
  mjetPU_data->SetLineColor(kBlack);
  mjetPU_data->SetMarkerColor(kBlack);
  mjetPU_data->SetMarkerStyle(8);
  mjetPU_data->SetMarkerSize(1);
  mjetPU->SetLineColor(kRed);
  mjetPU->SetFillColor(kRed);
  mjetPU->SetMarkerColor(kRed);
  mjetPU->SetMarkerStyle(8);
  mjetPU->SetMarkerSize(0);
  mjetPU->Draw("E2 SAME");
  mjetPU_data->Draw("P SAME");
  p1->RedrawAxis();

  p2->cd();
  mWPU_data->Draw("AP");
  // mWPU_data->GetHistogram()->SetMaximum(85);
  mWPU_data->GetXaxis()->SetRangeUser(0, 50);
  mWPU_data->GetYaxis()->SetRangeUser(60, 97);
  mWPU_data->GetXaxis()->SetTitle("Number of primary vertices");
  mWPU_data->GetYaxis()->SetTitle(" ");
  mWPU_data->SetTitle(" ");
  mWPU_data->GetXaxis()->SetNdivisions(505);
  mWPU_data->GetYaxis()->SetNdivisions(505);
  // mWPU_data->GetYaxis()->SetTitleOffset(1.1);
  mWPU_data->GetXaxis()->SetTitleOffset(0.8);
  // mWPU_data->GetYaxis()->SetTitleSize(0.05);
  mWPU_data->GetXaxis()->SetTitleSize(0.1);
  mWPU_data->GetXaxis()->SetLabelSize(0.08);
  mWPU_data->GetYaxis()->SetLabelSize(0.08);
  mWPU_data->SetLineColor(14);
  mWPU_data->SetMarkerColor(14);
  mWPU_data->SetMarkerStyle(8);
  mWPU_data->SetMarkerSize(1);
  mWPU->SetLineColor(kOrange+7);
  mWPU->SetFillColor(kOrange+7);
  mWPU->SetMarkerColor(kOrange+7);
  mWPU->SetMarkerStyle(8);
  mWPU->SetMarkerSize(0);
  mWPU->Draw("E2 SAME");
  mWPU_data->Draw("P SAME");
  p2->RedrawAxis();

  a->cd();
  CMSLabel(true, 0.2, 0.9);

  TPad *b = new TPad("b","b",0.0,0.46,1.0,0.54);
  b->SetBorderMode(0);
  b->Draw();
  b->cd();

  gStyle->SetLineStyleString(17,"4 10");
  TLine *line = new TLine(0.15,0.1,0.15,0.9); //vertical
  line->SetLineWidth(2);
  line->SetLineStyle(17);
  line->Draw();
  TLine *line2 = new TLine(0.9,0.1,0.9,0.9); //vertical
  line2->SetLineWidth(2);
  line2->SetLineStyle(17);
  line2->Draw();
  // line = new TLine(0.15,0.6,0.15,1); // vertical
  // line->Draw();
  // line = new TLine(0.13,0.1,0.17,0.7); //horizontal
  // line->Draw();
  // line = new TLine(0.13,0.3,0.17,0.9); // horizontal
  // line->Draw();

  a->cd();
  TLegend* leg = new TLegend(0.50, 0.35, 0.85, 0.65);
  leg->AddEntry(mjetPU_data,"t decay data","pe");
  leg->AddEntry(mjetPU,"t decay t#bar{t}","pf");
  leg->AddEntry(mWPU_data,"W decay data","pe");
  leg->AddEntry(mWPU,"W decay t#bar{t}","pf");
  leg->Draw();

  TString infotext = "137.1 fb^{-1} (13 TeV)";
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(0.9);
  text1->SetY(0.995);
  text1->SetTextSize(0.04);
  text1->Draw();

  TString filename = "/afs/desy.de/user/s/schwarzd/Plots/PileupVsMass_Run2/MjetVsPU.pdf";
  a->SaveAs(filename);
  delete a;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}





vector<vector<TH1F*>> get_hists(TFile* file, vector<int> binning){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  vector<TH1F*> h_mjet, h_mW;

  int Nbins = binning.size()-1;
  TH1F* dummyhist = new TH1F("hh", "hh", 60, 0, 300);

  for(unsigned int i=1; i<=Nbins; i++){
    TString name_mjet = "mjet"+to_string(i);
    TString name_mW = "mW"+to_string(i);
    h_mjet.push_back((TH1F*) dummyhist->Clone(name_mjet));
    h_mW.push_back((TH1F*) dummyhist->Clone(name_mW));
  }


  Double_t mjet, mW;
  Int_t npv;
  Bool_t passed_measurement_rec;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("passed_measurement_rec", &passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("Mass_Rec", &mjet);
  tree->SetBranchAddress("mW", &mW);
  tree->SetBranchAddress("NPV", &npv);
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    weight = rec_weight * gen_weight;
    if(passed_measurement_rec){
      for(int bin=0; bin<binning.size()-1; bin++){
        if(npv > binning[bin] && npv < binning[bin+1]){
          if(mjet > 120 && mjet < 240) h_mjet[bin]->Fill(mjet, weight);
          if(mW > 65 && mW < 95) h_mW[bin]->Fill(mW, weight);
        }
      }
    }
  }
  vector<vector<TH1F*>> hists;
  hists.push_back(h_mjet);
  hists.push_back(h_mW);
  return hists;
}
