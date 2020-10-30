#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){

  TString year, year_v;

  bool isCut = false;
  if(argc == 3 && strcmp(argv[2], "cut") == 0) isCut = true;
  if(argc > 1 && strcmp(argv[1], "2017") == 0){
    year = "2017";
    year_v = "2017v2";
  }
  if(argc > 1 && strcmp(argv[1], "2018") == 0){
    year = "2018";
    year_v = "2018";
  }

  vector<TString> variations = {"up_4", "down_4"};

  TString semilep                 = "/muon/test/uhh2.AnalysisModuleRunner.MC.TTbar_Mtt0000to0700_SemiLep_";
  vector<TString> hist_collection = {"Gen_Weights", "Gen_Weights_pass_gen", "Gen_Weights_Massbin_145", "Gen_Weights_Massbin_275"};
  vector<TString> hists           = {"gen_weight_big_1", "gen_weight_medium_1", "gen_weight_small_1"};
  vector<TString> hists_negetiv   = {"gen_weight_negativ_big_1", "gen_weight_negativ_medium_1", "gen_weight_negativ_small_1"};

  TFile* f;
  vector<TH1F*> hists_all, hists_pass_gen, hists_mass145, hists_mass275;
  vector<vector<TH1F*>> every_hist;
  TH1F* hist;

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  TLegend *leg;

  for(unsigned int i=0; i<variations.size(); i++){ // Variations
    f = new TFile(dir+year+semilep+year_v+"_"+variations[i]+".root");

    // for(unsigned int j=0; j<hist_collection.size(); j++){ // Hists Collection
    for(unsigned int k=0; k<hists.size(); k++){ // Single Hists
      hists_all.push_back((TH1F*)f->Get(hist_collection[0]+"/"+hists[k]));
      hists_pass_gen.push_back((TH1F*)f->Get(hist_collection[1]+"/"+hists[k]));
      hists_mass145.push_back((TH1F*)f->Get(hist_collection[2]+"/"+hists[k]));
      hists_mass275.push_back((TH1F*)f->Get(hist_collection[3]+"/"+hists[k]));
    }
    // }
    every_hist.push_back(hists_all);
    every_hist.push_back(hists_pass_gen);
    every_hist.push_back(hists_mass145);
    every_hist.push_back(hists_mass275);

    for(unsigned int j=0; j<every_hist.size(); j++){ // Go into vector with all Hists
      for(unsigned int k=0; k<hists.size(); k++){ // Single Hists
        if(isCut){
          every_hist[j].at(k)->SetTitle("");
          every_hist[j].at(k)->GetYaxis()->SetRangeUser(0, 10);
          every_hist[j].at(k)->GetXaxis()->SetNdivisions(505);
          every_hist[j].at(k)->GetYaxis()->SetNdivisions(505);
          every_hist[j].at(k)->GetXaxis()->SetTitleSize(0.05);
          every_hist[j].at(k)->GetYaxis()->SetTitleSize(0.04);
          every_hist[j].at(k)->GetXaxis()->SetTitleOffset(0.9);
          every_hist[j].at(k)->GetYaxis()->SetTitleOffset(1.5);
          every_hist[j].at(k)->GetXaxis()->SetTitle("weight");
          every_hist[j].at(k)->GetYaxis()->SetTitle("");
          every_hist[j].at(k)->SetLineWidth(2);
          every_hist[j].at(k)->SetLineColor(kRed);
          TCanvas *CUT = new TCanvas("Cut", "Cut", 600, 600);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.12);
          every_hist[j].at(k)->Draw("HIST");
          gPad->RedrawAxis();
          CUT->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Weights/"+year+"/cut/"+hist_collection[j]+"_"+hists[k]+"_FSR"+variations[i]+"_cut_on_y10.pdf");
          delete CUT;
        }
        else{
          every_hist[j].at(k)->SetTitle("");
          every_hist[j].at(k)->GetYaxis()->SetRangeUser(0, every_hist[j].at(k)->GetMaximum()*1.2);
          every_hist[j].at(k)->GetXaxis()->SetNdivisions(505);
          every_hist[j].at(k)->GetYaxis()->SetNdivisions(505);
          every_hist[j].at(k)->GetXaxis()->SetTitleSize(0.05);
          every_hist[j].at(k)->GetYaxis()->SetTitleSize(0.04);
          every_hist[j].at(k)->GetXaxis()->SetTitleOffset(0.9);
          every_hist[j].at(k)->GetYaxis()->SetTitleOffset(1.5);
          every_hist[j].at(k)->GetXaxis()->SetTitle("weight");
          every_hist[j].at(k)->GetYaxis()->SetTitle("");
          every_hist[j].at(k)->SetLineWidth(2);
          every_hist[j].at(k)->SetLineColor(kRed);
          TCanvas *A = new TCanvas("A", "A", 600, 600);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.12);
          every_hist[j].at(k)->Draw("HIST");
          gPad->RedrawAxis();
          A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/Weights/"+year+"/"+hist_collection[j]+"_"+hists[k]+"_FSR"+variations[i]+".pdf");
          delete A;
        }
      }
    }
    every_hist = {};
    hists_all = {};
    hists_pass_gen = {};
    hists_mass145 = {};
    hists_mass275 = {};
  }


}
