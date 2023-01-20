#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(){

  TString save_path = get_save_path();
  vector<TString> years = {"2016v3", "2017v2", "2018"};
  vector<TString> variations = {"HOME","JEC_up", "JEC_down", "COR_up", "COR_down", "JMS_upup", "JMS_downdown", "hdampup", "hdampdown", "ISRup_2", "ISRdown_2", "TuneUp", "TuneDown"};
  vector<TString> LegendEntries = {"nominal", "JEC up", "JEC down", "XCone up", "XCone down", "JMS upup", "JMS downdown", "hdamp up", "hdamp down", "isr up", "isr down", "Tune Up", "Tune Down"};
  vector<Color_t> colors = {kBlack, kAzure+7, kAzure+7, kOrange+1, kOrange+1, kRed-4, kRed-4, kGreen-3, kGreen-3, kMagenta-3, kMagenta-3, kBlue-6, kBlue-6};
  vector<int> styles = {1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};

  TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau32";
  int rebin = 5;

  for(auto year: years){
    vector<TH1F*> hists;
    for(auto variation: variations){
      TFile *f_tt_mu, *f_tt_el;
      if(variation == "HOME"){
        f_tt_mu = new TFile(dir+"/muon/"+prefix+"TTbar_"+year+".root");
        f_tt_el = new TFile(dir+"/elec/"+prefix+"TTbar_"+year+".root");
      }
      else if(!variation.Contains("_")){
        f_tt_mu = new TFile(dir+"/muon/"+prefix+"TTbar_"+variation+"_"+year+".root");
        f_tt_el = new TFile(dir+"/elec/"+prefix+"TTbar_"+variation+"_"+year+".root");
      }
      else if(variation.Contains("ISR")){
        if(year=="2016v3"){
          TString var = (variation.Contains("up"))?"isrup":"isrdown";
          f_tt_mu = new TFile(dir+"/muon/"+prefix+"TTbar_"+var+"_"+year+".root");
          f_tt_el = new TFile(dir+"/elec/"+prefix+"TTbar_"+var+"_"+year+".root");
        }
        else{
          f_tt_mu = new TFile(dir+"/muon/"+variation+"/"+prefix+"TTbar_"+year+".root");
          f_tt_el = new TFile(dir+"/elec/"+variation+"/"+prefix+"TTbar_"+year+".root");
        }
      }
      else{
        f_tt_mu = new TFile(dir+"/muon/"+variation+"/"+prefix+"TTbar_"+year+".root");
        f_tt_el = new TFile(dir+"/elec/"+variation+"/"+prefix+"TTbar_"+year+".root");
      }
      TH1F* h_mu = (TH1F*) f_tt_mu->Get(histname);
      TH1F* h_el = (TH1F*) f_tt_el->Get(histname);
      TH1F* h = (TH1F*) h_mu->Clone();
      h->Add(h_el);
      double integral = h->Integral();
      h->Scale(1/integral);
      h->Rebin(rebin);
      hists.push_back(h);
    }
    SetupGlobalStyle();
    TCanvas *c = new TCanvas("c", "c", 600, 600);
    gPad->SetLeftMargin(0.19);
    TLegend *leg = new TLegend(0.22, 0.6, 0.47, 0.85);
    leg->SetBorderSize(0);
    for(unsigned int i=0; i<hists.size(); i++){
      hists[i]->SetTitle("");
      hists[i]->GetXaxis()->SetTitle("#tau_{32}");
      hists[i]->GetYaxis()->SetTitle("a.u.");
      hists[i]->GetXaxis()->SetTitleSize();
      hists[i]->GetYaxis()->SetTitleSize();
      hists[i]->GetXaxis()->SetNdivisions(505);
      hists[i]->GetYaxis()->SetNdivisions(505);
      hists[i]->GetYaxis()->SetTitleOffset(1.8);
      hists[i]->SetLineWidth(3);
      hists[i]->SetLineColor(colors[i]);
      hists[i]->SetLineStyle(styles[i]);
      leg->AddEntry(hists[i], LegendEntries[i], "l");
    }
    hists[0]->GetYaxis()->SetRangeUser(0, hists[0]->GetMaximum()*1.2);
    hists[0]->Draw("HIST");
    for(auto h: hists) h->Draw("HIST SAME");
    leg->Draw();
    gPad->RedrawAxis();
    c->SaveAs(save_path+"/Tau32_JEC/Tau32_"+year+".pdf");
  }


  return 0;

}
