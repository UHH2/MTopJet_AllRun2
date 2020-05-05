#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc != 2){
    cout << "Usage: ./Radiation_SYS <year>" << endl;
    return 0;
  }

  bool debug = false;

  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  TString year = argv[1];
  TString tau32 = "AK8_tau_pass_rec_masscut_140/tau32_hadjet";
  TString save_path = "/afs/desy.de/user/p/paaschal/Plots/Radiation_SYS/";
  std::vector<int> color = {kOrange, kBlue, kGreen};

  // #############################################################################################
  // Get Background ##############################################################################
  if(debug) cout << "Background" << endl;

  TString other_path = "uhh2.AnalysisModuleRunner.MC.other.root";
  TString ST_path = "uhh2.AnalysisModuleRunner.MC.ST.root";
  TString WJets_path = "uhh2.AnalysisModuleRunner.MC.WJets.root";

  TFile *other_file = new TFile(dir+year+"/muon/"+other_path);
  TFile *ST_file = new TFile(dir+year+"/muon/"+ST_path);
  TFile *WJets_file = new TFile(dir+year+"/muon/"+WJets_path);

  TH1F *Background = (TH1F*)other_file->Get(tau32); // Add other hists to this one.
  TH1F *ST_hist = (TH1F*)ST_file->Get(tau32);
  TH1F *WJets_hist = (TH1F*)WJets_file->Get(tau32);

  Background->Add(ST_hist, 1);
  Background->Add(WJets_hist, 1);

  // #############################################################################################
  // Get Data ####################################################################################
  if(debug) cout << "Data" << endl;

  TString data_path = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file = new TFile(dir+year+"/muon/"+data_path);
  TH1F *data_hist = (TH1F*)data_file->Get(tau32);

  // #############################################################################################
  // Get TTbar and FSR ###########################################################################
  if(debug) cout << "TTbar" << endl;

  vector<TString> FSR = {"FSR_sqrt2", "FSR_2", "FSR_4"};
  vector<TString> FSRup = {"FSRup_sqrt2", "FSRup_2", "FSRup_4"};
  vector<TString> FSRdown = {"FSRdown_sqrt2", "FSRdown_2", "FSRdown_4"};
  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";

  // TTbar
  TFile *ttbar_file = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(tau32);
  ttbar_hist->Add(Background, 1);

  // FSR
  TH1F* sum_fsr_bkg;
  vector<TFile*> FSRup_file, FSRdown_file;
  vector<TH1F*> FSRup_hist, FSRdown_hist, sum_up_hist, sum_down_hist;

  // Fill FSR vectors
  for(unsigned int i=0; i<FSRup.size();i++){
    FSRup_file.push_back(new TFile(dir+year+"/muon/"+FSRup[i]+"/"+ttbar_path));
    FSRdown_file.push_back(new TFile(dir+year+"/muon/"+FSRdown[i]+"/"+ttbar_path));

    FSRup_hist.push_back((TH1F*)FSRup_file[i]->Get(tau32));
    FSRdown_hist.push_back((TH1F*)FSRdown_file[i]->Get(tau32));
  }

  // #############################################################################################
  // Add hists from FSR and background ###########################################################
  if(debug) cout << "FSR+Background" << endl;

  for(unsigned int i=0; i<FSRup_hist.size();i++){
    sum_fsr_bkg = (TH1F*)FSRup_hist[i]->Clone();
    sum_fsr_bkg->Add(Background, 1);
    sum_up_hist.push_back(sum_fsr_bkg);

    sum_fsr_bkg = (TH1F*)FSRdown_hist[i]->Clone();
    sum_fsr_bkg->Add(Background, 1);
    sum_down_hist.push_back(sum_fsr_bkg);
  }

  // #############################################################################################
  // Rebin 5 #######################################################################################
  if(debug) cout << "Rebin 5" << endl;

  TH1F *ttbar_rebin_5bins = rebin(ttbar_hist, 10);
  TH1F *data_rebin_5bins = rebin(data_hist, 10);

  vector<TH1F*> sum_up_rebin_5bins = rebin(sum_up_hist, 10);
  vector<TH1F*> sum_down_rebin_5bins = rebin(sum_down_hist, 10);

  // #############################################################################################
  // Rebin 10 #######################################################################################
  if(debug) cout << "Rebin 10" << endl;

  TH1F *ttbar_rebin_10bins = rebin(ttbar_hist, 5);
  TH1F *data_rebin_10bins = rebin(data_hist, 5);

  vector<TH1F*> sum_up_rebin_10bins = rebin(sum_up_hist, 5);
  vector<TH1F*> sum_down_rebin_10bins = rebin(sum_down_hist, 5);

  // #############################################################################################
  // Vector for compact Plots ####################################################################
  if(debug) cout << "Compact Plots" << endl;

  vector<TH1F*> ttbar_all = {ttbar_hist, ttbar_rebin_5bins, ttbar_rebin_10bins};
  vector<TH1F*> data_all = {data_hist, data_rebin_5bins, data_rebin_10bins};

  vector<vector<TH1F*>> sum_up_all = {sum_up_hist, sum_up_rebin_5bins, sum_up_rebin_10bins};
  vector<vector<TH1F*>> sum_down_all = {sum_down_hist, sum_down_rebin_5bins, sum_down_rebin_10bins};

  // #############################################################################################
  // Normalize Vectors ###########################################################################
  if(debug) cout << "Normalized Vector" << endl;
  // Vectors have the same order like the ones above
  vector<TH1F*> ttbar_all_norm = normalize(ttbar_all);
  vector<TH1F*> data_all_norm = normalize(data_all);

  vector<vector<TH1F*>> sum_up_all_norm = normalize(sum_up_all);
  vector<vector<TH1F*>> sum_down_all_norm = normalize(sum_down_all);

  // #############################################################################################
  // Combine Non normalized and normalized #######################################################
  if(debug) cout << "(Non-)Normalized Vector" << endl;

  vector<vector<TH1F*>> ttbar_general = {ttbar_all, ttbar_all_norm};
  vector<vector<TH1F*>> data_general = {data_all, data_all_norm};

  vector<vector<vector<TH1F*>>> sum_up_general = {sum_up_all, sum_up_all_norm};
  vector<vector<vector<TH1F*>>> sum_down_general = {sum_down_all, sum_down_all_norm};

  /*
  ██████  ██       ██████  ████████ ███████
  ██   ██ ██      ██    ██    ██    ██
  ██████  ██      ██    ██    ██    ███████
  ██      ██      ██    ██    ██         ██
  ██      ███████  ██████     ██    ███████
  */
  if(debug) cout << "Plots" << endl;

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);
  TLegend *leg;

  for(unsigned int a=0; a<ttbar_general.size();a++){ //norm and not norm
    for(unsigned int b=0; b<ttbar_general[a].size();b++){ // not rebin and rebin
      // ttbar hist style
      ttbar_general[a].at(b)->SetLineWidth(2);
      ttbar_general[a].at(b)->SetLineColor(kRed);
      // data hist style
      data_general[a].at(b)->SetMarkerStyle(8);
      data_general[a].at(b)->SetMarkerColor(kBlack);
    }
  }

  // #############################################################################################
  // Single Plots ################################################################################
  if(debug) cout << "Sinlge Plots" << endl;

  for(unsigned int a=0; a<sum_up_general.size();a++){ // not norm and norm (ttbar, data and variations)
    for(unsigned int b=0; b<sum_up_general.at(a).size();b++){ // not rebin and rebin (ttbar, data and variations)
      for(unsigned int c=0; c<sum_up_general.at(a).at(b).size();c++){ // variation factors (sqrt_2, 2, 4)
        sum_up_general[a].at(b).at(c)->SetTitle(FSR[c]);
        sum_up_general[a].at(b).at(c)->GetXaxis()->SetRangeUser(0, 400);
        sum_up_general[a].at(b).at(c)->GetYaxis()->SetRangeUser(0, sum_up_general[a].at(b).at(c)->GetMaximum()*1.2);
        sum_up_general[a].at(b).at(c)->GetXaxis()->SetNdivisions(505);
        sum_up_general[a].at(b).at(c)->GetYaxis()->SetNdivisions(505);
        sum_up_general[a].at(b).at(c)->GetXaxis()->SetTitleSize(0.05);
        sum_up_general[a].at(b).at(c)->GetYaxis()->SetTitleSize(0.04);
        sum_up_general[a].at(b).at(c)->GetXaxis()->SetTitleOffset(0.9);
        sum_up_general[a].at(b).at(c)->GetYaxis()->SetTitleOffset(1.5);
        sum_up_general[a].at(b).at(c)->GetXaxis()->SetTitle("#tau_{32}");
        sum_up_general[a].at(b).at(c)->GetYaxis()->SetTitle("");
        sum_up_general[a].at(b).at(c)->SetLineWidth(2);
        sum_up_general[a].at(b).at(c)->SetLineColor(color[c]);

        sum_down_general[a].at(b).at(c)->SetLineWidth(2);
        sum_down_general[a].at(b).at(c)->SetLineColor(color[c]);

        TCanvas *A = new TCanvas("A", "A", 600, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        sum_up_general[a].at(b).at(c)->Draw("HIST");
        sum_down_general[a].at(b).at(c)->Draw("SAME HIST");
        ttbar_general[a].at(b)->Draw("SAME HIST");
        data_general[a].at(b)->Draw("SAME P");
        if(a==0){ // not norm
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_"+year+".pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_Rebin5_"+year+".pdf");
          if(b==2) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_Rebin10_"+year+".pdf");
        }
        if(a==1){ // norm
          if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_norm_"+year+".pdf");
          if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_Rebin5_norm_"+year+".pdf");
          if(b==2) A->SaveAs(save_path+year+"/tau32_SYS_"+FSR[c]+"_Rebin10_norm_"+year+".pdf");
        }
        delete A;
      }
    }
  }

  // #############################################################################################
  // All in one Plot #############################################################################
  if(debug) cout << "All" << endl;

  for(unsigned int a=0; a<sum_up_general.size();a++){ // not norm and norm
    for(unsigned int b=0; b<sum_up_general[a].size();b++){ // not rebin and rebin
      ttbar_general[a].at(b)->SetTitle("");
      ttbar_general[a].at(b)->GetXaxis()->SetRangeUser(0, 400);
      ttbar_general[a].at(b)->GetYaxis()->SetRangeUser(0, get_highest_peak(sum_up_general[a].at(b))*1.2);
      ttbar_general[a].at(b)->GetXaxis()->SetNdivisions(505);
      ttbar_general[a].at(b)->GetYaxis()->SetNdivisions(505);
      ttbar_general[a].at(b)->GetXaxis()->SetTitleSize(0.05);
      ttbar_general[a].at(b)->GetYaxis()->SetTitleSize(0.04);
      ttbar_general[a].at(b)->GetXaxis()->SetTitleOffset(0.9);
      ttbar_general[a].at(b)->GetYaxis()->SetTitleOffset(1.5);
      ttbar_general[a].at(b)->GetXaxis()->SetTitle("#tau_{32}");
      ttbar_general[a].at(b)->GetYaxis()->SetTitle("");

      for(unsigned int c=0; c<sum_up_general[a].at(b).size(); c++){
        sum_up_general[a].at(b).at(c)->SetLineWidth(2);
        sum_up_general[a].at(b).at(c)->SetLineColor(color[c]);

        sum_down_general[a].at(b).at(c)->SetLineWidth(2);
        sum_down_general[a].at(b).at(c)->SetLineColor(color[c]);
      }

      TCanvas *A = new TCanvas("A", "A", 600, 600);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      ttbar_general[a].at(b)->Draw("HIST");
      data_general[a].at(b)->Draw("SAME P");
      for(unsigned int c=0; c<sum_up_general[a].at(b).size();c++){
        sum_up_general[a].at(b).at(c)->Draw("SAME HIST");
        sum_down_general[a].at(b).at(c)->Draw("SAME HIST");
      }
      leg = new TLegend(0.25,0.65,0.45,0.85);
      leg->SetTextSize(0.5);
      leg->AddEntry(ttbar_general[a].at(b),"nominal","l");
      for(unsigned int c=0; c<sum_up_general[a].at(b).size();c++) leg->AddEntry(sum_up_general[a].at(b).at(c),FSR[c],"l");
      leg->SetTextSize(0.05);
      leg->Draw();
      gPad->RedrawAxis();
      if(a==0){
        if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_all_"+year+".pdf");
        if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin5_"+year+".pdf");
        if(b==2) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin10_"+year+".pdf");
      }
      if(a==1){
        if(b==0) A->SaveAs(save_path+year+"/tau32_SYS_all_norm_"+year+".pdf");
        if(b==1) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin5_norm_"+year+".pdf");
        if(b==2) A->SaveAs(save_path+year+"/tau32_SYS_all_Rebin10_norm_"+year+".pdf");
      }
      leg->Clear();
      delete A;
    }
  }
}
