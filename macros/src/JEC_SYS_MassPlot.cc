#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

bool forTalk = true;
bool debug   = false;

using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2, bool isData);

// -------------------------------------------------------------------------------------------------------
// ---------------------------------------- FUNCTIONS ----------------------------------------------------
// -------------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------------
void main_plot_settings(TH1F* hist, int x_min, int x_max, int color, TString title, TString xAxis, TString yAxis)
{
  hist->SetTitle(title);
  hist->GetXaxis()->SetRangeUser(x_min, x_max);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.1);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(19);
  hist->GetXaxis()->SetTitle(xAxis);
  hist->GetYaxis()->SetTitle(yAxis);
  hist->SetLineWidth(2);
  hist->SetLineColor(color);
}

// -------------------------------------------------------------------------------------------------------
void data_ratio_settings(TH1F* data, int x_min, int x_max)
{
  data->GetXaxis()->SetTickLength(0.07);
  data->GetXaxis()->SetTitleSize(25);
  data->GetXaxis()->SetTitleFont(43);
  data->GetXaxis()->SetTitleOffset(4.0);
  data->GetXaxis()->SetLabelFont(43);
  data->GetXaxis()->SetLabelSize(21);
  data->GetXaxis()->SetLabelOffset(0.035);
  data->GetYaxis()->SetTitle("#frac{MC}{Data}");
  data->GetYaxis()->CenterTitle();
  data->GetYaxis()->SetTitleSize(20);
  data->GetYaxis()->SetTitleFont(43);
  data->GetYaxis()->SetTitleOffset(2.2);
  data->GetYaxis()->SetLabelFont(43);
  data->GetYaxis()->SetLabelSize(19);
  data->GetYaxis()->SetLabelOffset(0.009);
  data->GetYaxis()->SetNdivisions(505);
  data->SetTitle(" ");
  data->GetYaxis()->SetRangeUser(0.2, 1.8);
  data->GetXaxis()->SetRangeUser(x_min, x_max);

  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(8);
  data->SetMarkerSize(0.5);
}

// -------------------------------------------------------------------------------------------------------
TH1F* GetRatio(TH1F* h1, TH1F* h2, vector<double> data_bin_error, bool isData){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetNbinsX();
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      if(isData) ratio->SetBinContent(i, 0);
      else       ratio->SetBinContent(i, 1);

      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}

// -------------------------------------------------------------------------------------------------------
void draw_plot(vector<TH1F*> hists, vector<TString> path, vector<double> mean)
{
  // path    = {save_path, channel, var, addition, year, norm}
  // hists   = {nom, up, down, data}
  // rations = {nom/data, up/data, down/data, data/data}

  TString unique_path = path[0]+"/"+path[1]+"/Wjet_mass_"+path[2]+path[3]+"_"+path[4]+path[5]+".pdf";

  TLegend *leg = new TLegend(0.4,0.15,0.6,0.35);
  TCanvas *A = new TCanvas(unique_path, unique_path, 600, 600);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  hists[0]->Draw("HIST");
  hists[1]->Draw("SAME HIST");
  hists[2]->Draw("SAME HIST");
  hists[3]->Draw("SAME P");

  TString nom  = "Nominal  <#it{m}_{W}^{reco}>="+dtos(mean[0], 2)+" GeV";
  TString up   = path[2]+" up   <#it{m}_{W}^{reco}>="+dtos(mean[1], 2)+" GeV";
  TString down = path[2]+" down <#it{m}_{W}^{reco}>="+dtos(mean[2], 2)+" GeV";

  leg->AddEntry(hists[3],"Data","pl");
  leg->AddEntry(hists[0],nom ,"l");
  leg->AddEntry(hists[1],up  ,"l");
  leg->AddEntry(hists[2],down,"l");

  leg->SetTextSize(0.03);
  leg->Draw();

  gPad->RedrawAxis();

  A->SaveAs(unique_path);

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
void draw_plot_ratio(vector<TH1F*> hists, vector<TH1F*> ratios, vector<TString> path, vector<double> mean)
{
  //          path    = {save_path, channel, var, addition, year, norm}
  //          hists   = {nom, up, down, data}
  // plotAll; hists   = {nom, jup, jdown, xup, xdown, data}
  //          rations = {nom/data, up/data, down/data, data/data}

  // for TCanvas name; Also used as save_path
  TString unique_path = path[0]+"/"+path[1]+"/Wjet_mass_"+path[2]+path[3]+"_"+path[4]+path[5]+".pdf";
  bool diffLeg_ll = false;
  bool diffLeg_lh = false;
  bool plotAll    = false;
  double nHists = hists.size();
  if(nHists==6) plotAll = true;
  if(nHists!=ratios.size()) throw runtime_error("Not the same amount of histograms and ratios in 'draw_plot_ratio'");

  if(debug) cout << path[1] << " " << path[2]  << " " << path[3] << " " << path[4] << endl;
  if(strcmp(path[1], "combine")==0&&strcmp(path[3], "_ll")==0) diffLeg_ll = true;
  if(strcmp(path[1], "combine")==0&&strcmp(path[3], "_lh")==0) diffLeg_lh = true;
  if(debug) cout << "Different Legend (combine, ll): " << diffLeg_ll << endl;
  if(debug) cout << "Different Legend (combine, lh): " << diffLeg_lh << endl;

  TLegend *leg = new TLegend(0.40,0.15,0.65,0.35);
  if(plotAll)             leg = new TLegend(0.40,0.10,0.65,0.37);
  if(diffLeg_ll)          leg = new TLegend(0.25,0.6,0.5,0.8);
  if(diffLeg_lh)          leg = new TLegend(0.45,0.15,0.70,0.35);
  if(diffLeg_lh&&plotAll) leg = new TLegend(0.45,0.10,0.70,0.37);
  if(diffLeg_ll&&plotAll) leg = new TLegend(0.25,0.55,0.5,0.82);

  TCanvas *A = new TCanvas(unique_path, unique_path, 600, 600);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  hists[0]->GetXaxis()->SetLabelSize(0);
  hists[0]->Draw("HIST");
  for(unsigned int i=1; i<nHists-1; i++) hists[i]->Draw("SAME HIST");
  hists[nHists-1]->Draw("SAME P");

  TString jup, jdown, xup, xdown, up, down;
  TString nom   = "t#bar{t}    ";
  if(strcmp(path[2], "JEC")==0)
  {
    up   = "t#bar{t} f^{JEC}=+1 ";
    down = "t#bar{t} f^{JEC}=-1  ";
  }
  else
  {
    up   = "t#bar{t} f^{XCone}=+1 ";
    down = "t#bar{t} f^{XCone}=-1  ";
  }
  if(plotAll){
    jup   = "t#bar{t} f^{JEC}=+1   ";
    jdown = "t#bar{t} f^{JEC}=-1    ";
    xup   = "t#bar{t} f^{XCone}=+1 ";
    xdown = "t#bar{t} f^{XCone}=-1  ";
  }
  TString data  = "Data    ";

  leg->SetNColumns(2);

  leg->AddEntry(hists[nHists-1],data,"pl");
  leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[nHists-1], 2)+" GeV", "");
  leg->AddEntry(hists[0],nom ,"l");
  leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[0], 2)+" GeV", "");
  if(plotAll){
    leg->AddEntry(hists[1],jup  ,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[1], 2)+" GeV", "");
    leg->AddEntry(hists[2],jdown,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[2], 2)+" GeV", "");
    leg->AddEntry(hists[3],xup  ,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[3], 2)+" GeV", "");
    leg->AddEntry(hists[4],xdown,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[4], 2)+" GeV", "");
  }
  else{
    leg->AddEntry(hists[1],up  ,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[1], 2)+" GeV", "");
    leg->AddEntry(hists[2],down,"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{W}^{reco}>="+dtos(mean[2], 2)+" GeV", "");
  }

  leg->SetTextSize(0.03);
  leg->Draw();

  // Ratio Plot kommt in unteres pad
  A->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  ratios[nHists-1]->GetXaxis()->SetTitle("#it{m}_{W}^{reco} [GeV]");
  // ratios[nHists-1]->GetXaxis()->SetTitleOffset(4);

  ratios[nHists-1]->Draw("PE");
  ratios[0]->Draw("HIST SAME");
  for(unsigned int i=1; i<nHists-1; i++) ratios[i]->Draw("SAME HIST");

  gPad->RedrawAxis();
  A->SaveAs(unique_path);

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
// ------------------------------------------ MAIN -------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]){
  TString mTop_ = "mtop1695";           // CHANGE_MTOP

  // #################################################################################################
  // Settings ########################################################################################

  // Declare different variables used in the code ------------------------------
  print_seperater();

  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_SYS <rebin_number>\n";
    cout << "hist->Rebin(rebin_number)\n";
    return 0;
  }

  // Default -------------------------------------------------------------------
  cout.precision(6);
  Int_t   oldLevel     = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ...
  gErrorIgnoreLevel    = kWarning;          // suppress TCanvas output

  TString reconst      = "btag";            // match btag_cut btag_sel compare min_mass
  TString MC_uncert    = "central";         // central avg nominal
  bool    add_mTop     = false;             // !!!
  bool    usePeak_in   = false;
  bool    one_Ptbin    = false;

  // Input ---------------------------------------------------------------------
  int bin_width   = stoi(argv[1]);

  // Directories ---------------------------------------------------------------

  TString save_path = get_save_path();
  save_path += "/Plots/JEC_SYS";
  save_path = creat_folder_and_path(save_path, "MassPlots");
  // if reconst is not btag include another folder here
  save_path = creat_folder_and_path(save_path, "BinWidth_"+to_string(bin_width));
  // if add_mTop is true include another folder here
  cout << save_path << endl;
  creat_folder(save_path+"/muon");
  creat_folder(save_path+"/elec");
  creat_folder(save_path+"/combine");

  // Directories ---------------------------------------------------------------
  vector<int> index_year   = {0,1,2,3};
  vector<TString> years    = {"2016","2017","2018","combine"};
  vector<TString> channels = {"muon","elec","combine"};

  // -------------------------------------------------------------------------------------------------------
  // ------------------------------------------ MAIN -------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------

  for(int ptbin=0; ptbin<5; ptbin++){

    // Define ptbins -----------------------------------------------------------
    if((ptbin==4)) continue;
    TString addition="";
    if(ptbin==0) addition="_hh";
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="";
    string add = (string) addition;

    if(debug) cout << setw(12) << "\n|----------------- " << setw(9) << centered(add) << " -------------------|\n" << endl;

    cout << '\n';
    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
    // if(reconst=="btag")      w_mass = hist_class+"wmass_match";
    if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(reconst=="btag"&&ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(reconst=="btag"&&ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";
    if(reconst=="btag"&&ptbin==4) w_mass = hist_class+"wmass_match";
    // if(reconst=="btag")      w_mass = hist_class+"wjet_pt_match_S1divW_W";

    // #################################################################################################
    // Get Hists #######################################################################################
    if(debug) cout << "Get Hists" << endl;

    vector<TH1F*> h_data_muon          = get_all_hists(data_muon, w_mass);
    vector<TH1F*> h_data_elec          = get_all_hists(data_elec, w_mass);
    vector<TH1F*> h_data_combine       = combine_channels(h_data_muon, h_data_elec);

    vector<TH1F*> h_ttbar_muon         = get_all_hists(ttbar_muon, w_mass);
    vector<TH1F*> h_ttbar_elec         = get_all_hists(ttbar_elec, w_mass);
    vector<TH1F*> h_ttbar_combine      = combine_channels(h_ttbar_muon, h_ttbar_elec);

    vector<TH1F*> h_st_muon            = get_all_hists(st_muon, w_mass);
    vector<TH1F*> h_st_elec            = get_all_hists(st_elec, w_mass);
    vector<TH1F*> h_st_combine         = combine_channels(h_st_muon, h_st_elec);

    vector<TH1F*> h_wjets_muon         = get_all_hists(wjets_muon, w_mass);
    vector<TH1F*> h_wjets_elec         = get_all_hists(wjets_elec, w_mass);
    vector<TH1F*> h_wjets_combine      = combine_channels(h_wjets_muon, h_wjets_elec);

    vector<TH1F*> h_other_muon         = get_all_hists(other_muon, w_mass);
    vector<TH1F*> h_other_elec         = get_all_hists(other_elec, w_mass);
    vector<TH1F*> h_other_combine      = combine_channels(h_other_muon, h_other_elec);

    // jec -----------------------------------------------------------------------
    vector<TH1F*> h_jec_up_muon        = get_all_hists(jec_up_muon, w_mass);
    vector<TH1F*> h_jec_up_elec        = get_all_hists(jec_up_elec, w_mass);
    vector<TH1F*> h_jec_up_combine     = combine_channels(h_jec_up_muon, h_jec_up_elec);

    vector<TH1F*> h_jec_down_muon      = get_all_hists(jec_down_muon, w_mass);
    vector<TH1F*> h_jec_down_elec      = get_all_hists(jec_down_elec, w_mass);
    vector<TH1F*> h_jec_down_combine   = combine_channels(h_jec_down_muon, h_jec_down_elec);

    // cor -----------------------------------------------------------------------
    vector<TH1F*> h_cor_up_muon        = get_all_hists(cor_up_muon, w_mass);
    vector<TH1F*> h_cor_up_elec        = get_all_hists(cor_up_elec, w_mass);
    vector<TH1F*> h_cor_up_combine     = combine_channels(h_cor_up_muon, h_cor_up_elec);

    vector<TH1F*> h_cor_down_muon      = get_all_hists(cor_down_muon, w_mass);
    vector<TH1F*> h_cor_down_elec      = get_all_hists(cor_down_elec, w_mass);
    vector<TH1F*> h_cor_down_combine   = combine_channels(h_cor_down_muon, h_cor_down_elec);

    // #################################################################################################
    // Combine Hists ###################################################################################
    if(debug) cout << "Combine Hists" << endl;

    // Data ----------------------------------------------------------------------

    // dummy: in order to keep code clean and faster to write

    vector<TH1F*> all_data_muon              = h_data_muon;        // dummy
    vector<TH1F*> all_data_elec              = h_data_elec;        // dummy
    vector<TH1F*> all_data_combine           = h_data_combine;     // dummy

    vector<TH1F*> all_data_muon_rebin        = rebin(all_data_muon, bin_width);
    vector<TH1F*> all_data_elec_rebin        = rebin(all_data_elec, bin_width);
    vector<TH1F*> all_data_combine_rebin     = rebin(all_data_combine, bin_width);

    vector<TH1F*> all_data_muon_norm         = normalize(all_data_muon_rebin);
    vector<TH1F*> all_data_elec_norm         = normalize(all_data_elec_rebin);
    vector<TH1F*> all_data_combine_norm      = normalize(all_data_combine_rebin);

    // Bkg -----------------------------------------------------------------------

    vector<TH1F*> h_bkg_muon                 = AddHists(h_st_muon, h_wjets_muon, h_other_muon, 1);
    vector<TH1F*> h_bkg_elec                 = AddHists(h_st_elec, h_wjets_elec, h_other_elec, 1);
    vector<TH1F*> h_bkg_combine              = AddHists(h_st_combine, h_wjets_combine, h_other_combine, 1);

    // TTbar ---------------------------------------------------------------------
    if(debug) cout << "TTbar" << endl;

    vector<TH1F*> all_ttbar_muon             = AddHists(h_bkg_muon, h_ttbar_muon, 1);
    vector<TH1F*> all_ttbar_elec             = AddHists(h_bkg_elec, h_ttbar_elec, 1);
    vector<TH1F*> all_ttbar_combine          = AddHists(h_bkg_combine, h_ttbar_combine, 1);

    vector<TH1F*> all_ttbar_muon_rebin       = rebin(all_ttbar_muon, bin_width);
    vector<TH1F*> all_ttbar_elec_rebin       = rebin(all_ttbar_elec, bin_width);
    vector<TH1F*> all_ttbar_combine_rebin    = rebin(all_ttbar_combine, bin_width);

    vector<TH1F*> all_ttbar_muon_norm        = normalize(all_ttbar_muon_rebin);
    vector<TH1F*> all_ttbar_elec_norm        = normalize(all_ttbar_elec_rebin);
    vector<TH1F*> all_ttbar_combine_norm     = normalize(all_ttbar_combine_rebin);

    // JEC -----------------------------------------------------------------------
    if(debug) cout << "JEC" << endl;

    vector<TH1F*> all_jec_up_muon            = AddHists(h_bkg_muon, h_jec_up_muon, 1);
    vector<TH1F*> all_jec_up_elec            = AddHists(h_bkg_elec, h_jec_up_elec, 1);
    vector<TH1F*> all_jec_up_combine         = AddHists(h_bkg_combine, h_jec_up_combine, 1);

    vector<TH1F*> all_jec_down_muon          = AddHists(h_bkg_muon, h_jec_down_muon, 1);
    vector<TH1F*> all_jec_down_elec          = AddHists(h_bkg_elec, h_jec_down_elec, 1);
    vector<TH1F*> all_jec_down_combine       = AddHists(h_bkg_combine, h_jec_down_combine, 1);

    vector<TH1F*> all_jec_down_muon_rebin    = rebin(all_jec_down_muon, bin_width);
    vector<TH1F*> all_jec_down_elec_rebin    = rebin(all_jec_down_elec, bin_width);
    vector<TH1F*> all_jec_down_combine_rebin = rebin(all_jec_down_combine, bin_width);

    vector<TH1F*> all_jec_up_muon_rebin      = rebin(all_jec_up_muon, bin_width);
    vector<TH1F*> all_jec_up_elec_rebin      = rebin(all_jec_up_elec, bin_width);
    vector<TH1F*> all_jec_up_combine_rebin   = rebin(all_jec_up_combine, bin_width);

    vector<TH1F*> all_jec_down_muon_norm     = normalize(all_jec_down_muon_rebin);
    vector<TH1F*> all_jec_down_elec_norm     = normalize(all_jec_down_elec_rebin);
    vector<TH1F*> all_jec_down_combine_norm  = normalize(all_jec_down_combine_rebin);

    vector<TH1F*> all_jec_up_muon_norm       = normalize(all_jec_up_muon_rebin);
    vector<TH1F*> all_jec_up_elec_norm       = normalize(all_jec_up_elec_rebin);
    vector<TH1F*> all_jec_up_combine_norm    = normalize(all_jec_up_combine_rebin);

    // COR -----------------------------------------------------------------------
    if(debug) cout << "COR" << endl;

    vector<TH1F*> all_cor_up_muon            = AddHists(h_bkg_muon, h_cor_up_muon, 1);
    vector<TH1F*> all_cor_up_elec            = AddHists(h_bkg_elec, h_cor_up_elec, 1);
    vector<TH1F*> all_cor_up_combine         = AddHists(h_bkg_combine, h_cor_up_combine, 1);

    vector<TH1F*> all_cor_down_muon          = AddHists(h_bkg_muon, h_cor_down_muon, 1);
    vector<TH1F*> all_cor_down_elec          = AddHists(h_bkg_elec, h_cor_down_elec, 1);
    vector<TH1F*> all_cor_down_combine       = AddHists(h_bkg_combine, h_cor_down_combine, 1);

    vector<TH1F*> all_cor_up_muon_rebin      = rebin(all_cor_up_muon, bin_width);
    vector<TH1F*> all_cor_up_elec_rebin      = rebin(all_cor_up_elec, bin_width);
    vector<TH1F*> all_cor_up_combine_rebin   = rebin(all_cor_up_combine, bin_width);

    vector<TH1F*> all_cor_down_muon_rebin    = rebin(all_cor_down_muon, bin_width);
    vector<TH1F*> all_cor_down_elec_rebin    = rebin(all_cor_down_elec, bin_width);
    vector<TH1F*> all_cor_down_combine_rebin = rebin(all_cor_down_combine, bin_width);

    vector<TH1F*> all_cor_up_muon_norm       = normalize(all_cor_up_muon_rebin);
    vector<TH1F*> all_cor_up_elec_norm       = normalize(all_cor_up_elec_rebin);
    vector<TH1F*> all_cor_up_combine_norm    = normalize(all_cor_up_combine_rebin);

    vector<TH1F*> all_cor_down_muon_norm     = normalize(all_cor_down_muon_rebin);
    vector<TH1F*> all_cor_down_elec_norm     = normalize(all_cor_down_elec_rebin);
    vector<TH1F*> all_cor_down_combine_norm  = normalize(all_cor_down_combine_rebin);

    // -------------------------------------------------------------------------------------------------------
    // ------------------------------------------ Plots ------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------

    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetLegendBorderSize(0);

    TString xAxis      = "#it{m}_{W}^{reco} [GeV]";
    TString yAxis_norm = "#Delta N/N";

    for(int year=0; year<4; year++){

      string sYear = (string) years[year];
      if(debug) cout << setw(12) << "\n|------------ " << setw(8) << centered(sYear) << " --------------|\n" << endl;
      if(year!=3) continue; // skip single years

      // #################################################################################################
      // Settings ########################################################################################

      // Masspeack -------------------------------------------------------------
      if(debug) cout << "Masspeak Bins" << endl;

      // Get bins with bin-content>Limit
      vector<double> PeakLimit;
      vector<int> peak_bins_muon, peak_bins_elec, peak_bins_combine;
      peak_bins_muon    = bins_upper_limit(all_data_muon_rebin[year], 75);
      peak_bins_elec    = bins_upper_limit(all_data_elec_rebin[year], 75);
      peak_bins_combine = bins_upper_limit(all_data_combine_rebin[year], 100);

      // Bin Center
      vector<double> bin_centers_muon    = bin_center_upper_limit(all_data_muon_rebin[year],    peak_bins_muon,    bin_width);
      vector<double> bin_centers_elec    = bin_center_upper_limit(all_data_elec_rebin[year],    peak_bins_elec,    bin_width);
      vector<double> bin_centers_combine = bin_center_upper_limit(all_data_combine_rebin[year], peak_bins_combine, bin_width);

      // Path
      vector<TString> path_muon_jec         = {save_path, "muon",    "JEC", addition, years[year], ""};
      vector<TString> path_elec_jec         = {save_path, "elec",    "JEC", addition, years[year], ""};
      vector<TString> path_combine_jec      = {save_path, "combine", "JEC", addition, years[year], ""};
      vector<TString> path_muon_jec_norm    = {save_path, "muon",    "JEC", addition, years[year], "_norm"};
      vector<TString> path_elec_jec_norm    = {save_path, "elec",    "JEC", addition, years[year], "_norm"};
      vector<TString> path_combine_jec_norm = {save_path, "combine", "JEC", addition, years[year], "_norm"};

      vector<TString> path_muon_cor         = {save_path, "muon",    "COR", addition, years[year], ""};
      vector<TString> path_elec_cor         = {save_path, "elec",    "COR", addition, years[year], ""};
      vector<TString> path_combine_cor      = {save_path, "combine", "COR", addition, years[year], ""};
      vector<TString> path_muon_cor_norm    = {save_path, "muon",    "COR", addition, years[year], "_norm"};
      vector<TString> path_elec_cor_norm    = {save_path, "elec",    "COR", addition, years[year], "_norm"};
      vector<TString> path_combine_cor_norm = {save_path, "combine", "COR", addition, years[year], "_norm"};

      vector<TString> path_muon_all         = {save_path, "muon",    "all", addition, years[year], ""};
      vector<TString> path_elec_all         = {save_path, "elec",    "all", addition, years[year], ""};
      vector<TString> path_combine_all      = {save_path, "combine", "all", addition, years[year], ""};
      vector<TString> path_muon_all_norm    = {save_path, "muon",    "all", addition, years[year], "_norm"};
      vector<TString> path_elec_all_norm    = {save_path, "elec",    "all", addition, years[year], "_norm"};
      vector<TString> path_combine_all_norm = {save_path, "combine", "all", addition, years[year], "_norm"};

      // ttbar -----------------------------------------------------------------

      // void main_plot_settings(hist, x_min, x_max, color, title, TString xAxis, TString yAxis,save_path)
      main_plot_settings(all_ttbar_muon_rebin[year],    0, 180, kRed, "", xAxis, "Events");
      main_plot_settings(all_ttbar_elec_rebin[year],    0, 180, kRed, "", xAxis, "Events");
      main_plot_settings(all_ttbar_combine_rebin[year], 0, 180, kRed, "", xAxis, "Events");
      main_plot_settings(all_ttbar_muon_norm[year],     0, 180, kRed, "", xAxis,   "a.u.");
      main_plot_settings(all_ttbar_elec_norm[year],     0, 180, kRed, "", xAxis,   "a.u.");
      main_plot_settings(all_ttbar_combine_norm[year],  0, 180, kRed, "", xAxis,   "a.u.");

      // data ------------------------------------------------------------------

      // data_plot_settings(TH1F* hist)
      data_plot_settings(all_data_muon_rebin[year]);
      data_plot_settings(all_data_elec_rebin[year]);
      data_plot_settings(all_data_combine_rebin[year]);
      data_plot_settings(all_data_muon_norm[year]);
      data_plot_settings(all_data_elec_norm[year]);
      data_plot_settings(all_data_combine_norm[year]);

      // jms -------------------------------------------------------------------

      // add_plot_settings(TH1F* hist, int color=1, int style=kSolid, int width=2)
      add_plot_settings(all_jec_up_muon_rebin[year], kBlue);
      add_plot_settings(all_jec_up_elec_rebin[year], kBlue);
      add_plot_settings(all_jec_up_combine_rebin[year], kBlue);
      add_plot_settings(all_jec_up_muon_norm[year], kBlue);
      add_plot_settings(all_jec_up_elec_norm[year], kBlue);
      add_plot_settings(all_jec_up_combine_norm[year], kBlue);

      add_plot_settings(all_jec_down_muon_rebin[year], kBlue, 7);
      add_plot_settings(all_jec_down_elec_rebin[year], kBlue, 7);
      add_plot_settings(all_jec_down_combine_rebin[year], kBlue, 7);
      add_plot_settings(all_jec_down_muon_norm[year], kBlue, 7);
      add_plot_settings(all_jec_down_elec_norm[year], kBlue, 7);
      add_plot_settings(all_jec_down_combine_norm[year], kBlue, 7);

      add_plot_settings(all_cor_up_muon_rebin[year], kGreen+2);
      add_plot_settings(all_cor_up_elec_rebin[year], kGreen+2);
      add_plot_settings(all_cor_up_combine_rebin[year], kGreen+2);
      add_plot_settings(all_cor_up_muon_norm[year], kGreen+2);
      add_plot_settings(all_cor_up_elec_norm[year], kGreen+2);
      add_plot_settings(all_cor_up_combine_norm[year], kGreen+2);

      add_plot_settings(all_cor_down_muon[year], kGreen+2, 7);
      add_plot_settings(all_cor_down_elec[year], kGreen+2, 7);
      add_plot_settings(all_cor_down_combine[year], kGreen+2, 7);
      add_plot_settings(all_cor_down_muon_norm[year], kGreen+2, 7);
      add_plot_settings(all_cor_down_elec_norm[year], kGreen+2, 7);
      add_plot_settings(all_cor_down_combine_norm[year], kGreen+2, 7);

      // #################################################################################################
      // Mass Plots ######################################################################################
      if(debug) cout << "Mass Plots" << endl;

      // Mean ------------------------------------------------------------------

      double mean_ttbar_muon    = trunc_mean(all_ttbar_muon_norm[year],    bin_centers_muon[0], bin_centers_muon[1]);
      double mean_data_muon     = trunc_mean(all_data_muon_norm[year],     bin_centers_muon[0], bin_centers_muon[1]);
      double mean_jec_up_muon   = trunc_mean(all_jec_up_muon_norm[year],   bin_centers_muon[0], bin_centers_muon[1]);
      double mean_jec_down_muon = trunc_mean(all_jec_down_muon_norm[year], bin_centers_muon[0], bin_centers_muon[1]);
      double mean_cor_up_muon   = trunc_mean(all_cor_up_muon_norm[year],   bin_centers_muon[0], bin_centers_muon[1]);
      double mean_cor_down_muon = trunc_mean(all_cor_down_muon_norm[year], bin_centers_muon[0], bin_centers_muon[1]);
      vector<double> means_jec_muon  = {mean_ttbar_muon, mean_jec_up_muon, mean_jec_down_muon};
      vector<double> means_cor_muon  = {mean_ttbar_muon, mean_cor_up_muon, mean_cor_down_muon};
      vector<double> means_all_muon  = {mean_ttbar_muon, mean_jec_up_muon, mean_jec_down_muon, mean_cor_up_muon, mean_cor_down_muon, mean_data_muon};

      double mean_ttbar_elec    = trunc_mean(all_ttbar_elec_norm[year],    bin_centers_elec[0], bin_centers_elec[1]);
      double mean_data_elec     = trunc_mean(all_data_elec_norm[year],     bin_centers_elec[0], bin_centers_elec[1]);
      double mean_jec_up_elec   = trunc_mean(all_jec_up_elec_norm[year],   bin_centers_elec[0], bin_centers_elec[1]);
      double mean_jec_down_elec = trunc_mean(all_jec_down_elec_norm[year], bin_centers_elec[0], bin_centers_elec[1]);
      double mean_cor_up_elec   = trunc_mean(all_cor_up_elec_norm[year],   bin_centers_elec[0], bin_centers_elec[1]);
      double mean_cor_down_elec = trunc_mean(all_cor_down_elec_norm[year], bin_centers_elec[0], bin_centers_elec[1]);
      vector<double> means_jec_elec  = {mean_ttbar_elec, mean_jec_up_elec, mean_jec_down_elec};
      vector<double> means_cor_elec  = {mean_ttbar_elec, mean_cor_up_elec, mean_cor_down_elec};
      vector<double> means_all_elec  = {mean_ttbar_elec, mean_jec_up_elec, mean_jec_down_elec, mean_cor_up_elec, mean_cor_down_elec, mean_data_elec};

      double mean_ttbar_combine    = trunc_mean(all_ttbar_combine_norm[year],    bin_centers_combine[0], bin_centers_combine[1]);
      double mean_data_combine     = trunc_mean(all_data_combine_norm[year] ,    bin_centers_combine[0], bin_centers_combine[1]);
      double mean_jec_up_combine   = trunc_mean(all_jec_up_combine_norm[year],   bin_centers_combine[0], bin_centers_combine[1]);
      double mean_jec_down_combine = trunc_mean(all_jec_down_combine_norm[year], bin_centers_combine[0], bin_centers_combine[1]);
      double mean_cor_up_combine   = trunc_mean(all_cor_up_combine_norm[year],   bin_centers_combine[0], bin_centers_combine[1]);
      double mean_cor_down_combine = trunc_mean(all_cor_down_combine_norm[year], bin_centers_combine[0], bin_centers_combine[1]);
      vector<double> means_jec_combine  = {mean_ttbar_combine, mean_jec_up_combine, mean_jec_down_combine, mean_data_combine};
      vector<double> means_cor_combine  = {mean_ttbar_combine, mean_cor_up_combine, mean_cor_down_combine, mean_data_combine};
      vector<double> means_all_combine  = {mean_ttbar_combine, mean_jec_up_combine, mean_jec_down_combine, mean_cor_up_combine, mean_cor_down_combine, mean_data_combine};

      cout << left; // for setw()
      cout << "ttbar    - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_ttbar_muon       << " GeV" << endl;
      cout << "data     - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_data_muon        << " GeV" << endl;
      cout << "jec up   - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_jec_up_muon      << " GeV" << endl;
      cout << "jec down - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_jec_down_muon    << " GeV" << endl;
      cout << "cor up   - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_cor_up_muon      << " GeV" << endl;
      cout << "cor down - " << years[year] << setw(12) << " - muon: "    << setw(4) << mean_cor_down_muon    << " GeV" << endl;
      cout << "ttbar    - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_ttbar_elec       << " GeV" << endl;
      cout << "data     - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_data_elec        << " GeV" << endl;
      cout << "jec up   - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_jec_up_elec      << " GeV" << endl;
      cout << "jec down - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_jec_down_elec    << " GeV" << endl;
      cout << "cor up   - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_cor_up_elec      << " GeV" << endl;
      cout << "cor down - " << years[year] << setw(12) << " - elec: "    << setw(4) << mean_cor_down_elec    << " GeV" << endl;
      cout << "ttbar    - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_ttbar_combine    << " GeV" << endl;
      cout << "data     - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_data_combine     << " GeV" << endl;
      cout << "jec up   - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_jec_up_combine   << " GeV" << endl;
      cout << "jec down - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_jec_down_combine << " GeV" << endl;
      cout << "cor up   - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_cor_up_combine   << " GeV" << endl;
      cout << "cor down - " << years[year] << setw(12) << " - combine: " << setw(4) << mean_cor_down_combine << " GeV" << endl;
      cout << endl;
      cout << right;

      // Plot ------------------------------------------------------------------

      // draw_plot(vec<Hists> {t, up, do, d}, path, channel, year, mean={0,0,0}, norm="", addition)
      // draw_plot_cor(vec<Hists> {t, up, do, d}, path, channel, year, mean={0,0,0}, norm="", addition)

      // draw_plot(jec_muon, path_muon_jec, means_jec_muon);

      // Peak ------------------------------------------------------------------

      // muon ------------------------------------------------------------------
      vector<double> norm_err_muon = normalize_error(all_data_muon_rebin[year]);
      set_new_bin_error(all_data_muon_norm[year], norm_err_muon);

      vector<TH1F*> ratios_jec_muon;
      ratios_jec_muon.push_back(GetRatio(all_ttbar_muon_norm[year], all_data_muon_norm[year], norm_err_muon, false));
      ratios_jec_muon.push_back(GetRatio(all_jec_up_muon_norm[year], all_data_muon_norm[year], norm_err_muon, false));
      ratios_jec_muon.push_back(GetRatio(all_jec_down_muon_norm[year], all_data_muon_norm[year], norm_err_muon, false));
      ratios_jec_muon.push_back(GetRatio(all_data_muon_norm[year], all_data_muon_norm[year], norm_err_muon, true));

      // Ratio of ttbar had some errors. Copy the bin content solved those. Check again later !!!
      TH1F* saftey_hist_muon  = copy_bin_content(ratios_jec_muon[0], 0, 180, "copy_muon");
      ratios_jec_muon[0] = saftey_hist_muon;

      main_plot_settings(all_ttbar_muon_norm[year],  bin_centers_muon[0], bin_centers_muon[1], kRed, "", xAxis,   "a.u.");
      data_ratio_settings(ratios_jec_muon[3], bin_centers_muon[0], bin_centers_muon[1]);
      add_plot_settings(ratios_jec_muon[0], kRed,  1, 1);
      add_plot_settings(ratios_jec_muon[1], kBlue, 1, 1);
      add_plot_settings(ratios_jec_muon[2], kBlue, 2, 1);

      vector<TH1F*> jec_muon    = {all_ttbar_muon_norm[year],    all_jec_up_muon_norm[year],    all_jec_down_muon_norm[year],    all_data_muon_norm[year]};
      vector<TH1F*> cor_muon    = {all_ttbar_muon_norm[year],    all_cor_up_muon_norm[year],    all_cor_down_muon_norm[year],    all_data_muon_norm[year]};
      vector<TH1F*> all_muon    = {all_ttbar_muon_norm[year],    all_jec_up_muon_norm[year],    all_jec_down_muon_norm[year],    all_cor_up_muon_norm[year],    all_cor_down_muon_norm[year],    all_data_muon_norm[year]};

      // draw_plot_ratio(jec_muon, ratios_jec_muon, path_muon_jec_norm, means_jec_muon);

      // elec ------------------------------------------------------------------
      vector<double> norm_err_elec = normalize_error(all_data_elec_rebin[year]);
      set_new_bin_error(all_data_elec_norm[year], norm_err_elec);

      vector<TH1F*> ratios_jec_elec;
      ratios_jec_elec.push_back(GetRatio(all_ttbar_elec_norm[year], all_data_elec_norm[year], norm_err_elec, false));
      ratios_jec_elec.push_back(GetRatio(all_jec_up_elec_norm[year], all_data_elec_norm[year], norm_err_elec, false));
      ratios_jec_elec.push_back(GetRatio(all_jec_down_elec_norm[year], all_data_elec_norm[year], norm_err_elec, false));
      ratios_jec_elec.push_back(GetRatio(all_data_elec_norm[year], all_data_elec_norm[year], norm_err_elec, true));

      // Ratio of ttbar had some errors. Copy the bin content solved those. Check again later !!!
      TH1F* saftey_hist_elec  = copy_bin_content(ratios_jec_elec[0], 0, 180, "copy_elec");
      ratios_jec_elec[0] = saftey_hist_elec;

      main_plot_settings(all_ttbar_elec_norm[year],  bin_centers_elec[0], bin_centers_elec[1], kRed, "", xAxis,   "a.u.");
      data_ratio_settings(ratios_jec_elec[3], bin_centers_elec[0], bin_centers_elec[1]);
      add_plot_settings(ratios_jec_elec[0], kRed,  1, 1);
      add_plot_settings(ratios_jec_elec[1], kBlue, 1, 1);
      add_plot_settings(ratios_jec_elec[2], kBlue, 2, 1);

      vector<TH1F*> jec_elec    = {all_ttbar_elec_norm[year],    all_jec_up_elec_norm[year],    all_jec_down_elec_norm[year],    all_data_elec_norm[year]};
      vector<TH1F*> cor_elec    = {all_ttbar_elec_norm[year],    all_cor_up_elec_norm[year],    all_cor_down_elec_norm[year],    all_data_elec_norm[year]};
      vector<TH1F*> all_elec    = {all_ttbar_elec_norm[year],    all_jec_up_elec_norm[year],    all_jec_down_elec_norm[year],    all_cor_up_elec_norm[year],    all_cor_down_elec_norm[year],    all_data_elec_norm[year]};

      // draw_plot_ratio(jec_elec, ratios_jec_elec, path_elec_jec_norm, means_jec_elec);

      // combine ------------------------------------------------------------------
      vector<double> norm_err_combine = normalize_error(all_data_combine_rebin[year]);
      set_new_bin_error(all_data_combine_norm[year], norm_err_combine);

      TH1F* ratio_ttbar_combine    = GetRatio(all_ttbar_combine_norm[year],    all_data_combine_norm[year], norm_err_combine, false);
      TH1F* ratio_jec_up_combine   = GetRatio(all_jec_up_combine_norm[year],   all_data_combine_norm[year], norm_err_combine, false);
      TH1F* ratio_jec_down_combine = GetRatio(all_jec_down_combine_norm[year], all_data_combine_norm[year], norm_err_combine, false);
      TH1F* ratio_cor_up_combine   = GetRatio(all_cor_up_combine_norm[year],   all_data_combine_norm[year], norm_err_combine, false);
      TH1F* ratio_cor_down_combine = GetRatio(all_cor_down_combine_norm[year], all_data_combine_norm[year], norm_err_combine, false);
      TH1F* ratio_data_combine     = GetRatio(all_data_combine_norm[year],     all_data_combine_norm[year], norm_err_combine, true);

      // Ratio of ttbar had some errors and was not plotted. Copy the bin content solved this problem. Check again later !!!
      TH1F* saftey_hist_combine  = copy_bin_content(ratio_ttbar_combine, 0, 180, "copy_combine");
      ratio_ttbar_combine = saftey_hist_combine;

      main_plot_settings(all_ttbar_combine_norm[year],  bin_centers_combine[0], bin_centers_combine[1], kRed, "", xAxis,   "a.u.");
      data_ratio_settings(ratio_data_combine,   bin_centers_combine[0], bin_centers_combine[1]);
      add_plot_settings(ratio_ttbar_combine,    kRed,     1, 1);
      add_plot_settings(ratio_jec_up_combine,   kBlue,    1, 1);
      add_plot_settings(ratio_jec_down_combine, kBlue,    2, 1);
      add_plot_settings(ratio_cor_up_combine,   kGreen+2, 1, 1);
      add_plot_settings(ratio_cor_down_combine, kGreen+2, 2, 1);

      vector<TH1F*> ratios_jec_combine = {ratio_ttbar_combine, ratio_jec_up_combine, ratio_jec_down_combine, ratio_data_combine};
      vector<TH1F*> ratios_cor_combine = {ratio_ttbar_combine, ratio_cor_up_combine, ratio_cor_down_combine, ratio_data_combine};
      vector<TH1F*> ratios_all_combine = {ratio_ttbar_combine, ratio_jec_up_combine, ratio_jec_down_combine, ratio_cor_up_combine, ratio_cor_down_combine, ratio_data_combine};

      vector<TH1F*> jec_combine = {all_ttbar_combine_norm[year], all_jec_up_combine_norm[year], all_jec_down_combine_norm[year], all_data_combine_norm[year]};
      vector<TH1F*> cor_combine = {all_ttbar_combine_norm[year], all_cor_up_combine_norm[year], all_cor_down_combine_norm[year], all_data_combine_norm[year]};
      vector<TH1F*> all_combine = {all_ttbar_combine_norm[year], all_jec_up_combine_norm[year], all_jec_down_combine_norm[year], all_cor_up_combine_norm[year], all_cor_down_combine_norm[year], all_data_combine_norm[year]};

      draw_plot_ratio(jec_combine, ratios_jec_combine, path_combine_jec_norm, means_jec_combine);
      draw_plot_ratio(cor_combine, ratios_cor_combine, path_combine_cor_norm, means_cor_combine);
      draw_plot_ratio(all_combine, ratios_all_combine, path_combine_all_norm, means_all_combine);

    } // year
  } // pt bin
} // main







// if(forTalk){
//   TString cmstext = "CMS";
//   TLatex *text2 = new TLatex(3.5, 24, cmstext);
//   text2->SetNDC();
//   text2->SetTextAlign(13);
//   text2->SetX(0.6);
//   text2->SetTextFont(62);
//   text2->SetTextSize(0.05);
//   text2->SetY(0.84);
//   text2->Draw();
//
//   TString preltext = "Work in Progress";
//   TLatex *text3 = new TLatex(3.5, 24, preltext);
//   text3->SetNDC();
//   text3->SetTextAlign(13);
//   text3->SetX(0.6);
//   text3->SetTextFont(52);
//   text3->SetTextSize(0.035);
//   text3->SetY(0.78);
//   text3->Draw();
// }
