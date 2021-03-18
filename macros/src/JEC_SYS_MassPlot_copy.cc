#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

bool forTalk = true;

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
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.5);
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
  int Nbins = h1->GetSize() - 2;
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
      // if(isData) ratio->SetBinError(i, data_bin_error[i-1]);
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
  TLegend *leg = new TLegend(0.35,0.15,0.6,0.35);

  TCanvas *A = new TCanvas("A"+path[0]+path[4], "A", 600, 600);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  hists[0]->Draw("HIST");
  hists[1]->Draw("SAME HIST");
  hists[2]->Draw("SAME HIST");
  hists[3]->Draw("SAME P");

  TString nom  = "Nominal  <m_{jet}>="+dtos(mean[0], 2)+" GeV";
  TString up   = path[2]+" up   <m_{jet}>="+dtos(mean[1], 2)+" GeV";
  TString down = path[2]+" down <m_{jet}>="+dtos(mean[2], 2)+" GeV";

  leg->AddEntry(hists[0],nom ,"l");
  leg->AddEntry(hists[1],up  ,"l");
  leg->AddEntry(hists[2],down,"l");
  leg->AddEntry(hists[3],"Data","pl");

  leg->SetTextSize(0.03);
  leg->Draw();

  gPad->RedrawAxis();

  A->SaveAs(path[0]+"/"+path[1]+"/Wjet_mass_"+path[2]+"_"+path[3]+"_"+path[4]+path[5]+".pdf");

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
void draw_plot_ratio(vector<TH1F*> hists, vector<TH1F*> ratios, vector<TString> path, vector<double> mean)
{
  // path    = {save_path, channel, var, addition, year, norm}
  // hists   = {nom, up, down, data}
  // rations = {nom/data, up/data, down/data, data/data}

  TLegend *leg = new TLegend(0.35,0.15,0.6,0.35);
  TCanvas *A = new TCanvas("A"+path[0]+path[4], "A", 600, 600);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  hists[0]->Draw("HIST");
  hists[1]->Draw("SAME HIST");
  hists[2]->Draw("SAME HIST");
  hists[3]->Draw("SAME P");

  TString nom  = "Nominal  <m_{jet}>="+dtos(mean[0], 2)+" GeV";
  TString up   = path[2]+" up   <m_{jet}>="+dtos(mean[1], 2)+" GeV";
  TString down = path[2]+" down <m_{jet}>="+dtos(mean[2], 2)+" GeV";

  leg->AddEntry(hists[0],nom ,"l");
  leg->AddEntry(hists[1],up  ,"l");
  leg->AddEntry(hists[2],down,"l");
  leg->AddEntry(hists[3],"Data","pl");

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

  ratios[3]->Draw("P");
  ratios[0]->Draw("HIST SAME");
  ratios[1]->Draw("HIST SAME");
  ratios[2]->Draw("HIST SAME");


  gPad->RedrawAxis();

  A->SaveAs(path[0]+"/"+path[1]+"/Wjet_mass_"+path[2]+path[3]+"_"+path[4]+path[5]+".pdf");

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
// ------------------------------------------ MAIN -------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]){
  bool debug = true;
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

    vector<TH1F*> allyears_bkg,  allyears_data,  allyears_ttbar;
    vector<TH1F*> allyears_JECup, allyears_JECdown, allyears_XCup, allyears_XCdown;


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
      peak_bins_combine = bins_upper_limit(all_data_combine_rebin[year], 75);

      // Bin Center
      vector<double> bin_centers_muon = bin_center_upper_limit(all_data_muon_rebin[year], peak_bins_muon, bin_width);

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

      // ttbar -----------------------------------------------------------------

      TString xAxis      = "m_{jet} [GeV]";
      TString yAxis_norm = "#Delta N/N";

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

      add_plot_settings(all_cor_up_muon_rebin[year], kBlue);
      add_plot_settings(all_cor_up_elec_rebin[year], kBlue);
      add_plot_settings(all_cor_up_combine_rebin[year], kBlue);
      add_plot_settings(all_cor_up_muon_norm[year], kBlue);
      add_plot_settings(all_cor_up_elec_norm[year], kBlue);
      add_plot_settings(all_cor_up_combine_norm[year], kBlue);

      add_plot_settings(all_cor_down_muon[year], kBlue, 7);
      add_plot_settings(all_cor_down_elec[year], kBlue, 7);
      add_plot_settings(all_cor_down_combine[year], kBlue, 7);
      add_plot_settings(all_cor_down_muon_norm[year], kBlue, 7);
      add_plot_settings(all_cor_down_elec_norm[year], kBlue, 7);
      add_plot_settings(all_cor_down_combine_norm[year], kBlue, 7);

      // #################################################################################################
      // Mass Plots ######################################################################################
      if(debug) cout << "Mass Plots" << endl;

      // Mean ------------------------------------------------------------------

      double mean_ttbar_muon    = trunc_mean(all_ttbar_muon_norm[year],    bin_centers_muon[0], bin_centers_muon[1]);
      double mean_jec_up_muon   = trunc_mean(all_jec_up_muon_norm[year],   bin_centers_muon[0], bin_centers_muon[1]);
      double mean_jec_down_muon = trunc_mean(all_jec_down_muon_norm[year], bin_centers_muon[0], bin_centers_muon[1]);
      double mean_cor_up_muon   = trunc_mean(all_cor_up_muon_norm[year],   bin_centers_muon[0], bin_centers_muon[1]);
      double mean_cor_down_muon = trunc_mean(all_cor_down_muon_norm[year], bin_centers_muon[0], bin_centers_muon[1]);
      vector<double> means_jec  = {mean_ttbar_muon, mean_jec_up_muon, mean_jec_down_muon};
      vector<double> means_cor  = {mean_ttbar_muon, mean_cor_up_muon, mean_cor_down_muon};

      cout << left; // for setw()
      cout << "ttbar    - " << years[year] << setw(12) << " - muon: " << setw(8) << mean_ttbar_muon    << endl;
      cout << "jec up   - " << years[year] << setw(12) << " - muon: " << setw(8) << mean_jec_up_muon   << endl;
      cout << "jec down - " << years[year] << setw(12) << " - muon: " << setw(8) << mean_jec_down_muon << endl;
      cout << "cor up   - " << years[year] << setw(12) << " - muon: " << setw(8) << mean_cor_up_muon   << endl;
      cout << "cor down - " << years[year] << setw(12) << " - muon: " << setw(8) << mean_cor_down_muon << endl;
      cout << endl;
      cout << right;

      // Plot ------------------------------------------------------------------

      vector<TH1F*> jec_muon = {all_ttbar_muon_norm[year], all_jec_up_muon_norm[year], all_jec_down_muon_norm[year], all_data_muon_norm[year]};
      vector<TH1F*> cor_muon = {all_ttbar_muon_norm[year], all_cor_up_muon_norm[year], all_cor_down_muon_norm[year], all_data_muon_norm[year]};

      // draw_plot(vec<Hists> {t, up, do, d}, path, channel, year, mean={0,0,0}, norm="", addition)
      // draw_plot_cor(vec<Hists> {t, up, do, d}, path, channel, year, mean={0,0,0}, norm="", addition)

      // draw_plot(jec_muon, path_muon_jec, means_jec);

      // Peak ------------------------------------------------------------------

      vector<double> norm_err = normalize_error(all_data_muon_rebin[year]);

      vector<TH1F*> ratios_jec_muon;
      ratios_jec_muon.push_back(GetRatio(all_ttbar_muon_norm[year], all_data_muon_norm[year], norm_err, false));
      ratios_jec_muon.push_back(GetRatio(all_jec_up_muon_norm[year], all_data_muon_norm[year], norm_err, false));
      ratios_jec_muon.push_back(GetRatio(all_jec_down_muon_norm[year], all_data_muon_norm[year], norm_err, false));
      ratios_jec_muon.push_back(GetRatio(all_data_muon_norm[year], all_data_muon_norm[year], norm_err, true));

      // for(int i=1; i<ratios_jec_muon[0]->GetNbinsX()+1; i++)
      // {
      //   cout << ratios_jec_muon[0]->GetBinContent(i) << " ---------- " << i << endl;
      // }
      cout << ratios_jec_muon[0]->GetMaximum() << endl;
      cout << ratios_jec_muon[1]->GetMaximum() << endl;
      cout << ratios_jec_muon[2]->GetMaximum() << endl;
      cout << ratios_jec_muon[3]->GetMaximum() << endl;

      main_plot_settings(all_ttbar_muon_norm[year],  bin_centers_muon[0], bin_centers_muon[1], kRed, "", xAxis,   "a.u.");
      data_ratio_settings(ratios_jec_muon[3], bin_centers_muon[0], bin_centers_muon[1]);
      add_plot_settings(ratios_jec_muon[0], kRed,  kSolid,  1);
      add_plot_settings(ratios_jec_muon[1], kBlue, kSolid,  1);
      add_plot_settings(ratios_jec_muon[2], kBlue, kDashed, 1);

      // plot_single_histogram(TH1F* hist, TString title, TString xAxis, int x_max, int color, TString save_path){
      // TEST
      plot_single_histogram(ratios_jec_muon[0], "", "", 180, kBlack, save_path+"/test.pdf");
      plot_single_histogram(ratios_jec_muon[2], "", "", 180, kBlack, save_path+"/test2.pdf");

      draw_plot_ratio(jec_muon, ratios_jec_muon, path_muon_jec, means_jec);

      // TString norm_str="";
      // for(unsigned int norm=0; norm<ttbar_.size(); norm++){
      //   if(norm==1) norm_str="_norm";
      //   for(int correction=0; correction<2; correction++){
      //
      //     corrections_up_[norm][correction]->SetLineWidth(2);
      //     corrections_up_[norm][correction]->SetLineColor(kBlue);
      //
      //     corrections_down_[norm][correction]->SetLineStyle(2);
      //     corrections_down_[norm][correction]->SetLineWidth(2);
      //     corrections_down_[norm][correction]->SetLineColor(kBlue);
      //
      //     ttbar_[norm]->GetYaxis()->SetTitleOffset(1.3);
      //
      //     TCanvas *A = new TCanvas("A", "A", 600, 600);
      //     gPad->SetLeftMargin(0.15);
      //     gPad->SetBottomMargin(0.12);
      //     ttbar_[norm]->Draw("HIST");
      //     corrections_down_[norm][correction]->Draw("SAME HIST");
      //     corrections_up_[norm][correction]->Draw("SAME HIST");
      //     data_[norm]->Draw("SAME P");
      //     // if(norm==0 && bin_width==1 && isAll) line->Draw("SAME");
      //     leg = new TLegend(0.15,0.65,0.35,0.85);
      //     if(bin_width==1) leg = new TLegend(0.45,0.15,0.6,0.35);;
      //     leg->AddEntry(ttbar_[norm],"Nominal","l");
      //     if(correction==0){
      //       leg->AddEntry(corrections_up_[norm][correction],"JEC up","l");
      //       leg->AddEntry(corrections_down_[norm][correction],"JEC down","l");
      //     }
      //     if(correction==1){
      //       leg->AddEntry(corrections_up_[norm][correction],"XCone up","l");
      //       leg->AddEntry(corrections_down_[norm][correction],"XCone down","l");
      //     }
      //     leg->AddEntry(data_[norm],"Data","pl");
      //     leg->SetTextSize(0.02);
      //     leg->Draw();
      //     gPad->RedrawAxis();
      //
      //
      //
      //     if(correction==0) A->SaveAs(save_path+"/Wjet_mass_sensitivity_JEC"+addition+norm_str+".pdf");
      //     if(correction==1) A->SaveAs(save_path+"/Wjet_mass_sensitivity_XCone"+addition+norm_str+".pdf");
      //     delete A;
      //     leg->Clear();
      //   }
      // }

    } // year
  } // pt bin
} // main




// ttbar_[norm]->SetTitle("");
// if     (ptbin==1&&bin_width==5) ttbar_[norm]->GetXaxis()->SetRangeUser(20, 180);
// else if(ptbin==1&&bin_width==6) ttbar_[norm]->GetXaxis()->SetRangeUser(20, 180);
// else                            ttbar_[norm]->GetXaxis()->SetRangeUser(0, 180);

// TLine *line = new TLine(0, Limit, 180, Limit);
// line->SetLineColor(kGray);
// line->SetLineWidth(1);





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
