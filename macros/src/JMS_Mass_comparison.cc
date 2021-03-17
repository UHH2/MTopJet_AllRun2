#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/CreatHists.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

// -------------------------------------------------------------------------------------------------------
void main_plot_settings(TH1F* hist, int x_max, int color, TString title, TString xAxis, TString yAxis, TString save_path)
{
  hist->SetTitle(title);
  hist->GetXaxis()->SetRangeUser(100, x_max);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
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
void add_plot_settings(TH1F* hist, int color=1, int style=kSolid, int width=2)
{
  hist->SetLineWidth(width);
  hist->SetLineStyle(style);
  hist->SetLineColor(color);
}

// -------------------------------------------------------------------------------------------------------
void data_plot_settings(TH1F* hist)
{
  hist->SetMarkerStyle(8);  // data hist style
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);
}

// -------------------------------------------------------------------------------------------------------
void draw_plot_jms(TH1F* h1, TH1F* h2, TH1F* h3, const TString path, const TString year, vector<double> mean, TString norm="")
{
  TLegend *leg = new TLegend(0.55,0.65,0.8,0.85);
  TCanvas *A = new TCanvas(path+year+norm, "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  h1->Draw("HIST");
  h2->Draw("SAME HIST");
  h3->Draw("SAME HIST");

  TString best = "Best fit; mean="+dtos(mean[0], 2)+" GeV";
  TString uu   = "JMS uu; mean="+dtos(mean[1], 2)+" GeV";
  TString dd   = "JMS dd; mean="+dtos(mean[2], 2)+" GeV";

  leg->AddEntry(h1,best,"l");
  leg->AddEntry(h2,uu,"l");
  leg->AddEntry(h3,dd,"l");

  leg->SetTextSize(0.02);
  leg->Draw();

  gPad->RedrawAxis();
  A->SaveAs(path+"/jms_"+year+norm+".pdf");

}

// -------------------------------------------------------------------------------------------------------
void draw_plot_comparison(TH1F* h1, TH1F* h2, TH1F* data, const TString path, const TString year, vector<double> mean={0, 0})
{
  TLegend *leg = new TLegend(0.5,0.65,0.75,0.85);

  TCanvas *A = new TCanvas("A"+path+year, "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  h1->Draw("HIST");
  h2->Draw("SAME HIST");

  TString jms = "JMS; mean="+dtos(mean[0], 2)+" GeV";
  TString old = "Old; mean="+dtos(mean[1], 2)+" GeV";

  leg->AddEntry(h1,jms,"l");
  leg->AddEntry(h2,old,"l");

  leg->SetTextSize(0.03);
  leg->Draw();

  gPad->RedrawAxis();
  A->SaveAs(path+"/comparison_jms_old_"+year+".pdf");

  data->Draw("SAME P");
  leg->AddEntry(data,"Data","pl");
  A->SaveAs(path+"/comparison_jms_old_data_"+year+".pdf");

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
double trunc_mean(TH1F* hist)
{
  TH1F* hist_trunc = (TH1F*) hist->Clone();
  // cout << setw(10) << "old" << "   " << setw(10) << "new" << endl;
  for(int bin=1; bin<hist->GetNbinsX()+1; bin++)
  {
    bool bin_ = (hist->GetBinCenter(bin)<120||hist->GetBinCenter(bin)>240);
    if(bin_) hist_trunc->SetBinContent(bin, 0);
    else     hist_trunc->SetBinContent(bin, hist->GetBinContent(bin));
    // cout << setw(10) << hist->GetBinContent(bin) << "   " << setw(10) << hist_trunc->GetBinContent(bin) << endl;
  }
  double mean_trunc = hist_trunc->GetMean();
  double mean_old   = hist->GetMean();
  // cout << setw(10) << mean_old << "   " << setw(10) << mean_trunc << endl;
  return mean_trunc;
}

// -------------------------------------------------------------------------------------------------------
// ------------------------------------------ MAIN -------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  /*
  Show the behaviour af the different mass scale variations uu, ud, du an dd.
  Show the difference between the old Mass Jet distribution and the new one
  using the BestFit method.
  */
  bool debug = true;
  TString save_path = get_save_path();

  print_seperater();

  if(argc != 1){
    cout << "\n" << "Usage: ./JEC_Mass_comparison\n";
    cout << "All years and channels are created simultanously." << endl;
    return 0;
  }

  TString save_path_general = save_path+"/Plots"; // CHANGE_PT
  save_path_general = creat_folder_and_path(save_path_general, "JMS");
  creat_folder(save_path_general+"/muon");
  creat_folder(save_path_general+"/elec");
  creat_folder(save_path_general+"/combine");

  vector<int> index_year   = {0,1,2,3};
  vector<TString> years    = {"2016","2017","2018","combine"};
  vector<TString> channels = {"muon","elec","combine"};

  // -------------------------------------------------------------------------------------------------------
  // ------------------------------------------ MAIN -------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------
  cout << '\n';
  TString hist; TString hist_class = "JetMassScaleHists/";
  hist = hist_class+"hadjet_jms_mass";

  TString hist_old; TString hist_class_old = "XCone_cor/";
  hist_old = hist_class_old+"M_jet1";

  // #################################################################################################
  // Get Hists #######################################################################################
  if(debug) cout << "Get Hists" << endl;

  vector<TH1F*> h_data_muon          = get_all_hists(data_muon, hist);
  vector<TH1F*> h_data_muon_old      = get_all_hists(data_muon, hist_old);
  vector<TH1F*> h_data_elec          = get_all_hists(data_elec, hist);
  vector<TH1F*> h_data_elec_old      = get_all_hists(data_elec, hist_old);
  vector<TH1F*> h_data_combine       = combine_channels(h_data_muon, h_data_elec);
  vector<TH1F*> h_data_combine_old   = combine_channels(h_data_muon_old, h_data_elec_old);

  vector<TH1F*> h_ttbar_muon         = get_all_hists(ttbar_muon, hist);
  vector<TH1F*> h_ttbar_muon_old     = get_all_hists(ttbar_muon, hist_old);
  vector<TH1F*> h_ttbar_elec         = get_all_hists(ttbar_elec, hist);
  vector<TH1F*> h_ttbar_elec_old     = get_all_hists(ttbar_elec, hist_old);
  vector<TH1F*> h_ttbar_combine      = combine_channels(h_ttbar_muon, h_ttbar_elec);
  vector<TH1F*> h_ttbar_combine_old  = combine_channels(h_ttbar_muon_old, h_ttbar_elec_old);

  vector<TH1F*> h_st_muon            = get_all_hists(st_muon, hist);
  vector<TH1F*> h_st_muon_old        = get_all_hists(st_muon, hist_old);
  vector<TH1F*> h_st_elec            = get_all_hists(st_elec, hist);
  vector<TH1F*> h_st_elec_old        = get_all_hists(st_elec, hist_old);
  vector<TH1F*> h_st_combine         = combine_channels(h_st_muon, h_st_elec);
  vector<TH1F*> h_st_combine_old     = combine_channels(h_st_muon_old, h_st_elec_old);

  vector<TH1F*> h_wjets_muon         = get_all_hists(wjets_muon, hist);
  vector<TH1F*> h_wjets_muon_old     = get_all_hists(wjets_muon, hist_old);
  vector<TH1F*> h_wjets_elec         = get_all_hists(wjets_elec, hist);
  vector<TH1F*> h_wjets_elec_old     = get_all_hists(wjets_elec, hist_old);
  vector<TH1F*> h_wjets_combine      = combine_channels(h_wjets_muon, h_wjets_elec);
  vector<TH1F*> h_wjets_combine_old  = combine_channels(h_wjets_muon_old, h_wjets_elec_old);

  vector<TH1F*> h_other_muon         = get_all_hists(other_muon, hist);
  vector<TH1F*> h_other_muon_old     = get_all_hists(other_muon, hist_old);
  vector<TH1F*> h_other_elec         = get_all_hists(other_elec, hist);
  vector<TH1F*> h_other_elec_old     = get_all_hists(other_elec, hist_old);
  vector<TH1F*> h_other_combine      = combine_channels(h_other_muon, h_other_elec);
  vector<TH1F*> h_other_combine_old  = combine_channels(h_other_muon_old, h_other_elec_old);

  // JMS -----------------------------------------------------------------------
  // if(debug) cout << "Get Hists- JMS_uu" << endl; // WORKS
  vector<TH1F*> h_jms_uu_muon        = get_all_hists(jms_uu_muon, hist);
  vector<TH1F*> h_jms_uu_elec        = get_all_hists(jms_uu_elec, hist);
  vector<TH1F*> h_jms_uu_combine     = combine_channels(h_jms_uu_muon, h_jms_uu_elec);

  // if(debug) cout << "Get Hists- JMS_ud" << endl;
  vector<TH1F*> h_jms_ud_muon        = get_all_hists(jms_ud_muon, hist);
  vector<TH1F*> h_jms_ud_elec        = get_all_hists(jms_ud_elec, hist);
  vector<TH1F*> h_jms_ud_combine     = combine_channels(h_jms_ud_muon, h_jms_ud_elec);

  // if(debug) cout << "Get Hists- JMS_du" << endl; // WORKS
  vector<TH1F*> h_jms_du_muon        = get_all_hists(jms_du_muon, hist);
  vector<TH1F*> h_jms_du_elec        = get_all_hists(jms_du_elec, hist);
  vector<TH1F*> h_jms_du_combine     = combine_channels(h_jms_du_muon, h_jms_du_elec);

  // if(debug) cout << "Get Hists- JMS_dd" << endl;
  vector<TH1F*> h_jms_dd_muon        = get_all_hists(jms_dd_muon, hist);
  vector<TH1F*> h_jms_dd_elec        = get_all_hists(jms_dd_elec, hist);
  vector<TH1F*> h_jms_dd_combine     = combine_channels(h_jms_dd_muon, h_jms_dd_elec);

  // #################################################################################################
  // Combine Hists ###################################################################################
  if(debug) cout << "Combine Hists" << endl;

  // Data ----------------------------------------------------------------------

  // dummy: in order to keep code clean and faster to write

  vector<TH1F*> all_data_muon             = h_data_muon;        // dummy
  vector<TH1F*> all_data_elec             = h_data_elec;        // dummy
  vector<TH1F*> all_data_combine          = h_data_combine;     // dummy

  vector<TH1F*> all_data_muon_old         = h_data_muon_old;    // dummy
  vector<TH1F*> all_data_elec_old         = h_data_elec_old;    // dummy
  vector<TH1F*> all_data_combine_old      = h_data_combine_old; // dummy

  vector<TH1F*> all_data_muon_norm        = normalize(all_data_muon);
  vector<TH1F*> all_data_elec_norm        = normalize(all_data_elec);
  vector<TH1F*> all_data_combine_norm     = normalize(all_data_combine);

  vector<TH1F*> all_data_muon_norm_old    = normalize(all_data_muon_old);
  vector<TH1F*> all_data_elec_norm_old    = normalize(all_data_elec_old);
  vector<TH1F*> all_data_combine_norm_old = normalize(all_data_combine_old);

  // Bkg -----------------------------------------------------------------------

  vector<TH1F*> h_bkg_muon        = AddHists(h_st_muon, h_wjets_muon, h_other_muon, 1);
  vector<TH1F*> h_bkg_elec        = AddHists(h_st_elec, h_wjets_elec, h_other_elec, 1);
  vector<TH1F*> h_bkg_combine     = AddHists(h_st_combine, h_wjets_combine, h_other_combine, 1);

  vector<TH1F*> h_bkg_muon_old    = AddHists(h_st_muon_old, h_wjets_muon_old, h_other_muon_old, 1);
  vector<TH1F*> h_bkg_elec_old    = AddHists(h_st_elec_old, h_wjets_elec_old, h_other_elec_old, 1);
  vector<TH1F*> h_bkg_combine_old = AddHists(h_st_combine_old, h_wjets_combine_old, h_other_combine_old, 1);

  // TTbar ---------------------------------------------------------------------
  if(debug) cout << "TTbar" << endl;

  vector<TH1F*> all_ttbar_muon             = AddHists(h_bkg_muon, h_ttbar_muon, 1);
  vector<TH1F*> all_ttbar_elec             = AddHists(h_bkg_elec, h_ttbar_elec, 1);
  vector<TH1F*> all_ttbar_combine          = AddHists(h_bkg_combine, h_ttbar_combine, 1);
  vector<TH1F*> all_ttbar_muon_old         = AddHists(h_bkg_muon_old, h_ttbar_muon_old, 1);
  vector<TH1F*> all_ttbar_elec_old         = AddHists(h_bkg_elec_old, h_ttbar_elec_old, 1);
  vector<TH1F*> all_ttbar_combine_old      = AddHists(h_bkg_combine_old, h_ttbar_combine_old, 1);

  vector<TH1F*> all_ttbar_muon_norm        = normalize(all_ttbar_muon);
  vector<TH1F*> all_ttbar_elec_norm        = normalize(all_ttbar_elec);
  vector<TH1F*> all_ttbar_combine_norm     = normalize(all_ttbar_combine);
  vector<TH1F*> all_ttbar_muon_norm_old    = normalize(all_ttbar_muon_old);
  vector<TH1F*> all_ttbar_elec_norm_old    = normalize(all_ttbar_elec_old);
  vector<TH1F*> all_ttbar_combine_norm_old = normalize(all_ttbar_combine_old);

  // JMS -----------------------------------------------------------------------
  if(debug) cout << "JMS" << endl;

  vector<TH1F*> all_jms_uu_muon         = AddHists(h_bkg_muon, h_jms_uu_muon, 1);
  vector<TH1F*> all_jms_uu_elec         = AddHists(h_bkg_elec, h_jms_uu_elec, 1);
  vector<TH1F*> all_jms_uu_combine      = AddHists(h_bkg_combine, h_jms_uu_combine, 1);

  vector<TH1F*> all_jms_ud_muon         = AddHists(h_bkg_muon, h_jms_ud_muon, 1);
  vector<TH1F*> all_jms_ud_elec         = AddHists(h_bkg_elec, h_jms_ud_elec, 1);
  vector<TH1F*> all_jms_ud_combine      = AddHists(h_bkg_combine, h_jms_ud_combine, 1);

  vector<TH1F*> all_jms_du_muon         = AddHists(h_bkg_muon, h_jms_du_muon, 1);
  vector<TH1F*> all_jms_du_elec         = AddHists(h_bkg_elec, h_jms_du_elec, 1);
  vector<TH1F*> all_jms_du_combine      = AddHists(h_bkg_combine, h_jms_du_combine, 1);

  vector<TH1F*> all_jms_dd_muon         = AddHists(h_bkg_muon, h_jms_dd_muon, 1);
  vector<TH1F*> all_jms_dd_elec         = AddHists(h_bkg_elec, h_jms_dd_elec, 1);
  vector<TH1F*> all_jms_dd_combine      = AddHists(h_bkg_combine, h_jms_dd_combine, 1);

  vector<TH1F*> all_jms_uu_muon_norm    = normalize(all_jms_uu_muon);
  vector<TH1F*> all_jms_uu_elec_norm    = normalize(all_jms_uu_elec);
  vector<TH1F*> all_jms_uu_combine_norm = normalize(all_jms_uu_combine);

  vector<TH1F*> all_jms_ud_muon_norm    = normalize(all_jms_ud_muon);
  vector<TH1F*> all_jms_ud_elec_norm    = normalize(all_jms_ud_elec);
  vector<TH1F*> all_jms_ud_combine_norm = normalize(all_jms_ud_combine);

  vector<TH1F*> all_jms_du_muon_norm    = normalize(all_jms_du_muon);
  vector<TH1F*> all_jms_du_elec_norm    = normalize(all_jms_du_elec);
  vector<TH1F*> all_jms_du_combine_norm = normalize(all_jms_du_combine);

  vector<TH1F*> all_jms_dd_muon_norm    = normalize(all_jms_dd_muon);
  vector<TH1F*> all_jms_dd_elec_norm    = normalize(all_jms_dd_elec);
  vector<TH1F*> all_jms_dd_combine_norm = normalize(all_jms_dd_combine);

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

    // #################################################################################################
    // Settings ########################################################################################

    // ttbar -------------------------------------------------------------------

    // void main_plot_settings(TH1F* hist, int x_max, int color, TString title, TString xAxis, TString yAxis, TString save_path)
    TString xAxis      = "m_{jet} [GeV]";
    TString yAxis_norm = "#Delta N/N";
    main_plot_settings(all_ttbar_muon[year], 300, kRed, "", xAxis, "Events", save_path_general);
    main_plot_settings(all_ttbar_elec[year], 300, kRed, "", xAxis, "Events", save_path_general);
    main_plot_settings(all_ttbar_combine[year], 300, kRed, "", xAxis, "Events", save_path_general);
    main_plot_settings(all_ttbar_muon_norm[year], 300, kRed, "", xAxis, "a.u.", save_path_general);
    main_plot_settings(all_ttbar_elec_norm[year], 300, kRed, "", xAxis, "a.u.", save_path_general);
    main_plot_settings(all_ttbar_combine_norm[year], 300, kRed, "", xAxis, "a.u.", save_path_general);

    add_plot_settings(all_ttbar_muon_old[year], kBlue);
    add_plot_settings(all_ttbar_elec_old[year], kBlue);
    add_plot_settings(all_ttbar_combine_old[year], kBlue);
    add_plot_settings(all_ttbar_muon_norm_old[year], kBlue);
    add_plot_settings(all_ttbar_elec_norm_old[year], kBlue);
    add_plot_settings(all_ttbar_combine_norm_old[year], kBlue);

    // data --------------------------------------------------------------------

    // data_plot_settings(TH1F* hist)
    data_plot_settings(all_data_muon[year]);
    data_plot_settings(all_data_elec[year]);
    data_plot_settings(all_data_combine[year]);
    data_plot_settings(all_data_muon_norm[year]);
    data_plot_settings(all_data_elec_norm[year]);
    data_plot_settings(all_data_combine_norm[year]);

    data_plot_settings(all_data_muon_old[year]);
    data_plot_settings(all_data_elec_old[year]);
    data_plot_settings(all_data_combine_old[year]);
    data_plot_settings(all_data_muon_norm_old[year]);
    data_plot_settings(all_data_elec_norm_old[year]);
    data_plot_settings(all_data_combine_norm_old[year]);

    // jms ---------------------------------------------------------------------

    // add_plot_settings(TH1F* hist, int color=1, int style=kSolid, int width=2)
    add_plot_settings(all_jms_uu_muon[year], kBlue);
    add_plot_settings(all_jms_uu_elec[year], kBlue);
    add_plot_settings(all_jms_uu_combine[year], kBlue);
    add_plot_settings(all_jms_uu_muon_norm[year], kBlue);
    add_plot_settings(all_jms_uu_elec_norm[year], kBlue);
    add_plot_settings(all_jms_uu_combine_norm[year], kBlue);

    add_plot_settings(all_jms_dd_muon[year], kGreen+2);
    add_plot_settings(all_jms_dd_elec[year], kGreen+2);
    add_plot_settings(all_jms_dd_combine[year], kGreen+2);
    add_plot_settings(all_jms_dd_muon_norm[year], kGreen+2);
    add_plot_settings(all_jms_dd_elec_norm[year], kGreen+2);
    add_plot_settings(all_jms_dd_combine_norm[year], kGreen+2);

    add_plot_settings(all_jms_ud_muon[year], kYellow+2, 7);
    add_plot_settings(all_jms_ud_elec[year], kYellow+2, 7);
    add_plot_settings(all_jms_ud_combine[year], kYellow+2, 7);
    add_plot_settings(all_jms_ud_muon_norm[year], kYellow+2, 7);
    add_plot_settings(all_jms_ud_elec_norm[year], kYellow+2, 7);
    add_plot_settings(all_jms_ud_combine_norm[year], kYellow+2, 7);

    add_plot_settings(all_jms_du_muon[year], kGray+2, 7);
    add_plot_settings(all_jms_du_elec[year], kGray+2, 7);
    add_plot_settings(all_jms_du_combine[year], kGray+2, 7);
    add_plot_settings(all_jms_du_muon_norm[year], kGray+2, 7);
    add_plot_settings(all_jms_du_elec_norm[year], kGray+2, 7);
    add_plot_settings(all_jms_du_combine_norm[year], kGray+2, 7);

    // #################################################################################################
    // Plots ###########################################################################################
    if(debug) cout << "Mass Plots" << endl;

    TString path_muon    = save_path_general+"/muon";
    TString path_elec    = save_path_general+"/elec";
    TString path_combine = save_path_general+"/combine";

    // mean --------------------------------------------------------------------

    double mean_data_muon           = trunc_mean(all_data_muon[year]); // all_data_muon[year]->GetMean();
    double mean_ttbar_muon          = trunc_mean(all_ttbar_muon[year]); // all_ttbar_muon[year]->GetMean();
    double mean_uu_muon             = trunc_mean(all_jms_uu_muon[year]); // all_jms_uu_muon[year]->GetMean();
    double mean_dd_muon             = trunc_mean(all_jms_dd_muon[year]); // all_jms_dd_muon[year]->GetMean();
    double mean_ud_muon             = trunc_mean(all_jms_ud_muon[year]); // all_jms_ud_muon[year]->GetMean();
    double mean_du_muon             = trunc_mean(all_jms_du_muon[year]); // all_jms_du_muon[year]->GetMean();

    double mean_data_elec           = trunc_mean(all_data_elec[year]); // all_data_elec[year]->GetMean();
    double mean_ttbar_elec          = trunc_mean(all_ttbar_elec[year]); // all_ttbar_elec[year]->GetMean();
    double mean_uu_elec             = trunc_mean(all_jms_uu_elec[year]); // all_jms_uu_elec[year]->GetMean();
    double mean_dd_elec             = trunc_mean(all_jms_dd_elec[year]); // all_jms_dd_elec[year]->GetMean();
    double mean_ud_elec             = trunc_mean(all_jms_ud_elec[year]); // all_jms_ud_elec[year]->GetMean();
    double mean_du_elec             = trunc_mean(all_jms_du_elec[year]); // all_jms_du_elec[year]->GetMean();

    double mean_data_combine        = trunc_mean(all_data_combine[year]); // all_data_combine[year]->GetMean();
    double mean_ttbar_combine       = trunc_mean(all_ttbar_combine[year]); // all_ttbar_combine[year]->GetMean();
    double mean_uu_combine          = trunc_mean(all_jms_uu_combine[year]); // all_jms_uu_combine[year]->GetMean();
    double mean_dd_combine          = trunc_mean(all_jms_dd_combine[year]); // all_jms_dd_combine[year]->GetMean();
    double mean_ud_combine          = trunc_mean(all_jms_ud_combine[year]); // all_jms_ud_combine[year]->GetMean();
    double mean_du_combine          = trunc_mean(all_jms_du_combine[year]); // all_jms_du_combine[year]->GetMean();

    vector<double> mean_muon        = {mean_ttbar_muon, mean_uu_muon, mean_dd_muon};
    vector<double> mean_elec        = {mean_ttbar_elec, mean_uu_elec, mean_dd_elec};
    vector<double> mean_combine     = {mean_ttbar_combine, mean_uu_combine, mean_dd_combine};

    double mean_ttbar_muon_old      = trunc_mean(all_ttbar_muon_old[year]); // all_ttbar_muon_old[year]->GetMean();
    double mean_ttbar_elec_old      = trunc_mean(all_ttbar_elec_old[year]); // all_ttbar_elec_old[year]->GetMean();
    double mean_ttbar_combine_old   = trunc_mean(all_ttbar_combine_old[year]); // all_ttbar_combine_old[year]->GetMean();

    vector<double> mean_muon_old    = {mean_ttbar_muon, mean_ttbar_muon_old};
    vector<double> mean_elec_old    = {mean_ttbar_elec, mean_ttbar_elec_old};
    vector<double> mean_combine_old = {mean_ttbar_combine, mean_ttbar_combine_old};

    cout << "data  - " << years[year] << " - muon:    " << mean_data_muon  << endl;
    cout << "jmsuu - " << years[year] << " - muon:    " << mean_uu_muon    << endl;
    cout << "ttbar - " << years[year] << " - muon:    " << mean_ttbar_muon;
    cout << " | (old: " << mean_ttbar_muon_old << ")" << endl;
    cout << "jmsdd - " << years[year] << " - muon:    " << mean_dd_muon    << endl;
    cout << endl;
    cout << "data  - " << years[year] << " - elec:    " << mean_data_elec  << endl;
    cout << "jmsuu - " << years[year] << " - elec:    " << mean_uu_elec    << endl;
    cout << "ttbar - " << years[year] << " - elec:    " << mean_ttbar_elec;
    cout << " | (old: " << mean_ttbar_elec_old << ")" << endl;
    cout << "jmsdd - " << years[year] << " - elec:    " << mean_dd_elec    << endl;
    cout << endl;
    cout << "data  - " << years[year] << " - combine: " << mean_data_combine  << endl;
    cout << "jmsuu - " << years[year] << " - combine: " << mean_uu_combine    << endl;
    cout << "ttbar - " << years[year] << " - combine: " << mean_ttbar_combine;
    cout << " | (old: " << mean_ttbar_combine_old << ")" << endl;
    cout << "jmsdd - " << years[year] << " - combine: " << mean_dd_combine    << endl;
    cout << endl;

    // draw --------------------------------------------------------------------

    // draw_plot_jms(TH1F* h1, TH1F* h2, TH1F* h3, const TString path, const TString year, vector<double> mean, TString norm="")
    draw_plot_jms(all_ttbar_muon[year], all_jms_uu_muon[year], all_jms_dd_muon[year], path_muon, years[year], mean_muon);
    draw_plot_jms(all_ttbar_elec[year], all_jms_uu_elec[year], all_jms_dd_elec[year], path_elec, years[year], mean_elec);
    draw_plot_jms(all_ttbar_combine[year], all_jms_uu_combine[year], all_jms_dd_combine[year], path_combine, years[year], mean_combine);

    draw_plot_jms(all_ttbar_muon_norm[year], all_jms_uu_muon_norm[year], all_jms_dd_muon_norm[year], path_muon, years[year], mean_muon, "_norm");
    draw_plot_jms(all_ttbar_elec_norm[year], all_jms_uu_elec_norm[year], all_jms_dd_elec_norm[year], path_elec, years[year], mean_elec, "_norm");
    draw_plot_jms(all_ttbar_combine_norm[year], all_jms_uu_combine_norm[year], all_jms_dd_combine_norm[year], path_combine, years[year], mean_combine, "_norm");

    // draw_plot_comparison(TH1F* h1, TH1F* h2, TH1F* data, const TString path, const TString year, vector<double> mean={0, 0})
    draw_plot_comparison(all_ttbar_muon[year], all_ttbar_muon_old[year], all_data_muon[year], path_muon, years[year], mean_muon_old);
    draw_plot_comparison(all_ttbar_elec[year], all_ttbar_elec_old[year], all_data_elec[year], path_elec, years[year], mean_elec_old);
    draw_plot_comparison(all_ttbar_combine[year], all_ttbar_combine_old[year], all_data_combine[year], path_combine, years[year], mean_combine_old);

  }

}
