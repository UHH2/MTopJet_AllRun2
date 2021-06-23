#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/CreatHists.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

TH1F* GetRatio(TH1F* h1, TH1F* h2);

bool isCorr = true;
bool debug = false;
TString save_path;

// -------------------------------------------------------------------------------------------------------
void plot_settings(TH1F* hist, int color, int style, TString title, TString xAxis, TString yAxis)
{
  hist->SetTitle(title);
  hist->GetXaxis()->SetRangeUser(100, 300);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.1);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitle(xAxis);
  hist->GetYaxis()->SetTitle(yAxis);
  hist->SetLineWidth(1);
  hist->SetLineColor(color);
  hist->SetLineStyle(style);
}

// -------------------------------------------------------------------------------------------------------
void plot_settings_ratio(TH1F* hist, int color, int style, TString title, TString xAxis, TString yAxis)
{
  hist->GetXaxis()->SetTickLength(0.07);
  hist->GetXaxis()->SetTitleSize(25);
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleOffset(4.0);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(21);
  hist->GetXaxis()->SetLabelOffset(0.035);
  hist->GetYaxis()->SetTitle(yAxis);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(20);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleOffset(2.2);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(19);
  hist->GetYaxis()->SetLabelOffset(0.009);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetTitle(" ");
  hist->GetYaxis()->SetRangeUser(0.92, 1.08);
  hist->SetLineWidth(1);
  hist->SetLineColor(color);
  hist->SetLineStyle(style);
}


// -------------------------------------------------------------------------------------------------------
void draw_plot_comparison(vector<TH1F*> hists, const TString path, const TString year, const TString channel, vector<TString> leg1, vector<double> mean)
{
  if(hists.size()!=leg1.size()) cerr << "Check leg and hist vector size";

  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  TCanvas *A = new TCanvas("A"+path+year, "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  hists[0]->Draw("HIST");
  for(unsigned int i=1; i<hists.size(); i++) hists[i]->Draw("SAME HIST");

  leg->SetNColumns(2);
  for(unsigned int i=0; i<leg1.size(); i++){
    leg->AddEntry(hists[i],leg1[i],"l");
    leg->AddEntry((TObject*)0, "<#it{m}_{jet}^{reco}>="+dtos(mean[i], 2)+" GeV", "");
  }

  leg->SetBorderSize(0);
  leg->SetTextSize(0.02);
  leg->Draw();

  gPad->RedrawAxis();
  if(isCorr) A->SaveAs(path+"/"+channel+"/Compare_Uncertainty_Method_"+leg1[1]+"_"+year+"_corr.pdf");
  else       A->SaveAs(path+"/"+channel+"/Compare_Uncertainty_Method_"+leg1[1]+"_"+year+".pdf");

  delete A;
  leg->Clear();
}

// -------------------------------------------------------------------------------------------------------
void DrawHist(TH1F* hist, TLegend* leg, double mean, TString entry, TString option, bool isRatio=false)
{
  hist->Draw(option);
  leg->AddEntry(hist,entry,"l");
  leg->AddEntry((TObject*)0, "<#it{m}_{jet}^{reco}>="+dtos(mean, 2)+" GeV", "");
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

  SetupGlobalStyle();
  gErrorIgnoreLevel = kWarning;

  if(argc != 1){
    cout << "\n" << "Usage: ./JMS_compare_PoU\n";
    cout << "All years and channels are created simultanously." << endl;
    return 0;
  }

  if(debug) cout << "Create required folders ..." << endl;
  save_path = get_save_path();
  save_path = save_path+"/Plots"; // CHANGE_PT
  save_path = creat_folder_and_path(save_path, "JMS");
  creat_folder(save_path+"/muon");
  creat_folder(save_path+"/elec");
  creat_folder(save_path+"/combine");

  // vector<int> index_year   = {0,1,2,3};
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
  cout << "Get all histograms ..." << endl;
  // vector<TH1F*> h_data_muon          = get_all_hists(data_muon, hist);
  // vector<TH1F*> h_data_elec          = get_all_hists(data_elec, hist);
  // vector<TH1F*> h_data_combine       = combine_channels(h_data_muon, h_data_elec);
  //
  vector<TH1F*> h_ttbar_muon         = get_all_hists(ttbar_muon, hist);
  vector<TH1F*> h_ttbar_elec         = get_all_hists(ttbar_elec, hist);
  vector<TH1F*> h_ttbar_combine      = combine_channels(h_ttbar_muon, h_ttbar_elec);


  // JMS -----------------------------------------------------------------------
  if(debug) cout << "\t ... JMS_uu" << endl;
  vector<TH1F*> h_jms_uu_muon        = get_all_hists(jms_uu_muon, hist);
  vector<TH1F*> h_jms_uu_elec        = get_all_hists(jms_uu_elec, hist);
  vector<TH1F*> h_jms_uu_combine     = combine_channels(h_jms_uu_muon, h_jms_uu_elec);

  if(debug) cout << "\t ... JMS_dd" << endl;
  vector<TH1F*> h_jms_dd_muon        = get_all_hists(jms_dd_muon, hist);
  vector<TH1F*> h_jms_dd_elec        = get_all_hists(jms_dd_elec, hist);
  vector<TH1F*> h_jms_dd_combine     = combine_channels(h_jms_dd_muon, h_jms_dd_elec);

  if(debug) cout << "\t ... JMS_u" << endl;
  vector<TH1F*> h_jms_u_muon        = get_all_hists(jms_u_muon, hist);
  vector<TH1F*> h_jms_u_elec        = get_all_hists(jms_u_elec, hist);
  vector<TH1F*> h_jms_u_combine     = combine_channels(h_jms_u_muon, h_jms_u_elec);

  if(debug) cout << "\t ... JMS_d" << endl;
  vector<TH1F*> h_jms_d_muon        = get_all_hists(jms_d_muon, hist);
  vector<TH1F*> h_jms_d_elec        = get_all_hists(jms_d_elec, hist);
  vector<TH1F*> h_jms_d_combine     = combine_channels(h_jms_d_muon, h_jms_d_elec);

  if(debug) cout << "\t ... JMS_u_acorr" << endl;
  vector<TH1F*> h_jms_u_muon_corr        = get_all_hists(jms_u_muon_corr, hist);
  vector<TH1F*> h_jms_u_elec_corr        = get_all_hists(jms_u_elec_corr, hist);
  vector<TH1F*> h_jms_u_combine_corr     = combine_channels(h_jms_u_muon_corr, h_jms_u_elec_corr);

  if(debug) cout << "\t ... JMS_d_acorr" << endl;
  vector<TH1F*> h_jms_d_muon_corr        = get_all_hists(jms_d_muon_corr, hist);
  vector<TH1F*> h_jms_d_elec_corr        = get_all_hists(jms_d_elec_corr, hist);
  vector<TH1F*> h_jms_d_combine_corr     = combine_channels(h_jms_d_muon_corr, h_jms_d_elec_corr);

  int index_year = 0; // combine

  // plot_settings(h_ttbar_muon[index_Com], kGray, 1,  "", "m_{jet}", "Events");
  // plot_settings(h_jms_uu_muon[index_Com], kRed+2, 1,  "", "m_{jet}", "Events");
  // plot_settings(h_jms_dd_muon[index_Com], kRed+2, 1,  "", "m_{jet}", "Events");
  // plot_settings(h_jms_u_muon[index_Com], kBlue+2, 1, "", "m_{jet}", "Events");
  // plot_settings(h_jms_d_muon[index_Com], kBlue+2, 1, "", "m_{jet}", "Events");
  // plot_settings(h_jms_u_muon_corr[index_Com], kGreen+2, 1, "", "m_{jet}", "Events");
  // plot_settings(h_jms_d_muon_corr[index_Com], kGreen+2, 1, "", "m_{jet}", "Events");

  cout << "Set plot settings ..." << endl;
  plot_settings(h_ttbar_combine[index_Com], 14, 1,  "", "m_{jet}", "Events");
  plot_settings(h_jms_uu_combine[index_Com], 810, 1,  "", "m_{jet}", "Events");
  plot_settings(h_jms_dd_combine[index_Com], 810, 2,  "", "m_{jet}", "Events");
  plot_settings(h_jms_u_combine[index_Com], kAzure+7, 1, "", "m_{jet}", "Events");
  plot_settings(h_jms_d_combine[index_Com], kAzure+7, 2, "", "m_{jet}", "Events");
  plot_settings(h_jms_u_combine_corr[index_Com], kGreen+2, 1, "", "m_{jet}", "Events");
  plot_settings(h_jms_d_combine_corr[index_Com], kGreen+2, 2, "", "m_{jet}", "Events");

  // double mean_ttbar_muon = trunc_mean(h_ttbar_muon[index_Com], 120, 240); // all_jms_uu_muon[year]->GetMean();
  // double mean_uu_muon = trunc_mean(h_jms_uu_muon[index_Com], 120, 240); // all_jms_uu_muon[year]->GetMean();
  // double mean_dd_muon = trunc_mean(h_jms_dd_muon[index_Com], 120, 240); // all_jms_dd_muon[year]->GetMean();
  // double mean_u_muon = trunc_mean(h_jms_u_muon[index_Com], 120, 240); // all_jms_ud_muon[year]->GetMean();
  // double mean_d_muon = trunc_mean(h_jms_d_muon[index_Com], 120, 240); // all_jms_du_muon[year]->GetMean();
  // double mean_u_muon_corr = trunc_mean(h_jms_u_muon_corr[index_Com], 120, 240); // all_jms_ud_muon[year]->GetMean();
  // double mean_d_muon_corr = trunc_mean(h_jms_d_muon_corr[index_Com], 120, 240); // all_jms_du_muon[year]->GetMean();

  cout << "Get truncated mean ..." << endl;
  double mean_ttbar_combine = trunc_mean(h_ttbar_combine[index_Com], 120, 240); // all_jms_uu_muon[year]->GetMean();
  double mean_uu_combine = trunc_mean(h_jms_uu_combine[index_Com], 120, 240); // all_jms_uu_muon[year]->GetMean();
  double mean_dd_combine = trunc_mean(h_jms_dd_combine[index_Com], 120, 240); // all_jms_dd_muon[year]->GetMean();
  double mean_u_combine = trunc_mean(h_jms_u_combine[index_Com], 120, 240); // all_jms_ud_muon[year]->GetMean();
  double mean_d_combine = trunc_mean(h_jms_d_combine[index_Com], 120, 240); // all_jms_du_muon[year]->GetMean();
  double mean_u_combine_corr = trunc_mean(h_jms_u_combine_corr[index_Com], 120, 240); // all_jms_ud_muon[year]->GetMean();
  double mean_d_combine_corr = trunc_mean(h_jms_d_combine_corr[index_Com], 120, 240); // all_jms_du_muon[year]->GetMean();

  gStyle->SetOptStat(0);
  if(gErrorIgnoreLevel == 2000) cout << "Start plotting (kWarning)..." << endl;
  else if(gErrorIgnoreLevel == 3000) cout << "Start plotting (kError)..." << endl;
  else cout << "Start plotting ..." << endl;
  // vector<TH1F*> muon_up = {h_jms_uu_muon[index_Com], h_jms_u_muon[index_Com]};
  // vector<TH1F*> muon_down = {h_jms_dd_muon[index_Com], h_jms_d_muon[index_Com]};
  // vector<TH1F*> muon_all = {h_jms_uu_muon[index_Com], h_jms_u_muon[index_Com], h_jms_u_muon_corr[index_Com], h_ttbar_muon[index_Com], h_jms_dd_muon[index_Com], h_jms_d_muon[index_Com], h_jms_d_muon_corr[index_Com]};

  // if(debug) cout << "\t ... muon" << endl;
  // draw_plot_comparison(muon_up, save_path, years[index_Com], "muon", {"First", "PoU_up"}, {mean_uu_muon, mean_u_muon});
  // draw_plot_comparison(muon_down, save_path, years[index_Com], "muon", {"First", "PoU_down"}, {mean_dd_muon, mean_d_muon});
  // draw_plot_comparison(muon_all, save_path, years[index_Com], "muon", {"Ellipse up", "PoU up"}, {mean_dd_muon, mean_d_muon});

  // vector<TH1F*> combine_up = {h_jms_uu_combine[index_Com], h_jms_u_combine[index_Com]};
  // vector<TH1F*> combine_down = {h_jms_dd_combine[index_Com], h_jms_d_combine[index_Com]};
  // vector<TH1F*> combine_all = {h_jms_uu_combine[index_Com], h_jms_u_combine[index_Com], h_jms_u_muon_corr[index_Com], h_ttbar_combine[index_Com], h_jms_dd_combine[index_Com], h_jms_d_combine[index_Com], h_jms_d_muon_corr[index_Com]};

  // if(debug) cout << "\t ... combine" << endl;
  // draw_plot_comparison(combine_up, save_path, years[index_Com], "combine", {"First", "PoU_up"}, {mean_uu_combine, mean_u_combine});
  // draw_plot_comparison(combine_down, save_path, years[index_Com], "combine", {"First", "PoU_down"}, {mean_dd_combine, mean_d_combine});
  // draw_plot_comparison(combine_all, save_path, years[index_Com], "combine", {"UpUp", "PoU down", "PoU down fully acorr", "PoU_down"}, {mean_dd_combine, mean_d_combine});

  if(debug) cout << "\t ... rho -1 & calculated" << endl;
  TCanvas *A = new TCanvas("A_combine", "A_combine", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  TLegend *leg = new TLegend(0.52,0.55,0.85,0.85);
  leg->SetNColumns(2);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  // void DrawHist(TH1F* hist, TLegend* leg, double mean, TString leg,TString option="hist")
  TString rho = "-0.21";
  h_ttbar_combine[index_Com]->GetXaxis()->SetLabelSize(0);
  DrawHist(h_ttbar_combine[index_Com], leg, mean_ttbar_combine, "Nominal", "HIST");
  DrawHist(h_jms_uu_combine[index_Com], leg, mean_uu_combine, "UpUp", "HIST SAME");
  DrawHist(h_jms_dd_combine[index_Com], leg, mean_dd_combine, "DownDown", "HIST SAME");
  DrawHist(h_jms_u_combine[index_Com], leg, mean_u_combine, "Up (#rho="+rho+")", "HIST SAME");
  DrawHist(h_jms_d_combine[index_Com], leg, mean_d_combine, "Down (#rho="+rho+")", "HIST SAME");
  DrawHist(h_jms_u_combine_corr[index_Com], leg, mean_u_combine_corr, "Up (#rho=-1)", "HIST SAME");
  DrawHist(h_jms_d_combine_corr[index_Com], leg, mean_d_combine_corr, "Down (#rho=-1)", "HIST SAME");

  leg->Draw();
  A->cd();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  if(debug) cout << "\t ... rho -1 & calculated in ratio" << endl;
  TH1F* ratio_tt = GetRatio(h_ttbar_combine[index_Com], h_ttbar_combine[index_Com]);
  TH1F* ratio_uu = GetRatio(h_ttbar_combine[index_Com], h_jms_uu_combine[index_Com]);
  TH1F* ratio_dd = GetRatio(h_ttbar_combine[index_Com], h_jms_dd_combine[index_Com]);
  TH1F* ratio_u = GetRatio(h_ttbar_combine[index_Com], h_jms_u_combine[index_Com]);
  TH1F* ratio_d = GetRatio(h_ttbar_combine[index_Com], h_jms_d_combine[index_Com]);
  TH1F* ratio_uc = GetRatio(h_ttbar_combine[index_Com], h_jms_u_combine_corr[index_Com]);
  TH1F* ratio_dc = GetRatio(h_ttbar_combine[index_Com], h_jms_d_combine_corr[index_Com]);

  plot_settings_ratio(ratio_tt, 14, 1,  "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_uu, 810, 1,  "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_dd, 810, 2,  "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_u, kAzure+7, 1, "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_d, kAzure+7, 2, "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_uc, kGreen+2, 1, "", "m_{jet}", "#frac{JMS}{nom}");
  plot_settings_ratio(ratio_dc, kGreen+2, 2, "", "m_{jet}", "#frac{JMS}{nom}");

  // void DrawHist(TH1F* hist, TLegend* leg, double mean, TString leg,TString option="hist")
  ratio_tt->Draw("HIST");
  ratio_uu->Draw("HIST SAME");
  ratio_dd->Draw("HIST SAME");
  ratio_u->Draw("HIST SAME");
  ratio_d->Draw("HIST SAME");
  ratio_uc->Draw("HIST SAME");
  ratio_dc->Draw("HIST SAME");

  gPad->RedrawAxis();
  A->SaveAs(save_path+"/combine/Compare_Uncertainty_Method.pdf");
  delete A;
  leg->Clear();

  cout << endl;
}

// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------------
TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetNbinsX();
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
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
