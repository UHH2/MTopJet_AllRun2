#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  /*
  Show the behaviour af the different mass scale variations uu, ud, du an dd.
  Show the difference between the old Mass Jet distribution and the new one
  using the BestFit method.
  */
  bool debug = true;
  TString save_path = get_save_path();

  print_seperater();


  if(argc != 3){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin>\n";
    return 0;
  }

  // #################################################################################################
  // Only one fit for all bins #######################################################################
  if(debug) cout << "String into bool" << endl;

  // Default -------------------------------------------------------------------
  Int_t   oldLevel     = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ... - functions_explain
  gErrorIgnoreLevel    = kWarning;          // suppress TCanvas output

  // Input ---------------------------------------------------------------------
  int     bin_width   = stoi(argv[2]);
  TString year        = argv[1];

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Year" << endl;

  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0) is16   = true;
  else if(strcmp(year, "2017")==0) is17   = true;
  else if(strcmp(year, "2018")==0) is18   = true;
  else if(strcmp(year,  "combined")==0) isAll  = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or combined)");

  // Rebin -------------------------------------------------------------------------------------------
  if(debug) cout << "Set Number Bins" << endl;
  int     number_bins     = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders


  /*
  ██████  ██ ██████  ███████  ██████ ████████  ██████  ██████  ██ ███████ ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██      ██
  ██   ██ ██ ██████  █████   ██         ██    ██    ██ ██████  ██ █████   ███████
  ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██           ██
  ██████  ██ ██   ██ ███████  ██████    ██     ██████  ██   ██ ██ ███████ ███████
  */

  // #################################################################################################
  // creating subdirectories. Necessary, because different binning examined ##########################
  /*
  Explanation for mkdir(char **, mode_t mode) and mode_t: https://jameshfisher.com/2017/02/24/what-is-mode_t/
  example:    0002 is -------w-
  .           0003 is -------wx
  .           0004 is ------r--
  .           0005 is ------r-x
  .           0006 is ------rw-
  .           0007 is ------rwx
  positions:  xyzw - x: type (folder, etc.) not necessary here;- y: owner;- z: group;- w: other;
  */

  TString save_path_general = save_path+"/Plots"; // CHANGE_PT
  save_path_general = creat_folder_and_path(save_path_general, "JMS");


  /*
  .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██        ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
  .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
  */

  cout << '\n'; // CHANGE_PT
  TString hist; TString hist_class = "JetMassScaleHists/";
  hist = hist_class+"hadjet_jms_mass";

  TString hist_old; TString hist_class_old = "XCone_cor/";
  hist_old = hist_class_old+"M_jet1";

  // #################################################################################################
  // Get Background ##################################################################################
  if(debug) cout << "Background" << endl;

  vector<TFile*> file_bkg_v;
  vector<TH1F*> hists_bkg_v, hists_bkg_old_v;
  vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
  for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+year+"/muon/"+path_bkg_v[i]));
  for(unsigned int i=0; i<file_bkg_v.size(); i++){
    hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(hist));
    hists_bkg_old_v.push_back((TH1F*)file_bkg_v[i]->Get(hist_old));
  }
  TH1F *bkg = AddHists(hists_bkg_v, 1);
  TH1F *bkg_old = AddHists(hists_bkg_old_v, 1);

  // #################################################################################################
  // Get Data ########################################################################################
  if(debug) cout << "Data" << endl;

  TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  TFile *data_file      = new TFile(dir+year+"/muon/"+data_path);
  TH1F  *data           = (TH1F*)data_file->Get(hist);
  TH1F  *data_old       = (TH1F*)data_file->Get(hist_old);

  TH1F  *data_norm       = normalize(data);
  TH1F  *data_rebin      = rebin(data, bin_width);
  TH1F  *data_rebin_norm = normalize(data_rebin);


  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file  = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F  *ttbar       = (TH1F*)ttbar_file->Get(hist);
  TH1F  *ttbar_old   = (TH1F*)ttbar_file->Get(hist_old);

  ttbar->Add(bkg, 1);
  ttbar_old->Add(bkg_old, 1);

  TH1F* ttbar_norm              = normalize(ttbar);
  TH1F* ttbar_norm_old          = normalize(ttbar_old);
  TH1F* ttbar_rebin             = rebin(ttbar, bin_width);
  TH1F* ttbar_rebin_old         = rebin(ttbar_old, bin_width);
  TH1F* ttbar_rebin_norm        = normalize(ttbar_rebin);
  TH1F* ttbar_rebin_norm_old    = normalize(ttbar_rebin_old);


  // #################################################################################################
  // Get SYS #########################################################################################
  if(debug) cout << "JMS" << endl;

  TFile *JECuu_file   = new TFile(dir+year+"/muon/JMS_upup/"+ttbar_path);
  TFile *JECdd_file   = new TFile(dir+year+"/muon/JMS_downdown/"+ttbar_path);
  TFile *JECdu_file   = new TFile(dir+year+"/muon/JMS_downup/"+ttbar_path);
  TFile *JECud_file   = new TFile(dir+year+"/muon/JMS_updown/"+ttbar_path);
  TH1F  *JMSuu        = (TH1F*)JECuu_file->Get(hist);
  TH1F  *JMSdd        = (TH1F*)JECdd_file->Get(hist);
  TH1F  *JMSud        = (TH1F*)JECud_file->Get(hist);
  TH1F  *JMSdu        = (TH1F*)JECdu_file->Get(hist);

  // Add hists from JEC and background -----------------------------------------
  if(debug) cout << "JEC+Background" << endl;

  JMSuu->Add(bkg, 1);
  JMSdd->Add(bkg, 1);
  JMSud->Add(bkg, 1);
  JMSdu->Add(bkg, 1);

  TH1F* JMSuu_norm       = normalize(JMSuu);
  TH1F* JMSdd_norm       = normalize(JMSdd);
  TH1F* JMSud_norm       = normalize(JMSud);
  TH1F* JMSdu_norm       = normalize(JMSdu);


  /*
  ██████  ██       ██████  ████████ ███████
  ██   ██ ██      ██    ██    ██    ██
  ██████  ██      ██    ██    ██    ███████
  ██      ██      ██    ██    ██         ██
  ██      ███████  ██████     ██    ███████
  */

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  if(debug) cout << "Plots Sensitivity" << endl;

  // #################################################################################################
  // Settings ########################################################################################

  ttbar->SetTitle("");
  ttbar->GetXaxis()->SetRangeUser(100, 300);
  ttbar->GetYaxis()->SetRangeUser(0, ttbar->GetMaximum()*1.1);
  ttbar->GetXaxis()->SetNdivisions(505);
  ttbar->GetYaxis()->SetNdivisions(505);
  ttbar->GetXaxis()->SetTitleSize(0.05);
  ttbar->GetYaxis()->SetTitleSize(0.04);
  ttbar->GetXaxis()->SetTitleOffset(0.9);
  ttbar->GetYaxis()->SetTitleOffset(0.9);
  ttbar->GetXaxis()->SetTitle("m_{jet} [GeV]");
  ttbar->GetYaxis()->SetTitle("Events");
  ttbar->GetYaxis()->SetTitleOffset(1.5);
  ttbar->SetLineWidth(2);  // ttbar hist style
  ttbar->SetLineColor(kRed);

  ttbar_norm->SetTitle("");
  ttbar_norm->GetXaxis()->SetRangeUser(100, 300);
  ttbar_norm->GetYaxis()->SetRangeUser(0, ttbar_norm->GetMaximum()*1.1);
  ttbar_norm->GetXaxis()->SetNdivisions(505);
  ttbar_norm->GetYaxis()->SetNdivisions(505);
  ttbar_norm->GetXaxis()->SetTitleSize(0.05);
  ttbar_norm->GetYaxis()->SetTitleSize(0.04);
  ttbar_norm->GetXaxis()->SetTitleOffset(0.9);
  ttbar_norm->GetYaxis()->SetTitleOffset(0.9);
  ttbar_norm->GetXaxis()->SetTitle("m_{jet} [GeV]");
  ttbar_norm->GetYaxis()->SetTitle("#Delta N/N");
  ttbar_norm->GetYaxis()->SetTitleOffset(1.5);
  ttbar_norm->SetLineWidth(2);  // ttbar hist style
  ttbar_norm->SetLineColor(kRed);

  // Data --------------------------------------------------------------------------------------------
  data->SetMarkerStyle(8);  // data hist style
  data->SetMarkerColor(kBlack);
  data->SetLineColor(kBlack);
  data_norm->SetMarkerStyle(8);  // data hist style
  data_norm->SetMarkerColor(kBlack);
  data_norm->SetLineColor(kBlack);
  // upup --------------------------------------------------------------------------------------------
  JMSuu->SetLineWidth(2);
  JMSuu->SetLineColor(kBlue);
  JMSuu_norm->SetLineWidth(2);
  JMSuu_norm->SetLineColor(kBlue);
  // downdown ----------------------------------------------------------------------------------------
  JMSdd->SetLineWidth(2);
  JMSdd->SetLineColor(kGreen+2);
  JMSdd_norm->SetLineWidth(2);
  JMSdd_norm->SetLineColor(kGreen+2);
  // updown ------------------------------------------------------------------------------------------
  JMSud->SetLineWidth(2);
  JMSud->SetLineStyle(7);
  JMSud->SetLineColor(kGray+2);
  JMSud_norm->SetLineWidth(2);
  JMSud_norm->SetLineStyle(7);
  JMSud_norm->SetLineColor(kYellow+1);
  // downup ------------------------------------------------------------------------------------------
  JMSdu->SetLineWidth(2);
  JMSdu->SetLineStyle(7);
  JMSdu->SetLineColor(kMagenta+2);
  JMSdu_norm->SetLineWidth(2);
  JMSdu_norm->SetLineStyle(7);
  JMSdu_norm->SetLineColor(kMagenta+2);

  // Legend ------------------------------------------------------------------------------------------
  TLegend *leg;


  /*
  ███    ███  █████  ███████ ███████     ██████  ██       ██████  ████████ ███████
  ████  ████ ██   ██ ██      ██          ██   ██ ██      ██    ██    ██    ██
  ██ ████ ██ ███████ ███████ ███████     ██████  ██      ██    ██    ██    ███████
  ██  ██  ██ ██   ██      ██      ██     ██      ██      ██    ██    ██         ██
  ██      ██ ██   ██ ███████ ███████     ██      ███████  ██████     ██    ███████
  */
  if(debug) cout << "Mass Plots" << endl;

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar->Draw("HIST");
  JMSuu->Draw("SAME HIST");
  JMSdd->Draw("SAME HIST");
  // data->Draw("SAME P");
  leg = new TLegend(0.55,0.55,0.85,0.85);
  leg->AddEntry(ttbar,"Best fit","l");
  leg->AddEntry(JMSuu ,"JMS upup","l");
  leg->AddEntry(JMSdd ,"JMS downdown","l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path_general+"/mass_comparison_jms_"+year+".pdf");

  JMSud->Draw("SAME HIST");
  JMSdu->Draw("SAME HIST");
  leg->AddEntry(JMSud ,"JMS updown","l");
  leg->AddEntry(JMSdu ,"JMS downup","l");
  A->SaveAs(save_path_general+"/mass_comparison_jms_add_"+year+".pdf");
  delete A;
  leg->Clear();

  // #################################################################################################
  // Norm ############################################################################################
  if(debug) cout << "Mass Plots Norm" << endl;

  A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar_norm->Draw("HIST");
  // data_norm->Draw("SAME P");
  JMSuu_norm->Draw("SAME HIST");
  JMSdd_norm->Draw("SAME HIST");
  leg = new TLegend(0.55,0.55,0.85,0.85);
  leg->AddEntry(ttbar_norm,"Best fit","l");
  leg->AddEntry(data_norm,"Data","pl");
  leg->AddEntry(JMSuu_norm ,"JMS up","l");
  leg->AddEntry(JMSdd_norm ,"JMS down","l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path_general+"/mass_comparison_jms_norm_"+year+".pdf");

  JMSud_norm->Draw("SAME HIST");
  JMSdu_norm->Draw("SAME HIST");
  leg->AddEntry(JMSud_norm ,"JMS updown","l");
  leg->AddEntry(JMSdu_norm ,"JMS downup","l");
  A->SaveAs(save_path_general+"/mass_comparison_jms_add_norm_"+year+".pdf");
  delete A;
  leg->Clear();

  /*
  .██████  ██████  ███    ███ ██████   █████  ██████  ███████     ██████  ██       ██████  ████████ ███████
  ██      ██    ██ ████  ████ ██   ██ ██   ██ ██   ██ ██          ██   ██ ██      ██    ██    ██    ██
  ██      ██    ██ ██ ████ ██ ██████  ███████ ██████  █████       ██████  ██      ██    ██    ██    ███████
  ██      ██    ██ ██  ██  ██ ██      ██   ██ ██   ██ ██          ██      ██      ██    ██    ██         ██
  .██████  ██████  ██      ██ ██      ██   ██ ██   ██ ███████     ██      ███████  ██████     ██    ███████
  */
  if(debug) cout << "\nCompare Plots" << endl;

  // Old ttbar ---------------------------------------------------------------------------------------
  ttbar_old->SetLineWidth(2);
  ttbar_old->SetLineColor(kBlue);

  A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  ttbar->Draw("HIST");
  ttbar_old->Draw("SAME HIST");
  leg = new TLegend(0.60,0.65,0.75,0.85);
  leg->AddEntry(ttbar,"Best fit","l");
  leg->AddEntry(ttbar_old ,"Old","l");
  leg->SetTextSize(0.04);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path_general+"/method_comparison_jms_old_"+year+".pdf");
  data->Draw("SAME P");
  leg->AddEntry(data,"Data","pl");
  A->SaveAs(save_path_general+"/method_comparison_jms_old_data_"+year+".pdf");
  delete A;
  leg->Clear();

  cout << "\nIntegral New: " << ttbar->Integral()<<endl;
  cout << "Integral Old: "   << ttbar_old->Integral()<<endl;

}
