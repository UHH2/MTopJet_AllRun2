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
  print_seperater();


  if(argc != 2){
    cout << "\n" << "Usage: ./JEC_SYS <year>\n";
    return 0;
  }

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Year" << endl;

  TString year = argv[1];
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0) is16   = true;
  else if(strcmp(year, "2017")==0) is17   = true;
  else if(strcmp(year, "2018")==0) is18   = true;
  else if(strcmp(year,  "combined")==0) isAll  = true;
  else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or combined)");


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

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots"; // CHANGE_PT
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

  // #################################################################################################
  // Get TTbar #######################################################################################
  if(debug) cout << "TTbar" << endl;

  TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TFile *ttbar_file  = new TFile(dir+year+"/muon/"+ttbar_path);
  TH1F  *ttbar       = (TH1F*)ttbar_file->Get(hist);
  TH1F  *ttbar_old   = (TH1F*)ttbar_file->Get(hist_old);

  ttbar->Add(bkg, 1);
  ttbar_old->Add(bkg_old, 1);

  TH1F* ttbar_norm       = normalize(ttbar);
  TH1F* ttbar_norm_old   = normalize(ttbar_old);


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
  ttbar->GetYaxis()->SetTitleOffset(1.5);
  ttbar->GetXaxis()->SetTitle("m_{jet}");
  ttbar->GetYaxis()->SetTitle("");
  ttbar->SetLineWidth(2);  // ttbar hist style
  ttbar->SetLineColor(kRed);

  // Data --------------------------------------------------------------------------------------------
  data->SetMarkerStyle(8);  // data hist style
  data->SetMarkerColor(kBlack);

  // upup --------------------------------------------------------------------------------------------
  JMSuu->SetLineWidth(2);
  JMSuu->SetLineColor(kBlue);
  // downdown ----------------------------------------------------------------------------------------
  JMSdd->SetLineWidth(2);
  JMSdd->SetLineColor(kGreen);
  // updown ------------------------------------------------------------------------------------------
  JMSud->SetLineWidth(2);
  JMSud->SetLineStyle(7);
  JMSud->SetLineColor(kGray);
  // downup ------------------------------------------------------------------------------------------
  JMSdu->SetLineWidth(2);
  JMSdu->SetLineStyle(7);
  JMSdu->SetLineColor(kOrange);

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
  JMSud->Draw("SAME HIST");
  JMSdu->Draw("SAME HIST");
  leg = new TLegend(0.62,0.65,0.82,0.85);
  leg->AddEntry(ttbar,"nominal","l");
  leg->AddEntry(JMSuu ,"JMS_upup","l");
  leg->AddEntry(JMSdd ,"JMS_downdown","l");
  leg->AddEntry(JMSud ,"JMS_updown","l");
  leg->AddEntry(JMSdu ,"JMS_downup","l");
  leg->SetTextSize(0.02);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path_general+"/mass_comparison_jms_"+year+".pdf");
  data->Draw("SAME P");
  leg->AddEntry(data,"Data","pl");
  A->SaveAs(save_path_general+"/mass_comparison_jms_data_"+year+".pdf");
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
  leg = new TLegend(0.62,0.65,0.82,0.85);
  leg->AddEntry(ttbar,"best fit","l");
  leg->AddEntry(ttbar_old ,"old","l");
  leg->SetTextSize(0.02);
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path_general+"/method_comparison_jms_"+year+".pdf");
  data->Draw("SAME P");
  leg->AddEntry(data,"Data","pl");
  A->SaveAs(save_path_general+"/method_comparison_jms_data_"+year+".pdf");
  delete A;
  leg->Clear();

  cout << "\nIntegral New: " << ttbar->Integral()<<endl;
  cout << "Integral Old: " << ttbar_old->Integral()<<endl;

}
