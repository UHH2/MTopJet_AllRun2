#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

// #include <io.h>
using namespace std;

int main(int argc, char* argv[]){
  bool debug = true;

  // #################################################################################################
  // Modification ####################################################################################
  TString mTop_ = "mtop1695";           // CHANGE_MTOP
  // Changes of how many pt bins to use // CHANGE_PT
  // Modification ####################################################################################
  // #################################################################################################

  // #################################################################################################
  // Declare different variables used in the code ####################################################
  print_seperater();

  if(argc != 4){
    cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number> <one_Ptbin>\n";
    cout << "rebin_number:        hist->Rebin(rebin_number)\n";
    return 0;
  }

  /*
  ██████  ███████ ███████ ██ ███    ██ ███████
  ██   ██ ██      ██      ██ ████   ██ ██
  ██   ██ █████   █████   ██ ██ ██  ██ █████
  ██   ██ ██      ██      ██ ██  ██ ██ ██
  ██████  ███████ ██      ██ ██   ████ ███████
  */

  // Default -------------------------------------------------------------------
  cout.precision(6);
  Int_t   oldLevel     = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ... - functions_explain
  gErrorIgnoreLevel    = kWarning;          // suppress TCanvas output
  TString reconst      = "btag";    // match btag_cut btag_sel compare min_mass
  TString MC_uncert    = "central"; // central avg nominal
  bool    add_mTop     = false;
  bool    usePeak_in   = true;

  // Input ---------------------------------------------------------------------
  bool    one_Ptbin   = stob(argv[3]);
  int     bin_width   = stoi(argv[2]);
  TString year        = argv[1];

  // Rebin -------------------------------------------------------------------------------------------
  int number_bins = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders

  // Year --------------------------------------------------------------------------------------------
  if(debug) cout << "Getting Input Variables" << endl;

  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16   = true;
  else if(strcmp(year, "2017")==0)     is17   = true;
  else if(strcmp(year, "2018")==0)     is18   = true;
  else if(strcmp(year, "combined")==0) isAll  = true;
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

  TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS";
  save_path_general = creat_folder_and_path(save_path_general, "MassPlots");  // CHANGE_PT

  save_path_general = creat_folder_and_path(save_path_general, year);
  save_path_general = creat_folder_and_path(save_path_general, reconst);
  save_path_general = creat_folder_and_path(save_path_general, "rebin"+str_number_bins);
  if(add_mTop) save_path_general = creat_folder_and_path(save_path_general, "mtop_combi");

  /*
  ██████  ████████ ██████  ██ ███    ██ ███████
  ██   ██    ██    ██   ██ ██ ████   ██ ██
  ██████     ██    ██████  ██ ██ ██  ██ ███████
  ██         ██    ██   ██ ██ ██  ██ ██      ██
  ██         ██    ██████  ██ ██   ████ ███████
  */

  for(int ptbin=0; ptbin<5; ptbin++){ // CHANGE_PT
    if((ptbin==4)) continue;
    TString addition="";
    if(ptbin==0) addition="_hh"; // CHANGE_PT
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    if(ptbin==4) addition="";

    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */

    vector<TH1F*> allyears_bkg,  allyears_data,  allyears_ttbar;
    vector<TH1F*> allyears_JECup, allyears_JECdown, allyears_XCup, allyears_XCdown;

    cout << '\n'; // CHANGE_PT
    TString w_mass; TString hist_class = "comparison_topjet_xcone_pass_rec/";
    // if(reconst=="btag")      w_mass = hist_class+"wmass_match";                    // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet
    if(reconst=="btag"&&ptbin==0) w_mass = hist_class+"wmass_match_ptdiv_hh";
    if(reconst=="btag"&&ptbin==1) w_mass = hist_class+"wmass_match_ptdiv_hl";
    if(reconst=="btag"&&ptbin==2) w_mass = hist_class+"wmass_match_ptdiv_lh";
    if(reconst=="btag"&&ptbin==3) w_mass = hist_class+"wmass_match_ptdiv_ll";
    if(reconst=="btag"&&ptbin==4) w_mass = hist_class+"wmass_match";
    // if(reconst=="btag")      w_mass = hist_class+"wjet_pt_match_S1divW_W";              // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet

    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    vector<TFile*> file_bkg_v;
    vector<TH1F*>  hists_bkg_v;
    vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
    for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+year+"/muon/"+path_bkg_v[i]));
    for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
    TH1F *bkg = AddHists(hists_bkg_v, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    TFile *data_file      = new TFile(dir+year+"/muon/"+data_path);
    TH1F *data            = (TH1F*)data_file->Get(w_mass);
    TH1F *data_rebin      = rebin(data, bin_width);
    TH1F *data_norm       = normalize(data);
    TH1F *data_rebin_norm = normalize(data_rebin);

    vector<double> data_rebin_norm_err = normalize_error(data_rebin);
    for(int ipar=0; ipar<data_rebin_norm->GetNbinsX(); ipar++) data_rebin_norm->SetBinError(ipar+1, data_rebin_norm_err[ipar]);

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_path       = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
    TFile  *ttbar_file       = new TFile(dir+year+"/muon/"+ttbar_path);
    TH1F   *ttbar            = (TH1F*)ttbar_file->Get(w_mass);
    ttbar->Add(bkg, 1);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC" << endl;

    TFile *JECup_file   = new TFile(dir+year+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    TFile *JECdown_file = new TFile(dir+year+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    TH1F  *JECup        = (TH1F*)JECup_file->Get(w_mass);
    TH1F  *JECdown      = (TH1F*)JECdown_file->Get(w_mass);

    // Add hists from JEC and background ---------------------------------------------------------------
    if(debug) cout << "JEC+Background" << endl;

    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone" << endl;

    TFile *XCup_file   = new TFile(dir+year+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    TFile *XCdown_file = new TFile(dir+year+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    TH1F  *XConeup     = (TH1F*)XCup_file->Get(w_mass);
    TH1F  *XConedown   = (TH1F*)XCdown_file->Get(w_mass);

    // Add hists from FSR and background ---------------------------------------------------------------
    if(debug) cout << "XCone+Background" << endl;

    XConeup->Add(bkg, 1);
    XConedown->Add(bkg, 1);

    /*
    ███    ███ ████████  ██████  ██████
    ████  ████    ██    ██    ██ ██   ██
    ██ ████ ██    ██    ██    ██ ██████
    ██  ██  ██    ██    ██    ██ ██
    ██      ██    ██     ██████  ██
    */

    /*
    If mTop is used, the histograms need to be added to the nominal ones before the normalization
    is done. Therefore, the Code normalizes, rebins and calculates the error after the decision is made
    to include the mTop samples
    */

    vector<TH1F*> mTop_nominal_v, mTop_jecup_v, mTop_jecdown_v, mTop_xcup_v, mTop_xcdown_v;
    vector<TString> mtop_names = {"1665", "1695", "1715", "1735", "1755", "1785"};
    if(add_mTop&&isAll){
      for(int mtop=0; mtop<6; mtop++){
        if(mtop!=1&&mtop!=4) continue;
        // if(debug) cout << "mTop "+mtop_names[mtop] << endl;
        cout << "mTop "+mtop_names[mtop] << endl;

        TString mTop_path       = "uhh2.AnalysisModuleRunner.MC.TTbar.root"; // CHANGE_NORMAL
        TFile  *mTop_file       = new TFile(dir+year+"/muon/"+mTop_path);
        TH1F   *mTop            = (TH1F*)mTop_file->Get(w_mass);
        mTop->Add(bkg, 1);

        // #################################################################################################
        // Get SYS #########################################################################################
        if(debug) cout << "JEC" << endl;

        TFile *mTop_JECup_file   = new TFile(dir+year+"/muon/JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_names[mtop]+".root");
        TFile *mTop_JECdown_file = new TFile(dir+year+"/muon/JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_names[mtop]+".root");
        TH1F  *mTop_JECup        = (TH1F*)mTop_JECup_file->Get(w_mass);
        TH1F  *mTop_JECdown      = (TH1F*)mTop_JECdown_file->Get(w_mass);
        mTop_JECup->Add(bkg, 1);
        mTop_JECdown->Add(bkg, 1);

        // ------------------------------------------------------------------------------------------------
        if(debug) cout << "XCone" << endl;

        TFile *mTop_XCup_file   = new TFile(dir+year+"/muon/COR_up/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_names[mtop]+".root");
        TFile *mTop_XCdown_file = new TFile(dir+year+"/muon/COR_down/uhh2.AnalysisModuleRunner.MC.TTbar_mtop"+mtop_names[mtop]+".root");
        TH1F  *mTop_XCup        = (TH1F*)mTop_XCup_file->Get(w_mass);
        TH1F  *mTop_XCdown      = (TH1F*)mTop_XCdown_file->Get(w_mass);
        mTop_XCup->Add(bkg, 1);
        mTop_XCdown->Add(bkg, 1);

        // #################################################################################################
        // All in one ######################################################################################
        mTop_nominal_v.push_back(mTop);
        mTop_jecup_v.push_back(mTop_JECup);
        mTop_jecdown_v.push_back(mTop_JECdown);
        mTop_xcup_v.push_back(mTop_XCup);
        mTop_xcdown_v.push_back(mTop_XCdown);
      }

      TH1F* mTop_nominal = AddHists(mTop_nominal_v ,1);
      TH1F* mTop_jecup   = AddHists(mTop_jecup_v ,1);
      TH1F* mTop_jecdown = AddHists(mTop_jecdown_v ,1);
      TH1F* mTop_xcup    = AddHists(mTop_xcup_v ,1);
      TH1F* mTop_xcdown  = AddHists(mTop_xcdown_v ,1);

      // #################################################################################################
      // Nominal+samples #################################################################################

      ttbar->Add(mTop_nominal ,1);
      JECup->Add(mTop_jecup ,1);
      JECdown->Add(mTop_jecdown ,1);
      XConeup->Add(mTop_xcup ,1);
      XConedown->Add(mTop_xcdown ,1);
    }

    TH1F* ttbar_rebin          = rebin(ttbar, bin_width);
    TH1F* ttbar_norm           = normalize(ttbar);
    TH1F* ttbar_rebin_norm     = normalize(ttbar_rebin);

    TH1F* JECup_norm           = normalize(JECup);
    TH1F* JECdown_norm         = normalize(JECdown);
    TH1F* JECup_rebin          = rebin(JECup, bin_width);
    TH1F* JECdown_rebin        = rebin(JECdown, bin_width);
    TH1F* JECup_rebin_norm     = rebin(JECup_norm, bin_width);
    TH1F* JECdown_rebin_norm   = rebin(JECdown_norm, bin_width);

    TH1F* XConeup_norm         = normalize(XConeup);
    TH1F* XConedown_norm       = normalize(XConedown);
    TH1F* XConeup_rebin        = rebin(XConeup, bin_width);
    TH1F* XConedown_rebin      = rebin(XConedown, bin_width);
    TH1F* XConeup_rebin_norm   = rebin(XConeup_norm, bin_width);
    TH1F* XConedown_rebin_norm = rebin(XConedown_norm, bin_width);

    vector<double> ttbar_rebin_norm_err     = normalize_error(ttbar_rebin);
    vector<double> JECup_rebin_norm_err     = normalize_error(JECup_rebin);
    vector<double> JECdown_rebin_norm_err   = normalize_error(JECdown_rebin);
    vector<double> XConeup_rebin_norm_err   = normalize_error(XConeup_rebin);
    vector<double> XConedown_rebin_norm_err = normalize_error(XConedown_rebin);

    vector<TH1F*> corrections_up_norm = {JECup_rebin_norm, XConeup_rebin_norm};
    vector<TH1F*> corrections_down_norm = {JECdown_rebin_norm, XConedown_rebin_norm};

    // #################################################################################################
    // Masspeack #######################################################################################
    if(debug) cout << "Masspeak Bins" << endl;

    bool usePeak =false;
    vector<double> PeakLimit;
    if(isAll&&usePeak_in) usePeak = true;

    double Limit = 75; // 190; // CHANGE_PT
    vector<int> peak_bins_v;
    if(usePeak) peak_bins_v = bins_upper_limit(data_rebin, Limit); // Get bins withc bin-content>Limit

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
    /*
    Plotting the histograms included above.
    */
    if(debug) cout << "Plots Sensitivity" << endl;

    // #################################################################################################
    // Settings ########################################################################################
    /* The vectors are only used here */
    vector<TH1F*> ttbar_ = {ttbar_rebin, ttbar_rebin_norm};
    vector<TH1F*> data_ = {data_rebin, data_rebin_norm};
    vector<vector<TH1F*>> corrections_up_ = {{JECup_rebin, XConeup_rebin}, corrections_up_norm};
    vector<vector<TH1F*>> corrections_down_ = {{JECdown_rebin, XConedown_rebin}, corrections_down_norm};

    for(unsigned int norm=0; norm<ttbar_.size(); norm++){
      if(debug) cout << "Inside Norm - " << norm << endl;

      ttbar_[norm]->SetTitle("");
      if     (ptbin==1&&bin_width==5) ttbar_[norm]->GetXaxis()->SetRangeUser(20, 180);
      else if(ptbin==1&&bin_width==6) ttbar_[norm]->GetXaxis()->SetRangeUser(20, 180);
      else                            ttbar_[norm]->GetXaxis()->SetRangeUser(0, 180);

      if(bin_width==1)
      {
        ttbar_[norm]->GetXaxis()->SetRangeUser(peak_bins_v[0], peak_bins_v[peak_bins_v.size()-1]);
      }
      ttbar_[norm]->GetYaxis()->SetRangeUser(0, ttbar_[norm]->GetMaximum()*1.2);
      ttbar_[norm]->GetXaxis()->SetNdivisions(505);
      ttbar_[norm]->GetYaxis()->SetNdivisions(505);
      ttbar_[norm]->GetXaxis()->SetTitleSize(0.04);
      ttbar_[norm]->GetYaxis()->SetTitleSize(0.04);
      ttbar_[norm]->GetXaxis()->SetTitleOffset(0.9);
      ttbar_[norm]->GetYaxis()->SetTitleOffset(1.0);
      ttbar_[norm]->GetXaxis()->SetTitle("m_{Wjet} [GeV]");
      if(norm==0) ttbar_[norm]->GetYaxis()->SetTitle("Events");
      else        ttbar_[norm]->GetYaxis()->SetTitle("#DeltaN/N");
      ttbar_[norm]->GetYaxis()->SetTitleOffset(1.2);
      ttbar_[norm]->SetLineWidth(2);  // ttbar hist style
      ttbar_[norm]->SetLineColor(kRed);

      // Data --------------------------------------------------------------------------------------------
      data_[norm]->SetMarkerStyle(8);  // data hist style
      data_[norm]->SetMarkerColor(kBlack);
      data_[norm]->SetLineColor(kBlack);
    }

    // Legend ------------------------------------------------------------------------------------------
    TLegend *leg;

    /*
    ███    ███  █████  ███████ ███████     ██████  ██       ██████  ████████ ███████
    ████  ████ ██   ██ ██      ██          ██   ██ ██      ██    ██    ██    ██
    ██ ████ ██ ███████ ███████ ███████     ██████  ██      ██    ██    ██    ███████
    ██  ██  ██ ██   ██      ██      ██     ██      ██      ██    ██    ██         ██
    ██      ██ ██   ██ ███████ ███████     ██      ███████  ██████     ██    ███████
    */
    // TLine *line = new TLine(0, Limit, 180, Limit);
    // line->SetLineColor(kGray);
    // line->SetLineWidth(1);
    if(debug) cout << "Mass Plots" << endl;
    TString norm_str="";
    for(unsigned int norm=0; norm<ttbar_.size(); norm++){
      if(norm==1) norm_str="_norm";
      for(int correction=0; correction<2; correction++){
        gStyle->SetPadTickY(1);
        gStyle->SetPadTickX(1);
        gStyle->SetOptStat(kFALSE);
        gStyle->SetLegendBorderSize(0);
        /*
        Plotting the histograms included above.
        */
        if(debug) cout << "Mass Plots" << endl;

        corrections_up_[norm][correction]->SetLineWidth(2);
        corrections_up_[norm][correction]->SetLineColor(kBlue);

        corrections_down_[norm][correction]->SetLineStyle(2);
        corrections_down_[norm][correction]->SetLineWidth(2);
        corrections_down_[norm][correction]->SetLineColor(kBlue);

        ttbar_[norm]->GetYaxis()->SetTitleOffset(1.3);

        TCanvas *A = new TCanvas("A", "A", 600, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        ttbar_[norm]->Draw("HIST");
        corrections_down_[norm][correction]->Draw("SAME HIST");
        corrections_up_[norm][correction]->Draw("SAME HIST");
        data_[norm]->Draw("SAME P");
        // if(norm==0 && bin_width==1 && isAll) line->Draw("SAME");
        leg = new TLegend(0.6,0.65,0.8,0.85);
        if(bin_width==1) leg = new TLegend(0.45,0.15,0.6,0.35);;
        // leg->SetTextSize(0.02);
        leg->AddEntry(ttbar_[norm],"Nominal","l");
        if(correction==0){
          leg->AddEntry(corrections_up_[norm][correction],"JEC up","l");
          leg->AddEntry(corrections_down_[norm][correction],"JEC down","l");
        }
        if(correction==1){
          leg->AddEntry(corrections_up_[norm][correction],"XCone up","l");
          leg->AddEntry(corrections_down_[norm][correction],"XCone down","l");
        }
        leg->AddEntry(data_[norm],"Data","pl");
        leg->SetTextSize(0.03);
        leg->Draw();
        gPad->RedrawAxis();
        if(correction==0) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_JEC"+addition+norm_str+".pdf");
        if(correction==1) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_XCone"+addition+norm_str+".pdf");
        delete A;
        leg->Clear();
      }
    }

    if(bin_width==1){
      corrections_up_[1][0]->SetLineWidth(2);
      corrections_up_[1][0]->SetLineColor(kBlue);
      corrections_down_[1][0]->SetLineStyle(2);
      corrections_down_[1][0]->SetLineWidth(2);
      corrections_down_[1][0]->SetLineColor(kBlue);

      corrections_up_[1][1]->SetLineWidth(2);
      corrections_up_[1][1]->SetLineColor(kGreen+3);
      corrections_down_[1][1]->SetLineStyle(2);
      corrections_down_[1][1]->SetLineWidth(2);
      corrections_down_[1][1]->SetLineColor(kGreen+3);

      TLine *line_up  = new TLine(0.,0.,1.,1.);
      TLine *line_down= new TLine(0.,0.,1.,1.);
      line_up->SetLineStyle(1);
      line_up->SetLineWidth(2);
      line_up->SetLineColor(kGray+2);
      line_down->SetLineStyle(2);
      line_down->SetLineWidth(2);
      line_down->SetLineColor(kGray+2);

      ttbar_[1]->GetYaxis()->SetTitleOffset(1.2);

      TCanvas *A = new TCanvas("A", "A", 600, 600);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.12);
      ttbar_[1]->Draw("HIST");
      corrections_down_[1][0]->Draw("SAME HIST");
      corrections_up_[1][0]->Draw("SAME HIST");
      corrections_down_[1][1]->Draw("SAME HIST");
      corrections_up_[1][1]->Draw("SAME HIST");
      data_[1]->Draw("SAME P");
      // if(norm==0 && bin_width==1 && isAll) line->Draw("SAME");
      leg = new TLegend(0.45,0.15,0.6,0.35);
      leg->SetNColumns(2);
      // leg->SetTextSize(0.02);
      leg->AddEntry(ttbar_[0],"Nominal","l");
      leg->AddEntry(data_[0],"Data","pl");
      leg->AddEntry(corrections_up_[0][0],"JEC","l");
      leg->AddEntry(line_up,"Up","l");
      leg->AddEntry(corrections_up_[0][1],"XCone","l");
      leg->AddEntry(line_down,"Down","l");
      leg->SetTextSize(0.03);
      leg->Draw();
      gPad->RedrawAxis();
      A->SaveAs(save_path_general+"/Wjet_mass_sensitivity.pdf");
      delete A;
      leg->Clear();
    }
  }
}
