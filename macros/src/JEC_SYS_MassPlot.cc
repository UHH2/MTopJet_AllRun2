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
  Explanation:
  1) This script loops many times over bins. The bin number starts at 1 (e.g. for the function GetBinContent)
  On the other hand, the bin content or titles are stored in a vector. The vector starts at 0.
  Mostly, the loops start at bin=0. Therefor, using the function GetBinContent the bin needs to be bin+1
  */

  // Changes of how many pt bins to use - Search for CHANGE_PT
  // #################################################################################################
  // Declare different variables used in the code ####################################################
  for(int ptbin=0; ptbin<4; ptbin++){ // CHANGE_PT
    print_seperater();

    TString reconst   = "btag";    // match btag_cut btag_sel compare min_mass
    TString MC_uncert = "central"; // central avg nominal
    cout.precision(6);

    if(argc != 3){
      cout << "\n" << "Usage: ./JEC_SYS <year> <rebin_number>\n";
      cout << "rebin_number:        hist->Rebin(rebin_number)\n";
      return 0;
    }

    // #################################################################################################
    // Only one fit for all bins #######################################################################
    bool debug       = false;

    TString addition="";
    if(ptbin==0) addition="_hh"; // CHANGE_PT
    if(ptbin==1) addition="_hl";
    if(ptbin==2) addition="_lh";
    if(ptbin==3) addition="_ll";
    // #################################################################################################
    // cout setting ####################################################################################
    cout << "debug:   " << debug       << endl;

    // #################################################################################################
    // Print Options ###################################################################################
    Int_t oldLevel = gErrorIgnoreLevel;
    // Set by: gErrorIgnoreLevel = ...

    // Rebin -------------------------------------------------------------------------------------------
    /*
    Right now every Data bin will be excluded, if the BinContent is zero. This leads to a problem,
    since in the first bins (at a binning of 1 GeV) the bin Content is (e.g.) 000001000.
    */
    int bin_width = atoi(argv[2]); // atoi() c-string representing integer into int
    int number_bins = 180/bin_width;
    TString str_number_bins = to_string(number_bins); // For creating folders

    // Year --------------------------------------------------------------------------------------------
    if(debug) cout << "Getting Input Variables" << endl;

    TString year = argv[1];
    bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
    if     (strcmp(year, "2016")==0) is16   = true;
    else if(strcmp(year, "2017")==0) is17   = true;
    else if(strcmp(year, "2018")==0) is18   = true;
    else if(strcmp(year,  "all")==0) isAll  = true;
    else throw runtime_error("Give me the correct year please (2016, 2017, 2018 or all)");

    vector<TString> all_years;
    if(is16 || is17 || is18) all_years = {year};
    else if(isAll) all_years = {"2016", "2017", "2018"};
    else throw runtime_error("Give me the correct year please (2016, 2017, 2018, all or 1718)");

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

    TString save_path_general = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS"; // CHANGE_PT
    save_path_general = creat_folder_and_path(save_path_general, "pt_bins");

    // #################################################################################################
    // creat subdirectories ############################################################################

    save_path_general = creat_folder_and_path(save_path_general, year);
    save_path_general = creat_folder_and_path(save_path_general, reconst);
    save_path_general = creat_folder_and_path(save_path_general, "rebin"+str_number_bins);


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
    // if(reconst=="btag")      w_mass = hist_class+"wjet_pt_match_S1divW_W";              // mass from subjet matched with ak4 with highest_btag close to the xcone fathadjet



    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      vector<TFile*> file_bkg_v;
      vector<TH1F*> hists_bkg_v;
      vector<TString> path_bkg_v = {"uhh2.AnalysisModuleRunner.MC.other.root", "uhh2.AnalysisModuleRunner.MC.SingleTop.root", "uhh2.AnalysisModuleRunner.MC.WJets.root"};
      for(unsigned int i=0; i<path_bkg_v.size(); i++) file_bkg_v.push_back(new TFile(dir+all_years[iyear]+"/muon/"+path_bkg_v[i]));
      for(unsigned int i=0; i<file_bkg_v.size(); i++) hists_bkg_v.push_back((TH1F*)file_bkg_v[i]->Get(w_mass));
      TH1F *bkg = AddHists(hists_bkg_v, 1);
      allyears_bkg.push_back(bkg);
    }

    TH1F *bkg = AddHists(allyears_bkg, 1);

    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TString data_path     = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *data_file      = new TFile(dir+all_years[iyear]+"/muon/"+data_path);
      TH1F  *data           = (TH1F*)data_file->Get(w_mass);
      allyears_data.push_back(data);
    }

    TH1F *data            = AddHists(allyears_data, 1);
    TH1F *data_rebin      = rebin(data, bin_width);
    TH1F *data_norm       = normalize(data);
    TH1F *data_rebin_norm = normalize(data_rebin);

    vector<double> data_rebin_norm_err = normalize_error(data_rebin);
    cout<< data_rebin_norm->GetNbinsX()<<"   "<< data_rebin_norm_err.size()<< endl;
    for(int ipar=0; ipar<data_rebin_norm->GetNbinsX(); ipar++) data_rebin_norm->SetBinError(ipar+1, data_rebin_norm_err[ipar]);

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TString ttbar_path = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *ttbar_file  = new TFile(dir+all_years[iyear]+"/muon/"+ttbar_path);
      TH1F  *ttbar       = (TH1F*)ttbar_file->Get(w_mass);
      allyears_ttbar.push_back(ttbar);
    }

    TH1F  *ttbar           = AddHists(allyears_ttbar, 1);
    ttbar->Add(bkg, 1);
    TH1F* ttbar_rebin      = rebin(ttbar, bin_width);
    TH1F* ttbar_norm       = normalize(ttbar);
    TH1F* ttbar_rebin_norm = normalize(ttbar_rebin);

    vector<double> ttbar_rebin_norm_err = normalize_error(ttbar_rebin);

    // #################################################################################################
    // Get SYS #########################################################################################
    if(debug) cout << "JEC" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *JECup_file   = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECup.root");
      TFile *JECdown_file = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_JECdown.root");
      TH1F  *JECup        = (TH1F*)JECup_file->Get(w_mass);
      TH1F  *JECdown      = (TH1F*)JECdown_file->Get(w_mass);
      allyears_JECup.push_back(JECup);
      allyears_JECdown.push_back(JECdown);
    }

    TH1F *JECup   = AddHists(allyears_JECup, 1);
    TH1F *JECdown = AddHists(allyears_JECdown, 1);

    // Add hists from JEC and background ---------------------------------------------------------------
    if(debug) cout << "JEC+Background" << endl;

    JECup->Add(bkg, 1);
    JECdown->Add(bkg, 1);

    TH1F* JECup_norm         = normalize(JECup);
    TH1F* JECdown_norm       = normalize(JECdown);
    TH1F* JECup_rebin        = rebin(JECup, bin_width);
    TH1F* JECdown_rebin      = rebin(JECdown, bin_width);
    TH1F* JECup_rebin_norm   = rebin(JECup_norm, bin_width);
    TH1F* JECdown_rebin_norm = rebin(JECdown_norm, bin_width);

    vector<double> JECup_rebin_norm_err   = normalize_error(JECup_rebin);
    vector<double> JECdown_rebin_norm_err = normalize_error(JECdown_rebin);

    // ------------------------------------------------------------------------------------------------
    if(debug) cout << "XCone" << endl;

    for(unsigned int iyear=0; iyear<all_years.size(); iyear++){
      TFile *XCup_file   = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConeup.root");
      TFile *XCdown_file = new TFile(dir+all_years[iyear]+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar_XConedown.root");
      TH1F  *XCup        = (TH1F*)XCup_file->Get(w_mass);
      TH1F  *XCdown      = (TH1F*)XCdown_file->Get(w_mass);
      allyears_XCup.push_back(XCup);
      allyears_XCdown.push_back(XCdown);
    }

    TH1F *XConeup   = AddHists(allyears_XCup, 1);
    TH1F *XConedown = AddHists(allyears_XCdown, 1);

    // Add hists from FSR and background ---------------------------------------------------------------
    if(debug) cout << "XCone+Background" << endl;

    XConeup->Add(bkg, 1);
    XConedown->Add(bkg, 1);

    TH1F* XConeup_norm         = normalize(XConeup);
    TH1F* XConedown_norm       = normalize(XConedown);
    TH1F* XConeup_rebin        = rebin(XConeup, bin_width);
    TH1F* XConedown_rebin      = rebin(XConedown, bin_width);
    TH1F* XConeup_rebin_norm   = rebin(XConeup_norm, bin_width);
    TH1F* XConedown_rebin_norm = rebin(XConedown_norm, bin_width);

    vector<double> XConeup_rebin_norm_err   = normalize_error(XConeup_rebin);
    vector<double> XConedown_rebin_norm_err = normalize_error(XConedown_rebin);

    vector<TH1F*> corrections_up_norm = {JECup_rebin_norm, XConeup_rebin_norm};
    vector<TH1F*> corrections_down_norm = {JECdown_rebin_norm, XConedown_rebin_norm};

    // #################################################################################################
    // Empty Data Bins #################################################################################
    if(debug) cout << "Empty Bins Data" << endl;

    vector<int> empty_bins_v;
    for(int bin=0; bin < number_bins; bin++){
      if(!(abs(data_rebin_norm->GetBinContent(bin+1))>0)){
        empty_bins_v.push_back(bin+1);
        if(debug) cout << "Empty Bins: " << empty_bins_v[bin] << "\n";
      }
    }
    int number_empty_bins = empty_bins_v.size();
    if(debug) cout << "Number Empty Bins: " << number_empty_bins << "\n";

    // #################################################################################################
    // Masspeack #######################################################################################
    if(debug) cout << "Masspeak Bins" << endl;

    vector<double> PeakLimit;
    int Limit = 75;


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
      ttbar_[norm]->SetTitle("");
      ttbar_[norm]->GetXaxis()->SetRangeUser(0, 180);
      ttbar_[norm]->GetYaxis()->SetRangeUser(0, ttbar_[norm]->GetMaximum()*1.2);
      ttbar_[norm]->GetXaxis()->SetNdivisions(505);
      ttbar_[norm]->GetYaxis()->SetNdivisions(505);
      ttbar_[norm]->GetXaxis()->SetTitleSize(0.05);
      ttbar_[norm]->GetYaxis()->SetTitleSize(0.04);
      ttbar_[norm]->GetXaxis()->SetTitleOffset(0.9);
      ttbar_[norm]->GetYaxis()->SetTitleOffset(1.5);
      ttbar_[norm]->GetXaxis()->SetTitle("m_{Wjet}");
      ttbar_[norm]->GetYaxis()->SetTitle("");
      ttbar_[norm]->SetLineWidth(2);  // ttbar hist style
      ttbar_[norm]->SetLineColor(kRed);

      // Data --------------------------------------------------------------------------------------------
      data_[norm]->SetMarkerStyle(8);  // data hist style
      data_[norm]->SetMarkerColor(kBlack);
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
    cout << Limit << endl;
    TLine *line = new TLine(0, Limit, 180, Limit);
    line->SetLineColor(kGray);
    line->SetLineWidth(1);

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
        corrections_down_[norm][correction]->SetLineColor(kBlue);

        TCanvas *A = new TCanvas("A", "A", 1000, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);
        ttbar_[norm]->Draw("HIST");
        corrections_down_[norm][correction]->Draw("SAME HIST");
        corrections_up_[norm][correction]->Draw("SAME HIST");
        data_[norm]->Draw("SAME P");
        if(norm==0 && bin_width==1 && isAll) line->Draw("SAME");
        leg = new TLegend(0.6,0.65,0.8,0.85);
        leg->SetTextSize(0.2);
        leg->AddEntry(ttbar_[norm],"nominal","l");
        if(correction==0){
          leg->AddEntry(corrections_up_[norm][correction],"JECup","l");
          leg->AddEntry(corrections_down_[norm][correction],"JECdown","l");
        }
        if(correction==1){
          leg->AddEntry(corrections_up_[norm][correction],"XConeup","l");
          leg->AddEntry(corrections_down_[norm][correction],"XConedown","l");
        }
        leg->AddEntry(data_[norm],"Data","p");
        leg->SetTextSize(0.05);
        leg->Draw();
        gPad->RedrawAxis();
        if(correction==0) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_JEC"+addition+norm_str+".pdf");
        if(correction==1) A->SaveAs(save_path_general+"/Wjet_mass_sensitivity_XCone"+addition+norm_str+".pdf");
        delete A;
        leg->Clear();
      }
    }

    /*
    ██████  ██████
    .    ██ ██   ██
    .█████  ██   ██
    ██      ██   ██
    ███████ ██████
    */

    // TLine *horizontal = new TLine(0, 0.7, 1000, 0.7);
    // horizontal->SetLineColor(kRed);
    // horizontal->SetLineWidth(1);
    //
    // TLine *vertical   = new TLine(300, 0.5, 300, 1);
    // vertical->SetLineColor(kRed);
    // vertical->SetLineWidth(1);
    //
    // cout << ttbar->Integral() << endl;
    // cout << ttbar->Integral(1, 30, 1, 70) << endl; 7287.97
    // cout << ttbar->Integral(31, 101, 0, 70) << endl; 7357.99
    // cout << ttbar->Integral(1, 30, 71, 101) << endl; 6767.4
    // cout << ttbar->Integral(31, 101, 71, 101) << endl; 8885.25


    // ttbar->GetYaxis()->SetRangeUser(0.5, 1);
    // gStyle->SetPadTickY(1);
    // gStyle->SetPadTickX(1);
    // gStyle->SetOptStat(kFALSE);
    // gStyle->SetLegendBorderSize(0);
    // /*
    // Plotting the histograms included above.
    // */
    // TCanvas *A = new TCanvas("A", "A", 1000, 600);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetBottomMargin(0.12);
    // ttbar->Draw("Colz");
    // horizontal->Draw("same");
    // vertical->Draw("same");
    // A->SaveAs(save_path_general+"/S1divW.pdf");
    // delete A;

  }
}
