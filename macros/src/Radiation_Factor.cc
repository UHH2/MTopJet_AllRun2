#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

using namespace std;

int main(int argc, char* argv[]){

  TString save_path = get_save_path();

  if(argc != 2){
    cout << "Usage: ./Radiation_SYS <year>" << endl;
    return 0;
  }

  vector<vector<double>> f_fsr;

  vector<TString> masswindows = {"140to200", "140to160", "160to170", "170to180", "180to190", "190to200"};
  for(auto mass_string: masswindows){

    TFile* infile = new TFile("FSR_hists_mjet"+mass_string+".root");

    // #################################################################################################
    // Only one fit for all bins #######################################################################
    // Default -------------------------------------------------------------------
    Int_t   oldLevel  = gErrorIgnoreLevel; // Set by: gErrorIgnoreLevel = ... - functions_explain
    gErrorIgnoreLevel = kWarning;          // suppress TCanvas output
    int     bin_width = 1;                 // functions_explain
    bool    debug     = true;
    bool    rel_error = false;

    // Input ---------------------------------------------------------------------
    TString year      = argv[1];

    // #################################################################################################
    // cout settings ###################################################################################
    cout << "\n============================= General Settings\n";
    cout << "debug:      " << debug       << endl;
    cout << "rel. error: " << rel_error       << endl;
    cout << "\n============================= Calculation Settings\n";
    cout << "Bin Width:  " << bin_width   << endl << endl;

    // #################################################################################################
    // cout settings ###################################################################################
    if(debug) cout << "Getting Input Variables" << endl;

    bool is16=false; bool is17=false; bool is18=false; bool is1718=false;
    if(strcmp(year, "2016")==0)          is16   = true;
    else if(strcmp(year, "2017")==0)     is17   = true;
    else if(strcmp(year, "2018")==0)     is18   = true;
    else if(strcmp(year, "combined")==0) is1718 = true;
    else{
      cout << "" << endl;
      cout << "WHAT IS THE ************* YEAR?" << endl;
      cout << "... Sorry, had a long day, but can you please give me the correct year? (2016, 2017, 2018, combined -(17&18))" << endl;
      cout << "" << endl;
      throw runtime_error("By the way ... it is line 17");
    }

    /*
    ██████  ██ ██████  ███████  ██████ ████████  ██████  ██████  ██ ███████ ███████
    ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██      ██
    ██   ██ ██ ██████  █████   ██         ██    ██    ██ ██████  ██ █████   ███████
    ██   ██ ██ ██   ██ ██      ██         ██    ██    ██ ██   ██ ██ ██           ██
    ██████  ██ ██   ██ ███████  ██████    ██     ██████  ██   ██ ██ ███████ ███████
    */
    if(debug) cout << "Directiories" << endl;

    TString save_path = save_path+"/Plots/Radiation_SYS/mjet"+mass_string;
    save_path = creat_folder_and_path(save_path);
    save_path = creat_folder_and_path(save_path, year);
    save_path = creat_folder_and_path(save_path, "rebin"+to_string(bin_width));
    // #################################################################################################
    // creat subdirectories ############################################################################
    if(debug) cout << "Sub-Directiories" << endl;
    creat_folder(save_path, "single_bins");
    // creat_folder(save_path, "Addition_err");

    /*
    .██████  ███████ ████████     ██   ██ ██ ███████ ████████ ███████
    ██       ██         ██        ██   ██ ██ ██         ██    ██
    ██   ███ █████      ██        ███████ ██ ███████    ██    ███████
    ██    ██ ██         ██        ██   ██ ██      ██    ██         ██
    .██████  ███████    ██        ██   ██ ██ ███████    ██    ███████
    */

    TString tau32 = "comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau32";
    if(is1718) year = "2017"; // this will be changed to "combined" again

    // Difference between 2016 and 17&18
    std::vector<int> color;
    vector<TString> FSR_strings, FSRup_strings, FSRdown_strings;
    if(is16){
      color           = {kBlue};
      FSR_strings     = {"FSR_2"};
      FSRup_strings   = {"FSRup_2"};
      FSRdown_strings = {"FSRdown_2"};
    }
    if(is17 || is18 || is1718){
      color           = {kOrange-1, kBlue, kGreen+3};
      FSR_strings     = {"FSR_sqrt2", "FSR_2", "FSR_4"};
      FSRup_strings   = {"FSRup_sqrt2", "FSRup_2", "FSRup_4"};
      FSRdown_strings = {"FSRdown_sqrt2", "FSRdown_2", "FSRdown_4"};
    }

    // #################################################################################################
    // Get Background ##################################################################################
    if(debug) cout << "Background" << endl;

    TH1F *bkg;
    if(is16) bkg = (TH1F*) infile->Get("bgr_16");
    else if(is17) bkg = (TH1F*) infile->Get("bgr_17");
    else if(is18) bkg = (TH1F*) infile->Get("bgr_18");
    else if(is1718){
      bkg = (TH1F*) infile->Get("bgr_17");
      TH1F *bkg2 = (TH1F*) infile->Get("bgr_18");
      bkg->Add(bkg2);
    }
    // #################################################################################################
    // Get Data ########################################################################################
    if(debug) cout << "Data" << endl;

    TH1F *data;
    if(is16) data = (TH1F*) infile->Get("data_16");
    else if(is17) data = (TH1F*) infile->Get("data_17");
    else if(is18) data = (TH1F*) infile->Get("data_18");
    else if(is1718){
      data = (TH1F*) infile->Get("data_17");
      TH1F *data2 = (TH1F*) infile->Get("data_18");
      data->Add(data2);
    }

    // #################################################################################################
    // Get TTbar #######################################################################################
    if(debug) cout << "TTbar" << endl;

    TH1F *ttbar;
    if(is16) ttbar = (TH1F*) infile->Get("nominal_16");
    else if(is17) ttbar = (TH1F*) infile->Get("nominal_17");
    else if(is18) ttbar = (TH1F*) infile->Get("nominal_18");
    else if(is1718){
      ttbar = (TH1F*) infile->Get("nominal_17");
      TH1F *ttbar2 = (TH1F*) infile->Get("nominal_18");
      ttbar->Add(ttbar2);
    }

    // #################################################################################################
    // Get FSR #########################################################################################
    if(debug) cout << "FSR" << endl;

    vector<TH1F*> FSRup, FSRdown;
    if(is16){
      FSRup.push_back((TH1F*) infile->Get("fsrup2_16"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2_16"));
    }
    else if(is17){
      FSRup.push_back((TH1F*) infile->Get("fsrupsqrt2_17"));
      FSRup.push_back((TH1F*) infile->Get("fsrup2_17"));
      FSRup.push_back((TH1F*) infile->Get("fsrup4_17"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdownsqrt2_17"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2_17"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown4_17"));
    }
    else if(is18){
      FSRup.push_back((TH1F*) infile->Get("fsrupsqrt2_18"));
      FSRup.push_back((TH1F*) infile->Get("fsrup2_18"));
      FSRup.push_back((TH1F*) infile->Get("fsrup4_18"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdownsqrt2_18"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2_18"));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown4_18"));
    }
    else if(is1718){
      TH1F* h_upsqrt2 = (TH1F*) infile->Get("fsrupsqrt2_17");
      h_upsqrt2->Add((TH1F*) infile->Get("fsrupsqrt2_18"));
      FSRup.push_back(h_upsqrt2);
      TH1F* h_up2 = (TH1F*) infile->Get("fsrup2_17");
      h_up2->Add((TH1F*) infile->Get("fsrup2_18"));
      FSRup.push_back(h_up2);
      TH1F* h_up4 = (TH1F*) infile->Get("fsrup4_17");
      h_up4->Add((TH1F*) infile->Get("fsrup4_18"));
      FSRup.push_back(h_up4);
      TH1F* h_downsqrt2 = (TH1F*) infile->Get("fsrdownsqrt2_17");
      h_downsqrt2->Add((TH1F*) infile->Get("fsrdownsqrt2_18"));
      FSRdown.push_back(h_downsqrt2);
      TH1F* h_down2 = (TH1F*) infile->Get("fsrdown2_17");
      h_down2->Add((TH1F*) infile->Get("fsrdown2_18"));
      FSRdown.push_back(h_down2);
      TH1F* h_down4 = (TH1F*) infile->Get("fsrdown4_17");
      h_down4->Add((TH1F*) infile->Get("fsrdown4_18"));
      FSRdown.push_back(h_down4);
    }
    if(debug) cout << "Number FSR Hists: " << FSRup.size() << endl;

    // Add hists from FSR and background ---------------------------------------------------------------
    if(debug) cout << "FSR+Background" << endl;
    for(unsigned int i=0; i<FSRup.size();i++){
      if(debug) cout << "Number bins Bkg: " << bkg->GetNbinsX()      << endl;
      if(debug) cout << "Number bins FSR: " << FSRup[i]->GetNbinsX() << endl;
      FSRup[i]->Add(bkg, 1);
      FSRdown[i]->Add(bkg, 1);
    }
    // Add hists from ttbar and background -------------------------------------------------------------
    ttbar->Add(bkg);

    // #################################################################################################
    // Rebin ###########################################################################################

    // Rebin 10 ----------------------------------------------------------------------------------------
    if(debug) cout << "Rebin 10" << endl;
    TH1F *ttbar_rebin = rebin(ttbar, bin_width);
    TH1F *data_rebin  = rebin(data,  bin_width);

    vector<TH1F*> FSRup_rebin   = rebin(FSRup, bin_width);
    vector<TH1F*> FSRdown_rebin = rebin(FSRdown, bin_width);

    // For tau bin dependency
    vector<TH1F*> FSRup_rebin_norm   = normalize(FSRup_rebin);
    vector<TH1F*> FSRdown_rebin_norm = normalize(FSRdown_rebin);

    // #################################################################################################
    // Vector for compact Plots ########################################################################
    if(debug) cout << "Compact Plots" << endl;

    vector<TH1F*> ttbar_all = {ttbar, ttbar_rebin};
    vector<TH1F*> data_all  = {data, data_rebin};

    vector<vector<TH1F*>> FSRup_all = {FSRup, FSRup_rebin};
    vector<vector<TH1F*>> FSRdown_all = {FSRdown, FSRdown_rebin};

    // #################################################################################################
    // Normalize Vectors ###############################################################################
    if(debug) cout << "Normalized Vector" << endl;

    // Vectors have the same order like the ones above
    vector<TH1F*> ttbar_all_norm = normalize(ttbar_all);
    vector<TH1F*> data_all_norm  = normalize(data_all);

    vector<vector<TH1F*>> FSRup_all_norm   = normalize(FSRup_all);
    vector<vector<TH1F*>> FSRdown_all_norm = normalize(FSRdown_all);

    // #################################################################################################
    // Normalize Error ###############################################################################
    if(debug) cout << "Normalized Error" << endl;

    vector<double> ttbar_rebin_norm_err = normalize_error(ttbar_all[1]);
    vector<double> data_rebin_norm_err  = normalize_error(data_all[1]);

    vector<vector<double>> FSRup_rebin_norm_err   = normalize_error(FSRup_rebin);
    vector<vector<double>> FSRdown_rebin_norm_err = normalize_error(FSRdown_rebin);

    // #################################################################################################
    // Combine Non normalized and normalized ###########################################################
    if(debug) cout << "Norm and not norm Vectors" << endl;

    vector<vector<TH1F*>> ttbar_general = {ttbar_all, ttbar_all_norm};
    vector<vector<TH1F*>> data_general  = {data_all, data_all_norm};

    vector<vector<vector<TH1F*>>> FSRup_general   = {FSRup_all, FSRup_all_norm};
    vector<vector<vector<TH1F*>>> FSRdown_general = {FSRdown_all, FSRdown_all_norm};

    /*
    ██████  ██       ██████  ████████ ███████
    ██   ██ ██      ██    ██    ██    ██
    ██████  ██      ██    ██    ██    ███████
    ██      ██      ██    ██    ██         ██
    ██      ███████  ██████     ██    ███████
    */
    if(debug) cout << "Plots" << endl;
    if(is1718) year = "combined";


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
        data_general[a].at(b)->SetLineColor(kBlack);
      }
    }

    // #############################################################################################
    // Single Plots ################################################################################
    if(debug) cout << "Sinlge Plots" << endl;

    for(unsigned int a=0; a<FSRup_general.size();a++){ // not norm and norm (ttbar, data and variations)
      for(unsigned int b=0; b<FSRup_general.at(a).size();b++){ // not rebin and rebin (ttbar, data and variations)
        if(b==0) continue;
        for(unsigned int c=0; c<FSRup_general.at(a).at(b).size();c++){ // variation factors (sqrt_2, 2, 4)
          // FSRup_general[a].at(b).at(c)->SetTitle(FSR_strings[c]);
          FSRup_general[a].at(b).at(c)->SetTitle("");
          FSRup_general[a].at(b).at(c)->GetXaxis()->SetRangeUser(0, 1.1);
          if(a==0) FSRup_general[a].at(b).at(c)->GetYaxis()->SetRangeUser(0, FSRup_general[a].at(b).at(c)->GetMaximum()*1.2);
          else     FSRup_general[a].at(b).at(c)->GetYaxis()->SetRangeUser(0, 0.27);
          FSRup_general[a].at(b).at(c)->GetXaxis()->SetNdivisions(505);
          FSRup_general[a].at(b).at(c)->GetYaxis()->SetNdivisions(505);
          FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitleSize(0.05);
          FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitleSize(0.05);
          FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitleOffset(0.9);
          FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitleOffset(0.9);
          FSRup_general[a].at(b).at(c)->GetXaxis()->SetTitle("#tau_{32}");
          if(a==0)
          {
            FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitle("Events");
            FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitleOffset(0.9);
          }
          else     FSRup_general[a].at(b).at(c)->GetYaxis()->SetTitle("#DeltaN/N");
          FSRup_general[a].at(b).at(c)->SetLineWidth(2);
          FSRup_general[a].at(b).at(c)->SetLineColor(color[c]);

          FSRdown_general[a].at(b).at(c)->SetLineWidth(2);
          FSRdown_general[a].at(b).at(c)->SetLineStyle(2);
          FSRdown_general[a].at(b).at(c)->SetLineColor(color[c]);

          TCanvas *A = new TCanvas("A", "A", 600, 600);
          gPad->SetLeftMargin(0.15);
          gPad->SetBottomMargin(0.12);
          FSRup_general[a].at(b).at(c)->Draw("HIST");
          FSRdown_general[a].at(b).at(c)->Draw("SAME HIST");
          ttbar_general[a].at(b)->Draw("SAME HIST");
          data_general[a].at(b)->Draw("SAME P");
          leg = new TLegend(0.15,0.65,0.35,0.85);
          leg->AddEntry(ttbar_general[a].at(b),"Nominal","l");
          leg->AddEntry(FSRup_general[a].at(b).at(c) ,"FSR up","l");
          leg->AddEntry(FSRdown_general[a].at(b).at(c) ,"FSR down","l");
          leg->AddEntry(data_general[a].at(b) ,"Data","pl");
          leg->SetTextSize(0.03);
          leg->Draw();
          gPad->RedrawAxis();
          if(a==0){ // not norm
            if(b==0) A->SaveAs(save_path+"/tau32_SYS_"+FSR_strings[c]+".pdf");
            if(b==1) A->SaveAs(save_path+"/tau32_SYS_"+FSR_strings[c]+"_Rebin.pdf");
            if(is16) A->SaveAs(save_path+"/tau32_SYS_all.pdf");
          }
          if(a==1){ // norm
            if(b==0) A->SaveAs(save_path+"/tau32_SYS_"+FSR_strings[c]+"_norm.pdf");
            if(b==1) A->SaveAs(save_path+"/tau32_SYS_"+FSR_strings[c]+"_Rebin_norm.pdf");
            if(is16) A->SaveAs(save_path+"/tau32_SYS_all_norm.pdf");
          }
          delete A;
        }
      }
    }

    // #############################################################################################
    // All in one Plot #############################################################################
    cout << "" << endl;

    for(auto & hist: FSRup_general[1][0]) cout << "blubb: " << hist->GetMaximum() << endl;

    if(is17 || is18 || is1718){
      if(debug) cout << "All Plots in one" << endl;

      for(unsigned int a=0; a<FSRup_general.size();a++){ // not norm and norm
        // if(a==0) continue;
        for(unsigned int b=0; b<FSRup_general[a].size();b++){ // not rebin and rebin
          ttbar_general[a].at(b)->SetTitle("");
          ttbar_general[a].at(b)->GetXaxis()->SetRangeUser(0, 1.1);
          if(a==0) ttbar_general[a].at(b)->GetYaxis()->SetRangeUser(0, get_highest_peak(FSRup_general[a].at(b))*1.2);
          else     ttbar_general[a].at(b)->GetYaxis()->SetRangeUser(0, 0.27);
          ttbar_general[a].at(b)->GetXaxis()->SetNdivisions(505);
          ttbar_general[a].at(b)->GetYaxis()->SetNdivisions(505);
          ttbar_general[a].at(b)->GetXaxis()->SetTitleSize(0.05);
          ttbar_general[a].at(b)->GetYaxis()->SetTitleSize(0.05);
          ttbar_general[a].at(b)->GetXaxis()->SetTitleOffset(0.9);
          ttbar_general[a].at(b)->GetYaxis()->SetTitleOffset(0.9);
          ttbar_general[a].at(b)->GetXaxis()->SetTitle("#tau_{32}");
          ttbar_general[a].at(b)->GetYaxis()->SetTitle("");
          cout << "Highest Peak: " << get_highest_peak(FSRup_general[a].at(b)) << endl;

          for(unsigned int c=0; c<FSRup_general[a].at(b).size(); c++){
            if(c==0) continue;
            FSRup_general[a].at(b).at(c)->SetLineWidth(2);
            FSRup_general[a].at(b).at(c)->SetLineColor(color[c]);

            FSRdown_general[a].at(b).at(c)->SetLineWidth(2);
            FSRdown_general[a].at(b).at(c)->SetLineColor(color[c]);
          }

          if(a==0)
          {
            ttbar_general[a].at(b)->GetYaxis()->SetTitle("Events");
            ttbar_general[a].at(b)->GetYaxis()->SetTitleOffset(1.5);
          }
          else ttbar_general[a].at(b)->GetYaxis()->SetTitle("#DeltaN/N");

          // Dummy for Legend ----------------------------------------------------
          TLine *line_up  = new TLine(0.,0.,1.,1.);
          TLine *line_down= new TLine(0.,0.,1.,1.);
          line_up->SetLineStyle(1);
          line_up->SetLineWidth(2);
          line_up->SetLineColor(kGray+2);
          line_down->SetLineStyle(2);
          line_down->SetLineWidth(2);
          line_down->SetLineColor(kGray+2);
          // ---------------------------------------------------------------------

          TCanvas *A = new TCanvas("A", "A", 600, 600);
          gPad->SetLeftMargin(0.12);
          gPad->SetBottomMargin(0.12);
          ttbar_general[a].at(b)->Draw("HIST");
          data_general[a].at(b)->Draw("SAME P");
          for(unsigned int c=0; c<FSRup_general[a].at(b).size();c++){
            if(c==0) continue;
            FSRup_general[a].at(b).at(c)->Draw("SAME HIST");
            FSRdown_general[a].at(b).at(c)->Draw("SAME HIST");
          }
          leg = new TLegend(0.15,0.65,0.35,0.85);
          leg->SetNColumns(2);
          leg->AddEntry(ttbar_general[a].at(b),"Nominal","l");
          leg->AddEntry((TObject*)0,"","");
          for(unsigned int c=0; c<FSRup_general[a].at(b).size();c++){
            if(c==0) continue;
            // leg->AddEntry(FSRup_general[a].at(b).at(c),FSRup_strings[c],"l");
            else if(c==1){
              leg->AddEntry(FSRup_general[a].at(b).at(c),"FSR 2","l");
              leg->AddEntry(line_up,"Up","l");
            } else if(c==2){
              leg->AddEntry(FSRup_general[a].at(b).at(c),"FSR 4","l");
              leg->AddEntry(line_down,"Down","l");
            }
            // leg->AddEntry(FSRdown_general[a].at(b).at(c),FSRdown_strings[c],"l");
          }
          leg->SetTextSize(0.03);
          leg->Draw();
          gPad->RedrawAxis();
          if(a==0&&b==1) A->SaveAs(save_path+"/tau32_SYS_all.pdf");
          if(a==1&&b==1) A->SaveAs(save_path+"/tau32_SYS_all_norm.pdf");
          leg->Clear();
          delete A;
        }
      }
    }


    /*
    ███████ ██ ████████
    ██      ██    ██
    █████   ██    ██
    ██      ██    ██
    ██      ██    ██
    */


    /*
    **********************************************************************************************
    **************************** FIT  ************************************************************
    **********************************************************************************************
    */

    // set first fit function
    TF1 *fit = new TF1("fit", "[0] + [1]*log(x)");
    TF1 *fit_alt = new TF1("fit_alt", "[0] + [1]*log(x)");


    /*
    OVERVIEW ---------------------------------------------------------------------
    xxx - ttbar; data; FSRxx;
    FSRxx - FSRup; FSRdown;
    xxx_general = {xxx_all, xxx_all_norm}
    xxx_all     = {xxx, xxx_rebin}
    FSRxx       = {FSRxx_sqrt2, FSRxx_2, FSRxx_4}
    */

    if(debug) cout << "-----------------------------------------------------------" << endl;
    if(debug) cout << "Start: Dependency in tau32 bins" << endl;

    // At first only for Masscut_140 and Rebin10
    int n_factors;
    vector<TH1F*> ordered_fsr, ordered_fsr_norm;

    vector<vector<double>> ordered_fsr_norm_err;
    vector<vector<double>> fit_bin_values, fit_bin_values_norm; // chi2, parameters and errors --- {a, a_err, b, b_err, chi2}
    vector<double>         a, a_err, b, b_err, chi2, fit_error_estimate;
    vector<TF1*>           fit_functions, fit_functions_norm;
    int fit_value_size = 5;

    // first index: factor (order: 1/4, 1/2, 1/sqrt2, 1, sqrt2, 2, 4) --- second index: bin
    vector<double> mu, mu2, events, events_norm;
    vector<double> data_norm_bin_content, data_norm_bin_error, ttbar_norm_bin_content;
    vector<double> mu_err, events_err, events_norm_err;
    vector<double> events_rel_err, events_norm_rel_err;

    if(is16){
      n_factors = 3;
      ordered_fsr = {FSRdown_rebin[0], ttbar_all[1], FSRup_rebin[0]};
      ordered_fsr_norm = {FSRdown_rebin_norm[0], ttbar_all_norm[1], FSRup_rebin_norm[0]};
      mu = {0.5, 1, 2};
      mu_err = {0, 0, 0};
    }

    if(is17 || is18 || is1718){
      n_factors = 7;
      ordered_fsr = {FSRdown_rebin[2], FSRdown_rebin[1], FSRdown_rebin[0], ttbar_all[1], FSRup_rebin[0], FSRup_rebin[1], FSRup_rebin[2]};
      ordered_fsr_norm = {FSRdown_rebin_norm[2], FSRdown_rebin_norm[1], FSRdown_rebin_norm[0], ttbar_all_norm[1], FSRup_rebin_norm[0], FSRup_rebin_norm[1], FSRup_rebin_norm[2]};
      mu = {0.25, 0.5, 1/sqrt(2), 1, sqrt(2), 2, 4};
      mu_err = {0, 0, 0, 0, 0, 0, 0};
    }

    mu2 = square_vector(mu);
    vector<TString> tau_bin = {"0 < #tau_{32} < 0.1", "0.1 < #tau_{32} < 0.2", "0.2 < #tau_{32} < 0.3", "0.3 < #tau_{32} < 0.4", "0.4 < #tau_{32} < 0.5", "0.5 < #tau_{32} < 0.6", "0.6 < #tau_{32} < 0.7", "0.7 < #tau_{32} < 0.8", "0.8 < #tau_{32} < 0.9", "0.9 < #tau_{32} < 1"};
    vector<TString> number_bin = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    for(unsigned int factor=0; factor<n_factors; factor++) ordered_fsr_norm_err.push_back(normalize_error(ordered_fsr[factor]));

    for(unsigned int isnorm=0; isnorm<2; isnorm++){ // isnorm=0: normal --- isnorm=1: norm
      if(isnorm==0) continue;
      if(debug) cout << "Star: filling axis" << endl;

      cout << "" << endl;
      cout << "**********************************************************************************************************" << endl;
      if(isnorm==0) cout << "*********** yAxis: normal --- xAxis: log *****************************************************************" << endl;
      if(isnorm==1) cout << "*********** yAxis: norm   --- xAxis: log *****************************************************************" << endl;
      cout << "**********************************************************************************************************" << endl;
      cout << "" << endl;

      TCanvas *E = new TCanvas(); // All Graphs in one
      E->Divide(3,3);
      E->SetCanvasSize(800, 1000);
      E->SetWindowSize(800, 1000);

      // START FOR LOOP FOR FITS -------------------------------------------------
      for(unsigned int bin=3; bin < 11; bin++){ // The first two Bins are empty (for Data), therefor they are excluded
        cout << "------------------------------------------------------------------------------------------------ " << bin << endl;
        if(isnorm == 1){ // data_norm_rebin - for the chi2 Method
          data_norm_bin_content.push_back(data_all_norm[1]->GetBinContent(bin));
          data_norm_bin_error.push_back(data_rebin_norm_err[bin-1]); // vector!!!
          ttbar_norm_bin_content.push_back(ttbar_all_norm[1]->GetBinContent(bin));
        }
        for(int factor=0; factor<n_factors; factor++){
          events.push_back(ordered_fsr[factor]->GetBinContent(bin));
          events_err.push_back(ordered_fsr[factor]->GetBinError(bin));
          events_norm.push_back(ordered_fsr_norm[factor]->GetBinContent(bin));
          events_norm_err.push_back(ordered_fsr_norm_err[factor].at(bin-1));
        }

        // &vector[0] gives the vector as array (How to do it with vector<vector<xxx>> ?)
        if(debug) cout << "Start: TGraphError - Single - bin " << bin << endl;
        TGraphErrors *single_bin;

        TCanvas *A = new TCanvas("A","A", 700, 700); // each Graph seperatly

        if(isnorm==0) single_bin = new TGraphErrors(n_factors, &mu2[0], &events[0]     , &mu_err[0], &events_err[0]);
        if(isnorm==1) single_bin = new TGraphErrors(n_factors, &mu2[0], &events_norm[0], &mu_err[0], &events_norm_err[0]);

        single_bin->SetTitle(tau_bin[bin-1]);
        gPad->SetLogx();
        single_bin->GetXaxis()->SetTitle("log(x_{FSR})");
        single_bin->GetYaxis()->SetTitle("#DeltaN/N");
        if(is16){
          single_bin->GetXaxis()->SetTitleOffset(1);
          single_bin->GetYaxis()->SetTitleOffset(1.5);
        } else if(is1718){
          single_bin->GetXaxis()->SetTitleOffset(1.3);
          single_bin->GetYaxis()->SetTitleOffset(1.55);
        }
        single_bin->SetMarkerStyle(2);

        single_bin->Fit("fit", "Q");

        TGraphErrors* fit_err_68;
        if(isnorm==0) fit_err_68 = new TGraphErrors(n_factors, &mu2[0], &events[0], &mu_err[0], &events_err[0]);
        if(isnorm==1) fit_err_68 = new TGraphErrors(n_factors, &mu2[0], &events_norm[0], &mu_err[0], &events_norm_err[0]);
        TFitter* fitter = (TFitter*) TVirtualFitter::GetFitter();
        fitter->GetConfidenceIntervals(fit_err_68, 0.68);
        fit_err_68->SetFillColorAlpha(kOrange-4,0.5); // FillColor transparency

        vector<double> error_const, fit_values;
        double errm4, errm2, errms2, err0, errps2, errp2, errp4, average;
        double fxm4, fym4, fxm2, fym2, fxms2, fyms2;
        double fx0, fy0;
        double fxps2, fyps2, fxp2, fyp2, fxp4, fyp4;
        if(is16){
          fit_err_68->GetPoint(0, fxm2, fym2);
          fit_err_68->GetPoint(1, fx0, fy0);
          fit_err_68->GetPoint(2, fxp2, fyp2);

          errm2   = fit_err_68->GetErrorY(0);
          err0    = fit_err_68->GetErrorY(1);
          errp2   = fit_err_68->GetErrorY(2);

          fit_values = {fym2, fy0, fyp2};

          average = (errm2+err0+errp2)/3.;
          if(debug) cout << "Estimated Error: " << average << endl;
          error_const = {average, average, average};
          fit_error_estimate.push_back(average);
        } else{
          fit_err_68->GetPoint(0, fxm4, fym4);
          fit_err_68->GetPoint(1, fxm2, fym2);
          fit_err_68->GetPoint(2, fxms2, fyms2);
          fit_err_68->GetPoint(3, fx0, fy0);
          fit_err_68->GetPoint(4, fxps2, fyps2);
          fit_err_68->GetPoint(5, fxp2, fyp2);
          fit_err_68->GetPoint(6, fxp4, fyp4);

          errm4  = fit_err_68->GetErrorY(0);
          errm2  = fit_err_68->GetErrorY(1);
          errms2 = fit_err_68->GetErrorY(2);
          err0   = fit_err_68->GetErrorY(3);
          errps2 = fit_err_68->GetErrorY(4);
          errp2  = fit_err_68->GetErrorY(5);
          errp4  = fit_err_68->GetErrorY(6);
          if(debug){
            cout << "m4 "  << errm4 << endl;
            cout << "m2 "  << errm2 << endl;
            cout << "ms2 " << errms2 << endl;
            cout << "0 "   << err0 << endl;
            cout << "ps2 " << errps2 << endl;
            cout << "p2 "  << errp2 << endl;
            cout << "p4 "  << errp4 << endl;
          }

          fit_values  = {fym4, fym2, fyms2, fy0, fyps2, fyp2, fyp4};
          average = (errm4+errm2+errms2+err0+errps2+errp2+errp4)/7.;
          if(debug) cout << "Estimated Error: " << average << endl;
          error_const = {average, average, average, average, average, average, average};
          fit_error_estimate.push_back(average);
        }
        if(debug) for(unsigned int i =0; i<fit_values.size(); i++) cout << i << "  " << fit_values[i] << endl;
        TGraphErrors* fit_estimate;
        fit_estimate = new TGraphErrors(n_factors, &mu2[0], &fit_values[0], &mu_err[0], &error_const[0]);
        fit_estimate->SetFillColorAlpha(kGreen-7,0.5); // FillColor transparency

        // Stretch xAxis of 17, 18 and combined ----------------------------------
        // vector<double> stretch_x     = {1/20, 20};
        // vector<double> stretch_e     = {events_norm[0], events_norm[6]};
        // vector<double> stretch_x_err = {0, 0};
        // vector<double> stretch_e_err = {events_norm_err[0], events_norm_err[6]};
        // TGraphErrors* graph_stretch  = new TGraphErrors(2, &stretch_x[0], &stretch_e[0], &stretch_x_err[0], &stretch_e_err[0]);
        // graph_stretch->SetMarkerStyle(1);
        // graph_stretch->SetMarkerColor(0);
        // graph_stretch->GetHistogram()->GetXaxis()->SetLimits(1/20, 20);
        //
        // graph_stretch->SetTitle(tau_bin[bin-1]);
        // graph_stretch->GetXaxis()->SetTitle("log(#mu_{FSR})");
        // graph_stretch->GetYaxis()->SetTitle("#DeltaN/N");
        // if(is16){
        //   graph_stretch->GetXaxis()->SetTitleOffset(1);
        //   graph_stretch->GetYaxis()->SetTitleOffset(1.5);
        // } else if(is1718){
        //   graph_stretch->GetXaxis()->SetTitleOffset(1.25);
        //   graph_stretch->GetYaxis()->SetTitleOffset(1.55);
        // }

        // single_bin->SetPoint(7, 1/20, events_norm[0]);
        if(is17||is18||is1718){
          single_bin->SetPoint(7, 30, events_norm[6]);
          single_bin->SetPoint(8, 0.01, events_norm[0]);
        }

        // Draw ------------------------------------------------------------------
        if(bin==3 || bin==4 || bin==5 || bin==6)  leg = new TLegend(0.15,0.15,0.5,0.3);
        if(bin==7 || bin==8 || bin==9 || bin==10) leg = new TLegend(0.15,0.7,0.5,0.85);
        leg->AddEntry(fit_err_68,"68% CL","f");
        leg->AddEntry(fit_estimate,"Estimated","f");
        leg->SetTextSize(0.04);
        single_bin->SetMarkerStyle(2);
        // if(is1718){
        //   graph_stretch->Draw("APE");
        //   single_bin->Draw("PE same");
        // } else{
        //   single_bin->Draw("APE");
        // }
        // graph_stretch->Draw("APE");
        single_bin->Draw("APE");
        if(is17||is18||is1718) single_bin->GetXaxis()->SetLimits(0.05,21);
        // if(is16)               single_bin->GetXaxis()->SetLimits(0.2,5);

        // if(is17||is18){
        //   single_bin->GetXaxis()->SetLimits(1/20,20);
        //   single_bin->Draw("APE same");
        //   A->Update();
        // }
        fit_estimate->Draw("L3 same");
        fit_err_68->Draw("L3 same");
        single_bin->Draw("PE same");

        // -----------------------------------------------------------------------
        // Set each Plotrange individually ---------------------------------------
        if(bin==6){
          single_bin->SetMaximum(0.23); // One has to redefine the axis after drawing
          single_bin->SetMinimum(0.19);
        }else if(bin==5&&is16){
          single_bin->SetMaximum(0.215); // One has to redefine the axis after drawing
          single_bin->SetMinimum(0.164);
        }else if(bin==8&&is16){
          single_bin->SetMaximum(0.204); // One has to redefine the axis after drawing
          single_bin->SetMinimum(0.142);
        }

        leg->Draw();
        if(isnorm==0) A->SaveAs(save_path+"/single_bins/tau32_bin_"+number_bin[bin-1]+".pdf");
        if(isnorm==1) A->SaveAs(save_path+"/single_bins/tau32_bin_"+number_bin[bin-1]+"_norm.pdf");
        delete A;

        // Now plot all graphs into one file
        E->cd(bin-1);
        if(debug) cout << "Start: TGraphError - all in one - bin " << bin << endl;

        single_bin->SetTitle(tau_bin[bin-1]);
        // single_bin->SetMarkerStyle(2);
        gPad->SetLogx();
        single_bin->GetXaxis()->SetTitle("log(x_{FSR})");
        single_bin->Draw("APE");
        if(is17||is18){
          single_bin->GetXaxis()->SetLimits(1/20,20);
          single_bin->Draw("APE");
          E->Update();
        }
        fit_estimate->Draw("L3 same");
        fit_err_68->Draw("L3 same");
        single_bin->Draw("PE");
        leg->Draw();

        events = {}; events_err = {}; events_norm = {}; events_norm_err = {}; // clear vectors

        if(isnorm==1){
          /*
          Get values of the parameter of the fit
          Values du not depend on scale of ((not-)log) xAxis, but only on the scale of the yAxis ((not-)norm)
          Therefore the vectors for the Value are filled once
          */
          TF1 *fit_copy = single_bin->GetFunction("fit");
          a.push_back(fit_copy->GetParameter(0));
          a_err.push_back(fit_copy->GetParError(0));
          b.push_back(fit_copy->GetParameter(1));
          b_err.push_back(fit_copy->GetParError(1));
          chi2.push_back(fit_copy->GetChisquare());

        }

        cout << "" << endl;
      }

      if(isnorm==0) E->SaveAs(save_path+"/tau_bin_all.pdf");
      if(isnorm==1) E->SaveAs(save_path+"/tau_bin_all_norm.pdf");

    }
    if(debug) cout << "Size: " << fit_error_estimate.size() << endl;
    if(debug) for(auto number: fit_error_estimate) cout << number << endl;

    cout.precision(6);
    cout << "" << '\n' << fixed;
    cout << "|====================================================================|" << endl;
    cout << "|---------------------------- Fit Values ----------------------------|" << endl;
    cout << "|--------------------------------------------------------------------|" << endl;
    cout << "|   Bin  |     a     |   a_err   |     b     |    b_err  |    chi2   |" << endl;
    cout << "|--------|-----------|-----------|-----------|-----------|-----------|" << endl;

    int bin;
    double a_n, a_n_err, b_n, b_n_err, chi2_n;

    for(unsigned int i=0; i<8; i++){
      bin     = i+3;
      a_n     = a.at(i);
      a_n_err = a_err.at(i);
      b_n     = b.at(i);
      b_n_err = b_err.at(i);
      chi2_n  = chi2.at(i);
      if( i % 2 == 0) cout << DGRAY;
      else            cout << RESET;
      cout<<"|"<<setw(8) <<centered(to_string(bin))    <<"|"<<setw(11)<<(to_string(a_n));
      cout<<"|"<<setw(11)<<(to_string(a_n_err))<<"|"<<setw(11)<<(to_string(b_n));
      cout<<"|"<<setw(11)<<(to_string(b_n_err))<<"|"<<setw(11)<<(to_string(chi2_n))<<"|"<<endl;
    }
    cout << "|--------------------------------------------------------------------|" << endl;
    cout << "" << endl;

    /*
    .██████ ██   ██ ██ ██████
    ██      ██   ██ ██      ██
    ██      ███████ ██  █████
    ██      ██   ██ ██ ██
    .██████ ██   ██ ██ ███████
    */
    if(debug) cout << "Chi2 - Parameters" << endl;
    // At first the chi2 fit is done for the normilized case.
    // [0]=data_bin_content | [1]=a | [2]=b | [3]=error (at first only data_bin_error)

    int NParams = 4;
    vector<vector<double>> StartParameters;
    // TString chi2_formula = "(([0]-[1]-[2]*log(x))/[3])^2";

    // for 2017, 2018
    // Chi2 term for all bins
    TF1 *b3  = new TF1("b3", "(([0]-[1]-[2]*log(x))^2/[3])");
    TF1 *b4  = new TF1("b4", "(([4]-[5]-[6]*log(x))^2/[7])");
    TF1 *b5  = new TF1("b5", "(([8]-[9]-[10]*log(x))^2/[11])");
    TF1 *b6  = new TF1("b6", "(([12]-[13]-[14]*log(x))^2/[15])");
    TF1 *b7  = new TF1("b7", "(([16]-[17]-[18]*log(x))^2/[19])");
    TF1 *b8  = new TF1("b8", "(([20]-[21]-[22]*log(x))^2/[23])");
    TF1 *b9  = new TF1("b9", "(([24]-[25]-[26]*log(x))^2/[27])");
    TF1 *b10 = new TF1("b10", "(([28]-[29]-[30]*log(x))^2/[31])");

    TF1 *chi2_function;
    if(is16) chi2_function   = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 0.1, 10);
    if(is17) chi2_function   = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 10, 200);
    if(is18) chi2_function   = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 1, 50);
    if(is1718) chi2_function = new TF1("chi2_function", "b3+b4+b5+b6+b7+b8+b9+b10", 1, 50);

    // get parameters for terms
    vector<double> parameters;
    for(int bin=0; bin<8; bin++){
      // double sigma_total_square = pow(data_rebin_norm_err[bin+2], 2)+pow(ttbar_rebin_norm_err[bin+2], 2);
      if(debug) cout << "Fit estimated error: " << fit_error_estimate[bin] << endl;
      double sigma_total_square = pow(data_rebin_norm_err[bin+2], 2)+pow(fit_error_estimate[bin], 2);
      parameters.push_back(data_norm_bin_content[bin]);
      parameters.push_back(a[bin]);
      parameters.push_back(b[bin]);
      parameters.push_back(sigma_total_square);
    }

    // Set all parameters - befor 36 ... but why ?
    for(int ipar=0; ipar<32; ipar++){
      // cout << parameters[ipar] << endl;
      chi2_function->SetParameter(ipar, parameters[ipar]);
    }
    // chi2_function->SetTitle("#chi^{2} of "+year);
    chi2_function->SetTitle("");
    chi2_function->GetHistogram()->GetXaxis()->SetTitle("log(x_{FSR})");
    chi2_function->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    chi2_function->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
    chi2_function->GetHistogram()->GetYaxis()->SetTitle("#chi^{2}");
    chi2_function->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    chi2_function->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);


    if(debug) cout << "Chi2 - Function" << endl;
    TCanvas *B = new TCanvas("B","B", 600, 600);
    chi2_function->Draw();
    gPad->SetLogx();
    B->SaveAs(save_path+"/chi2.pdf");
    // Get Minimum of chi2
    if(debug) cout << "Chi2 - Minimum" << endl;
    double minY         = chi2_function->GetMinimum();
    double minX2        = chi2_function->GetX(minY, 0.01, 1000);
    double minX         = sqrt(minX2);
    double X2_sigmaup   = chi2_function->GetX(minY+1, minX2, 2000);
    double X2_sigmadown = chi2_function->GetX(minY+1, 0.01, minX2);

    double X_sigmaup    = sqrt(X2_sigmaup);
    double sigmaup      = X_sigmaup - minX;

    double X_sigmadown  = sqrt(X2_sigmadown);
    double sigmadown    = minX - X_sigmadown;

    // \u03BC = Mu
    cout << "" << endl;
    cout << " \u03BC = " << minX << " +" << sigmaup << " -" << sigmadown << endl;
    print_seperater();
    f_fsr.push_back({minX, sigmaup, sigmadown});

    fstream fsr_txt;
    fsr_txt.open(save_path+"/fsr_factor.txt", ios::out);
    fsr_txt << "ymin     : " << minY << "\nxmin     : " << minX << "\nxmin_up  :  " << sigmaup << "\nxmin_down:  " << sigmadown  << endl;
    fsr_txt.close();


    TGraph *extreme_points = new TGraph();
    extreme_points->SetPoint(0, minX2       , minY);
    extreme_points->SetPoint(1, X2_sigmaup  , minY+1);
    extreme_points->SetPoint(2, X2_sigmadown, minY+1);
    extreme_points->SetMarkerStyle(8);
    extreme_points->SetMarkerSize(0.2);
    extreme_points->SetMarkerColor(kBlue);
    // #################################################################################################
    // unified x axis interval  ########################################################################
    TF1 *chi2_function_shifted;
    if(debug) cout << "Chi2 - shifted Function" << endl;
    double upper, lower;
    if(is16){
      // upper = chi2_function->GetX(minY+200, minX2, 200);
      // lower = chi2_function->GetX(minY+200, 0.01, minX2);
      upper = minX2/10;
      lower = minX2*10;
    } else if(is17){
      upper = chi2_function->GetX(minY+10, minX2, 200);
      lower = chi2_function->GetX(minY+10, 0.01, minX2);
    } else if(is18){
      upper = chi2_function->GetX(minY+20, minX2, 200);
      lower = chi2_function->GetX(minY+20, 0.01, minX2);
    } else{
      // upper = chi2_function->GetX(minY+20, minX2, 200);
      // lower = chi2_function->GetX(minY+20, 0.01, minX2);
      upper = minX2/10;
      lower = minX2*10;
    }
    cout << lower << endl;
    cout << upper << endl;

    chi2_function_shifted = new TF1("chi2_function_shifted", "b3+b4+b5+b6+b7+b8+b9+b10", upper, lower);
    chi2_function_shifted->SetTitle("");
    // if(is16)      chi2_function_shifted = new TF1("chi2_function_shifted", "b3+b4+b5+b6+b7+b8+b9+b10", minX2-1, minX2+1);
    // else if(is17) chi2_function_shifted = new TF1("chi2_function_shifted", "b3+b4+b5+b6+b7+b8+b9+b10", minX2-40, minX2+40);
    // else if(is18) chi2_function_shifted = new TF1("chi2_function_shifted", "b3+b4+b5+b6+b7+b8+b9+b10", minX2-10, minX2+10);
    // else          chi2_function_shifted = new TF1("chi2_function_shifted", "b3+b4+b5+b6+b7+b8+b9+b10", minX2-10, minX2+10);


    // Set all parameters ------------------------------------------------------------------------------
    for(int ipar=0; ipar<32; ipar++) chi2_function_shifted->SetParameter(ipar, parameters[ipar]);
    // chi2_function_shifted->SetTitle("#chi^{2} of "+year);
    chi2_function_shifted->SetTitle("");
    chi2_function_shifted->GetHistogram()->GetXaxis()->SetTitle("log(x_{FSR})");
    chi2_function_shifted->GetHistogram()->GetXaxis()->SetTitleOffset(1.0);
    chi2_function_shifted->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
    chi2_function_shifted->GetHistogram()->GetYaxis()->SetTitle("#chi^{2}");
    chi2_function_shifted->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);
    chi2_function_shifted->GetHistogram()->GetYaxis()->SetTitleSize(0.05);

    TCanvas *T = new TCanvas("T","T", 600, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    // chi2_function_shifted->GetXaxis()->SetRangeUser(shifted_limit_down, shifted_limit_up);
    chi2_function_shifted->Draw();
    // extreme_points->Draw("P same");
    gPad->SetLogx();
    T->SaveAs(save_path+"/chi2_shifted.pdf");
    T->Clear();

    /*
    .██████  ██████  ██████  ██████  ███████ ████████ ███████ ██████      ██   ██ ██ ███████ ████████
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██ ██         ██
    ██      ██    ██ ██████  ██████  █████      ██    █████   ██   ██     ███████ ██ ███████    ██
    ██      ██    ██ ██   ██ ██   ██ ██         ██    ██      ██   ██     ██   ██ ██      ██    ██
    .██████  ██████  ██   ██ ██   ██ ███████    ██    ███████ ██████      ██   ██ ██ ███████    ██
    */
    if(debug) cout << "New Histogram with correct Minimum from Chi2" << endl;
    cout << '\n';
    /*
    Data is described mainly through the FSRup(_4) variation. Therefore, the calculated variation
    is only used for the up variation. The error for the calculated factor is therefore used for the up variation.
    */

    vector<double> FSRup_min = {0, 0}; vector<double> FSRup_min_up_err = {0, 0}; vector<double> FSRup_min_down_err = {0, 0};
    for(unsigned int bin=0; bin<8;bin++){ // Starting at Bin 3. For Completeness the first two bins are set to 0.
      double a_par = a[bin];
      double b_par = b[bin];
      double sigmaup_value   = a_par + b_par * log(X2_sigmaup); // pow(minX*sqrt(2),2)
      double sigmadown_value = a_par + b_par * log(X2_sigmadown); // pow(minX*(1/sqrt(2)),2)
      double nominal_value   = a_par + b_par * log(minX2);
      FSRup_min.push_back(nominal_value);
      if(sigmaup_value<nominal_value && nominal_value<sigmadown_value){
        if(debug) cout << "Bin " << bin << ": up(down) is down(up)" << endl;
        FSRup_min_up_err.push_back(abs(FSRup_min[bin+2]-sigmadown_value)); // Here and next line: (FSRup_min_up & _down)
        FSRup_min_down_err.push_back(abs(FSRup_min[bin+2]-sigmaup_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else if(sigmadown_value<nominal_value && nominal_value<sigmaup_value){
        if(debug) cout << "Bin " << bin << ": up(down) is up(down)" << endl;
        FSRup_min_up_err.push_back(abs(FSRup_min[bin+2]-sigmaup_value)); // Here and next line: (FSRup_min_up & _down)
        FSRup_min_down_err.push_back(abs(FSRup_min[bin+2]-sigmadown_value)); // take minus to get error for TGraphAsymmErrors (below)
      }
      else throw runtime_error("Something is wrong with the bin value for the FSR variation using the calculated minimum. both variations are smaller or greater than the nominal value");
    }


    TH1F *h_FSRup_min = new TH1F("FSRup_min", "", 10, 0, 1);
    for(unsigned int bin=1; bin<11;bin++) h_FSRup_min->SetBinContent(bin, FSRup_min[bin-1]);

    // Integral-----------------------------------------------------------------------------------------------
    // double fsr_min_integral = h_FSRup_min->Integral();
    // cout << "Integral of calculated \u03C4_{32} distribution: " << fsr_min_integral << endl;
    cout << h_FSRup_min->GetMaximum()     << endl;
    cout << h_FSRup_min->GetMaximum()*1.5 << endl;

    h_FSRup_min->SetTitle("");
    h_FSRup_min->GetXaxis()->SetRangeUser(0, 1);
    h_FSRup_min->GetYaxis()->SetRangeUser(0, h_FSRup_min->GetMaximum()*1.2);
    h_FSRup_min->GetXaxis()->SetNdivisions(505);
    h_FSRup_min->GetYaxis()->SetNdivisions(505);
    h_FSRup_min->GetXaxis()->SetTitleSize(0.05);
    h_FSRup_min->GetYaxis()->SetTitleSize(0.05);
    h_FSRup_min->GetXaxis()->SetTitleOffset(0.9);
    h_FSRup_min->GetYaxis()->SetTitleOffset(0.9);
    h_FSRup_min->GetXaxis()->SetTitle("#tau_{32}");
    h_FSRup_min->GetYaxis()->SetTitle("#DeltaN/N");
    h_FSRup_min->SetLineWidth(1);
    h_FSRup_min->SetLineColor(kBlack);

    if(is17 || is18 || is1718) FSRdown_general[1].at(1).at(2)->SetLineStyle(2);
    if(is16) FSRdown_general[1].at(1).at(0)->SetLineStyle(2);

    vector<double> xbins     = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
    vector<double> xbins_err = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    auto errors = new TGraphAsymmErrors(10, &xbins[0], &FSRup_min[0], &xbins_err[0], &xbins_err[0], &FSRup_min_down_err[0], &FSRup_min_up_err[0]);
    errors->SetFillColor(kGray);
    errors->SetFillStyle(1001); // 1001: solid - same as nothing
    errors->SetTitle(year);

    if(debug) cout << "Draw New Histogram with correct Minimum from Chi2" << endl;
    TCanvas *C = new TCanvas("C", "C", 600, 600);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    h_FSRup_min->Draw("HIST");
    errors->Draw("SAME 2");
    h_FSRup_min->Draw("SAME HIST"); // Draw again, because it will be over over otherwise
    ttbar_all_norm[1]->Draw("SAME HIST");
    data_all_norm[1]->Draw("SAME P");

    leg = new TLegend(0.15,0.65,0.35,0.85);
    leg->AddEntry(ttbar_all_norm[1],"Nominal","l");
    leg->AddEntry(h_FSRup_min, "FSR best fit value","l");
    leg->AddEntry(data_all_norm[1], "Data","pl");
    leg->SetTextSize(0.03);
    leg->Draw();
    gPad->RedrawAxis();
    C->SaveAs(save_path+"/tau32_minimum_factor.pdf");
    // C->SaveAs(save_path+year+"/Addition_err/tau32_minimum_factor_err_factorsqrt2.pdf");
    if(is17 || is18 || is1718) {
      FSRup_general[1].at(1).at(1)->Draw("SAME HIST");
      // FSRdown_general[1].at(1).at(2)->Draw("SAME HIST");
      leg->AddEntry(FSRup_general[1].at(1).at(1), "FSRup 2","l");
      // leg->AddEntry(FSRdown_general[1].at(1).at(2), "FSRdown 4","l");
    }
    if(is16) {
      FSRup_general[1].at(1).at(0)->Draw("SAME HIST");
      // FSRdown_general[1].at(1).at(0)->Draw("SAME HIST");
      leg->AddEntry(FSRup_general[1].at(1).at(0), "FSRup 2","l");
      // leg->AddEntry(FSRdown_general[1].at(1).at(0), "FSRdown 2","l");
    }
    C->SaveAs(save_path+"/tau32_minimum_factor_compare.pdf");
    // C->SaveAs(save_path+year+"/Addition_err/tau32_minimum_factor_compare_err_factorsqrt2.pdf");
    leg->Clear();
    delete C;
  }

  vector<double> xvalues, xup, xdown, central, up, down;
  for(unsigned int i=0; i<f_fsr.size(); i++){
    central.push_back(f_fsr[i][0]);
    up.push_back(f_fsr[i][1]);
    down.push_back(f_fsr[i][2]);
    xvalues.push_back(i+1);
    xup.push_back(0.);
    xdown.push_back(0.);
  }

  TString year = argv[1];
  double ymin = 0.0;
  double ymax = 2.0;
  if(year != "2016"){
    ymin = 0.0;
    ymax = 10.0;
  }

  TH1F* dummy = new TH1F("dummy", "dummy", f_fsr.size(), 0.5, 0.5+f_fsr.size());
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  for(unsigned int i=0; i<masswindows.size(); i++) dummy->GetXaxis()->SetBinLabel(i+1, masswindows[i]);
  dummy->SetTitle("");
  dummy->GetYaxis()->SetTitle("#it{f}^{FSR}");
  dummy->GetYaxis()->SetRangeUser(ymin, ymax);
  dummy->Draw();
  TGraphAsymmErrors * results = new TGraphAsymmErrors(f_fsr.size(), &xvalues[0], &central[0], &xdown[0], &xup[0], &down[0], &up[0]);
  results->SetMarkerStyle(8);
  results->Draw("P");
  c->SaveAs(save_path+"/Plots/Radiation_SYS/Results_"+year+".pdf");

  return 0;

}
