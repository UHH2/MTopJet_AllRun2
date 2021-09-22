#include "../include/CentralInclude.h"

using namespace std;

void PlotRatio(TGraphErrors* data_, TF1* fit_, TGraph* uncert_, int etabin, TString year_);


// class to transform fit parameters into uncorrelated ones
class TransformParameters{
public:
  TransformParameters(TFitResultPtr fitfunction, TString oldformula, int n_param);
  TString getNewFunction();
private:
  TString newFormula;
};

// class to vary correction
class Variation{
public:
  Variation(vector<double> params, vector<double> errors_up, vector<double> errors_down, TString function);
  vector<TF1*> getVariedFits() { return VariedFits; };
  TGraph *CalculateEnvelope(TF1* central, string updown);

private:
  vector<TF1*> VariedFits;
};
////

vector<double> CalculateParameterUncertainty(TGraphErrors* data_, TString function, TFitResultPtr fitresult, vector<double> parameters, TString updown);
TGraph *UncertaintyFromDiff(TF1* central, TF1* pol1, TF1* pol2, TString mode);

int main(int argc, char* argv[]){

  TString save_path = get_save_path();

  /*
  ██████  ███████ ███████ ██ ███    ██ ███████     ███████ ██ ████████     ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██
  ██   ██ ██      ██      ██ ████   ██ ██          ██      ██    ██        ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██
  ██   ██ █████   █████   ██ ██ ██  ██ █████       █████   ██    ██        █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██
  ██   ██ ██      ██      ██ ██  ██ ██ ██          ██      ██    ██        ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██
  ██████  ███████ ██      ██ ██   ████ ███████     ██      ██    ██        ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████
  */

  bool use_median = false;
  bool show_fit = true;

  int Nparams = 3; // REMEMBER: pol3 has 4 parameters!
  TString fitformula;
  if( Nparams == 1)      fitformula = "[0]";
  else if( Nparams == 2) fitformula = "[0] + [1]*x";
  else if( Nparams == 3) fitformula = "[0] + [1]*log(x) + [2]*log(x)*log(x)";
  else if( Nparams == 4) fitformula = "[0] + [1]*x + [2]*x*x + [3]*x*x*x";
  else if( Nparams == 5) fitformula = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x";
  else if( Nparams == 6) fitformula = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x";
  bool ConstrainParams = true;
  vector<vector<double>> StartParameters;
  StartParameters.push_back({0.9, 4.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 4.0E-4, -8.0E-7, 8.0E-13, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 3.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 8.5E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 5.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 6.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 2.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 3.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.9, 3.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});
  StartParameters.push_back({0.7, 5.0E-3, -5.0E-5, 3.0E-7, -6.0E-10, -1.7E-14});
  StartParameters.push_back({0.9, 5.0E-4, -8.0E-7,  8.0E-8, 1.0E-13, -1.7E-14});
  // StartParameters.push_back({0.9, 4.0E-4, -8.0E-7, 8.0E-10, 1.0E-13, -1.7E-14});


  // define in what range the parameters can be varied (upper limit = parameter + parameter * RangeFraction)
  // there is one value for avery eta bin
  vector<double> RangeFraction {2,2,2,2,2,2};
  ////


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  bool year_16 = false;
  bool year_17 = false;
  bool year_18 = false;
  TString filename, year;
  if(argc > 1 && strcmp(argv[1], "2016") == 0) {
    year = "2016";
    filename = dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2016v3.root";
  }
  if(argc > 1 && strcmp(argv[1], "2017") == 0){
    year = "2017";
    filename = dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2017v2.root";
  }
  if(argc > 1 && strcmp(argv[1], "2018") == 0){
    year = "2018";
    filename = dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2018.root";
  }

  string year_s = (string) year;

  TFile *file = new TFile(filename);
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get plots ------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  Float_t xbins[] = {0, 50, 70, 90, 120, 150, 180, 250, 350, 500};
  Float_t ybins[] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139};

  int no_ptbins = 9; // rows
  int no_etabins = 11; // columns

  TH2F * h_ratio_mean = new TH2F("h_ratio_mean","Mean",no_ptbins, xbins, no_etabins, ybins);
  TH2F * h_ratio_mean_err = new TH2F("h_ratio_mean_err","Mean Error",no_ptbins, xbins, no_etabins, ybins);
  TH2F * h_ptrec_mean = new TH2F("h_ptrec_mean","Mean",no_ptbins, xbins, no_etabins, ybins);
  TH2F * h_ptrec_mean_err = new TH2F("h_ptrec_mean_err","Mean Error",no_ptbins, xbins, no_etabins, ybins);
  TH2F * h_count = new TH2F("h_count","Event Count",no_ptbins, xbins, no_etabins, ybins);




  // get histograms from root file
  std::string histdir = "CorrectionHists/";
  std::string name_ratio = "PtReso_";
  std::string name_ptrec = "PtRec_";
  std::string name_count = "Count_";


  std::string filename_ratio, filename_ptrec, filename_count;
  TH1F * ratio[no_ptbins][no_etabins];
  TH1F * ptrec[no_ptbins][no_etabins];
  TH1F * count[no_ptbins][no_etabins];

  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      filename_ratio = histdir + name_ratio + std::to_string(pt_bin) + std::to_string(eta_bin) ;
      filename_ptrec = histdir + name_ptrec + std::to_string(pt_bin) + std::to_string(eta_bin) ;
      filename_count = histdir + name_count + std::to_string(pt_bin) + std::to_string(eta_bin) ;

      ratio[pt_bin][eta_bin] = (TH1F*)file->Get(filename_ratio.c_str());
      ptrec[pt_bin][eta_bin] = (TH1F*)file->Get(filename_ptrec.c_str());
      count[pt_bin][eta_bin] = (TH1F*)file->Get(filename_count.c_str());
    }
  }
  ////

  // get mean values and write to array
  double ratio_mean[no_ptbins][no_etabins], ratio_mean_err[no_ptbins][no_etabins];
  double ratio_rms[no_ptbins][no_etabins], ratio_rms_err[no_ptbins][no_etabins];
  double ptrec_mean[no_ptbins][no_etabins], ptrec_mean_err[no_ptbins][no_etabins];
  double ptrec_rms[no_ptbins][no_etabins], ptrec_rms_err[no_ptbins][no_etabins];
  double factor[no_ptbins][no_etabins], factor_err[no_ptbins][no_etabins];
  double n_events[no_ptbins][no_etabins];
  TF1 *ratio_fit[no_ptbins][no_etabins];
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){

      // get mean (ptrec/ptgen) with gaussian fit
      cout << "Get ptrec/ptgen bin pt" << pt_bin << " eta" << eta_bin << endl;
      double upper = ratio[pt_bin][eta_bin]->GetMean() + ratio[pt_bin][eta_bin]->GetRMS();
      double lower = ratio[pt_bin][eta_bin]->GetMean() - ratio[pt_bin][eta_bin]->GetRMS();
      TF1*ratio_func = new TF1("ratio_func","gaus",lower,upper);
      ratio[pt_bin][eta_bin]->Fit("ratio_func","RQ");
      ratio_fit[pt_bin][eta_bin] = ratio[pt_bin][eta_bin]->GetFunction("ratio_func");
      double gaussmean = ratio_fit[pt_bin][eta_bin]->GetParameter(1);
      // GET MEDIAN
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      ratio_fit[pt_bin][eta_bin]->GetQuantiles(1, values, quantile);
      double median = values[0];
      ////
      if(!use_median) ratio_mean[pt_bin][eta_bin] = gaussmean;
      else            ratio_mean[pt_bin][eta_bin] = median;
      ratio_mean_err[pt_bin][eta_bin] = ratio_fit[pt_bin][eta_bin]->GetParError(1);

      if(use_median){
        cout << "------------------------------------------------------------" << endl;
        cout << "PT bin: " << pt_bin << "  , ETA bin: " << eta_bin << endl;
        cout << "Gaussian Mean     = " << gaussmean << endl;
        cout << "Median            = " << median << endl;
        cout << "Relative Diff [%] = " << (median-gaussmean)/median * 100<< endl;
      }

      // get correction factor (f = 1/mean(ptrec/ptgen))
      factor[pt_bin][eta_bin] = 1/ratio_mean[pt_bin][eta_bin];
      factor_err[pt_bin][eta_bin] = abs(-ratio_mean_err[pt_bin][eta_bin]/(ratio_mean[pt_bin][eta_bin] * ratio_mean[pt_bin][eta_bin]));

      // get mean from ptrec with RootMean (because ptrec is only needed to get bin center in fit right)
      cout << "Get ptrec bin pt" << pt_bin << " eta" << eta_bin << endl;
      double upper_ = ptrec[pt_bin][eta_bin]->GetMean() + 2*ptrec[pt_bin][eta_bin]->GetRMS();
      double lower_ = ptrec[pt_bin][eta_bin]->GetMean() - 2*ptrec[pt_bin][eta_bin]->GetRMS();
      cout << "Mean before = " << ptrec[pt_bin][eta_bin]->GetMean() << endl;
      cout << "lower = " << lower_ << endl;
      cout << "upper = " << upper_ << endl;
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetRangeUser(lower_, upper_);
      ptrec_mean[pt_bin][eta_bin] = ptrec[pt_bin][eta_bin]->GetMean();
      ptrec_mean_err[pt_bin][eta_bin] = ptrec[pt_bin][eta_bin]->GetMeanError();
      cout << "Mean after = " << ptrec_mean[pt_bin][eta_bin] << endl;
      cout << "-------------------------------" << endl;

      // get event number in each bin
      n_events[pt_bin][eta_bin] = count[pt_bin][eta_bin]->Integral();
    }
  }
  ////

  // fill Hists with values from arrays
  TAxis *xaxis_ratio = h_ratio_mean->GetXaxis();
  TAxis *yaxis_ratio = h_ratio_mean->GetYaxis();
  TAxis *xaxis_ratio_err = h_ratio_mean_err->GetXaxis();
  TAxis *yaxis_ratio_err = h_ratio_mean_err->GetYaxis();
  TAxis *xaxis_ptrec = h_ptrec_mean->GetXaxis();
  TAxis *yaxis_ptrec = h_ptrec_mean->GetYaxis();
  TAxis *xaxis_ptrec_err = h_ptrec_mean_err->GetXaxis();
  TAxis *yaxis_ptrec_err = h_ptrec_mean_err->GetYaxis();
  TAxis *xaxis_count = h_count->GetXaxis();
  TAxis *yaxis_count = h_count->GetYaxis();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      h_ratio_mean->Fill(xaxis_ratio->GetBinCenter(pt_bin+1), yaxis_ratio->GetBinCenter(eta_bin+1), ratio_mean[pt_bin][eta_bin]);
      h_ratio_mean_err->Fill(xaxis_ratio->GetBinCenter(pt_bin+1), yaxis_ratio->GetBinCenter(eta_bin+1), ratio_mean_err[pt_bin][eta_bin]);
      h_ptrec_mean->Fill(xaxis_ptrec->GetBinCenter(pt_bin+1), yaxis_ptrec->GetBinCenter(eta_bin+1), ptrec_mean[pt_bin][eta_bin]);
      h_ptrec_mean_err->Fill(xaxis_ptrec->GetBinCenter(pt_bin+1), yaxis_ptrec->GetBinCenter(eta_bin+1), ptrec_mean_err[pt_bin][eta_bin]);
      h_count->Fill(xaxis_count->GetBinCenter(pt_bin+1), yaxis_count->GetBinCenter(eta_bin+1), n_events[pt_bin][eta_bin]);
    }
  }
  ////

  // now create graph Correction_factor vs <ptrec> in every eta bin
  Double_t g_factor[no_ptbins], g_factor_err[no_ptbins];
  Double_t g_ptrec[no_ptbins], g_ptrec_err[no_ptbins];
  TGraphErrors* factor_pt[no_etabins];
  TFitResultPtr OldFitResult[no_etabins];  // for first fit
  TFitResultPtr TempFitResult[no_etabins]; // for fit with uncorrelated Parameters
  TFitResultPtr NewFitResult[no_etabins]; // for fit with uncorrelated Parameters
  TFitResultPtr HighPol3FitResult[no_etabins];
  TFitResultPtr HighPol4FitResult[no_etabins];
  TF1* TempFit[no_etabins];
  TF1* NewFit[no_etabins];
  TF1* HighPol3Fit[no_etabins];
  TF1* HighPol4Fit[no_etabins];
  TString TempFunction[no_etabins];
  TString NewFunction[no_etabins];
  vector< vector<double> > params, errors, errors_up, errors_down;
  vector<double> chi2; // store chi2 value for every fit
  vector<int> ndf;     // store number of degrees of freedom for every fit
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    // cout << endl << endl;
    // cout << "=============================" << endl;
    if(eta_bin < 10) cout << "=== start with eta bin  " << eta_bin << " ===" << endl;
    else             cout << "=== start with eta bin " << eta_bin << " ===" << endl;
    // cout << "=============================" << endl << endl;
    for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
      g_factor[pt_bin] = factor[pt_bin][eta_bin];
      g_factor_err[pt_bin] = factor_err[pt_bin][eta_bin];
      g_ptrec[pt_bin] = ptrec_mean[pt_bin][eta_bin];
      g_ptrec_err[pt_bin] = ptrec_mean_err[pt_bin][eta_bin];
    }
    factor_pt[eta_bin] = new TGraphErrors(no_ptbins,g_ptrec,g_factor, g_ptrec_err, g_factor_err );

    // set first fit function
    TF1 *f1 = new TF1("f1",fitformula,0,500);

    // constrain parameters to make every fit converge
    if(ConstrainParams){
      for(unsigned int i=0; i<Nparams; i++){
        f1->SetParameter(i, StartParameters[eta_bin][i]);
        if(StartParameters[eta_bin][i] != 0){
          double up = StartParameters[eta_bin][i] + StartParameters[eta_bin][i] * RangeFraction[i];
          double down = StartParameters[eta_bin][i] - StartParameters[eta_bin][i] * RangeFraction[i];
          //f1->SetParLimits(i, up, down);
        }
      }
    }

    // perform the fit
    //cout << "=====================================================" << endl;
    OldFitResult[eta_bin] = factor_pt[eta_bin]->Fit(f1,"SR0Q");
    TString oldformula = fitformula;

    // now transform parameters to get uncorrelated ones
    TransformParameters Trafo(OldFitResult[eta_bin], oldformula, Nparams);
    TempFunction[eta_bin] = Trafo.getNewFunction();
    TempFit[eta_bin] = new TF1("fit", TempFunction[eta_bin],0,500);
    // Use Parameters from previous transformation as start value
    for(unsigned int i=0; i<Nparams; i++){
      TempFit[eta_bin]->SetParameter(i, OldFitResult[eta_bin]->Parameter(i));
    }
    TempFitResult[eta_bin] = factor_pt[eta_bin]->Fit(TempFit[eta_bin],"SR0Q");
    // cout << TempFunction[eta_bin] << endl;
    // TempFitResult[eta_bin]->Print("V");


    // Do Transformation Again
    // this has to be done because the numerical method does not deliver a
    // diagonal covariance matrix after one try
    TransformParameters Trafo2(TempFitResult[eta_bin], TempFunction[eta_bin], Nparams);
    NewFunction[eta_bin] = Trafo2.getNewFunction();
    NewFit[eta_bin] = new TF1("fit", NewFunction[eta_bin],0,500);
    // Use Parameters from previous transformation as start value
    for(unsigned int i=0; i<Nparams; i++){
      NewFit[eta_bin]->SetParameter(i, TempFitResult[eta_bin]->Parameter(i));
    }
    NewFitResult[eta_bin] = factor_pt[eta_bin]->Fit(NewFit[eta_bin],"SR0Q");
    chi2.push_back(NewFitResult[eta_bin]->Chi2());
    ndf.push_back(NewFitResult[eta_bin]->Ndf());
    // cout << NewFunction[eta_bin] << endl;
    // NewFitResult[eta_bin]->Print("V");

    // write new parameters and its errors
    vector<double> dummy;
    params.push_back(dummy);
    errors.push_back(dummy);
    errors_up.push_back(dummy);
    errors_down.push_back(dummy);
    for(unsigned int i=0; i<Nparams; i++){
      params[eta_bin].push_back(NewFitResult[eta_bin]->Parameter(i));
      errors[eta_bin].push_back(NewFitResult[eta_bin]->ParError(i));
    }
    //errors_up[eta_bin]   = CalculateParameterUncertainty(factor_pt[eta_bin], NewFunction[eta_bin], NewFitResult[eta_bin], params[eta_bin], "up");
    //errors_down[eta_bin] = CalculateParameterUncertainty(factor_pt[eta_bin], NewFunction[eta_bin], NewFitResult[eta_bin], params[eta_bin], "down");

    // perform fit with other polynomials and treat difference to original fit as uncertainty
    HighPol3Fit[eta_bin] = new TF1("fit2", "pol3",30,450);
    for(unsigned int i=0; i<4; i++){
      HighPol3Fit[eta_bin]->SetParameter(i, StartParameters[eta_bin][i]);
    }
    HighPol3FitResult[eta_bin] = factor_pt[eta_bin]->Fit(HighPol3Fit[eta_bin],"SR0Q");

    HighPol4Fit[eta_bin] = new TF1("fit2", "pol4",30,450);
    for(unsigned int i=0; i<6; i++){
      HighPol4Fit[eta_bin]->SetParameter(i, StartParameters[eta_bin][i]);
    }
    HighPol4FitResult[eta_bin] = factor_pt[eta_bin]->Fit(HighPol4Fit[eta_bin],"SR0Q");
  }

  // this is check comparing the method with increasing chi2 by to find uncertainty on fit parameters
  // to just getting the uncertainty from the root fit
  // for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
  //   for(unsigned int i=0; i<Nparams; i++){
  //     double percent_off1 = 100 * (abs(errors[eta_bin][i]) - abs(errors_up[eta_bin][i]))/abs(errors_up[eta_bin][i]);
  //     double percent_off2 = 100 * (abs(errors[eta_bin][i]) - abs(errors_down[eta_bin][i]))/abs(errors_up[eta_bin][i]);
  //     cout << "Error Diff in percent: " << abs(percent_off1) << " or " << abs(percent_off1) << endl;
  //   }
  // }


  // now compute variations with old parameters
  vector< vector<TF1*> > VariedFits;
  vector<TGraph*> UpVariation;
  vector<TGraph*> DownVariation;
  vector<TGraph*> Area;
  vector<TGraph*> AreaFromPol;

  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    // Variation Variations(params[eta_bin], errors_up[eta_bin], errors_down[eta_bin], NewFunction[eta_bin]);
    // UpVariation.push_back(Variations.CalculateEnvelope(NewFit[eta_bin], "up"));
    // DownVariation.push_back(Variations.CalculateEnvelope(NewFit[eta_bin], "down"));
    // Area.push_back(Variations.CalculateEnvelope(NewFit[eta_bin], "area"));
    AreaFromPol.push_back(UncertaintyFromDiff(NewFit[eta_bin], HighPol3Fit[eta_bin], HighPol4Fit[eta_bin], "area"));
    UpVariation.push_back(UncertaintyFromDiff(NewFit[eta_bin], HighPol3Fit[eta_bin], HighPol4Fit[eta_bin], "up"));
    DownVariation.push_back(UncertaintyFromDiff(NewFit[eta_bin], HighPol3Fit[eta_bin], HighPol4Fit[eta_bin], "down"));
  }

  // write up, down, central in root file
  TString rootname = "../CorrectionFile/Correction_allHad_"+year+".root";
  TFile * Correction_file = new TFile(rootname,"RECREATE");;
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    stringstream ss;
    ss << eta_bin;
    TString eta = ss.str();
    Correction_file->mkdir(eta);
    Correction_file->cd(eta);
    UpVariation[eta_bin]->Write("Up");
    DownVariation[eta_bin]->Write("Down");
    NewFit[eta_bin]->Write("Central");
  }
  Correction_file->Close();

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  TCanvas *A = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  h_ratio_mean->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ratio_mean->GetYaxis()->SetTitle("#eta^{rec}");
  h_ratio_mean->GetZaxis()->SetTitle("MEAN(p_{T}^{rec} / p_{T}^{gen})");
  h_ratio_mean->GetXaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetYaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetZaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetXaxis()->SetTitleOffset(0.9);
  h_ratio_mean->GetYaxis()->SetTitleOffset(0.8);
  h_ratio_mean->GetZaxis()->SetTitleOffset(0.9);
  h_ratio_mean->Draw("COLZ");
  h_ratio_mean->Draw("text:same");
  A->SaveAs(save_path+"/Correction_allHad/"+year+"/Ratio_Mean_"+year+".pdf");

  TCanvas *B = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  h_ptrec_mean->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ptrec_mean->GetYaxis()->SetTitle("#eta^{rec}");
  h_ptrec_mean->GetZaxis()->SetTitle("MEAN(p_{T}^{rec})");
  h_ptrec_mean->GetXaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetYaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetZaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetXaxis()->SetTitleOffset(0.9);
  h_ptrec_mean->GetYaxis()->SetTitleOffset(0.8);
  h_ptrec_mean->GetZaxis()->SetTitleOffset(0.9);
  h_ptrec_mean->Draw("COLZ");
  h_ptrec_mean->Draw("text:same");
  B->SaveAs(save_path+"/Correction_allHad/"+year+"/Ptrec_Mean_"+year+".pdf");

  TCanvas *C = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  h_ratio_mean_err->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ratio_mean_err->GetYaxis()->SetTitle("#eta^{rec}");
  h_ratio_mean_err->GetZaxis()->SetTitle("MEAN ERR(p_{T}^{rec} / p_{T}^{gen})");
  h_ratio_mean_err->GetXaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetYaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetZaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetXaxis()->SetTitleOffset(0.9);
  h_ratio_mean_err->GetYaxis()->SetTitleOffset(0.8);
  h_ratio_mean_err->GetZaxis()->SetTitleOffset(0.9);
  h_ratio_mean_err->Draw("COLZ");
  h_ratio_mean_err->Draw("text:same");
  C->SaveAs(save_path+"/Correction_allHad/"+year+"/Ratio_MEAN_ERR_"+year+".pdf");

  TCanvas *D = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  h_ptrec_mean_err->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ptrec_mean_err->GetYaxis()->SetTitle("#eta^{rec}");
  h_ptrec_mean_err->GetZaxis()->SetTitle("MEAN ERR(p_{T}^{rec})");
  h_ptrec_mean_err->GetXaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetYaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetZaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetXaxis()->SetTitleOffset(0.9);
  h_ptrec_mean_err->GetYaxis()->SetTitleOffset(0.8);
  h_ptrec_mean_err->GetZaxis()->SetTitleOffset(0.9);
  h_ptrec_mean_err->Draw("COLZ");
  h_ptrec_mean_err->Draw("text:same");
  D->SaveAs(save_path+"/Correction_allHad/"+year+"/Ptrec_MEAN_ERR_"+year+".pdf");


  TCanvas *E = new TCanvas();
  E->Divide(3,4);
  E->SetCanvasSize(800, 1200);
  E->SetWindowSize(800, 1200);
  TLegend *leg[no_etabins];
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    E->cd(eta_bin+1);
    gPad->SetLeftMargin(0.15);
    leg[eta_bin] = new TLegend(0.55,0.15,0.87,0.45);
    leg[eta_bin]->AddEntry(factor_pt[eta_bin],"correction factor","pl");
    if(show_fit){
      leg[eta_bin]->AddEntry(NewFit[eta_bin],"pol2 fit","l");
      leg[eta_bin]->AddEntry(HighPol3Fit[eta_bin],"pol3 fit","l");
      leg[eta_bin]->AddEntry(HighPol4Fit[eta_bin],"pol4 fit","l");
      leg[eta_bin]->AddEntry(AreaFromPol[eta_bin],"fit uncertainty","f");
    }
    std::stringstream title;
    title << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1];
    factor_pt[eta_bin]->GetYaxis()->SetTitle("correction factor");
    factor_pt[eta_bin]->GetXaxis()->SetRangeUser(30, 430);
    factor_pt[eta_bin]->GetYaxis()->SetRangeUser(0.8, 1.05);
    factor_pt[eta_bin]->SetTitle(title.str().c_str());
    factor_pt[eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}");
    factor_pt[eta_bin]->GetYaxis()->SetTitle("correction factor");
    factor_pt[eta_bin]->GetXaxis()->SetTitleSize(0.05);
    factor_pt[eta_bin]->GetXaxis()->SetTitleOffset(0.9);
    factor_pt[eta_bin]->GetXaxis()->SetNdivisions(505);
    factor_pt[eta_bin]->GetYaxis()->SetTitleSize(0.06);
    factor_pt[eta_bin]->GetYaxis()->SetTitleOffset(1.1);
    factor_pt[eta_bin]->GetYaxis()->SetNdivisions(505);
    factor_pt[eta_bin]->SetMarkerStyle(20);
    factor_pt[eta_bin]->SetMarkerSize(0.8);
    factor_pt[eta_bin]->SetLineColor(1);
    factor_pt[eta_bin]->Draw("AP");

    // Area[eta_bin]->SetFillColor(16);
    // Area[eta_bin]->Draw("f");

    AreaFromPol[eta_bin]->SetFillColor(16);
    if(show_fit) AreaFromPol[eta_bin]->Draw("f");

    factor_pt[eta_bin]->Draw("P");

    HighPol3Fit[eta_bin]->SetLineColor(kBlue);
    HighPol4Fit[eta_bin]->SetLineColor(kAzure+7);
    if(show_fit){
      HighPol3Fit[eta_bin]->Draw("SAME");
      HighPol4Fit[eta_bin]->Draw("SAME");
      NewFit[eta_bin]->Draw("SAME");
    }
    factor_pt[eta_bin]->Draw("P SAME");

    leg[eta_bin]->Draw();

    TLatex text;
    text.SetNDC(kTRUE);
    text.SetTextFont(43);
    text.SetTextSize(12);
    // format double precision for text box
    char chi2text[32];
    if(chi2[eta_bin] > 1) sprintf(chi2text, "%.4g", chi2[eta_bin]);
    else                  sprintf(chi2text, "%.3g", chi2[eta_bin]);
    TString chi2term = "#chi^{2}/ndf = ";
    chi2term += chi2text;
    chi2term += "/";
    chi2term += ndf[eta_bin];
    text.DrawLatex(.2,.2, chi2term);

    gPad->RedrawAxis();
  }
  E->SaveAs(save_path+"/Correction_allHad/"+year+"/Fits_"+year+".pdf");


  TCanvas *F = new TCanvas();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    F->Clear();
    F->Divide(3,4);
    F->SetCanvasSize(800, 1200);
    F->SetWindowSize(800, 1200);
    TLegend *leg_ratio[no_etabins];
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      F->cd(eta_bin+1);
      gPad->SetLeftMargin(0.15);
      TLine *ratio_line = new TLine(ratio_mean[pt_bin][eta_bin], 0, ratio_mean[pt_bin][eta_bin], ratio[pt_bin][eta_bin]->GetMaximum());
      ratio_line->SetLineColor(kRed);
      ratio_line->SetLineWidth(3);
      ratio_line->SetLineStyle(7);
      leg_ratio[eta_bin] = new TLegend(0.2,0.6,0.45,0.75);
      leg_ratio[eta_bin]->AddEntry(ratio[pt_bin][eta_bin],"p_{T}^{rec}/p_{T}^{gen}","f");
      leg_ratio[eta_bin]->AddEntry(ratio_fit[pt_bin][eta_bin],"fit","l");
      if(use_median)leg_ratio[eta_bin]->AddEntry(ratio_line,"median","l");
      else          leg_ratio[eta_bin]->AddEntry(ratio_line,"gaussian mean","l");
      std::stringstream title2;
      title2 << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1] << " | " << xbins[pt_bin] << " < p_{T}^{gen} < " << xbins[pt_bin+1];
      ratio[pt_bin][eta_bin]->SetTitle(title2.str().c_str());
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}/p_{T}^{gen}");
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitle("events");
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitleSize(0.05);
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitleOffset(0.9);
      ratio[pt_bin][eta_bin]->GetXaxis()->SetNdivisions(505);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitleSize(0.06);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitleOffset(1.1);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetNdivisions(505);
      ratio[pt_bin][eta_bin]->SetFillColor(kGray);
      ratio[pt_bin][eta_bin]->SetLineColor(1);
      ratio[pt_bin][eta_bin]->Draw("HIST");
      ratio_fit[pt_bin][eta_bin]->Draw("SAME");
      leg_ratio[eta_bin]->Draw();
      ratio_line->Draw("SAME");
      gPad->RedrawAxis();
    }
    TString name = save_path+"/Correction_allHad/"+year+"/Fits_ratio_ptbin";
    TString ending = std::to_string(pt_bin) + "_"+year+".pdf";
    F->SaveAs(name+ending);
  }


  TCanvas *F2 = new TCanvas();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    F2->Clear();
    F2->Divide(3,4);
    F2->SetCanvasSize(800, 1200);
    F2->SetWindowSize(800, 1200);
    TLegend *leg_ptrec[no_etabins];
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      F2->cd(eta_bin+1);
      gPad->SetLeftMargin(0.15);
      TLine *ptrec_line = new TLine(ptrec_mean[pt_bin][eta_bin], 0, ptrec_mean[pt_bin][eta_bin], ptrec[pt_bin][eta_bin]->GetMaximum());
      ptrec_line->SetLineColor(kRed);
      ptrec_line->SetLineWidth(3);
      ptrec_line->SetLineStyle(7);
      std::stringstream title3;
      title3 << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1] << " | " << xbins[pt_bin] << " < p_{T}^{gen} < " << xbins[pt_bin+1];
      ptrec[pt_bin][eta_bin]->SetTitle(title3.str().c_str());
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}");
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitle("events");
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetRangeUser(0, 1.1 * ptrec[pt_bin][eta_bin]->GetMaximum());
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitleSize(0.05);
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitleOffset(0.9);
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetNdivisions(505);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitleSize(0.06);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitleOffset(1.1);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetNdivisions(505);
      ptrec[pt_bin][eta_bin]->SetFillColor(kGray);
      ptrec[pt_bin][eta_bin]->SetLineColor(1);
      ptrec[pt_bin][eta_bin]->Draw("HIST");
      ptrec_line->Draw("SAME");
      gPad->RedrawAxis();
    }
    TString name = save_path+"/Correction_allHad/"+year+"/ptrec_mean_ptbin";
    TString ending = std::to_string(pt_bin) + "_"+year+".pdf";
    F2->SaveAs(name+ending);
  }

  ////----
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    PlotRatio(factor_pt[eta_bin], NewFit[eta_bin], AreaFromPol[eta_bin], eta_bin, year);
    TCanvas *G = new TCanvas("G", "G", 600, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.105);
    factor_pt[eta_bin]->Draw("AP");
    if(show_fit) AreaFromPol[eta_bin]->Draw("f");
    factor_pt[eta_bin]->Draw("P");
    if(show_fit){
      HighPol3Fit[eta_bin]->Draw("SAME");
      HighPol4Fit[eta_bin]->Draw("SAME");
      NewFit[eta_bin]->Draw("SAME");
    }
    leg[eta_bin]->Draw("SAME");
    factor_pt[eta_bin]->Draw("P SAME");

    TLatex text;
    text.SetNDC(kTRUE);
    text.SetTextFont(43);
    text.SetTextSize(18);
    // format double precision for text box
    char chi2text[32];
    if(chi2[eta_bin] > 1) sprintf(chi2text, "%.4g", chi2[eta_bin]);
    else                  sprintf(chi2text, "%.3g", chi2[eta_bin]);
    TString chi2term = "#chi^{2}/ndf = ";
    chi2term += chi2text;
    chi2term += "/";
    chi2term += ndf[eta_bin];
    text.DrawLatex(.2,.2, chi2term);

    gPad->RedrawAxis();

    TString filename = save_path+"/Correction_allHad/"+year+"/Fits_";
    filename += eta_bin;
    filename += "_"+year_s+".pdf";
    G->SaveAs(filename);
    delete G;
  }

  ////----

  TCanvas *G2 = new TCanvas("G2", "G2", 600, 600);
  gPad->SetLeftMargin(0.15);
  TLine *example_line = new TLine(ptrec_mean[2][3], 0, ptrec_mean[2][3], ptrec[2][3]->GetMaximum());
  example_line->SetLineColor(kRed);
  example_line->SetLineWidth(3);
  example_line->SetLineStyle(7);
  ptrec[2][3]->Draw("HIST");
  example_line->Draw("SAME");
  gPad->RedrawAxis();
  G2->SaveAs(save_path+"/Correction_allHad/"+year+"/ptrec_mean_example_"+year+".pdf");

  TCanvas *G3 = new TCanvas("G3", "G3", 600, 600);
  gPad->SetLeftMargin(0.15);
  ratio[2][3]->Draw("HIST");
  ratio_fit[2][3]->Draw("SAME");
  gPad->RedrawAxis();
  G3->SaveAs(save_path+"/Correction_allHad/"+year+"/ratio_mean_example_"+year+".pdf");

  TCanvas *H = new TCanvas("H", "H", 1800, 600);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  h_count->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_count->GetYaxis()->SetTitle("#eta^{rec}");
  h_count->GetZaxis()->SetTitle("Events");
  h_count->GetXaxis()->SetTitleSize(0.05);
  h_count->GetYaxis()->SetTitleSize(0.05);
  h_count->GetZaxis()->SetTitleSize(0.05);
  h_count->GetXaxis()->SetTitleOffset(0.9);
  h_count->GetYaxis()->SetTitleOffset(0.8);
  h_count->GetZaxis()->SetTitleOffset(0.9);
  h_count->Draw("COLZ");
  h_count->Draw("text:same");
  H->SaveAs(save_path+"/Correction_allHad/"+year+"/Event_Count_"+year+".pdf");

  return 0;
}

/*
████████     ██████   █████  ██████   █████  ███    ███
   ██        ██   ██ ██   ██ ██   ██ ██   ██ ████  ████
   ██        ██████  ███████ ██████  ███████ ██ ████ ██
   ██        ██      ██   ██ ██   ██ ██   ██ ██  ██  ██
   ██ ██     ██      ██   ██ ██   ██ ██   ██ ██      ██
*/


//
// Class to transform fit parameters into uncorrelated parameters
//

TransformParameters::TransformParameters(TFitResultPtr fitresult, TString oldformula, int n_param){
  //
  // To get uncorrelated Parameters from a fit, a matrix equation has to be solved.
  // From the regular fit, one gets a set of Parameters P = (a,b,c,...).
  // Now you have to find a Transformation O^T Cov(P) O = D, where D is a Diagonal Matrix with Eigenvalues
  // of Cov(P) on its diagonal.
  // The new (uncorrelated) Parameters are then P' = O^T P
  //

  // get Covariance Matrix
  TMatrixDSym oldCov = fitresult->GetCovarianceMatrix();
  // extract Eigenvectors
  const TMatrixDSymEigen oldEigen(oldCov);
  // get Transformation Matrix O and O^T
  const TMatrixD TrafoMatrix = oldEigen.GetEigenVectors();
  TMatrixD TranspMatrix(TrafoMatrix);
  TranspMatrix.T();

  // get old parameters
  vector<double> oldParams;
  for(int i=0; i<n_param; i++){
    oldParams.emplace_back(fitresult->Parameter(i));
  }
  const TVectorD V_oldParams(n_param,&oldParams[0]);

  // calculate new parameters in new basis
  TDecompSVD ParameterSVD(TrafoMatrix);
  Bool_t ok;
  const TVectorD V_newParams= ParameterSVD.Solve(V_oldParams, ok);

  // transform back new parameters into old ones
  // a = a (a',b',c',...)
  TDecompSVD ParameterSVD2(TranspMatrix);
  Bool_t ok2;
  const TVectorD V_old_from_newParams= ParameterSVD2.Solve(V_newParams, ok2);

  //parametrize old parameters by using new ones with Trafo Matrix
  TMatrixD M_Dummy(TrafoMatrix);

  // Fill vector of vectors: 1 vector for each old parameter
  // containing the values of the new parameters expressing the old one.
  vector<vector<double>> TransformationValues;
  for(int i=0; i<n_param; i++){
    vector<double> singleParam;
    for(int j=0; j<n_param; j++){
      singleParam.emplace_back(TMatrixDRow(M_Dummy,i)[j]);
    }
    TransformationValues.emplace_back(singleParam);
    for(unsigned int j=0; j<singleParam.size(); j++){
      singleParam.pop_back();
    }
  }

  // now write new fit function with new parameters
  oldformula;
  for(int i=0; i<n_param; i++){
    TString idx = "[";
    idx += i;
    idx += "]";
    TString new_expr = "( [0] *(";
    new_expr += TransformationValues[i][0];
    new_expr += ")";
    for(int j=1; j<n_param; j++){
      TString newidx = "";
      newidx += j;
      new_expr += " + ["+newidx+"] *(";
      new_expr += TransformationValues[i][j];
      new_expr += ")";
    }
    new_expr += ")";
    oldformula.ReplaceAll(idx,new_expr);
  }

  newFormula = oldformula;
}

/*
 ██████  ███████ ████████     ███████ ██    ██ ███    ██  ██████
██       ██         ██        ██      ██    ██ ████   ██ ██
██   ███ █████      ██        █████   ██    ██ ██ ██  ██ ██
██    ██ ██         ██        ██      ██    ██ ██  ██ ██ ██
 ██████  ███████    ██        ██       ██████  ██   ████  ██████
*/


// Some Get Functions
TString TransformParameters::getNewFunction(){
  return newFormula;
}

/*
██    ██  █████  ██████  ██  █████  ████████ ██  ██████  ███    ██
██    ██ ██   ██ ██   ██ ██ ██   ██    ██    ██ ██    ██ ████   ██
██    ██ ███████ ██████  ██ ███████    ██    ██ ██    ██ ██ ██  ██
 ██  ██  ██   ██ ██   ██ ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
  ████   ██   ██ ██   ██ ██ ██   ██    ██    ██  ██████  ██   ████
*/


Variation::Variation(vector<double> params, vector<double> errors_up, vector<double> errors_down, TString function){
  vector<double> pup;
  vector<double> pdown;

  int Nparams = params.size();

  for(unsigned int i=0; i<Nparams; i++){
    pup.push_back(params[i] + sqrt(errors_up[i]*errors_up[i]));
    pdown.push_back(params[i] - sqrt(errors_down[i]*errors_down[i]));
  }

  // Only one parameter of the function is either set to 'up' or 'down'.
  // Therefore, for 4 parameters there are 4 functions for 'up' and 4 functions for 'down' variations.
  int Nvariations = Nparams*2;

  // now get a function for every possible variation
  // first, iterate through 'up' variations (first half of variations)
  // then, iterate through 'down' variations (second half)
  // the i==j statement makes sure to iterate the variations through every parameter
  for(unsigned int i=0; i<Nvariations; i++){
    TF1* f = new TF1("f",function,0,500);
    for(unsigned int j=0; j<Nparams; j++){
      if(i<Nparams){
        if(i==j) f->SetParameter(j, pup[j]);
        else     f->SetParameter(j, params[j]);
      }
      else{
        if(i==j) f->SetParameter(j, pdown[j]);
        else     f->SetParameter(j, params[j]);
      }
    }
    VariedFits.push_back(f);
  }
}

/*
███████ ███    ██ ██    ██ ███████ ██       ██████  ██████  ███████
██      ████   ██ ██    ██ ██      ██      ██    ██ ██   ██ ██
█████   ██ ██  ██ ██    ██ █████   ██      ██    ██ ██████  █████
██      ██  ██ ██  ██  ██  ██      ██      ██    ██ ██      ██
███████ ██   ████   ████   ███████ ███████  ██████  ██      ███████
*/


TGraph *Variation::CalculateEnvelope(TF1* central, string updown){

  double c;
  double v;
  double delta;
  const int steps = 100;
  double stepsize = 500/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  for(int i=0; i<steps; i++){
    x = stepsize * i;
    c = central->Eval(x);
    double down2 = 0;
    double up2 = 0;
    for(unsigned int j=0; j<VariedFits.size(); j++){
      v = VariedFits[j]->Eval(x);
      delta = c-v;
      if(delta < 0) delta *= -1;
      if(v < c) down2 += pow(delta,2);
      if(v >= c) up2 +=  pow(delta,2);
    }
    up[i] = c + sqrt(up2);
    down[i] = c - sqrt(down2);
    xvalues[i] = x;
  }

  TGraph *upgraph = new TGraph(xvalues, up);
  TGraph *downgraph = new TGraph(xvalues, down);

  TGraph *area = new TGraph(2*steps);
  for (int i=0;i<steps;i++) {
    area->SetPoint(i,xvalues[i],up[i]);
    area->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
  }

  if(updown == "up") return upgraph;
  if(updown == "down") return downgraph;
  else return area;

}

/*
██    ██ ███    ██  ██████        ███████ ██████   ██████  ███    ███     ██████  ██ ███████ ███████
██    ██ ████   ██ ██             ██      ██   ██ ██    ██ ████  ████     ██   ██ ██ ██      ██
██    ██ ██ ██  ██ ██             █████   ██████  ██    ██ ██ ████ ██     ██   ██ ██ █████   █████
██    ██ ██  ██ ██ ██             ██      ██   ██ ██    ██ ██  ██  ██     ██   ██ ██ ██      ██
 ██████  ██   ████  ██████ ██     ██      ██   ██  ██████  ██      ██     ██████  ██ ██      ██
*/


TGraph *UncertaintyFromDiff(TF1* central, TF1* pol1, TF1* pol2, TString mode){
  int steps = 100;
  double stepsize = 500/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  for(int i=0; i<steps; i++){
    x = stepsize * i;
    double c = central->Eval(x);
    double err1 = sqrt( (c - pol1->Eval(x)) * (c - pol1->Eval(x)) );
    double err2 = sqrt( (c - pol2->Eval(x)) * (c - pol2->Eval(x)) );
    double err;
    if(err1 > err2) err = err1;
    else            err = err2;
    up[i] = c + err;
    down[i] = c - err;
    xvalues[i] = x;
  }
  if(mode == "area"){
    TGraph *area = new TGraph(2*steps);
    for (int i=0;i<steps;i++) {
      area->SetPoint(i,xvalues[i],up[i]);
      area->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
    }
    return area;
  }
  else if(mode == "up"){
    TGraph *up_var = new TGraph(steps);
    for (int i=0;i<steps;i++) up_var->SetPoint(i,xvalues[i],up[i]);
    return up_var;
  }
  else if(mode == "down"){
    TGraph *down_var = new TGraph(steps);
    for (int i=0;i<steps;i++) down_var->SetPoint(i,xvalues[i],down[i]);
    return down_var;
  }
}

/*
██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████          ██    ██ ███    ██  ██████
██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██         ██    ██ ████   ██ ██
██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████          ██    ██ ██ ██  ██ ██
██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██         ██    ██ ██  ██ ██ ██
██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██          ██████  ██   ████  ██████ ██
*/


// vary one parameter until the original chi2 grows by 1
// the delta to the original parameter is considered the 1sigma uncertainty
vector<double> CalculateParameterUncertainty(TGraphErrors* data_, TString function, TFitResultPtr fitresult, vector<double> parameters, TString updown){
  TGraphErrors* data = (TGraphErrors*) data_->Clone();
  TF1 *fit = new TF1("fit", function, 0, 500);
  double originalChi2 = fitresult->Chi2();
  vector<double> uncertainties;
  for(unsigned int i=0; i<parameters.size(); i++){
    int i_step = 0;
    double stepsize = 0.000003;
    double delta = 0;
    double diff = 0;
    bool foundUncert = false;
    double newChi2 = originalChi2;
    while(!foundUncert){
      i_step++;
      if(updown == "up")        delta =  i_step*sqrt(parameters[i]*parameters[i])*stepsize;
      else if(updown == "down") delta = -i_step*sqrt(parameters[i]*parameters[i])*stepsize;
      else throw std::invalid_argument("Select up or down in CalculateParameterUncertainty");
      for(unsigned int j=0; j<parameters.size(); j++){
        if(i != j) fit->FixParameter(j, parameters[j]);
        else       fit->FixParameter(j, parameters[j]+delta);
      }
      TFitResultPtr newresult = data->Fit(fit,"QSR");
      newChi2 = newresult->Chi2();
      diff = newChi2 - originalChi2;
      if(diff >= 1) foundUncert = true;
      if(i_step > 50000){
        cout << "had to break parameter uncertainty calculation after 50000 steps!!!" << endl;
        foundUncert = true;
      }
    }
    cout << "It took " << i_step << " steps to find uncertainty of parameter " << i
    << " resulting in a Chi2 difference of " << diff << endl;
    // cout << "Chi2 before = " << originalChi2 << " Chi2 after = " << newChi2 << endl;
    // cout << "Parameter before = " << parameters[i] << " Parameter after = " << parameters[i]+delta << endl;
    uncertainties.push_back(delta);
  }
  return uncertainties;
}

/*
██████  ██       ██████  ████████     ██████   █████  ████████ ██  ██████
██   ██ ██      ██    ██    ██        ██   ██ ██   ██    ██    ██ ██    ██
██████  ██      ██    ██    ██        ██████  ███████    ██    ██ ██    ██
██      ██      ██    ██    ██        ██   ██ ██   ██    ██    ██ ██    ██
██      ███████  ██████     ██        ██   ██ ██   ██    ██    ██  ██████
*/


void PlotRatio(TGraphErrors* data_, TF1* fit_, TGraph* uncert_, int etabin, TString year_){
  TString save_path = get_save_path();

  TGraphErrors* data = (TGraphErrors*) data_->Clone();
  TF1* fit = (TF1*) fit_->Clone();
  TGraph* uncert = (TGraphErrors*) uncert_->Clone();
  TString year = year_;

  // now make data ratio to fit
  int Npoints = data->GetN();
  Double_t *x = data->GetX();
  Double_t *y = data->GetY();
  Double_t *ex = data->GetEX();
  Double_t *ey = data->GetEY();

  for(int i=0; i<Npoints; i++){
    double ynew = y[i] / fit->Eval(x[i]);
    double eynew = ey[i] / fit->Eval(x[i]);
    data->SetPoint(i, x[i], ynew);
    data->SetPointError(i, ex[i], eynew);
  }

  // now make uncert ratio to fit
  int Npoints_u = uncert->GetN();
  Double_t *x_u = uncert->GetX();
  Double_t *y_u = uncert->GetY();

  for(int i=0; i<Npoints_u; i++){
    double ynew = y_u[i] / fit->Eval(x_u[i]);
    uncert->SetPoint(i, x_u[i], ynew);
  }

  TLine *line = new TLine(30,1,430,1);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);

  TCanvas *C = new TCanvas("C", "C", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.105);

  uncert->SetFillColor(16);

  data->Draw("AP");
  uncert->Draw("f");
  line->Draw("SAME");
  data->Draw("P SAME");


  data->GetXaxis()->SetTitle("p_{T}^{rec}");
  data->GetYaxis()->SetTitle("#frac{correction factor}{fit}");
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetXaxis()->SetTitleOffset(0.9);
  data->GetXaxis()->SetNdivisions(505);
  data->GetYaxis()->SetTitleSize(0.06);
  data->GetYaxis()->SetTitleOffset(1.1);
  data->GetYaxis()->SetNdivisions(505);

  data->GetXaxis()->SetRangeUser(30, 430);
  data->GetYaxis()->SetRangeUser(0.95, 1.05);
  gPad->RedrawAxis();

  TLegend* leg = new TLegend(0.55,0.15,0.87,0.40);
  leg->AddEntry(line,"fit","l");
  leg->AddEntry(uncert,"fit uncertainty","f");
  leg->Draw();

  TString filename = save_path+"/Correction_allHad/"+year_+"/FitRatio_";
  filename += etabin;
  filename += "_"+year_+".pdf";
  C->SaveAs(filename);
  delete C;
  return;
}
