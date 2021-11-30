#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"
#include "../include/normalise.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

#include "TDecompLU.h"

using namespace std;

void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown);
void PlotFit(TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin);
void PlotChi2(TF1* chi2function, vector<double> FSRvalues);
void PlotError(TH1F* hist, TString sysname);
void PlotResults(vector<vector<double>> Values, vector<TString> masswindows);
void printCovMatrix(TH2* m, TString name, double factor=1, int prec=2, bool diag=false);
void printCovMatrix(TMatrixD m, TString name, double factor=1, int prec=2, bool diag=false);
void DrawCov(TH2* m, TString text);
TString ConstructChi2(TMatrixD m, vector<TString> vecS, VecD vecD, int nbins);
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname);
TH2D* GetCovMatrix(TH1F* hist, double width, TString text, vector<bool> skip);
TH2D* NormCovMatrix(TH1* hist_, TH2* old_cov, int lower_, int upper_, bool width_);
vector<double> ExtractFSRValues(TF1* chi2function);

TString year, channel;
TString save_path;

bool debug = false;
bool forTalk = true;

int precision = 5; // For double to string function - dtos()

int main(int argc, char* argv[]){

  if(argc != 3){
    cout << "Usage: ./FSRuncertainty <year> <channel>" << endl;
    cout << "  year = 2016/2017/2018/combine" << endl;
    cout << "  channel = muon/elec/combine" << endl;
    return 0;
  }
  year = argv[1];
  TString ichannel = argv[2];
  if(ichannel.EqualTo("muon")) channel = "_muon";
  if(ichannel.EqualTo("elec")) channel = "_elec";
  if(ichannel.EqualTo("combine")) channel = "";

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();

  cout << "\n-------------" << endl;
  cout << "year:    " << year << endl;
  cout << "channel: " << ichannel << endl;
  cout << "-------------" << endl;

  // vector<TString> masswindows = {"140toInf", "140to160", "160to170", "170to180", "180to190", "190to200", "200toInf"};
  vector<TString> masswindows = {"140toInf"};
  vector<TString> ptwindows = {"400toInf", "400to450","450to500", "500toInf"};
  vector<vector<double>> f_fsr;

  for(auto mass_string: masswindows){
    cout << "Start with mass " << mass_string << " ... " << endl;
    if(debug) cout << "\t ... Create save_path" << endl;
    // save_path = get_save_path()+"/FSRuncertainty/";
    save_path = get_save_path();
    save_path = creat_folder_and_path(save_path, "FSRuncertainty_TEST");
    save_path = creat_folder_and_path(save_path, year);
    save_path = creat_folder_and_path(save_path, ichannel);
    save_path = creat_folder_and_path(save_path, mass_string);

    TFile* infile = new TFile("FSR_hists_mjet"+mass_string+".root");

    // -------------------------------------------------------------------------------------------------
    // Get Hists

    if(debug) cout << "\t ... Get Hists" << endl;
    TH1F *data = (TH1F*) infile->Get("data"+channel+"_"+year);
    TH1F *ttbar = (TH1F*) infile->Get("nominal"+channel+"_"+year);
    vector<TH1F*> bkg;
    bkg.push_back((TH1F*) infile->Get("bgr_st"+channel+"_"+year));
    bkg.push_back((TH1F*) infile->Get("bgr_wj"+channel+"_"+year));
    bkg.push_back((TH1F*) infile->Get("bgr_ot"+channel+"_"+year));
    vector<double> bkgsys = {0.23, 0.19, 1.0};

    if(debug) cout << "\t ... Get Hists for Jet" << endl;
    // TH1F *jecup = (TH1F*) infile->Get("jecup"+channel+"_"+year);
    // TH1F *jecdown = (TH1F*) infile->Get("jecdown"+channel+"_"+year);
    // TH1F *corup = (TH1F*) infile->Get("corup"+channel+"_"+year);
    // TH1F *cordown = (TH1F*) infile->Get("cordown"+channel+"_"+year);
    // TH1F *jmsup = (TH1F*) infile->Get("jmsup"+channel+"_"+year);
    // TH1F *jmsdown = (TH1F*) infile->Get("jmsdown"+channel+"_"+year);

    if(debug) cout << "\t ... Get Hists for FSR" << endl;
    vector<TH1F*> FSRup, FSRdown;
    vector<double> FSRvalues;
    if(year == "2016"){
      // FSRvalues = {1./2, 1., 2.};
      FSRvalues = {0.25, 1., 4.}; // squared
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+year));
    }
    else{
      // FSRvalues = {1./4, 1./2, 1./sqrt(2), 1., sqrt(2), 2., 4.};
      FSRvalues = {0.0625, 0.25, 0.5, 1., 2., 4., 16.}; // squared
      FSRdown.push_back((TH1F*) infile->Get("fsrdown4"+channel+"_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdownsqrt2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrupsqrt2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup4"+channel+"_"+year));
    }

    // -------------------------------------------------------------------------------------------------
    // Subtract backgournd

    if(debug) cout << "\t ... Subtract Bkg" << endl;
    TH1F* data_sub = SubtractBackgrounds(data, bkg, bkgsys);

    // -------------------------------------------------------------------------------------------------
    // Rebin

    if(debug) cout << "\t ... Rebin" << endl;
    int rebin = 5;
    double width = 0.01*rebin;
    if(rebin > 1){
      data_sub->Rebin(rebin);
      ttbar->Rebin(rebin);
      // jecup->Rebin(rebin);
      // jecdown->Rebin(rebin);
      // corup->Rebin(rebin);
      // cordown->Rebin(rebin);
      // jmsup->Rebin(rebin);
      // jmsdown->Rebin(rebin);
      for(auto fsr: FSRup) fsr->Rebin(rebin);
      for(auto fsr: FSRdown) fsr->Rebin(rebin);
      cout << "Rebin with a factor of " << rebin << " - resulting in a bin width of " <<  data_sub->GetBinWidth(1) << endl;
    }

    // -------------------------------------------------------------------------------------------------
    // Skip Bins

    // Assumption: 0 entries (on diag) in CovMatrix lead to problems Inversion
    // Remove all bins where we cannot perform a fit -> data=0 or not enough MC points

    vector<bool> skipBin;
    int freePara = (year=="2016")?2:3; // Hard-coded
    int numberBins = data_sub->GetNbinsX(); // Equal for all
    for(unsigned int i=1; i<=numberBins; i++){

      // negative bin content possible since data-bkg(MC)
      bool skipData = data_sub->GetBinContent(i)<=0;

      int MCs = 0;
      if(ttbar->GetBinContent(i)>0) MCs++;
      for(auto fsr: FSRup){if(fsr->GetBinContent(i)>0) MCs++;}
      for(auto fsr: FSRdown){if(fsr->GetBinContent(i)>0) MCs++;}
      bool skipMC = MCs<=freePara;

      bool skip = (skipData||skipMC);
      // cout << i << "\t Data " << skipData << "\t MC " << skipMC << "\t SkipBin: " << skip << "\tContent " << data_sub->GetBinContent(i) << endl;
      skipBin.push_back(skip);
    }

    // -------------------------------------------------------------------------------------------------
    // Create covariance matrix

    int shift = 0; // Skip first bins which are zero
    for(auto s: skipBin){if(s) shift++;}
    int upper_ = data_sub->GetNbinsX(); // full dimension
    if(debug) cout << "Lower " << shift << "\tUpper " << upper_ << endl;

    TH2D* covData = GetCovMatrix(data_sub, width, "data_sub", skipBin);
    TH2D* covttbar = GetCovMatrix(ttbar, width, "ttbar", skipBin);

    if(debug) cout << "Dimension Data " << covData->GetNbinsX() << endl;
    if(debug) cout << "Dimension ttbar " << covttbar->GetNbinsX() << endl;

    TH2D* covData_norm = NormCovMatrix(data_sub, covData, shift, upper_, false);
    TH2D* covttbar_norm = NormCovMatrix(ttbar, covttbar, shift, upper_, false);

    DrawCov(covData, "covData_nominal");   DrawCov(covData_norm, "covData_norm");
    DrawCov(covttbar, "covttbar_nominal"); DrawCov(covttbar_norm, "covttbar_norm");

    vector<TH2D*> covFSRup, covFSRup_norm;
    vector<TH2D*> covFSRdown, covFSRdown_norm;

    for(unsigned int i=0; i<FSRup.size(); i++){
      TString bin = to_string(i);
      covFSRup.push_back(GetCovMatrix(FSRup[i], width, "FSRup"+bin, skipBin));       covFSRup_norm.push_back(NormCovMatrix(FSRup[i], covFSRup[i], shift, upper_, false));
      covFSRdown.push_back(GetCovMatrix(FSRdown[i], width, "FSRdown"+bin, skipBin)); covFSRdown_norm.push_back(NormCovMatrix(FSRdown[i], covFSRdown[i], shift, upper_, false));

      DrawCov(covFSRup[i], "covFSRup"+bin+"_nominal");     DrawCov(covFSRup_norm[i], "covFSRup"+bin+"_norm");
      DrawCov(covFSRdown[i], "covFSRdown"+bin+"_nominal"); DrawCov(covFSRdown_norm[i], "covFSRdown"+bin+"_norm");
    }

    if(debug) cout << "Dimension FSRup " << covFSRup[0]->GetNbinsX() << endl;
    if(debug) cout << "Dimension FSRdown " << covFSRdown[0]->GetNbinsX() << endl;

    // printCovMatrix(covData, "covData");
    // printCovMatrix(covttbar, "covttbar");
    // printCovMatrix(covFSRup[0], "covFSRup");
    // printCovMatrix(covFSRdown[0], "covFSRdown");

    // ------------ Create covariance matrix for fit

    // Remove Zero-Bins from dimension
    double mDim = covData->GetNbinsX();
    if(mDim!=covData->GetNbinsY()) cerr << __LINE__ << ": Matrix is NOT symmetric" << endl;

    vector<double> diagValues;
    TH2D* covLargestUnc = new TH2D("covLargestUnc", "", mDim, shift*width, 1, covData->GetNbinsX(), shift*width, 1);
    for(int x=1; x<=mDim; x++){
      diagValues.push_back(covttbar_norm->GetBinContent(x,x));
      for(auto cov: covFSRup_norm) diagValues.push_back(cov->GetBinContent(x,x));
      for(auto cov: covFSRdown_norm) diagValues.push_back(cov->GetBinContent(x,x));
      double max = *max_element(diagValues.begin(), diagValues.end());

      if(debug){ for(auto v: diagValues) cout << v << "\t"; }
      if(debug){ cout<<"Max value: "<<max<<endl; }

      covLargestUnc->SetBinContent(x,x,max);
      diagValues = {};
    }

    // printCovMatrix(covttbar_norm, 100000);
    // printCovMatrix(covFSRup_norm[0], 100000);
    // printCovMatrix(covFSRdown_norm[0], 100000);
    // printCovMatrix(covLargestUnc, 100000);

    // ------------ Final covariance matrix

    TMatrixD covDataStat = TMatrix(mDim, mDim);
    TMatrixD covFit = TMatrix(mDim, mDim);
    for(int x = 0; x<mDim; x++){
      covFit[x][x] = covLargestUnc->GetBinContent(x+1,x+1);
      for(int y = 0; y<mDim; y++){
        covDataStat[x][y] = covData_norm->GetBinContent(x+1,y+1);

      }
    }
    TMatrixD covTotalError = TMatrix(covDataStat);
    covTotalError += covFit;

    // printCovMatrix(covTotalError, "covTotalError", 100000); // covTotalError.Print();
    // printCovMatrix(covDataStat, "covDataStat", 100000);  // covDataStat.Print();
    // printCovMatrix(covFit, "covFit", 100000);       // covFit.Print();

    // -------------------------------------------------------------------------------------------------
    // Normalize

    if(debug) cout << "\t ... Normalize" << endl;
    TH1F *data_norm = Normalize(data_sub);
    TH1F *ttbar_norm = Normalize(ttbar);
    // TH1F *jecup_norm = Normalize(jecup);
    // TH1F *jecdown_norm = Normalize(jecdown);
    // TH1F *corup_norm = Normalize(corup);
    // TH1F *cordown_norm = Normalize(cordown);
    // TH1F *jmsup_norm = Normalize(jmsup);
    // TH1F *jmsdown_norm = Normalize(jmsdown);
    vector<TH1F*> FSRup_norm, FSRdown_norm;
    for(auto fsr: FSRup) FSRup_norm.push_back(Normalize(fsr));
    for(auto fsr: FSRdown) FSRdown_norm.push_back(Normalize(fsr));
    PlotTau32(data_norm, ttbar_norm, FSRup_norm, FSRdown_norm);
    PlotError(data_norm, "DATA");
    PlotError(ttbar_norm, "TTbar");

    // -------------------------------------------------------------------------------------------------
    // Get additional uncertainties
    if(debug) cout << "\t ... Plot systematics (jec,cor,jms)" << endl;
    vector<TH1F*> sys;
    // sys.push_back(GetSYS(ttbar_norm, jecup_norm, jecdown_norm, "jec"));
    // sys.push_back(GetSYS(ttbar_norm, corup_norm, cordown_norm, "cor"));
    // sys.push_back(GetSYS(ttbar_norm, jmsup_norm, jmsdown_norm, "jms"));



    // -------------------------------------------------------------------------------------------------
    // Fit per bin
    if(debug) cout << "\t ... Fit" << endl;
    vector<vector<double>> fitparameters;
    vector<int> validbins;
    vector<TString> BinsFit;
    vector<double> BinsData;
    int nbins = ttbar_norm->GetXaxis()->GetNbins();
    TString fitformula = "[0] + [1]*log(x) + [2]*x";
    if(year == "2016") fitformula = "[0] + [1]*log(x)";
    int Npar = 1;
    if(fitformula.Contains("[1]")) Npar = 2;
    if(fitformula.Contains("[2]")) Npar = 3;
    if(fitformula.Contains("[3]")) Npar = 4;
    if(fitformula.Contains("[4]")) Npar = 5;
    if(fitformula.Contains("[5]")) Npar = 6;

    vector<TGraphErrors*> graphs, bands1, bands2;
    vector<TF1*> fits;
    for(int bin=1; bin<=nbins; bin++){
      // Ignore empty bins
      double mincontent = 0.0001;
      bool ignorebin = false;
      if(ttbar_norm->GetBinContent(bin) < mincontent) ignorebin = true;
      if(data_norm->GetBinContent(bin) < mincontent) ignorebin = true;
      for(auto fsr: FSRup_norm){
        if(fsr->GetBinContent(bin) < mincontent) ignorebin = true;
      }
      for(auto fsr: FSRdown_norm){
        if(fsr->GetBinContent(bin) < mincontent) ignorebin = true;
      }
      if(ignorebin){
        cout << "Too few entries, ignore bin " << bin << "!" << endl;
        continue;
      }

      validbins.push_back(bin);
      // Write bin contents and error in vectors in correct order
      vector<double> bincontents;
      vector<double> binerrors;
      vector<double> binerrors_large;
      vector<double> xerrors;

      for(auto fsr: FSRdown_norm){
        bincontents.push_back(fsr->GetBinContent(bin));
        binerrors.push_back(fsr->GetBinError(bin));
        binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
        xerrors.push_back(0.0);
      }
      bincontents.push_back(ttbar_norm->GetBinContent(bin));
      binerrors.push_back(ttbar_norm->GetBinError(bin));
      binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
      xerrors.push_back(0.0);
      for(auto fsr: FSRup_norm){
        bincontents.push_back(fsr->GetBinContent(bin));
        binerrors.push_back(fsr->GetBinError(bin));
        binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
        xerrors.push_back(0.0);
      }

      // Create TGraph
      int npoints = FSRvalues.size();
      TGraphErrors* graph = new TGraphErrors(npoints, &FSRvalues[0], &bincontents[0], &xerrors[0], &binerrors[0]);
      TF1 *fit = new TF1("fit", fitformula);
      graph->Fit("fit", "Q");
      vector<double> params;
      for(int i=0; i<Npar; i++) params.push_back(fit->GetParameter(i));

      // Here the fit functions are change to contain the values of the parameters and not the numbering e.g. [1]
      // stored in Vector for matrix multiplication
      TString FitWithParameters = fitformula;
      for(int i=0; i<Npar; i++) FitWithParameters.ReplaceAll("["+to_string(i)+"]", dtos(fit->GetParameter(i), precision));
      if(debug) cout << "\t" << FitWithParameters << endl;
      BinsFit.push_back(FitWithParameters);
      BinsData.push_back(data_norm->GetBinContent(bin));

      // Set Fit uncertainty
      TGraphErrors* band1 = (TGraphErrors*) graph->Clone();
      TFitter* fitter = (TFitter*) TVirtualFitter::GetFitter();
      fitter->GetConfidenceIntervals(band1, 0.68);
      TGraphErrors* band2 = new TGraphErrors(npoints, &FSRvalues[0], &bincontents[0], &xerrors[0], &binerrors_large[0]);


      fitparameters.push_back(params);
      graphs.push_back(graph);
      fits.push_back(fit);
      bands1.push_back(band1);
      bands2.push_back(band2);
    }

    if(debug) cout <<  "Size Vector fit function: " << BinsFit.size() << endl;
    if(debug) cout <<  "Size Vector data content: " << BinsData.size() << endl;

    if(debug) for(TString s: BinsFit) cout << s << "\t";
    if(debug) cout << endl;
    if(debug) for(double s: BinsData) cout << s << "\t";
    if(debug) cout << endl;

    for(int i=0; i<validbins.size();i++){
      PlotFit(graphs[i], fits[i], bands1[i], bands2[i], validbins[i], validbins[0], validbins[validbins.size()-1]);
    }

    // -------------------------------------------------------------------------------------------------
    // Now calculate the Chi2

    if(debug) cout << "\t ... Chi2" << endl;

    TString chi2formula;
    for(int i=0; i<validbins.size();i++){
      double dataentry = data_norm->GetBinContent(validbins[i]);
      double dataerror = data_norm->GetBinError(validbins[i]);
      // double fiterror = ttbar_norm->GetBinError(validbins[i]);
      double fiterror = FSRup_norm[FSRup_norm.size()-1]->GetBinError(validbins[i]);
      double error2 = pow(dataerror,2) + pow(fiterror,2);
      if(debug){
        cout << "----------------------------" << endl;
        cout << "Bin " << validbins[i] << endl;
        cout << "Data Error " << dataerror << endl;
        cout << "Fit  Error " << fiterror << endl;
      }

      for(auto s: sys) error2 += pow(s->GetBinContent(validbins[i]), 2);

      double error = sqrt(error2);
      TString modifiedformula = fitformula;
      for(int par=0; par<Npar; par++){
        TString paramstring = "["+to_string(par)+"]";
        modifiedformula.ReplaceAll(paramstring, to_string(fitparameters[i][par]));
      }

      TString formula = "(("+to_string(dataentry)+"-(" + modifiedformula +"))/"+to_string(error)+")^2";
      if(debug) cout << formula << endl;
      if(i==0) chi2formula = formula;
      else{
        chi2formula += "+";
        chi2formula += formula;
      }
    }
    double min = 0.01;
    double max = 100.;
    if(year == "2016"){
      min = 0.25;
      max = 4;
    }
    TF1* chi2function = new TF1("chi2", chi2formula, min, max);
    vector<double> Values = ExtractFSRValues(chi2function);
    f_fsr.push_back(Values);
    PlotChi2(chi2function, Values);

    cout << "--------------------------------------------------------" << endl;
    cout << "f(FSR) = " << Values[0] << " + " << Values[1] << " - " << Values[2] << ", chi2min = " << Values[5] << endl;
    cout << "--------------------------------------------------------" << endl;

    TString chi2formulaTEST = ConstructChi2(covTotalError, BinsFit, BinsData, mDim);
    TF1* chi2functionTEST = new TF1("chi2TEST", chi2formulaTEST, min, max);

    vector<double> ValuesTEST = ExtractFSRValues(chi2functionTEST);
    f_fsr.push_back(ValuesTEST);
    // PlotChi2(chi2functionTEST, ValuesTEST);

    cout << "--------------------------------------------------------" << endl;
    cout << "f(FSR) = " << ValuesTEST[0] << " + " << ValuesTEST[1] << " - " << ValuesTEST[2] << ", chi2min = " << ValuesTEST[5] << endl;
    cout << "--------------------------------------------------------" << endl;

    fstream fsr_txt;
    fsr_txt.open(save_path+"/fsr_factor"+channel+".txt", ios::out);
    fsr_txt << "f_FSR  = " << Values[0] << " + " << Values[1] << " - " << Values[2] << endl;
    fsr_txt << "f_up   = " << Values[3] << endl;
    fsr_txt << "f_down = " << Values[4] << endl;
    fsr_txt.close();

  }

  PlotResults(f_fsr, masswindows);
  return 0;

}

// #################################################################################################
// #################################################################################################
// #################################################################################################

void PlotFit(TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin){
  TString binnr = to_string(bin);
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graph->SetMarkerStyle(8);
  graph->SetTitle("Bin "+binnr);
  graph->GetXaxis()->SetTitle("(f^{FSR})^{2}");
  graph->GetYaxis()->SetTitle("a.u.");
  graph->Draw("AP");
  //fit->SetLineColor(kBlack);

  band2->SetLineColorAlpha(kOrange-4,0.5);
  band2->SetFillColorAlpha(kOrange-4,0.5);
  band1->SetLineColorAlpha(kAzure-4,0.8);
  band1->SetFillColorAlpha(kAzure-4,0.8);
  band2->Draw("L3 SAME");
  // band1->Draw("L3 SAME");

  graph->Draw("P SAME");
  fit->Draw("SAME");

  if(forTalk){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.6);
    text2->SetTextFont(62);
    text2->SetTextSize(0.05);
    text2->SetY(0.84);
    text2->Draw();

    TString preltext = "Work in Progress";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.6);
    text3->SetTextFont(52);
    text3->SetTextSize(0.035);
    text3->SetY(0.78);
    text3->Draw();
  }

  c->SaveAs(save_path+"/Bin"+binnr+".pdf");
  if(bin == firstbin)     c->Print(save_path+"/AllBins.pdf(","pdf");
  else if(bin == lastbin) c->Print(save_path+"/AllBins.pdf)","pdf");
  else                    c->Print(save_path+"/AllBins.pdf","pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotChi2(TF1* chi2function, vector<double> FSRvalues){
  double chimin = FSRvalues[5];
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  // gPad->SetLogx();
  chi2function->SetTitle(" ");
  if(year=="2016"){
    chi2function->GetXaxis()->SetRangeUser(0.25, 4);
  }
  else{
    chi2function->GetXaxis()->SetRangeUser(0.1, 25);
  }
  chi2function->GetYaxis()->SetRangeUser(chimin-1, chimin+9);
  chi2function->GetXaxis()->SetTitle("(f^{FSR})^{2}");
  chi2function->GetYaxis()->SetTitle("#chi^{2}");
  chi2function->GetYaxis()->SetTitleOffset(1.2);
  chi2function->Draw();
  double xmin = pow(FSRvalues[0]-FSRvalues[2], 2);
  double xmax = pow(FSRvalues[0]+FSRvalues[1], 2);

  TLine * l1 = new TLine(xmin, chimin-1.5, xmin, chi2function->Eval(xmin));
  TLine * l2 = new TLine(xmax, chimin-1.5, xmax, chi2function->Eval(xmax));
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l1->Draw("SAME");
  l2->Draw("SAME");

  if(forTalk){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.4);
    text2->SetTextFont(62);
    text2->SetTextSize(0.05);
    text2->SetY(0.84);
    text2->Draw();

    TString preltext = "Work in Progress";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.4);
    text3->SetTextFont(52);
    text3->SetTextSize(0.035);
    text3->SetY(0.78);
    text3->Draw();
  }

  c->SaveAs(save_path+"/chi2.pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

vector<double> ExtractFSRValues(TF1* chi2function){
  double minChi2      = chi2function->GetMinimum();
  double f_FSR2       = chi2function->GetX(minChi2, 0.001, 1000);
  double f_FSR        = sqrt(f_FSR2);
  double f_up2        = chi2function->GetX(minChi2+1, f_FSR2, 1000);
  double f_down2      = chi2function->GetX(minChi2+1, 0.001, f_FSR2);
  double f_up          = sqrt(f_up2);
  double f_down        = sqrt(f_down2);

  double sigmaup      = f_up - f_FSR;
  double sigmadown    = f_FSR - f_down;
  vector<double> values = {f_FSR, sigmaup, sigmadown, f_up, f_down, minChi2};
  return values;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotResults(vector<vector<double>> f_fsr, vector<TString> masswindows){
  vector<double> xvalues, xup, xdown, central, up, down;
  for(unsigned int i=0; i<f_fsr.size(); i++){
    central.push_back(f_fsr[i][0]);
    up.push_back(f_fsr[i][1]);
    down.push_back(f_fsr[i][2]);
    xvalues.push_back(i+1);
    xup.push_back(0.);
    xdown.push_back(0.);
  }

  double ymin = 0.0;
  double ymax = 5.0;

  TH1F* dummy = new TH1F("dummy", "dummy", f_fsr.size(), 0.5, 0.5+f_fsr.size());
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  for(unsigned int i=0; i<masswindows.size(); i++) dummy->GetXaxis()->SetBinLabel(i+1, masswindows[i]);
  dummy->SetTitle("");
  dummy->GetYaxis()->SetTitle("#it{f}^{FSR}");
  dummy->GetXaxis()->SetTitle("#it{m}_{jet}");
  dummy->GetYaxis()->SetRangeUser(ymin, ymax);
  dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->Draw();
  TGraphAsymmErrors * results = new TGraphAsymmErrors(f_fsr.size(), &xvalues[0], &central[0], &xdown[0], &xup[0], &down[0], &up[0]);
  results->SetMarkerStyle(8);
  results->Draw("P");
  TLine *errbandup   = new TLine(0.5, central[0]+up[0], 0.5+f_fsr.size(), central[0]+up[0]);
  TLine *errbanddown = new TLine(0.5, central[0]-down[0], 0.5+f_fsr.size(), central[0]-down[0]);
  errbandup->SetLineColor(13);
  errbanddown->SetLineColor(13);
  errbandup->SetLineStyle(2);
  errbanddown->SetLineStyle(2);
  errbandup->Draw("SAME");
  errbanddown->Draw("SAME");
  results->Draw("P SAME");

  double x = .7;
  for(int i=0; i<f_fsr.size();i++){
    // round min chi2 and store in TString
    std::ostringstream out;
    out.precision(2);
    out << std::fixed << f_fsr[i][5];
    TString label = "#chi_{min}^{2} = " + out.str();
    ////
    TLatex latex;
    latex.SetTextFont(62);
    latex.SetTextSize(0.02);
    latex.DrawLatex(x, .5, label);
    x += 1.0;
  }

  c->SaveAs(save_path+"/Results_"+year+".pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown){
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(.2);
  data->SetMarkerStyle(8);
  data->SetMarkerColor(1);
  data->SetLineColor(1);
  ttbar->SetLineWidth(2);
  ttbar->SetLineColor(14);
  ttbar->SetTitle(" ");
  ttbar->GetXaxis()->SetTitle("#tau_{32}");
  ttbar->GetYaxis()->SetTitle("a.u.");
  ttbar->GetYaxis()->SetRangeUser(0., 1.4*ttbar->GetMaximum());

  ttbar->Draw("HIST");
  Int_t color[] = {kAzure+7, kGreen-4, 810};
  int j = FSRdown.size()-1;
  for(int i=0; i<FSRup.size(); i++){
    FSRup[i]->SetLineWidth(2);
    FSRdown[j]->SetLineWidth(2);
    FSRup[i]->SetLineStyle(1);
    FSRdown[j]->SetLineStyle(2);
    FSRup[i]->SetLineColor(color[i]);
    FSRdown[j]->SetLineColor(color[i]);
    if(FSRup.size() == 1){
      FSRup[i]->Draw("HIST SAME");
      FSRdown[j]->Draw("HIST SAME");
    }
    if(i == 2){
      FSRup[i]->Draw("HIST SAME");
      FSRdown[j]->Draw("HIST SAME");
    }
    j--;
  }
  data->Draw("EP SAME");
  TLegend * leg = new TLegend(.25, .65, .50, .85);
  leg->AddEntry(data, "Data - background", "pe");
  leg->AddEntry(ttbar, "t#bar{t} nominal", "l");
  vector<TString> legnamesdown = {"t#bar{t} f^{FSR} = #frac{1}{4}", "t#bar{t} f^{FSR} = #frac{1}{2}", "t#bar{t} f^{FSR} = #frac{1}{#sqrt{2}}"};
  vector<TString> legnamesup = {"t#bar{t} f^{FSR} = #sqrt{2}", "t#bar{t} f^{FSR} = 2", "t#bar{t} f^{FSR} = 4"};
  if( year == "2016" ){
    legnamesdown = {"t#bar{t} f^{FSR} = #frac{1}{2}"};
    legnamesup = {"t#bar{t} f^{FSR} = 2"};
  }
  for(int i=0; i<FSRdown.size(); i++){
    if(i == 0) leg->AddEntry(FSRdown[i], legnamesdown[i], "l");
  }
  for(int i=0; i<FSRup.size(); i++){
    if(i == FSRup.size()-1) leg->AddEntry(FSRup[i], legnamesup[i], "l");
  }
  leg->Draw();
  gPad->RedrawAxis();

  if(forTalk){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.6);
    text2->SetTextFont(62);
    text2->SetTextSize(0.05);
    text2->SetY(0.84);
    text2->Draw();

    TString preltext = "Work in Progress";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.6);
    text3->SetTextFont(52);
    text3->SetTextSize(0.035);
    text3->SetY(0.78);
    text3->Draw();
  }

  c->SaveAs(save_path+"/Tau32.pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize){
  int nbins = data->GetXaxis()->GetNbins();
  TH1F* result = (TH1F*) data->Clone();
  for(unsigned int i=0; i<bgr.size(); i++){
    result->Add(bgr[i], -1);
  }
  for(int bin=1; bin<=nbins; bin++){
    double syserror2 = 0;
    for(unsigned int i=0; i<bgr.size(); i++){
      syserror2 += pow(bgr[i]->GetBinContent(bin) * syssize[i], 2);
    }
    double olderror2 = pow(result->GetBinError(bin), 2);
    result->SetBinError(bin, sqrt(syserror2+olderror2));
  }
  return result;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double diff1 = fabs( hist->GetBinContent(bin) - up->GetBinContent(bin) );
    double diff2 = fabs( hist->GetBinContent(bin) - down->GetBinContent(bin) );
    if(diff1 > diff2) sys->SetBinContent(bin, diff1);
    else              sys->SetBinContent(bin, diff2);
  }
  TCanvas *c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("#tau_{32}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(0., 0.005);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST");
  c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  delete c;
  return sys;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotError(TH1F* hist, TString sysname){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double error = hist->GetBinError(bin);
    sys->SetBinContent(bin, error);
  }
  TCanvas *c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("#tau_{32}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(0., 0.005);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST");
  c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TH2D* GetCovMatrix(TH1F* hist, double width, TString text, vector<bool> skip){

  // Int_t nbins = hist->GetNbinsX();
  Int_t nbins = 0;
  double xmin = 0.;

  for(bool s: skip){
    // Here we assume only the first bins are empty.
    // Check tau32 distribution if unsure
    if(!s) nbins++;
    else xmin += width;
  }

  // To skip bins in histogram
  Int_t shift = (Int_t) skip.size() - nbins;
  if(debug) cout << text << " used bins " << nbins << " and start at " << xmin << " and shift bin by " << shift << endl;

  TH2D* cov = new TH2D(text, "", nbins, xmin, 1., nbins, xmin, 1.);
  for(int x=1; x<=nbins; x++){
    cov->SetBinContent(x, x, pow(hist->GetBinError(x+shift),2));
  }
  return cov;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TH2D* NormCovMatrix(TH1* hist_, TH2* old_cov, int lower_, int upper_, bool width_){

  int mDim = old_cov->GetNbinsX();

  TH2D* new_cov = (TH2D*) old_cov->Clone();
  new_cov->Reset();

  double integral = hist_->Integral(lower_, upper_);

  for(int i=1; i <= mDim; i++){
    for(int j=1; j <= mDim; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= mDim; k++){
        for(int l=1; l <= mDim; l++){
          old_entry = old_cov->GetBinContent(k,l);
          double binwidth_i = hist_->GetBinWidth(i);
          double binwidth_j = hist_->GetBinWidth(j);
          if(!width_){
            binwidth_i = 1.0;
            binwidth_j = 1.0;
          }
          if(i==k) derivation_i = (integral - hist_->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (hist_->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - hist_->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (hist_->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      new_cov->SetBinContent(i, j, sum);
    }
  }
  return new_cov;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void DrawCov(TH2* m, TString text){
  TCanvas *B = new TCanvas(text, "", 800, 800); // name used to avoid Warning

  B->SetRightMargin(0.09);
  B->SetLeftMargin(0.15);
  B->SetRightMargin(0.20);

  m->Draw("COLZ");
  B->SaveAs("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/"+text+".pdf");
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TString ConstructChi2(TMatrixD m, vector<TString> vecS, VecD vecD, int nbins){
  TString term = "(di - gi) * (dj - gj) * Vij-1";
  TString function = "";
  TString temp_term = term;

  TMatrixD minvert = TMatrixD(m);
  minvert.Invert();

  TMatrixD mId = TMatrixD(m);
  mId *= minvert;

  // printCovMatrix(m, "Nominal Matrix", 100000);
  // printCovMatrix(minvert, "Inverted Matrix");
  // printCovMatrix(mId, "Identity Test");

  cout << __LINE__ << endl;
  for(unsigned int i = 0; i<nbins; i++){
    for(unsigned int j = 0; j<nbins; j++){
      if(abs(minvert[i][j])<=10e-30) continue;
      temp_term = term;
      temp_term.ReplaceAll("di", dtos(vecD[i], precision));
      temp_term.ReplaceAll("gi", vecS[i]);
      temp_term.ReplaceAll("dj", dtos(vecD[j], precision));
      temp_term.ReplaceAll("gj", vecS[j]);
      temp_term.ReplaceAll("Vij-1", dtos(minvert[i][j], precision));
      // cout << j << "\t" << term << "\t" << temp_term << endl;
      function += " + "+temp_term;
    }
  }
  // cout << function << endl;
  return function;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void printCovMatrix(TH2* m, TString name, double factor, int prec, bool diag){

  cout << endl << name << "Number of bins x: " << m->GetNbinsX() << "\t y: " << m->GetNbinsY() << "\t | factor " << factor << endl;
  cout << "y\\x \t\t";
  for(unsigned int x=1; x<=m->GetNbinsX(); x++){
    if(diag&&m->GetBinContent(x,x) == 0) continue;
    cout << x << "\t";
  }
  cout << endl << endl;

  for(unsigned int x=1; x<=m->GetNbinsY(); x++){
    if(diag&&m->GetBinContent(x,x) == 0) continue;
    cout << x << "\t\t";
    for(unsigned int y=1; y<=m->GetNbinsX(); y++){
      if(diag&&m->GetBinContent(y,y) == 0) continue;
      cout << Round(m->GetBinContent(x,y)*factor,prec) << "\t";
    }
    cout << endl;
  }
  cout << endl;
}


void printCovMatrix(TMatrixD m, TString name, double factor, int prec, bool diag){

  cout << endl << name << " - Number of bins x: " << m.GetNcols() << "\t y: " << m.GetNrows() << "\t | factor " << factor << endl;
  cout << "y\\x \t\t";
  if(m.GetNcols()!=m.GetNrows()) cerr << "Not a symmetric Matrix (" << name << ")" << endl;
  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << x << "\t";
  }
  cout << endl << endl;

  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << x << "\t\t";
    for(unsigned int y=0; y<m.GetNcols(); y++){
      if(diag&&m[y][y] == 0) continue;
      cout << Round(m[x][y]*factor,prec) << "\t";
    }
    cout << endl;
  }
  cout << endl;
}
