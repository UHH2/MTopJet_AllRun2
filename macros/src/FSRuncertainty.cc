#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"
#include "../include/normalise.h"
// #include "../include/CovarianzMatrix.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

#include "TDecompLU.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

using namespace std;

void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown);
void PlotControl(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown, TString var);
void PlotFit(double data, TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin);
void PlotChi2(TF1* chi2function, vector<double> FSRvalues, TString name="chi2");
void PlotError(TH1F* hist, TString sysname);
void PlotResults(vector<vector<double>> Values, vector<TString> windows);
void printCovMatrix(TMatrixD m, TString name, double factor=1, int prec=2, bool diag=false);
void DrawCov(TH2D* h, TString text, TString add = "");
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TH1F* AddBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname, TString add = "");
TH2D* TMatrixDtoTH2D(TMatrixD m, double nbins, double xmin, double xmax, TString name);
TString ConstructChi2(TMatrixD m, vector<TString> vecS, VecD vecD, int nbins, vector<bool> skip);
TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel=0, TString add = "");
TMatrixD GetCovMatrix(TH1F* hist, vector<bool> skip);
TMatrixD NormCovMatrix(TH1* hist_, TMatrixD old_cov, int lower_=0, int upper_=1, bool width_=false);
TMatrixD TrimMatrix(TMatrixD matrix, vector<bool> skip, bool removeBin = false);
vector<bool> TrimSkipBin(TMatrixD matrix, vector<bool> vec, bool onlyZero=false, bool removeBin = false);
vector<double> ExtractFSRValues(TF1* chi2function);

TString year, channel, ichannel, var;
TString save_path, save_general, save_nfs;

bool debug = false;
bool forTalk = true;

bool saveNFS = false;
// bool weights;

bool onlyDataUnc = false;
bool onlyFitUnc = false;

int factor = 1e+6;
int precision = 5; // For double to string function - dtos()
int binTOremove = -1;

vector<TString> bins = {
  "0.00 < #tau_{32} < 0.05", "0.05 < #tau_{32} < 0.10",
  "0.10 < #tau_{32} < 0.15", "0.15 < #tau_{32} < 0.20",
  "0.20 < #tau_{32} < 0.25", "0.25 < #tau_{32} < 0.30",
  "0.30 < #tau_{32} < 0.35", "0.35 < #tau_{32} < 0.40",
  "0.40 < #tau_{32} < 0.45", "0.45 < #tau_{32} < 0.50",
  "0.50 < #tau_{32} < 0.55", "0.55 < #tau_{32} < 0.60",
  "0.60 < #tau_{32} < 0.65", "0.65 < #tau_{32} < 0.70",
  "0.70 < #tau_{32} < 0.75", "0.75 < #tau_{32} < 0.80",
  "0.80 < #tau_{32} < 0.85", "0.85 < #tau_{32} < 0.90",
  "0.90 < #tau_{32} < 0.95", "0.95 < #tau_{32} < 1.00"
};

TCanvas *c_sys = new TCanvas("sys", "sys", 600, 600); // to Plot all systematics in on pdf

int main(int argc, char* argv[]){

  if(onlyDataUnc&&onlyFitUnc) throw runtime_error(to_string(__LINE__)+": You cant't look at only Data and only Fit at the same time.");

  gErrorIgnoreLevel = kError;
  gStyle->SetOptTitle(0);

  if(argc != 4){
    cout << "Usage: ./FSRuncertainty <year> <channel> <binstudy>" << endl;
    cout << "  year = 2016/2017/2018/combine" << endl;
    cout << "  channel  = muon/elec/combine" << endl;
    cout << "  binstudy = pt/mjet" << endl;
    return 0;
  }
  year = argv[1];
  ichannel = argv[2];
  if(ichannel.EqualTo("muon")) channel = "_muon";
  if(ichannel.EqualTo("elec")) channel = "_elec";
  if(ichannel.EqualTo("combine")) channel = "";
  var = argv[3];
  // weights = stob(argv[4]);
  binTOremove = 19; // stoi(argv[4]);
  // cout << binTOremove << endl;

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();

  cout << "\n-------------" << endl;
  cout << "year:    " << year << endl;
  cout << "channel: " << ichannel << endl;
  cout << "-------------" << endl;

  vector<TString> masswindows = {"140toInf", "140to160", "160to170", "170to180", "180to190", "190to200", "200toInf"}; // default
  // vector<TString> masswindows = {"140toInf"};
  // vector<TString> ptwindows = {"0toInf", "400to450", "450to500", "500to600", "600toInf"}; // default
  vector<TString> ptwindows = {
      "0toInf", "420to450",
      "400to450", "450to500", "500to600", "600toInf",
      "400to480", "480to580", "580toInf",
      "400to440", "440to480", "480to550", "550to610", "610to700", "700toInf"
  }; // full
  vector<vector<double>> f_fsr;

  vector<TString> ERROR = {"ERROR"};
  vector<TString> windows = var.EqualTo("mjet")?masswindows:var.EqualTo("pt")?ptwindows:ERROR;

  save_path = get_save_path()+"/FSRuncertainty/"+year+"/"+ichannel;
  // if(!weights) save_path += "/noWeights";
  save_general = save_path; // For mass/pt windows

  save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/FSR/";
  // if(!weights) save_nfs += "noWeights/";
  CreateSavePath((string) (save_nfs+"/"+year));

  TString option = gSystem->AccessPathName("files/Tau32.root")?"recreate":"update";
  TFile *f_tau32 = new TFile("files/Tau32.root", option);

  for(auto win: windows){
    cout << "Start with " << var << " " << win << " ... " << endl;
    saveNFS = ( (win.EqualTo("140toInf") || win.EqualTo("0toInf")) && ichannel.EqualTo("combine") )?true:false;
    if(saveNFS) cout << "\t ... Save main run on nfs" << endl;

    save_path  = save_general;
    save_path += "/"+var+"_"+win;
    CreateSavePath((string) save_path);
    if(debug) cout << save_path << endl;

    TFile* infile = new TFile("files/FSR_hists_"+var+win+".root");
    // if(weights) infile = new TFile("files/FSR_hists_"+var+win+"_weights.root");
    // else        infile = new TFile("files/FSR_hists_"+var+win+"_noweights.root");
    // TFile* infile = new TFile("files/FSR_hists_"+var+win+".root");

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
    TH1F *jecup = (TH1F*) infile->Get("jecup"+channel+"_"+year);
    TH1F *jecdown = (TH1F*) infile->Get("jecdown"+channel+"_"+year);
    TH1F *corup = (TH1F*) infile->Get("corup"+channel+"_"+year);
    TH1F *cordown = (TH1F*) infile->Get("cordown"+channel+"_"+year);
    TH1F *jmsup = (TH1F*) infile->Get("jmsup"+channel+"_"+year);
    TH1F *jmsdown = (TH1F*) infile->Get("jmsdown"+channel+"_"+year);
    TH1F *hdampup = (TH1F*) infile->Get("hdampup"+channel+"_"+year);
    TH1F *hdampdown = (TH1F*) infile->Get("hdampdown"+channel+"_"+year);
    TH1F *isrup = (TH1F*) infile->Get("isrup2"+channel+"_"+year);
    TH1F *isrdown = (TH1F*) infile->Get("isrdown2"+channel+"_"+year);
    TH1F *tuneup = (TH1F*) infile->Get("tuneup"+channel+"_"+year);
    TH1F *tunedown = (TH1F*) infile->Get("tunedown"+channel+"_"+year);
    TString cr_name = "qcdbased"; // argv[4]; gluonmove
    TH1F *cr = (TH1F*) infile->Get(cr_name+channel+"_"+year);

    if(debug) cout << "\t ... Get Hists for FSR" << endl;
    vector<TH1F*> FSRup, FSRdown;
    vector<double> FSRvalues;
    vector<TString> FSRnameUp, FSRnameDown;
    if(year == "2016"){
      // FSRvalues = {1./2, 1., 2.};
      FSRvalues = {0.25, 1., 4.}; // squared
      FSRnameUp = {"2"};
      FSRnameDown = {"2"};
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+year));
    }
    else{
      // FSRvalues = {1./4, 1./2, 1./sqrt(2), 1., sqrt(2), 2., 4.};
      FSRvalues = {0.0625, 0.25, 0.5, 1., 2., 4., 16.}; // squared
      // FSRvalues = {0.25, 0.5, 1., 2., 4.}; // squared
      FSRnameUp = {"sqrt2", "2", "4"};
      FSRnameDown = {"4", "2", "sqrt2"};
      FSRdown.push_back((TH1F*) infile->Get("fsrdown4"+channel+"_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdownsqrt2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrupsqrt2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup4"+channel+"_"+year));

    }

    if(debug) cout << "\t ... Get control Hists" << endl;
    TH1F *data_var = (TH1F*) infile->Get("data"+channel+"_"+var+"_"+year);
    TH1F *ttbar_var = (TH1F*) infile->Get("nominal"+channel+"_"+var+"_"+year);
    vector<TH1F*> bkg_var;
    bkg_var.push_back((TH1F*) infile->Get("bgr_st"+channel+"_"+var+"_"+year));
    bkg_var.push_back((TH1F*) infile->Get("bgr_wj"+channel+"_"+var+"_"+year));
    bkg_var.push_back((TH1F*) infile->Get("bgr_ot"+channel+"_"+var+"_"+year));

    vector<TH1F*> FSRup_var, FSRdown_var;
    if(year == "2016"){
      FSRdown_var.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+var+"_"+year));
      FSRup_var.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+var+"_"+year));
    }
    else{
      FSRdown_var.push_back((TH1F*) infile->Get("fsrdown4"+channel+"_"+var+"_"+year));
      FSRdown_var.push_back((TH1F*) infile->Get("fsrdown2"+channel+"_"+var+"_"+year));
      FSRdown_var.push_back((TH1F*) infile->Get("fsrdownsqrt2"+channel+"_"+var+"_"+year));
      FSRup_var.push_back((TH1F*) infile->Get("fsrupsqrt2"+channel+"_"+var+"_"+year));
      FSRup_var.push_back((TH1F*) infile->Get("fsrup2"+channel+"_"+var+"_"+year));
      FSRup_var.push_back((TH1F*) infile->Get("fsrup4"+channel+"_"+var+"_"+year));
    }

    // -------------------------------------------------------------------------------------------------
    // Subtract backgournd

    if(debug) cout << "\t ... Subtract Bkg" << endl;
    TH1F* data_sub = SubtractBackgrounds(data, bkg, bkgsys);
    if(debug) cout << "\t ... Subtract Bkg 2" << endl;
    TH1F* data_sub_var = SubtractBackgrounds(data_var, bkg_var, bkgsys);

    // Add bkg for Control
    if(debug) cout << "\t ... Subtract Bkg 3" << endl;
    TH1F* ttbar_bkg_var = AddBackgrounds(ttbar_var, bkg_var, bkgsys);

    // -------------------------------------------------------------------------------------------------
    // Rebin

    if(debug) cout << "\t ... Rebin" << endl;
    int rebin = 5;
    double width = 0.01*rebin;
    if(rebin > 1){
      data_sub->Rebin(rebin);
      ttbar->Rebin(rebin);
      jecup->Rebin(rebin);
      jecdown->Rebin(rebin);
      corup->Rebin(rebin);
      cordown->Rebin(rebin);
      jmsup->Rebin(rebin);
      jmsdown->Rebin(rebin);
      hdampup->Rebin(rebin);
      hdampdown->Rebin(rebin);
      isrup->Rebin(rebin);
      isrdown->Rebin(rebin);
      tuneup->Rebin(rebin);
      tunedown->Rebin(rebin);
      cr->Rebin(rebin);
      for(auto fsr: FSRup) fsr->Rebin(rebin);
      for(auto fsr: FSRdown) fsr->Rebin(rebin);

      data_sub_var->Rebin(2*rebin);
      ttbar_var->Rebin(2*rebin);
      for(auto fsr: FSRup_var) fsr->Rebin(2*rebin);
      for(auto fsr: FSRdown_var) fsr->Rebin(2*rebin);
      if(debug) cout << "Rebin with a factor of " << rebin << " - resulting in a bin width of " <<  data_sub->GetBinWidth(1) << endl;
    }

    // PlotControl(data_sub_var, ttbar_var, FSRup_var, FSRdown_var, var);

    // -------------------------------------------------------------------------------------------------
    // Skip Bins

    // Assumption: 0 entries (on diag) in CovMatrix lead to problems Inversion
    // Remove all bins where we cannot perform a fit -> data=0 or not enough MC points

    if(debug) cout << "Fill SkipBin ... " << endl;
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
      if(debug) cout << i << "\t Data " << skipData << "\t MC " << skipMC << " (" << MCs << ") " << "\t SkipBin: " << skip << "\tContent " << data_sub->GetBinContent(i) << endl;
      skipBin.push_back(skip);
    }

    // -------------------------------------------------------------------------------------------------
    // Create covariance matrix

    if(debug) cout << "Construct covariance matricies ... " << endl;
    TMatrixD m_covData = GetCovMatrix(data_sub, skipBin);
    TMatrixD m_covData_norm = NormCovMatrix(data_sub, m_covData);

    TMatrixD m_covttbar = GetCovMatrix(ttbar, skipBin);
    TMatrixD m_covttbar_norm = NormCovMatrix(ttbar, m_covttbar);

    if(debug) printCovMatrix(m_covData, "m_covData");
    if(debug) printCovMatrix(m_covData_norm, "m_covData_norm", 1.e6);

    if(debug) printCovMatrix(m_covttbar, "m_covttbar");
    if(debug) printCovMatrix(m_covttbar_norm, "m_covttbar_norm", 1.e6, 2);

    vector<TMatrixD> m_covFSRup, m_covFSRup_norm, m_covFSRup_norm_trim;
    vector<TMatrixD> m_covFSRdown, m_covFSRdown_norm, m_covFSRdown_norm_trim;

    for(unsigned int i=0; i<FSRup.size(); i++){

      m_covFSRup.push_back(GetCovMatrix(FSRup[i], skipBin));
      m_covFSRup_norm.push_back(NormCovMatrix(FSRup[i], m_covFSRup[i]));

      if(debug) printCovMatrix(m_covFSRup[i], ((TString) "m_covFSRup "+to_string(i)));
      if(debug) printCovMatrix(m_covFSRup_norm[i], ((TString) "m_covFSRup_norm "+to_string(i)), 1.e6);

      m_covFSRdown.push_back(GetCovMatrix(FSRdown[i], skipBin));
      m_covFSRdown_norm.push_back(NormCovMatrix(FSRdown[i], m_covFSRdown[i]));

      if(debug) printCovMatrix(m_covFSRup[i], ((TString) "m_covFSRup "+to_string(i)));
      if(debug) printCovMatrix(m_covFSRdown_norm[i], ((TString) "m_covFSRdown_norm "+to_string(i)), 1.e6);

    }

    // -----------------------------------------------------
    // ------------ Create covariance matrix for fit

    // Remove Zero-Bins from dimension
    if(debug) cout << "Construct covariance matrix for fit ... " << endl;
    double mDim = m_covData_norm.GetNcols();
    if(mDim!=m_covData_norm.GetNrows()) cerr << __LINE__ << ": Matrix is NOT symmetric" << endl;

    vector<double> diagValues;
    TMatrix m_covFit = TMatrixD(mDim, mDim);
    for(int x=0; x<mDim; x++){
      diagValues.push_back(m_covttbar_norm[x][x]);
      for(auto cov: m_covFSRdown_norm) diagValues.push_back(cov[x][x]);
      for(auto cov: m_covFSRup_norm) diagValues.push_back(cov[x][x]);
      // diagValues.pop_back();
      double max = *max_element(diagValues.begin(), diagValues.end());

      if(debug){ for(auto v: diagValues) cout << v << "\t"; }
      if(debug){ cout<<"Max value: "<<max<<endl; }

      m_covFit[x][x] = max;
      diagValues = {};
    }

    // -----------------------------------------------------
    // Create covariance matrix for systematics

    if(debug) cout << "Construct covariance matricies for systematics ... " << endl;
    vector<TMatrix> sys_cor_up, sys_cor_down;
    c_sys->Divide(3,2);
    sys_cor_up.push_back(GetSYSCov(ttbar, jecup,   "jec_up",   1, win));
    sys_cor_up.push_back(GetSYSCov(ttbar, corup,   "cor_up",   2, win));
    sys_cor_up.push_back(GetSYSCov(ttbar, jmsup,   "jms_up",   3, win));
    sys_cor_up.push_back(GetSYSCov(ttbar, hdampup, "hdamp_up", 4, win)); // ---
    sys_cor_up.push_back(GetSYSCov(ttbar, isrup,   "isr_up",   5, win));
    sys_cor_up.push_back(GetSYSCov(ttbar, tuneup,  "tune_up",  6, win));
    sys_cor_up.push_back(GetSYSCov(ttbar, cr,      cr_name,    6, win));
    c_sys->Print(save_path+"/SYS_all.pdf(");

    sys_cor_down.push_back(GetSYSCov(ttbar, jecdown,   "jec_down",   1, win));
    sys_cor_down.push_back(GetSYSCov(ttbar, cordown,   "cor_down",   2, win));
    sys_cor_down.push_back(GetSYSCov(ttbar, jmsdown,   "jms_down",   3, win));
    sys_cor_down.push_back(GetSYSCov(ttbar, hdampdown, "hdamp_down", 4, win)); // ---
    sys_cor_down.push_back(GetSYSCov(ttbar, isrdown,   "isr_down",   5, win));
    sys_cor_down.push_back(GetSYSCov(ttbar, tunedown,  "tune_down",  6, win));
    c_sys->Print(save_path+"/SYS_all.pdf)");

    vector<TMatrix> sys_cor_up_norm, sys_cor_down_norm;
    // vector<TH1F*> h_syst = {jecup, corup, jmsup, hdampup, isrup, tuneup};
    vector<TH1F*> h_syst = {jecup, corup, jmsup, hdampup, isrup, tuneup, cr};
    for(unsigned int i=0; i<sys_cor_up.size(); i++){
      sys_cor_up_norm.push_back(NormCovMatrix(h_syst[i], sys_cor_up[i]));
    }

    // -----------------------------------------------------
    // Create total covariance matrix

    if(debug) cout << "Construct total covariance matrix ... " << endl;
    TMatrix m_covTotal = TMatrixD(m_covData_norm);
    m_covTotal += m_covFit;
    for(auto m: sys_cor_up_norm) m_covTotal += m;

    // -----------------------------------------------------
    // Print and Plot covariance matricies

    if(debug) printCovMatrix(m_covTotal,               "Both", 1e5, 3);
    if(debug) printCovMatrix(m_covData_norm, "m_covData_norm", 1e6, 3);
    if(debug) printCovMatrix(m_covFit,         "m_covFit_all", 1e6, 3);

    TH2D* h_covData_norm  = TMatrixDtoTH2D(m_covData_norm,  20, 0, 1, "covData_norm_"+year+win);
    TH2D* h_covttbar_norm = TMatrixDtoTH2D(m_covttbar_norm, 20, 0, 1, "covttbar_norm_"+year+win);

    vector<TH2D*> h_covFSRup_norm, h_covFSRdown_norm;
    for(unsigned int i=0; i<FSRup.size(); i++){
      h_covFSRup_norm.push_back(TMatrixDtoTH2D(m_covFSRup_norm[i], 20, 0, 1,  "covFSRup_norm_"+year+win+to_string(i)));
      h_covFSRdown_norm.push_back(TMatrixDtoTH2D(m_covFSRdown_norm[i], 20, 0, 1,  "covFSRdown_norm_"+year+win+to_string(i)));
    }

    if(debug) cout << "\t ... Draw matricies" << endl;
    DrawCov(h_covData_norm, "covData_norm_"+year, win);
    DrawCov(h_covttbar_norm, "covttbar_norm_"+year, win);
    for(unsigned int i=0; i<FSRup.size(); i++){
      DrawCov(h_covFSRup_norm[i], "covFSRup_norm"+FSRnameUp[i]+"_"+year, win);
      DrawCov(h_covFSRdown_norm[i], "covFSRdown_norm"+FSRnameDown[i]+"_"+year, win);
    }

    // -------------------------------------------------------------------------------------------------
    // Normalize

    if(debug) cout << "Normalize hists ... " << endl;
    if(debug) cout << "\t ... Normalize" << endl;
    TH1F *data_norm = Normalize(data_sub);
    TH1F *ttbar_norm = Normalize(ttbar);
    TH1F *jecup_norm = Normalize(jecup);
    TH1F *jecdown_norm = Normalize(jecdown);
    TH1F *corup_norm = Normalize(corup);
    TH1F *cordown_norm = Normalize(cordown);
    TH1F *jmsup_norm = Normalize(jmsup);
    TH1F *jmsdown_norm = Normalize(jmsdown);
    TH1F *hdampup_norm = Normalize(hdampup);
    TH1F *hdampdown_norm = Normalize(hdampdown);
    TH1F *isrup_norm = Normalize(isrup);
    TH1F *isrdown_norm = Normalize(isrdown);
    TH1F *tuneup_norm = Normalize(tuneup);
    TH1F *tunedown_norm = Normalize(tunedown);
    TH1F *cr_norm = Normalize(cr);
    vector<TH1F*> FSRup_norm, FSRdown_norm;
    for(auto fsr: FSRup) FSRup_norm.push_back(Normalize(fsr));
    for(auto fsr: FSRdown) FSRdown_norm.push_back(Normalize(fsr));

    // -------------------------------------------------------------------------------------------------
    // Set bin error from covarianz matrix
    if(debug) cout << "\t ... Set correlated bin error" << endl;

    double dim = data_norm->GetNbinsX();
    for(unsigned int i=1; i<=dim; i++){
      data_norm->SetBinError(i, sqrt(m_covData_norm[i-1][i-1]));
      ttbar_norm->SetBinError(i, sqrt(m_covttbar_norm[i-1][i-1]));
      for(unsigned int j=0; j<FSRup_norm.size(); j++){
        FSRup_norm[j]->SetBinError(i, sqrt(m_covFSRup_norm.at(j)[i-1][i-1]));
        FSRdown_norm[j]->SetBinError(i, sqrt(m_covFSRdown_norm.at(j)[i-1][i-1]));
      }
    }

    // -------------------------------------------------------------------------------------------------
    // Save Histogram in Root file for PaperPlots

    if(win.EqualTo("140toInf") || win.EqualTo("0toInf")){
      if(channel.EqualTo("")){
        if(debug) cout << "\t ... Save hist in root file" << endl;
        TString option = gSystem->AccessPathName("files/PaperPlots.root")?"recreate":"update";
        TFile *f_out = new TFile("files/PaperPlots.root", option);
        f_out->cd();
        data_norm->Write("tau32__DATA__"+year,   TObject::kOverwrite);
        ttbar_norm->Write("tau32__TTbar__"+year, TObject::kOverwrite);
        for(unsigned int j=0;  j<FSRup_norm.size(); j++){
          FSRup_norm[j]->Write("tau32__TTbarFSRup"+FSRnameUp[j]+"__"+year, TObject::kOverwrite);
          FSRdown_norm[j]->Write("tau32__TTbarFSRdown"+FSRnameDown[j]+"__"+year, TObject::kOverwrite);
        }
        f_out->Close();
      }
    }

    // -------------------------------------------------------------------------------------------------
    // Plot
    if(debug) cout << "\t ... Draw tau32" << endl;

    PlotTau32(data_sub, ttbar, FSRup, FSRdown);
    PlotTau32(data_norm, ttbar_norm, FSRup_norm, FSRdown_norm);
    PlotError(data_norm, "DATA");
    PlotError(ttbar_norm, "TTbar");

    f_tau32->cd();
    data_norm->Write(var+"_"+win+"__DATA__norm__"+year, TObject::kOverwrite);
    data_sub->Write(var+"_"+win+"__DATA__"+year, TObject::kOverwrite);

    ttbar_norm->Write(var+"_"+win+"__TTbar__norm__"+year, TObject::kOverwrite);
    ttbar->Write(var+"_"+win+"__TTbar__"+year, TObject::kOverwrite);

    for(unsigned int j=0;  j<FSRup_norm.size(); j++){
      FSRup_norm[j]->Write(var+"_"+win+"__TTbarFSRup"+FSRnameUp[j]+"__norm__"+year, TObject::kOverwrite);
      FSRup[j]->Write(var+"_"+win+"__TTbarFSRup"+FSRnameUp[j]+"__"+year, TObject::kOverwrite);
      FSRdown_norm[j]->Write(var+"_"+win+"__TTbarFSRdown"+FSRnameDown[j]+"__norm__"+year, TObject::kOverwrite);
      FSRdown[j]->Write(var+"_"+win+"__TTbarFSRdown"+FSRnameDown[j]+"__"+year, TObject::kOverwrite);
    }

    // -------------------------------------------------------------------------------------------------
    // Get additional uncertainties
    if(debug) cout << "\t ... Plot systematics (jec,cor,jms)" << endl;
    vector<TH1F*> sys;
    sys.push_back(GetSYS(ttbar_norm, jecup_norm,   jecdown_norm,   "jec",   win));
    sys.push_back(GetSYS(ttbar_norm, corup_norm,   cordown_norm,   "cor",   win));
    sys.push_back(GetSYS(ttbar_norm, jmsup_norm,   jmsdown_norm,   "jms",   win));
    sys.push_back(GetSYS(ttbar_norm, hdampup_norm, hdampdown_norm, "hdamp", win)); // ----
    sys.push_back(GetSYS(ttbar_norm, isrup_norm,   isrdown_norm,   "isr",   win));
    sys.push_back(GetSYS(ttbar_norm, tuneup_norm,  tunedown_norm,  "tune",  win));

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
      bool ignorebin = skipBin[bin-1];
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
      if(bin!=binTOremove){
        BinsFit.push_back(FitWithParameters);
        BinsData.push_back(data_norm->GetBinContent(bin));
      }

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
      PlotFit(BinsData[i], graphs[i], fits[i], bands1[i], bands2[i], validbins[i], validbins[0], validbins[validbins.size()-1]);
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
    // f_fsr.push_back(Values);
    PlotChi2(chi2function, Values);

    double value_paper = 1/Values[0];
    double up_paper    = TMath::Sqrt(1/pow(Values[0],4)*pow(Values[1],2));
    double down_paper  = TMath::Sqrt(1/pow(Values[0],4)*pow(Values[2],2));
    cout << "----------------------------------------------------------------------------" << endl;
    cout << "1/f(FSR) = " << Values[0] << " + " << Values[1] << " - " << Values[2] << ", chi2min = " << Values[5] << " (unCOR)" << endl;
    cout << "  f(FSR) = " << value_paper << " + " << up_paper << " - " << down_paper << endl;
    cout << "----------------------------------------------------------------------------" << endl;

    TString chi2formulaCOR = ConstructChi2(m_covTotal, BinsFit, BinsData, mDim, skipBin);
    TF1* chi2functionCOR = new TF1("chi2COR", chi2formulaCOR, min, max);
    chi2functionCOR->Draw();

    vector<double> ValuesCOR = ExtractFSRValues(chi2functionCOR);
    f_fsr.push_back(ValuesCOR);
    PlotChi2(chi2functionCOR, ValuesCOR, "chi2Cor");

    value_paper = 1/ValuesCOR[0];
    up_paper    = TMath::Sqrt(1/pow(ValuesCOR[0],4)*pow(ValuesCOR[1],2));
    down_paper  = TMath::Sqrt(1/pow(ValuesCOR[0],4)*pow(ValuesCOR[2],2));
    cout << "----------------------------------------------------------------------------" << endl;
    cout << "1/f(FSR) = " << ValuesCOR[0] << " + " << ValuesCOR[1] << " - " << ValuesCOR[2] << ", chi2min = " << ValuesCOR[5] << " (COR)" << endl;
    cout << "  f(FSR) = " << value_paper << " + " << up_paper << " - " << down_paper << endl;
    cout << "----------------------------------------------------------------------------" << endl;

    fstream fsr_txt;
    fsr_txt.open(save_path+"/fsr_factor"+channel+".txt", ios::out);
    fsr_txt << "f_FSR  = " << Values[0] << " + " << Values[1] << " - " << Values[2] << endl;
    fsr_txt << "f_up   = " << Values[3] << endl;
    fsr_txt << "f_down = " << Values[4] << endl;
    fsr_txt.close();

    fstream fsr_txt_cor;
    fsr_txt_cor.open(save_path+"/fsr_factor_cor"+channel+".txt", ios::out);
    fsr_txt_cor << "f_FSR  = " << ValuesCOR[0] << " + " << ValuesCOR[1] << " - " << ValuesCOR[2] << endl;
    fsr_txt_cor << "f_up   = " << ValuesCOR[3] << endl;
    fsr_txt_cor << "f_down = " << ValuesCOR[4] << endl;
    fsr_txt_cor.close();
  }

  if(windows.size()>1){
    save_path=save_general;
    PlotResults(f_fsr, windows);
  }
  f_tau32->Close();
  return 0;

}

// #################################################################################################
// #################################################################################################
// #################################################################################################

void PlotFit(double data, TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin){
  TString binnr = to_string(bin);
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  graph->SetMarkerStyle(8);
  graph->SetTitle(bins[bin-1]);
  // graph->SetTitle("Bin "+binnr);
  graph->GetXaxis()->SetTitle("1/(#it{f}^{FSR})^{2}");
  graph->GetXaxis()->SetTitleOffset(1.1);
  graph->GetYaxis()->SetTitleOffset(1.9);
  graph->GetYaxis()->SetTitle("a.u.");
  graph->GetFunction("fit")->SetLineColor(kRed+2);
  graph->Draw("AP");
  // graph->Draw("AP");

  // double range = year=="2016"?4.4:17.5;
  double range = graph->GetXaxis()->GetXmax();
  TLine* mark_data = new TLine(0, data, range, data);
  mark_data->SetLineColor(kGray+2);
  mark_data->SetLineStyle(2);
  // mark_data->Draw("SAME");

  band2->SetLineColorAlpha(kOrange-4,0.5);
  band2->SetFillColorAlpha(kOrange-4,0.5);
  band1->SetLineColorAlpha(kAzure-4,0.8);
  band1->SetFillColorAlpha(kAzure-4,0.8);
  band2->Draw("L3 SAME");
  // band1->Draw("L3 SAME");

  graph->Draw("P SAME");
  fit->SetLineColor(kRed+2);
  // fit->Draw("SAME");

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

    // TString preltext = "Preliminary";
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
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/Bin"+binnr+".pdf");

  if(bin == firstbin)     c->Print(save_path+"/AllBins.pdf(","pdf");
  else if(bin == lastbin) c->Print(save_path+"/AllBins.pdf)","pdf");
  else                    c->Print(save_path+"/AllBins.pdf","pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotChi2(TF1* chi2function, vector<double> FSRvalues, TString name){
  double chimin = FSRvalues[5];
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  // gPad->SetLogx();
  chi2function->SetTitle("");
  if(year=="2016"){
    chi2function->GetXaxis()->SetRangeUser(0.25, 4);
    chi2function->GetYaxis()->SetRangeUser(chimin-1, chimin+9);
  }
  else{
    chi2function->GetXaxis()->SetRangeUser(0.1, 25);
    chi2function->GetYaxis()->SetRangeUser(chimin-3, chimin+15);
    // chi2function->GetXaxis()->SetRangeUser(-0.1, 0.1);
  }
  // chi2function->GetYaxis()->SetRangeUser(20, 40);
  chi2function->GetXaxis()->SetTitle("1/(#it{f}^{FSR})^{2}");
  chi2function->GetXaxis()->SetTitleOffset(1.1);
  chi2function->GetYaxis()->SetTitle("#chi^{2}");
  chi2function->GetYaxis()->SetTitleOffset(1.3);
  chi2function->SetLineColor(kRed+2);
  chi2function->Draw();
  double xmin = pow(FSRvalues[0]-FSRvalues[2], 2);
  double xmax = pow(FSRvalues[0]+FSRvalues[1], 2);

  double ymin = (year=="2016")?1.5:3.8;
  TLine * l1 = new TLine(xmin, chimin-ymin, xmin, chi2function->Eval(xmin));
  TLine * l2 = new TLine(xmax, chimin-ymin, xmax, chi2function->Eval(xmax));
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l1->Draw("SAME");
  l2->Draw("SAME");

  if(forTalk){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.6);
    text2->SetTextFont(62);
    text2->SetTextSize(0.05);
    text2->SetY(0.36);
    text2->Draw();

    // TString preltext = "Preliminary";
    TString preltext = "Work in Progress";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.6);
    text3->SetTextFont(52);
    text3->SetTextSize(0.035);
    text3->SetY(0.3);
    text3->Draw();
  }

  c->SaveAs(save_path+"/"+name+".pdf");
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/"+name+".pdf");
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

// vector<TString> ptwindows = {
//     "0toInf", "420to450",
//     "400to450", "450to500", "500to600", "600toInf",
//     "400to480", "480to580", "580toInf",
//     "400to440", "440to480", "480to550", "550to610", "610to700", "700toInf"
// };

void PlotResults(vector<vector<double>> f_fsr, vector<TString> windows){
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
  TCanvas *c = new TCanvas("c", "c", 1500, 600);
  for(unsigned int i=0; i<windows.size(); i++) dummy->GetXaxis()->SetBinLabel(i+1, windows[i]);
  dummy->SetTitle("");
  dummy->GetYaxis()->SetTitle("1/#it{f}^{FSR}");
  TString xtitle = (var=="pt")?"#it{p}_{T}":"#it{m}_{jet}";
  dummy->GetXaxis()->SetTitle(xtitle);
  dummy->GetYaxis()->SetRangeUser(ymin, ymax);
  dummy->GetYaxis()->SetTitleOffset(1.3);
  dummy->Draw();
  TGraphAsymmErrors * results = new TGraphAsymmErrors(f_fsr.size(), &xvalues[0], &central[0], &xdown[0], &xup[0], &down[0], &up[0]);
  results->SetMarkerStyle(8);
  results->Draw("P");
  TLine *errbandup   = new TLine(0.5, central[0]+up[0], 0.5+f_fsr.size(), central[0]+up[0]);
  TLine *errbanddown = new TLine(0.5, central[0]-down[0], 0.5+f_fsr.size(), central[0]-down[0]);
  for(int i: {0,1,5,8}){
    TLine *line = new TLine(xvalues[i]+0.5, 0, xvalues[i]+0.5, 5);
    line->SetLineColor(13);
    line->SetLineStyle(1);
    line->Draw("SAME");
  }
  errbandup->SetLineColor(13);
  errbanddown->SetLineColor(13);
  errbandup->SetLineStyle(2);
  errbanddown->SetLineStyle(2);
  errbandup->Draw("SAME");
  errbanddown->Draw("SAME");
  results->Draw("P SAME");
  gPad->RedrawAxis();

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

  c->SaveAs(save_path+"/Results_"+var+"_"+year+".pdf");
  c->SaveAs(save_nfs+"/Results_"+var+"_"+year+".pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown){
  TString xtitle = data->GetMaximum()<1?"a.u.":"Events";
  TString add = data->GetMaximum()<1?"_norm":"";
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(.2);
  data->SetMarkerStyle(8);
  data->SetMarkerColor(1);
  data->SetLineColor(1);
  ttbar->SetLineWidth(2);
  ttbar->SetLineColor(14);
  ttbar->SetTitle(" ");
  ttbar->GetXaxis()->SetTitle("#tau_{32}");
  ttbar->GetYaxis()->SetTitle(xtitle);
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
  vector<TString> legnamesup = {"t#bar{t} f^{FSR} = #frac{1}{#sqrt{2}}", "t#bar{t} f^{FSR} = #frac{1}{2}", "t#bar{t} f^{FSR} = #frac{1}{4}"};
  vector<TString> legnamesdown = {"t#bar{t} f^{FSR} = 4", "t#bar{t} f^{FSR} = 2", "t#bar{t} f^{FSR} = #sqrt{2}"};
  if( year == "2016" ){
    legnamesup = {"t#bar{t} f^{FSR} = #frac{1}{2}"};
    legnamesdown = {"t#bar{t} f^{FSR} = 2"};
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

    // TString preltext = "Preliminary";
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

  c->SaveAs(save_path+"/Tau32"+add+".pdf");
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/Tau32"+add+".pdf");
  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void PlotControl(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown, TString var){

  TString xtitle = (var=="pt")?"#p_{T}":"#m_{jet}";

  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(.2);
  data->SetMarkerStyle(8);
  data->SetMarkerColor(1);
  data->SetLineColor(1);
  ttbar->SetLineWidth(2);
  ttbar->SetLineColor(14);
  ttbar->SetTitle(" ");
  ttbar->GetXaxis()->SetTitle(xtitle);
  ttbar->GetYaxis()->SetTitle("Events");
  ttbar->GetYaxis()->SetRangeUser(0., 1.2*ttbar->GetMaximum());

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

  TLegend * leg = new TLegend(.65, .65, .85, .85);
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

  c->SaveAs(save_path+"/"+var+".pdf");
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/"+var+".pdf");
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
  for(int bin=1; bin<=nbins; bin++){
    if(result->GetBinContent(bin)<0){
      result->SetBinContent(bin, 0);
      result->SetBinError(bin, 0);
    }
  }

  return result;
}

TH1F* AddBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize){
  int nbins = data->GetXaxis()->GetNbins();
  TH1F* result = (TH1F*) data->Clone();
  for(unsigned int i=0; i<bgr.size(); i++){
    result->Add(bgr[i], 1);
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

TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname, TString add){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double diff1 = fabs( hist->GetBinContent(bin) - up->GetBinContent(bin) );
    double diff2 = fabs( hist->GetBinContent(bin) - down->GetBinContent(bin) );
    if(diff1 > diff2) sys->SetBinContent(bin, diff1);
    else              sys->SetBinContent(bin, diff2);
  }
  TCanvas *c = new TCanvas(sysname+add, "", 600, 600);
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
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/SYS_"+sysname+".pdf");
  delete c;
  return sys;
}

TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel, TString add){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  vector<double> binUnc;
  for(int i=1; i<=nbins; i++){
    double diff = hist->GetBinContent(i) - var->GetBinContent(i);
    sys->SetBinContent(i, diff);
    binUnc.push_back(diff);
  }

  TMatrix matrix = TMatrixD(nbins, nbins);
  for(int i=0; i<nbins; i++){
    for(int j=0; j<nbins; j++){
      matrix[i][j] = binUnc[i]*binUnc[j];
    }
  }

  TH2D* h_norm = TMatrixDtoTH2D(matrix, 20, 0, 1, "cov_"+sysname+"_"+year+to_string(panel)+add);
  DrawCov(h_norm, "cov_"+sysname+"_"+year, to_string(panel)+add);

  TCanvas *c = new TCanvas(sysname+add, "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("#tau_{32}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(-1000, 1000);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST ][");
  c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/SYS_"+sysname+".pdf");

  c_sys->cd(panel);
  sys->Draw("HIST ][");

  delete c;
  return matrix;
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
  if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/SYS_"+sysname+".pdf");

  delete c;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TMatrixD GetCovMatrix(TH1F* hist, vector<bool> skip){

  Int_t nbins = hist->GetNbinsX();
  Int_t mDim = 0;
  std::vector<int> indexEmptyBins;
  for(unsigned int i=1; i<=nbins; i++){
    if(hist->GetBinContent(i)>0) mDim++;
    else indexEmptyBins.push_back(i);
  }

  if(debug) for(unsigned int i=0; i<indexEmptyBins.size(); i++) cout << i << "\t";
  if(debug) cout << endl;

  TMatrixD matrix = TMatrixD(nbins, nbins);
  for(int x=1; x<=nbins; x++){
    matrix[x-1][x-1] = (hist->GetBinContent(x)>0)?pow(hist->GetBinError(x),2):0.;
  }
  return matrix;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TMatrixD NormCovMatrix(TH1* hist_, TMatrixD old_cov, int lower_, int upper_, bool width_){

  int mDim = old_cov.GetNrows();
  if(mDim!=20) throw runtime_error("WRONG DIM");

  TMatrixD new_cov = TMatrix(mDim, mDim);

  double integral = hist_->Integral();

  for(int i=1; i <= mDim; i++){
    for(int j=1; j <= mDim; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= mDim; k++){
        for(int l=1; l <= mDim; l++){
          old_entry = old_cov[k-1][l-1];
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
      new_cov[i-1][j-1] = sum;
    }
  }
  return new_cov;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TMatrixD TrimMatrix(TMatrixD matrix, vector<bool> skip, bool removeBin){

  Int_t nbins = matrix.GetNrows();
  vector<Int_t> usedEntries;
  for(int x=0; x<nbins; x++){
    // cout << (matrix[x][x]!=0) << "\t" << !skip[x] << endl;
    if(removeBin){
      if(matrix[x][x]!=0&&!skip[x]&&(x!=binTOremove-1)) usedEntries.push_back(x);
    }
    else{
      if(matrix[x][x]!=0&&!skip[x]) usedEntries.push_back(x);
    }
  }

  Int_t dim = usedEntries.size();
  if(debug) cout << "New Dim is " << dim << endl;
  TMatrixD trimmed = TMatrixD(dim, dim);
  for(int x=0; x<dim; x++){
    Int_t i = usedEntries[x];
    for(int y=0; y<dim; y++){
      Int_t j = usedEntries[y];
      trimmed[x][y] = matrix[i][j];
    }
  }

  return trimmed;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

vector<bool> TrimSkipBin(TMatrixD matrix, vector<bool> vec, bool onlyZero, bool removeBin){

  Int_t nbins = matrix.GetNrows();
  vector<bool> vecTrim = vec;
  for(int x=nbins-1; x>=0; x--){
    if(onlyZero){
      if(removeBin){
        if(matrix[x][x]==0||x==binTOremove-1) vecTrim.erase(vecTrim.begin()+x);
      }
      else{
        if(matrix[x][x]==0) vecTrim.erase(vecTrim.begin()+x);
      }
    }
    else{
      if(removeBin){
        if(matrix[x][x]==0||vec[x]||x==binTOremove) vecTrim.erase(vecTrim.begin()+x);
      }
      else{
        if(matrix[x][x]==0||vec[x]) vecTrim.erase(vecTrim.begin()+x);
      }
    }
  }

  if(debug){
    for(int i: vec) cout << i << "\t"; cout << endl;
    for(int i: vecTrim) cout << i << "\t"; cout << endl;

    vector<int> test;
    for(int x=0; x<nbins; x++) test.push_back(x);
    for(int i: test) cout << i << "\t"; cout << endl;
    for(int x=nbins-1; x>=0; x--){
      if(onlyZero){
        if(matrix[x][x]==0) test.erase(test.begin()+x);
      }
      else{
        if(matrix[x][x]==0||vec[x]) test.erase(test.begin()+x);
      }
    }
    for(int i: test) cout << i << "\t"; cout << endl;
  }

  return vecTrim;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void DrawCov(TH2D* h, TString text, TString add){

  TCanvas *B = new TCanvas(text+add, "", 800, 800); // name used to avoid Warning

  B->SetRightMargin(0.09);
  B->SetLeftMargin(0.15);
  B->SetRightMargin(0.20);

  h->GetXaxis()->SetTitle("#tau_{32}");
  h->GetYaxis()->SetTitle("#tau_{32}");
  TString norm = (text.Contains("norm")?"normalized ":"");
  TString ztitle = norm + "uncertainty";
  h->GetZaxis()->SetTitle(ztitle);

  h->Draw("COLZ");
  B->SaveAs(save_path+"/"+text+".pdf");
  if(saveNFS) B->SaveAs(save_nfs+"/"+year+"/"+text+".pdf");

  h->Reset();
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TString ConstructChi2(TMatrixD m, vector<TString> vecS, VecD vecD, int nbins, vector<bool> skipBin){
  TString term = "(di - (gi)) * (dj - (gj)) * Vij-1";
  TString function = "";
  TString temp_term = term;

  vector<bool> dummy;
  for(unsigned int i=0; i<nbins; i++){dummy.push_back(false);}
  TMatrixD m_trim = TrimMatrix(m, skipBin, true); // Cut of diagonals with 0 for .Invert() and bins without data
  TString sizeM = to_string(m_trim.GetNcols());
  // cout << m_trim.GetNcols() << "   " << m_trim.GetNrows() << endl;
  if(debug) printCovMatrix(m, "standard", factor, 3);
  if(debug) printCovMatrix(m_trim, "trimmed", factor, 3);

  TMatrixF minvert = TMatrixF(m_trim);
  minvert.SetTol(1.e-30);
  minvert.Invert();

  TMatrixF mId = TMatrixF(m_trim);
  mId *= minvert;

  if(debug) printCovMatrix(m_trim, "Nominal Matrix", 1.e8, 3);
  if(debug) printCovMatrix(minvert, "Inverted Matrix", factor);
  if(debug) printCovMatrix(mId, "Identity Test", factor);

  TMatrixD m_final = TMatrixD(minvert);
  if(debug) cout << term << endl;
  nbins=m_final.GetNcols();
  if(debug) cout << "Final Check --- Matrix: " << nbins << "\t | data bin content: " << vecD.size() << "\t | terms: " << vecS.size() << endl;
  for(unsigned int i = 0; i<nbins; i++){
    for(unsigned int j = 0; j<nbins; j++){
      if(abs(m_final[i][j])<=10e-30) continue;
      temp_term = term;
      temp_term.ReplaceAll("di", dtos(vecD[i], precision));
      temp_term.ReplaceAll("gi", vecS[i]);
      temp_term.ReplaceAll("dj", dtos(vecD[j], precision));
      temp_term.ReplaceAll("gj", vecS[j]);
      temp_term.ReplaceAll("Vij-1", dtos(m_final[i][j], precision));
      if(debug) cout << j << "\t" << temp_term << endl;
      function += " + "+temp_term;
    }
  }
  return function;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

void printCovMatrix(TMatrixD m, TString name, double factor, int prec, bool diag){

  Int_t space = prec+8;
  cout << endl << name << " - Number of bins x: " << setw(space) << m.GetNcols() << " y: " << setw(space) << m.GetNrows() << " | factor " << setw(space) << factor << endl << endl;
  cout << setw(space) << "y\\x" << " | ";
  if(m.GetNcols()!=m.GetNrows()) cerr << "Not a symmetric Matrix (" << name << ")" << endl;
  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << setw(space) << x << " | ";
  }
  cout << endl;
  for(unsigned int x=0; x<=m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    TString s = "_";
    for(int i=0; i<space; i++) cout << "-";
    cout << " | ";
  }
  cout << endl;

  for(unsigned int x=0; x<m.GetNcols(); x++){
    if(diag&&m[x][x] == 0) continue;
    cout << setw(space) << x << " | ";
    for(unsigned int y=0; y<m.GetNrows(); y++){
      if(x==y){
        cout << RED<< setw(space) << Round(m[x][y]*factor, prec) << RESET << " | ";
      }
      else cout << setw(space) << Round(m[x][y]*factor, prec) << " | ";
    }
    cout << endl;
  }
  cout << endl;
}

// -------------------------------------------------------------------------------------------------
// ---                                                                                           ---
// -------------------------------------------------------------------------------------------------

TH2D* TMatrixDtoTH2D(TMatrixD m, double nbins, double xmin, double xmax, TString name){
  TH2D* h = new TH2D(name, "", nbins, xmin, xmax, nbins, xmin, xmax);
  for(unsigned int c=0; c<m.GetNcols(); c++){
    for(unsigned int r=0; r<m.GetNcols(); r++){
      h->SetBinContent(c+1, r+1, m[c][r]);
    }
  }
  return h;
}
