#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"
#include "../include/tdrstyle_all.h"

#include "TRandom2.h"
#include "TF3.h"
#include "TError.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"

// #include "../CovMatrices/JMS/CollectCovHeaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>
#include "TSystem.h"

using namespace std;

typedef map<int, TGraph2DErrors*> MapI2DGe;
typedef map<TString, MapI2DGe> MapSI2DGe;

typedef map<TString, TMatrixD> MapM;
typedef map<TString, MapM> MapMM;
typedef map<TString, MapMM> MapMMM;

typedef map<TString, TF2*> MapF2;
typedef map<TString, MapF2> MapFF2;
typedef map<int, TF2*> MapIF2;
typedef map<TString, MapIF2> MapSIF2;

// ---------------------------------------------------------------------------
// Declare functions

// MapHHH SetUpMap(VecTS hists);
void printVector(VecD& vector, TString info);
void printVector(VecI& vector, TString info);
void Plot2DGraph(TGraph2DErrors* bin_fit, int bin, TString title, VecD content, VecD error);
void Draw2DChi2(TF2* chi2, VecDD& points, VecDD& extrema, VecD& minimum, TString hist);
void drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange);
void PoUinTXT(double& nom, double& up, double& down, TString variation);
void PrintValues(VecD& nominal, VecD& uu, VecD& ud, VecD& du, VecD& dd, TString name, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms);
void SubtractBackgroundsMap(MapHHH& m_data, MapHHH& m_st, MapHHH& m_wjets, MapHHH& m_other);
void AddTwoMaps(MapHHH& map1, MapHHH& map2, int option);
bool sortcolx( const VecD& v1, const VecD& v2 );
TString GetChi2Cor(TMatrixD &m, vector<TString> &vecS, VecD &vecD, int nbins, vector<bool> &skipBin, TString wbin="");
MapF2 GetFits(MapSI2DGe& map, MapVI& peak, TString& hists);
MapHH GetHistograms(TString process, TString h_name);
MapH RebinAndNormalize(MapHHH& process, int& width);
MapH RebinMap(MapHHH& process, int& width);
MapH NormalizeMap(MapH& process);
MapM GetCovMatrixMap(MapH& map, const TString& save, const TString& process);
MapM GetCovMatrixNormMap(MapM& map, MapH& hists, const TString& save, const TString& process);
TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel=0, TString add = "");
MapM GetSYSCovMap(MapH& hists, MapH& vars, TString sysname, TString bin, int panel=0);
MapSI2DGe Creat2DGraph(MapHH& map, MapVI& peak, VecTS& hists);
void Creat3DFit(MapHH& map, MapVI& peak, VecTS& hists);
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TF1* GetSplitFunction(VecDD& ellipse);
TF2* GetChi2(MapF2& fits, MapHH& nhists, MapVI& peak, TString hist);
TEllipse* ApproxEllipse(const VecD& mid, const VecD& uu, const VecD& dd, const VecD& ud, const VecD& du);
TPolyLine3D* SetPolyLine(double x1, double x2, double y1, double y2, double z1, double z2, int color);
VecD ExtractJMSValues(const TF1* chi2function, const double& var);
VecD GetMinimumChi2(TF2* chi2, TString hist);
VecD StoreBinContentInVector(MapH& map, MapVI& peak, TString& hist);
VecDD GetSigmaEllipse(TF2* chi2, VecD& minimum, TString hist, double runs, double acc, double at=2.3);
VecDD GetSigmaEllipse_alt(TF2* chi2, VecD &minimum, TString hist, double runs, double acc, double at=2.3);
VecDD GetExtreme(VecDD& ellipseJMS, VecD& nominal_JMS);
VecTS GetFitFunctions(MapF2& fits, MapVI& peak, TString& hist);
vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist);

vector<vector<double>> FindXY_alt(TF2 *function, double zfix, double xmin, double xmax, double ymin, double ymax, int steps = 1000, double accuracy = 0, bool count=false);

TCanvas *c_sys = new TCanvas("sys", "sys", 600, 600); // to Plot all systematics in on pdf

bool debug = false;
bool fast = true; // Reads Cov Matrix for data and SYS from root file
bool NormToMC = false;
bool onlyLin = true;
bool addSyst = false;
TString save, year, channel;
TString save_path, save_general, save_nfs;
TString schannel = "";
TString addition = "";
TString fsr_sys = "2";
VecTS hists = {
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hl",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_lh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_ll"
};
TString cut = "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_";
int space;
int precision = 10; // For double to string function - dtos()
int cout_precision = 6;

double ymax_global = 0; // TODO

inline TString wbin(TString h){return h.ReplaceAll(cut, "");}

int main(int argc, char* argv[]){

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  // Input ---------------------------------------------------------------------
  year = argv[1];
  channel = argv[2];
  addSyst = stob(argv[3]);
  fast = stob(argv[4]);
  NormToMC = stob(argv[5]);

  cout << year << "   " << channel << "    syst " << addSyst << "    fast " << fast << "    norm to MC " << NormToMC << endl;
  // if(addSyst) addition = "syst_noJER";
  // if(addSyst) addition = "syst";
  if(addSyst){
    addition = NormToMC?"syst_ttbar":"syst_data";
  }
  if(channel.EqualTo("muon")||channel.EqualTo("elec")) schannel = "_"+channel;
  // TString sub = argv[3]; bool isSUB = false;
  // if(sub.EqualTo("1")) isSUB = true;
  // cout << isSUB << endl;
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptTitle(0);
  SetupGlobalStyle();

  // year ---------------------------------------------------------------------
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16  = true;
  else if(strcmp(year, "2017")==0)     is17  = true;
  else if(strcmp(year, "2018")==0)     is18  = true;
  else if(strcmp(year, "combine")==0)  isAll = true;
  // else cerr("Give me the correct year please (2016, 2017, 2018 or combine)");

  // Bins ----------------------------------------------------------------------
  int bin_width = 1;
  int number_bins = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders

  // Limit for bin content -----------------------------------------------------
  int limit = (isAll)?100:40;
  TString slimit = to_string(limit);


  // Create Directories --------------------------------------------------------
  if(debug) cout << "Create Directiories ..." << endl;

  save = get_save_path();
  save += "/JetCorrections/fit/"+year+"/"+channel+"/"+slimit+"/"+"BinWidth_"+to_string(bin_width)+"/"+addition+"/";
  TString temp_path = save+"{all_bins,projection}"; // Creats Subdirs x and y for {x,y}
  CreateSavePath((string) temp_path);
  cout << save << endl;

  save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/"+year+"/"+channel+"/"+slimit+"/"+addition+"/";
  temp_path = save_nfs+"SYS/"+year;
  CreateSavePath((string) save_nfs);
  CreateSavePath((string) temp_path);

  cout << "Start collecting Histograms ..." << endl;
  // ---------------------------------------------------------------------------
  // I load all histograms; This code was created after studies ended. Only one
  // option is available. If necessary include other options as well (no bins etc.)

  TString w_mass_hh = hists[0];
  TString w_mass_hl = hists[1];
  TString w_mass_lh = hists[2];
  TString w_mass_ll = hists[3];

  MapHHHH histograms;
  MapHHH m_data, m_ttbar, m_st, m_wjets, m_other, m_JECup, m_JECdown, m_CORup, m_CORdown, m_JERup, m_JERdown;
  MapHHH m_hdampup, m_tuneup, m_hdampdown, m_tunedown, m_CRqb, m_CRgm;
  MapHHH m_fsrup_sqrt2, m_isrup_sqrt2, m_fsrdown_sqrt2, m_isrdown_sqrt2;
  MapHHH m_fsrup_2, m_isrup_2, m_fsrdown_2, m_isrdown_2;
  MapHHH m_fsrup_4, m_isrup_4, m_fsrdown_4, m_isrdown_4;

  // MapHHH GetHistograms(channel, muon, elec, combine, nhist)
  if(debug) cout << "\t ... Get Histograms" << endl;
  for(TString hist: hists)
  {
    m_data[hist] = GetHistograms("data", hist);
    m_ttbar[hist] = GetHistograms("ttbar", hist);
    m_st[hist] = GetHistograms("st", hist);
    m_wjets[hist] = GetHistograms("wjets", hist);
    m_other[hist] = GetHistograms("other", hist);
    m_JECup[hist] = GetHistograms("JECup", hist);
    m_JECdown[hist] = GetHistograms("JECdown", hist);
    m_CORup[hist] = GetHistograms("CORup", hist);
    m_CORdown[hist] = GetHistograms("CORdown", hist);
    m_JERup[hist] = GetHistograms("JERup", hist);
    m_JERdown[hist] = GetHistograms("JERdown", hist);
    m_tuneup[hist] = GetHistograms("tuneup", hist);
    m_tunedown[hist] = GetHistograms("tunedown", hist);
    m_hdampup[hist] = GetHistograms("hdampup", hist);
    m_hdampdown[hist] = GetHistograms("hdampdown", hist);
    m_CRqb[hist] = GetHistograms("CRqb", hist);
    m_CRgm[hist] = GetHistograms("CRgm", hist);
    m_fsrup_sqrt2[hist] = GetHistograms("fsrup_sqrt2", hist);
    m_fsrdown_sqrt2[hist] = GetHistograms("fsrdown_sqrt2", hist);
    m_fsrup_2[hist] = GetHistograms("fsrup_2", hist);
    m_fsrdown_2[hist] = GetHistograms("fsrdown_2", hist);
    m_fsrup_4[hist] = GetHistograms("fsrup_4", hist);
    m_fsrdown_4[hist] = GetHistograms("fsrdown_4", hist);
    m_isrup_sqrt2[hist] = GetHistograms("isrup_sqrt2", hist);
    m_isrdown_sqrt2[hist] = GetHistograms("isrdown_sqrt2", hist);
    m_isrup_2[hist] = GetHistograms("isrup_2", hist);
    m_isrdown_2[hist] = GetHistograms("isrdown_2", hist);
    m_isrup_4[hist] = GetHistograms("isrup_4", hist);
    m_isrdown_4[hist] = GetHistograms("isrdown_4", hist);
  }

  MapHHH m_fsrup, m_isrup, m_fsrdown, m_isrdown;
  if(fsr_sys.EqualTo("sqrt2")){
    m_fsrup = m_fsrup_sqrt2;
    m_isrup = m_isrup_sqrt2;
    m_fsrdown = m_fsrdown_sqrt2;
    m_isrdown = m_isrdown_sqrt2;
  }
  else if(fsr_sys.EqualTo("2")){
    m_fsrup = m_fsrup_2;
    m_isrup = m_isrup_2;
    m_fsrdown = m_fsrdown_2;
    m_isrdown = m_isrdown_2;
  }
  else if(fsr_sys.EqualTo("4")){
    m_fsrup = m_fsrup_4;
    m_isrup = m_isrup_4;
    m_fsrdown = m_fsrdown_4;
    m_isrdown = m_isrdown_4;
  }
  else{
    if(addSyst) throw runtime_error("<E> You are adding systematics to fit but not defined a correct FSR/ISR variation (sqrt2,2,4)");
  }

  if(debug) cout << "\t ... Get empty/peak bins" << endl;
  // TODO: adjust limits for non combined years
  // int limit = (channel.EqualTo("muon")||channel.EqualTo("elec"))?50:100;
  // int limit = 100;

  MapVI empty_bins, peak_bins;
  for(TString h: hists){
    TString c = channel; TString y = year;
    empty_bins[h] = bins_empty(m_data[h][c][y]);
    peak_bins[h] = bins_upper_limit((TH1F*) m_data[h][c][y]->Rebin(bin_width), limit);
    if(debug){
      printVector(peak_bins[h], "PEAK:("+h+","+c+","+y+")");
      printVector(empty_bins[h], "EMPTY:("+h+","+c+","+y+")");
    }
  }

  SubtractBackgroundsMap(m_data, m_st, m_wjets, m_other); // Subtract Bkg from data
  if(debug) cout << "\t ... Rebin" << endl;
  MapH m_data_rebin = RebinMap(m_data, bin_width);
  MapH m_ttbar_rebin = RebinMap(m_ttbar, bin_width);
  MapH m_JECup_rebin = RebinMap(m_JECup, bin_width);
  MapH m_JECdown_rebin = RebinMap(m_JECdown, bin_width);
  MapH m_CORup_rebin = RebinMap(m_CORup, bin_width);
  MapH m_CORdown_rebin = RebinMap(m_CORdown, bin_width);
  MapH m_JERup_rebin = RebinMap(m_JERup, bin_width);
  MapH m_JERdown_rebin = RebinMap(m_JERdown, bin_width);
  MapH m_tuneup_rebin = RebinMap(m_tuneup, bin_width);
  MapH m_tunedown_rebin = RebinMap(m_tunedown, bin_width);
  MapH m_hdampup_rebin = RebinMap(m_hdampup, bin_width);
  MapH m_hdampdown_rebin = RebinMap(m_hdampdown, bin_width);
  MapH m_CRqb_rebin = RebinMap(m_CRqb, bin_width);
  MapH m_CRgm_rebin = RebinMap(m_CRgm, bin_width);
  MapH m_fsrup_rebin = RebinMap(m_fsrup, bin_width);
  MapH m_fsrdown_rebin = RebinMap(m_fsrdown, bin_width);
  MapH m_isrup_rebin = RebinMap(m_isrup, bin_width);
  MapH m_isrdown_rebin = RebinMap(m_isrdown, bin_width);

  int mDim = m_data_rebin[w_mass_hh]->GetNbinsX();
  if(mDim!=m_data_rebin[w_mass_hl]->GetNbinsX()||
  mDim!=m_data_rebin[w_mass_hl]->GetNbinsX()||
  mDim!=m_data_rebin[w_mass_hl]->GetNbinsX()) throw runtime_error("Number of bins for jet mass are not equal");

  if(debug) cout << "\t ... Get Cov Matrix" << endl;
  MapM m_cov_ttbar = GetCovMatrixMap(m_ttbar_rebin, save, "ttbar");
  MapM m_cov_JECup = GetCovMatrixMap(m_JECup_rebin, save, "JECup");
  MapM m_cov_JECdown = GetCovMatrixMap(m_JECdown_rebin, save, "JECdown");
  MapM m_cov_CORup = GetCovMatrixMap(m_CORup_rebin, save, "CORup");
  MapM m_cov_CORdown = GetCovMatrixMap(m_CORdown_rebin, save, "CORdown");
  MapM m_cov_JERup = GetCovMatrixMap(m_JERup_rebin, save, "JERup");
  MapM m_cov_JERdown = GetCovMatrixMap(m_JERdown_rebin, save, "JERdown");

  MapH m_norm = NormToMC?m_ttbar_rebin:m_data_rebin;
  MapM m_cov_norm_ttbar = GetCovMatrixNormMap(m_cov_ttbar, m_norm, save, "ttbar");
  MapM m_cov_norm_JECup = GetCovMatrixNormMap(m_cov_JECup, m_norm, save, "JECup");
  MapM m_cov_norm_JECdown = GetCovMatrixNormMap(m_cov_JECdown, m_norm, save, "JECdown");
  MapM m_cov_norm_CORup = GetCovMatrixNormMap(m_cov_CORup, m_norm, save, "CORup");
  MapM m_cov_norm_CORdown = GetCovMatrixNormMap(m_cov_CORdown, m_norm, save, "CORdown");
  MapM m_cov_norm_JERup = GetCovMatrixNormMap(m_cov_JERup, m_norm, save, "JERup");
  MapM m_cov_norm_JERdown = GetCovMatrixNormMap(m_cov_JERdown, m_norm, save, "JERdown");

  MapM m_cov_data = GetCovMatrixMap(m_data_rebin, save, "data");
  MapM m_cov_norm_data; // skip long runtime
  TString add2="";

  if(!addition.EqualTo("") && !addition.EqualTo("lin") && !addition.Contains("syst") && !addition.EqualTo("check")) add2 = "_"+addition;
  if(debug) cout << "add: " << add2 << endl;

  TString sys_file = NormToMC?"_ttbar":"_data";
  TFile *f = new TFile("files/CovMatrices_JMS"+sys_file+".root");
  f->cd();


  if(fast){
    cout << "Matricies for data " << RED << "from ROOT file" << RESET << endl;
    m_cov_norm_data[w_mass_hh].ResizeTo(mDim, mDim); m_cov_norm_data[w_mass_hh] = TMatrixD(*((TMatrixD*) f->Get("covData_hh_norm_"+year+"_"+channel+add2)));
    m_cov_norm_data[w_mass_hl].ResizeTo(mDim, mDim); m_cov_norm_data[w_mass_hl] = TMatrixD(*((TMatrixD*) f->Get("covData_hl_norm_"+year+"_"+channel+add2)));
    m_cov_norm_data[w_mass_lh].ResizeTo(mDim, mDim); m_cov_norm_data[w_mass_lh] = TMatrixD(*((TMatrixD*) f->Get("covData_lh_norm_"+year+"_"+channel+add2)));
    m_cov_norm_data[w_mass_ll].ResizeTo(mDim, mDim); m_cov_norm_data[w_mass_ll] = TMatrixD(*((TMatrixD*) f->Get("covData_ll_norm_"+year+"_"+channel+add2)));
  }
  else{
    m_cov_norm_data = GetCovMatrixNormMap(m_cov_data, m_norm, save, "data");
    m_cov_norm_data[w_mass_hh].Write("covData_hh_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
    m_cov_norm_data[w_mass_hl].Write("covData_hl_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
    m_cov_norm_data[w_mass_lh].Write("covData_lh_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
    m_cov_norm_data[w_mass_ll].Write("covData_ll_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
  }

  f->Close();

  // ==========================================================================================================
  // === Create covariance matrix for systematics

  if(debug) cout << "Construct covariance matricies for systematics ... " << endl;

  // vector<MapM> sys_cor_up, sys_cor_down;
  map<TString, vector<TMatrixD>>  sys_cor_up, sys_cor_down;
  map<TString, vector<TMatrixD>>  sys_cor_up_norm, sys_cor_down_norm;

  TString option = gSystem->AccessPathName("files/CovMatrices_JMS_SYS"+sys_file+".root")?"recreate":"update";
  cout << "opening files/CovMatrices_JMS_SYS"+sys_file+".root" << endl;
  TFile *f_sys = new TFile("files/CovMatrices_JMS_SYS"+sys_file+".root", option);
  f_sys->cd();

  if(addSyst){
    vector<TString> sys = {"tune","hdamp","fsr","isr","cr","jer"};
    for(TString b: hists){
      TString wbin = b;
      wbin.ReplaceAll(cut, "");
      c_sys->Divide(3,2);
      sys_cor_up[wbin] = {};
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_tuneup_rebin[b],  "tune_up", 1, wbin));
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_hdampup_rebin[b], "hdamp_up", 2, wbin));
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_fsrup_rebin[b],   "fsr_up", 3, wbin));
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_isrup_rebin[b],   "isr_up", 3, wbin));
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_CRqb_rebin[b],    "cr", 4, wbin)); // ---
      sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_JERup_rebin[b],   "jer_up", 1, wbin));
      // sys_cor_up[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_CRgm_rebin[b],    "cr", 4, wbin)); // ---
      c_sys->Print(save_nfs+"SYS/SYS_"+wbin+"_all.pdf(");

      sys_cor_down[wbin] = {};
      sys_cor_down[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_tunedown_rebin[b],  "tune_down", 1, wbin));
      sys_cor_down[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_hdampdown_rebin[b], "hdamp_down", 2, wbin));
      sys_cor_down[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_fsrdown_rebin[b],   "fsr_down", 3, wbin));
      sys_cor_down[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_isrdown_rebin[b],   "isr_down", 4, wbin));
      sys_cor_down[wbin].push_back(GetSYSCov(m_ttbar_rebin[b], m_JERdown_rebin[b],   "jer_down", 4, wbin));
      c_sys->Print(save_nfs+"SYS/SYS_"+wbin+"_all.pdf)");

      // // vector<TH1F*> h_syst = {jecup, corup, jmsup, hdampup, isrup, tuneup};
      vector<TH1F*> h_syst = {m_tuneup_rebin[b], m_hdampup_rebin[b], m_fsrup_rebin[b], m_isrup_rebin[b], m_CRqb_rebin[b], m_JERup_rebin[b]};
      sys_cor_up_norm[wbin] = {};

      if(fast){
        for(unsigned int i=0; i<sys_cor_up[wbin].size(); i++){
          TString name = "covSYS_"+wbin+"_norm_"+year+"_"+channel+add2+"_";
          // TString name = "covSYS_"+wbin+"_norm_"+year+"_"+channel+add2+"_"+sys[i];
          if(sys[i].Contains("sr")) name += "_"+fsr_sys;
          cout << "Get SYS cov matrix " << name << endl;
          sys_cor_up_norm[wbin].push_back(TMatrixD(*((TMatrixD*) f_sys->Get(name))));
        }
      }
      else{
        for(unsigned int i=0; i<sys_cor_up[wbin].size(); i++){
          auto start = high_resolution_clock::now(); // Calculation time - start
          // sys_cor_up_norm[wbin].push_back(NormCovMatrixAlt(h_syst[i], sys_cor_up[wbin][i], false, false));
          TH1F* h_norm = NormToMC?m_ttbar_rebin[b]:m_data_rebin[b];
          sys_cor_up_norm[wbin].push_back(NormCovMatrixAlt(h_norm, sys_cor_up[wbin][i], false, false));
          auto stop = high_resolution_clock::now();  // Calculation time - stop
          auto duration = duration_cast<seconds>(stop - start);
          TString name = "covSYS_"+wbin+"_norm_"+year+"_"+channel+add2+"_"+sys[i];
          if(sys[i].Contains("sr")) name += "_"+fsr_sys;
          cout << b << " " << sys[i] << " took " << RED << duration.count() << "s" << RESET << " | stored here: " << name << endl;
          sys_cor_up_norm[wbin][i].Write(name, TObject::kOverwrite);
        }
      }
    }
  }

  // ==========================================================================================================
  // === Construct CovMatrix and then normalize.
  // === Should not make a difference to construct covTotal from normalized matrices as before.

  MapM m_covFit_nom;
  for(TString h: hists){
    vector<double> diagValues;
    int mDim = m_data_rebin[h]->GetNbinsX();
    if(debug) cout << "\t ... Dimension of cov matricies is " << mDim << endl;
    TMatrixD covFit = TMatrixD(mDim, mDim);

    TMatrixD covTTbar   = m_cov_ttbar[h];
    TMatrixD covJECup   = m_cov_JECup[h];
    TMatrixD covJECdown = m_cov_JECdown[h];
    TMatrixD covCORup   = m_cov_CORup[h];
    TMatrixD covCORdown = m_cov_CORdown[h];

    if(debug) cout << "\t ... Fill covFit for " << h << endl;
    for(int x=0; x<mDim; x++){
      diagValues.push_back(covTTbar[x][x]);
      diagValues.push_back(covJECup[x][x]);
      diagValues.push_back(covJECdown[x][x]);
      diagValues.push_back(covCORup[x][x]);
      diagValues.push_back(covCORdown[x][x]);

      double max = *max_element(diagValues.begin(), diagValues.end());

      if(debug){ for(auto v: diagValues) cout << v << "\t"; }
      if(debug){ cout<<"Max value: "<<max<<endl; }

      covFit[x][x] = max;
      diagValues = {};
    }
    m_covFit_nom[h].ResizeTo(covFit); m_covFit_nom[h] = covFit;
    TString wbin = h; wbin.ReplaceAll(cut, "");
    if(debug) cout << "\t ... Draw covFit " << h << endl;
    TH2D* temp = TMatrixDtoTH2D(covFit, covFit.GetNcols(), 0, 180, h+channel+year+"Fit");
    DrawCov(temp, save+"Fit_"+wbin+"_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  }

  TMatrixD m_covTotal_nom_hh = TMatrixD(m_cov_data[w_mass_hh]); m_covTotal_nom_hh += m_covFit_nom[w_mass_hh];
  TMatrixD m_covTotal_nom_hl = TMatrixD(m_cov_data[w_mass_hl]); m_covTotal_nom_hl += m_covFit_nom[w_mass_hl];
  TMatrixD m_covTotal_nom_lh = TMatrixD(m_cov_data[w_mass_lh]); m_covTotal_nom_lh += m_covFit_nom[w_mass_lh];
  TMatrixD m_covTotal_nom_ll = TMatrixD(m_cov_data[w_mass_ll]); m_covTotal_nom_ll += m_covFit_nom[w_mass_ll];

  TMatrixD m_covTotal_norm_hh = TMatrixD(180, 180);
  TMatrixD m_covTotal_norm_hl = TMatrixD(180, 180);
  TMatrixD m_covTotal_norm_lh = TMatrixD(180, 180);
  TMatrixD m_covTotal_norm_ll = TMatrixD(180, 180);

  if(addSyst){
    for(auto m: sys_cor_up["hh"]) m_covTotal_nom_hh += m;
    for(auto m: sys_cor_up["hl"]) m_covTotal_nom_hl += m;
    for(auto m: sys_cor_up["lh"]) m_covTotal_nom_lh += m;
    for(auto m: sys_cor_up["ll"]) m_covTotal_nom_ll += m;

    if(fast){
      m_covTotal_norm_hh = TMatrixD(*((TMatrixD*) f_sys->Get("covTotal_hh_norm_"+year+"_"+channel)));
      m_covTotal_norm_hl = TMatrixD(*((TMatrixD*) f_sys->Get("covTotal_hl_norm_"+year+"_"+channel)));
      m_covTotal_norm_lh = TMatrixD(*((TMatrixD*) f_sys->Get("covTotal_lh_norm_"+year+"_"+channel)));
      m_covTotal_norm_ll = TMatrixD(*((TMatrixD*) f_sys->Get("covTotal_ll_norm_"+year+"_"+channel)));
    }
    else{
      auto start = high_resolution_clock::now(); // Calculation time - start
      TMatrixD m_covTotal_norm_hh = NormCovMatrixAlt(m_data_rebin[w_mass_hh], m_covTotal_nom_hh, false, false);
      m_covTotal_norm_hh.Write("covTotal_hh_norm_"+year+"_"+channel, TObject::kOverwrite);
      auto stop = high_resolution_clock::now();  // Calculation time - stop
      auto duration = duration_cast<seconds>(stop - start);
      cout << "Totel hh took " << RED << duration.count() << "s" << RESET << " | stored here: " << "covTotal_hh_norm_"+year+"_"+channel << endl;

      start = high_resolution_clock::now();
      TMatrixD m_covTotal_norm_hl = NormCovMatrixAlt(m_data_rebin[w_mass_hl], m_covTotal_nom_hl, false, false);
      m_covTotal_norm_hh.Write("covTotal_hl_norm_"+year+"_"+channel, TObject::kOverwrite);
      stop = high_resolution_clock::now();  // Calculation time - stop
      duration = duration_cast<seconds>(stop - start);
      cout << "Totel hl took " << RED << duration.count() << "s" << RESET << " | stored here: " << "covTotal_hl_norm_"+year+"_"+channel << endl;

      start = high_resolution_clock::now();
      TMatrixD m_covTotal_norm_lh = NormCovMatrixAlt(m_data_rebin[w_mass_lh], m_covTotal_nom_lh, false, false);
      m_covTotal_norm_hh.Write("covTotal_lh_norm_"+year+"_"+channel, TObject::kOverwrite);
      stop = high_resolution_clock::now();  // Calculation time - stop
      duration = duration_cast<seconds>(stop - start);
      cout << "Totel lh took " << RED << duration.count() << "s" << RESET << " | stored here: " << "covTotal_lh_norm_"+year+"_"+channel << endl;

      start = high_resolution_clock::now();
      TMatrixD m_covTotal_norm_ll = NormCovMatrixAlt(m_data_rebin[w_mass_ll], m_covTotal_nom_ll, false, false);
      m_covTotal_norm_hh.Write("covTotal_ll_norm_"+year+"_"+channel, TObject::kOverwrite);
      stop = high_resolution_clock::now();  // Calculation time - stop
      duration = duration_cast<seconds>(stop - start);
      cout << "Totel ll took " << RED << duration.count() << "s" << RESET << " | stored here: " << "covTotal_ll_norm_"+year+"_"+channel << endl;

    }
  }
  f_sys->Close();

  // ==========================================================================================================
  // === Normalize Map

  if(debug) cout << "\t ... Normalize" << endl;
  MapH m_data_norm = NormalizeMap(m_data_rebin);
  MapH m_ttbar_norm = NormalizeMap(m_ttbar_rebin);
  MapH m_JECup_norm = NormalizeMap(m_JECup_rebin);
  MapH m_JECdown_norm = NormalizeMap(m_JECdown_rebin);
  MapH m_CORup_norm = NormalizeMap(m_CORup_rebin);
  MapH m_CORdown_norm = NormalizeMap(m_CORdown_rebin);
  MapH m_JERup_norm = NormalizeMap(m_JERup_rebin);
  MapH m_JERdown_norm = NormalizeMap(m_JERdown_rebin);

  if(debug) cout << "Create TGraph2DErrors ... " << endl;
  MapHH m_all_norm;
  m_all_norm["data"] = m_data_norm;
  m_all_norm["ttbar"] = m_ttbar_norm;
  m_all_norm["JECup"] = m_JECup_norm;
  m_all_norm["JECdown"] = m_JECdown_norm;
  m_all_norm["CORup"] = m_CORup_norm;
  m_all_norm["CORdown"] = m_CORdown_norm;

  MapHH m_all_norm_jer;
  m_all_norm_jer["data"] = m_data_norm;
  m_all_norm_jer["ttbar"] = m_ttbar_norm;
  m_all_norm_jer["JECup"] = m_JECup_norm;
  m_all_norm_jer["JECdown"] = m_JECdown_norm;
  m_all_norm_jer["CORup"] = m_CORup_norm;
  m_all_norm_jer["CORdown"] = m_CORdown_norm;
  m_all_norm_jer["JERup"] = m_JERup_norm;
  m_all_norm_jer["JERdown"] = m_JERdown_norm;

  // MapSI2DGe Creat2DGraph(MapHHHH map, MapVIII peak, VecTS hists, TString channel)
  MapSI2DGe m_combine = Creat2DGraph(m_all_norm, peak_bins, hists);
  // Creat3DFit(m_all_norm_jer, peak_bins, hists);
  // MapSI2DGe m_combine = Creat3DGraph(m_all_norm_jer, peak_bins, hists);

  // =====================================================================================
  // === Creat chi2 function for pT bins                                               ===
  // =====================================================================================

  if(debug) cout << "Get single Chi2 ... " << endl;

  if(debug) cout << "\t ... Get Fit functions" << endl;
  MapF2 m_fit_combine_hh = GetFits(m_combine, peak_bins, w_mass_hh);
  MapF2 m_fit_combine_hl = GetFits(m_combine, peak_bins, w_mass_hl);
  MapF2 m_fit_combine_lh = GetFits(m_combine, peak_bins, w_mass_lh);
  MapF2 m_fit_combine_ll = GetFits(m_combine, peak_bins, w_mass_ll);
  cout << endl;

  if(debug) cout << "\t ... Construct Chi2" << endl;
  MapFF2 m_fit_combine;
  m_fit_combine[w_mass_hh] = m_fit_combine_hh;
  m_fit_combine[w_mass_hl] = m_fit_combine_hl;
  m_fit_combine[w_mass_lh] = m_fit_combine_lh;
  m_fit_combine[w_mass_ll] = m_fit_combine_ll;

  // TF2* GetChi2(MapF2 fits, MapVIII peak, TString hist)
  TF2* chi2_hh = GetChi2(m_fit_combine_hh, m_all_norm, peak_bins, w_mass_hh);
  TF2* chi2_hl = GetChi2(m_fit_combine_hl, m_all_norm, peak_bins, w_mass_hl);
  TF2* chi2_lh = GetChi2(m_fit_combine_lh, m_all_norm, peak_bins, w_mass_lh);
  TF2* chi2_ll = GetChi2(m_fit_combine_ll, m_all_norm, peak_bins, w_mass_ll);

  if(debug) cout << "\t ... Get Minimum" << endl;
  // vec_minimums = {x, y, z}
  VecD minimum_hh = GetMinimumChi2(chi2_hh, w_mass_hh);
  VecD minimum_hl = GetMinimumChi2(chi2_hl, w_mass_hl);
  VecD minimum_lh = GetMinimumChi2(chi2_lh, w_mass_lh);
  VecD minimum_ll = GetMinimumChi2(chi2_ll, w_mass_ll);

  if(debug) cout << "\t ... Get 1\u03C3 ellipse" << endl;
  VecDD ellipse_hh = GetSigmaEllipse(chi2_hh, minimum_hh, w_mass_hh, 1, 0.1);
  VecDD ellipse_hl = GetSigmaEllipse(chi2_hl, minimum_hl, w_mass_hl, 1, 0.1);
  VecDD ellipse_lh = GetSigmaEllipse(chi2_lh, minimum_lh, w_mass_lh, 1, 0.1);
  VecDD ellipse_ll = GetSigmaEllipse(chi2_ll, minimum_ll, w_mass_ll, 1, 0.1);

  if(debug) cout << "\t ... Draw Chi2" << endl;
  gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
  TColor::InvertPalette();

  VecDD dummy_vecDD = {};
  Draw2DChi2(chi2_hh, ellipse_hh, dummy_vecDD, minimum_hh, w_mass_hh);
  Draw2DChi2(chi2_hl, ellipse_hl, dummy_vecDD, minimum_hl, w_mass_hl);
  Draw2DChi2(chi2_lh, ellipse_lh, dummy_vecDD, minimum_lh, w_mass_lh);
  Draw2DChi2(chi2_ll, ellipse_ll, dummy_vecDD, minimum_ll, w_mass_ll);

  // =====================================================================================
  // === Get JMS factors for uncorrelated bin uncertainty                              ===
  // =====================================================================================

  // precision ---------------------------------------------------------------------
  space = cout_precision+4;
  cout << fixed << setprecision(cout_precision) << endl;
  cout << "Prepare results for JMS ..." << endl;

  TString function_hh = chi2_hh->GetFormula()->GetParName(0);
  TString function_hl = chi2_hl->GetFormula()->GetParName(0);
  TString function_lh = chi2_lh->GetFormula()->GetParName(0);
  TString function_ll = chi2_ll->GetFormula()->GetParName(0);
  TString function_JMS = function_hh+function_hl+function_lh+function_ll;
  TF2* chi2_JMS = new TF2("JMS", function_JMS, -3, 3, -3, 3);

  VecD nominal_JMS = GetMinimumChi2(chi2_JMS, "JMS");
  cout << nominal_JMS.size() << "\t" << nominal_JMS[0] << "\t" << nominal_JMS[1] << "\t" << nominal_JMS[2] << endl;

  VecDD ellipse_JMS = GetSigmaEllipse(chi2_JMS, nominal_JMS, "JMS", 1000, 0.1);
  cout << "Number of points for uncorrelated 1\u03C3 ellipse: " << ellipse_JMS.size() << endl;

  VecDD dummy = {{0,0},{0,0},{0,0},{0,0}};
  VecDD AllExtrema_JMS = ellipse_JMS.size()==0?dummy:GetExtreme(ellipse_JMS, nominal_JMS);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << "Extract factors from 1\u03C3 ellipse ..." << endl;
  TF1* fsplit = GetSplitFunction(ellipse_JMS);

  VecDD upper, lower;
  for(unsigned int i=0;i<ellipse_JMS.size();i++){
    double xp = ellipse_JMS[i][0]; double yp = ellipse_JMS[i][1];
    double limit = fsplit->Eval(xp);
    if(yp<=limit) lower.push_back(ellipse_JMS[i]);
    else upper.push_back(ellipse_JMS[i]);
  }

  sort(upper.begin(), upper.end(), sortcolx); // Precaution
  sort(lower.begin(), lower.end(), sortcolx); // Precaution
  reverse(lower.begin(), lower.end()); // Prepare to attach both vectors

  if(debug) cout << "\t ... ellipse properties" << endl;
  VecD dist_upper, dist_lower;
  for(auto point: lower) dist_lower.push_back(sqrt(pow(nominal_JMS[0]-point[0], 2)+pow(nominal_JMS[1]-point[1], 2)));
  for(auto point: upper) dist_upper.push_back(sqrt(pow(nominal_JMS[0]-point[0], 2)+pow(nominal_JMS[1]-point[1], 2)));

  double lmin = *min_element(dist_lower.begin(),dist_lower.end());
  double lmax = *max_element(dist_lower.begin(),dist_lower.end());
  int ilmin = find(dist_lower.begin(), dist_lower.end(), lmin)-dist_lower.begin();
  int ilmax = find(dist_lower.begin(), dist_lower.end(), lmax)-dist_lower.begin();
  if(debug) cout << lmax << " ("<< ilmax << ")" << "\t" << lmin << " ("<< ilmin << ")" << endl;

  double umin = *min_element(dist_upper.begin(),dist_upper.end());
  double umax = *max_element(dist_upper.begin(),dist_upper.end());
  int iumin = find(dist_upper.begin(), dist_upper.end(), umin)-dist_upper.begin();
  int iumax = find(dist_upper.begin(), dist_upper.end(), umax)-dist_upper.begin();
  if(debug) cout << umax << " ("<< iumax << ")" << "\t" << umin << " ("<< iumin << ")" << endl;

  // =====================================================================================
  // === extract JMS factors - uncorrelated

  if(debug) cout << "\t ... extract JMS factors" << endl;
  if(debug) cout << "\t\t ... from ellipse" << endl;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  VecD uu = {AllExtrema_JMS[0]};
  VecD dd = {AllExtrema_JMS[1]};
  VecD du = {AllExtrema_JMS[2]};
  VecD ud = {AllExtrema_JMS[3]};

  if(debug) cout << "\t\t ... for PoU" << endl;

  TString full_chi2_jec = function_JMS;
  TString full_chi2_nan = function_JMS;
  TString full_chi2_cor = function_JMS;

  full_chi2_jec.ReplaceAll("y", "("+to_string(nominal_JMS[1])+")");
  full_chi2_cor.ReplaceAll("x", to_string(nominal_JMS[0]));
  full_chi2_cor.ReplaceAll("y", "x");

  TF1 *full_chi2_function_jec = new TF1("chi2_function_jec", full_chi2_jec, -3, 3);
  TF1 *full_chi2_function_cor = new TF1("chi2_function_cor", full_chi2_cor, -3, 3);

  // returns: VecD {f_JMS, sigmaup, sigmadown, f_up, f_down, minChi2}
  VecD PoU_JEC = ExtractJMSValues(full_chi2_function_jec, 1);
  VecD PoU_COR = ExtractJMSValues(full_chi2_function_cor, 1);

  if(debug){
    printVector(PoU_JEC, "PoU for JEC");
    printVector(PoU_COR, "PoU for COR");
  }

  // drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange)
  drawChi2Projection(full_chi2_function_jec, "JEC", {nominal_JMS[0]-0.4, nominal_JMS[0]+0.4}, {215, 235});
  drawChi2Projection(full_chi2_function_cor, "COR", {nominal_JMS[1]-1, nominal_JMS[1]+1}, {215, 235});

  // PoUinTXT(double nom, double up, double down, TString variation)
  PoUinTXT(PoU_JEC[0], PoU_JEC[1], PoU_JEC[2], "JEC");
  PoUinTXT(PoU_COR[0], PoU_COR[1], PoU_COR[2], "COR");

  double sigma_jec = (PoU_JEC[1]+PoU_JEC[2])/2;
  double sigma_cor = (PoU_COR[1]+PoU_COR[2])/2;
  double sigma_jms = sqrt(pow(sigma_jec/PoU_JEC[0],2)+pow(sigma_cor/PoU_COR[0],2));

  if(debug){
    cout << "sigma/min (cor) = " << pow(sigma_cor/PoU_COR[0],2) << endl;
    cout << "sigma/min (jec) = " << pow(sigma_jec/PoU_JEC[0],2) << endl;
    cout << "Error f_{JMS}: " << sigma_jms << endl;
  }

  if(debug) cout << "\t ... calculate correlation factor \u03C1" << endl;
  VecD projection_cor_twosig = ExtractJMSValues(full_chi2_function_cor, 2.3);
  double sigma_cor_center = projection_cor_twosig[1];
  double sigma_cor_full = ymax_global; // just to avoid comfusion
  double rho = sqrt(1-sigma_cor_center/sigma_cor_full);

  // void PrintValues(VecD& nominal, VecD& uu, VecD& ud, VecD& du, VecD& dd, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms, TString name)
  PrintValues(nominal_JMS, uu, ud, du, dd, "JMS.txt", sigma_cor_center, rho, sigma_cor, sigma_jec, PoU_COR, PoU_JEC, sigma_jms);

  // VecDD ordered_ellipse = upper;
  // ordered_ellipse.insert(ordered_ellipse.end(), lower.begin(), lower.end());
  // TEllipse* ellipse = ApproxEllipse(nominal_JMS, upper[iumin], lower[ilmin], upper[iumax], lower[ilmax]);
  VecDD extrema = {uu, dd, ud, du};
  Draw2DChi2(chi2_JMS, ellipse_JMS, extrema, nominal_JMS, "JMS");

  // =====================================================================================
  // === Get JMS factors for correlated bin uncertainty                                ===
  // =====================================================================================

  VecDD ellipse_JMS_alt = GetSigmaEllipse_alt(chi2_JMS, nominal_JMS, "JMS", 1000, 0.1);
  cout << "Number of points for uncorrelated 1\u03C3 ellipse: " << ellipse_JMS_alt.size() << endl;

  cout << "Extract factors from 1\u03C3 ellipse ..." << endl;
  VecDD AllExtrema_JMS_alt = GetExtreme(ellipse_JMS_alt, nominal_JMS);
  VecD uu_alt = {AllExtrema_JMS_alt[0]};
  VecD dd_alt = {AllExtrema_JMS_alt[1]};
  VecD du_alt = {AllExtrema_JMS_alt[2]};
  VecD ud_alt = {AllExtrema_JMS_alt[3]};

  double zero = 0; VecD vec0 = {0};
  // void PrintValues(VecD& nominal, VecD& uu_alt, VecD& ud_alt, VecD& du_alt, VecD& dd_alt, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms, TString name)
  PrintValues(nominal_JMS, uu_alt, ud_alt, du_alt, dd_alt, "JMSalternative.txt", zero, zero, zero, zero, vec0, vec0, zero);

  // VecDD ordered_ellipse_alt = upper_alt;
  // ordered_ellipse_alt.insert(ordered_ellipse_alt.end(), lower_alt.begin(), lower_alt.end());
  // TEllipse* ellipse_alt = ApproxEllipse(nominal_JMS, upper_alt[iumin_alt], lower_alt[ilmin_alt], upper_alt[iumax_alt], lower_alt[ilmax_alt]);
  VecDD extrema_alt = {uu_alt, dd_alt, ud_alt, du_alt};
  Draw2DChi2(chi2_JMS, ellipse_JMS_alt, extrema_alt, nominal_JMS, "JMSalternative");

  // =====================================================================================
  // === Get JMS factors for correlated bin uncertainty                                ===
  // =====================================================================================

  if(debug) cout << "\t ... Construct FitFunctions" << endl;
  VecTS fits_hh = GetFitFunctions(m_fit_combine_hh, peak_bins, w_mass_hh);
  VecTS fits_hl = GetFitFunctions(m_fit_combine_hl, peak_bins, w_mass_hl);
  VecTS fits_lh = GetFitFunctions(m_fit_combine_lh, peak_bins, w_mass_lh);
  VecTS fits_ll = GetFitFunctions(m_fit_combine_ll, peak_bins, w_mass_ll);

  if(debug) cout << "\t ... Store Data bin Content" << endl;
  VecD content_data_hh = StoreBinContentInVector(m_data_norm, peak_bins, w_mass_hh);
  VecD content_data_hl = StoreBinContentInVector(m_data_norm, peak_bins, w_mass_hl);
  VecD content_data_lh = StoreBinContentInVector(m_data_norm, peak_bins, w_mass_lh);
  VecD content_data_ll = StoreBinContentInVector(m_data_norm, peak_bins, w_mass_ll);

  vector<bool> skipBins_hh = GetSkipBins(m_data, peak_bins, w_mass_hh);
  vector<bool> skipBins_hl = GetSkipBins(m_data, peak_bins, w_mass_hl);
  vector<bool> skipBins_lh = GetSkipBins(m_data, peak_bins, w_mass_lh);
  vector<bool> skipBins_ll = GetSkipBins(m_data, peak_bins, w_mass_ll);

  // ------------ Create covariance matrix for fit
  MapM m_covFit;
  for(TString h: hists){
    vector<double> diagValues;
    int mDim = m_data_rebin[h]->GetNbinsX();
    if(debug) cout << "\t ... Dimension of cov matricies is " << mDim << endl;
    TMatrixD covFit = TMatrixD(mDim, mDim);

    TMatrixD covTTbar   = m_cov_norm_ttbar[h];
    TMatrixD covJECup   = m_cov_norm_JECup[h];
    TMatrixD covJECdown = m_cov_norm_JECdown[h];
    TMatrixD covCORup   = m_cov_norm_CORup[h];
    TMatrixD covCORdown = m_cov_norm_CORdown[h];

    if(debug) cout << "\t ... Fill covFit for " << h << endl;
    for(int x=0; x<mDim; x++){
      diagValues.push_back(covTTbar[x][x]);
      diagValues.push_back(covJECup[x][x]);
      diagValues.push_back(covJECdown[x][x]);
      diagValues.push_back(covCORup[x][x]);
      diagValues.push_back(covCORdown[x][x]);

      // diagValues.pop_back();
      double max = *max_element(diagValues.begin(), diagValues.end());

      if(debug){ for(auto v: diagValues) cout << v*1e7 << "\t"; }
      if(debug){ cout<<"Max value: "<<max*1e7<<endl; }

      covFit[x][x] = max;
      diagValues = {};
    }
    m_covFit[h].ResizeTo(covFit); m_covFit[h] = covFit;
    TString wbin = h; wbin.ReplaceAll(cut, "");
    if(debug) cout << "\t ... Draw covFit " << h << endl;
    TH2D* temp = TMatrixDtoTH2D(covFit, covFit.GetNcols(), 0, 180, h+channel+year+"Fit");
    DrawCov(temp, save+"Fit_"+wbin+"_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  }

  TMatrixD m_covTotal_hh = TMatrixD(m_cov_norm_data[w_mass_hh]); m_covTotal_hh += m_covFit[w_mass_hh];
  TMatrixD m_covTotal_hl = TMatrixD(m_cov_norm_data[w_mass_hl]); m_covTotal_hl += m_covFit[w_mass_hl];
  TMatrixD m_covTotal_lh = TMatrixD(m_cov_norm_data[w_mass_lh]); m_covTotal_lh += m_covFit[w_mass_lh];
  TMatrixD m_covTotal_ll = TMatrixD(m_cov_norm_data[w_mass_ll]); m_covTotal_ll += m_covFit[w_mass_ll];

  if(addSyst){
    for(auto m: sys_cor_up_norm["hh"]) m_covTotal_hh += m;
    for(auto m: sys_cor_up_norm["hl"]) m_covTotal_hl += m;
    for(auto m: sys_cor_up_norm["lh"]) m_covTotal_lh += m;
    for(auto m: sys_cor_up_norm["ll"]) m_covTotal_ll += m;
  }

  if(debug){
    cout << "data hh - cols: " << m_cov_norm_data[w_mass_hh].GetNcols() << "\t rows: " << m_cov_norm_data[w_mass_hh].GetNrows() << "\t 90.90: " << m_cov_norm_data[w_mass_hh][90][90]*1e6 << "\t 90.95: " << m_cov_norm_data[w_mass_hh][90][95]*1e6 << endl;
    cout << "data hl - cols: " << m_cov_norm_data[w_mass_hl].GetNcols() << "\t rows: " << m_cov_norm_data[w_mass_hl].GetNrows() << "\t 90.90: " << m_cov_norm_data[w_mass_hl][90][90]*1e6 << "\t 90.95: " << m_cov_norm_data[w_mass_hl][90][95]*1e6 << endl;
    cout << "data lh - cols: " << m_cov_norm_data[w_mass_lh].GetNcols() << "\t rows: " << m_cov_norm_data[w_mass_lh].GetNrows() << "\t 90.90: " << m_cov_norm_data[w_mass_lh][90][90]*1e6 << "\t 90.95: " << m_cov_norm_data[w_mass_lh][90][95]*1e6 << endl;
    cout << "data ll - cols: " << m_cov_norm_data[w_mass_ll].GetNcols() << "\t rows: " << m_cov_norm_data[w_mass_ll].GetNrows() << "\t 90.90: " << m_cov_norm_data[w_mass_ll][90][90]*1e6 << "\t 90.95: " << m_cov_norm_data[w_mass_ll][90][95]*1e6 << endl;

    cout << "fit hh - cols: " << m_covFit[w_mass_hh].GetNcols() << "\t rows: " << m_covFit[w_mass_hh].GetNrows() << "\t 90.90: " << m_covFit[w_mass_hh][90][90]*1e6 << "\t 90.95: " << m_covFit[w_mass_hh][90][95]*1e6 << endl;
    cout << "fit hl - cols: " << m_covFit[w_mass_hl].GetNcols() << "\t rows: " << m_covFit[w_mass_hl].GetNrows() << "\t 90.90: " << m_covFit[w_mass_hl][90][90]*1e6 << "\t 90.95: " << m_covFit[w_mass_hl][90][95]*1e6 << endl;
    cout << "fit lh - cols: " << m_covFit[w_mass_lh].GetNcols() << "\t rows: " << m_covFit[w_mass_lh].GetNrows() << "\t 90.90: " << m_covFit[w_mass_lh][90][90]*1e6 << "\t 90.95: " << m_covFit[w_mass_lh][90][95]*1e6 << endl;
    cout << "fit ll - cols: " << m_covFit[w_mass_ll].GetNcols() << "\t rows: " << m_covFit[w_mass_ll].GetNrows() << "\t 90.90: " << m_covFit[w_mass_ll][90][90]*1e6 << "\t 90.95: " << m_covFit[w_mass_ll][90][95]*1e6 << endl;

    cout << "total hh - cols: " << m_covTotal_hh.GetNcols() << "\t rows: " << m_covTotal_hh.GetNrows() << "\t 90.90: " << m_covTotal_hh[90][90]*1e6 << "\t 90.95: " << m_covTotal_hh[90][95]*1e6 << endl;
    cout << "total hl - cols: " << m_covTotal_hl.GetNcols() << "\t rows: " << m_covTotal_hl.GetNrows() << "\t 90.90: " << m_covTotal_hl[90][90]*1e6 << "\t 90.95: " << m_covTotal_hl[90][95]*1e6 << endl;
    cout << "total lh - cols: " << m_covTotal_lh.GetNcols() << "\t rows: " << m_covTotal_lh.GetNrows() << "\t 90.90: " << m_covTotal_lh[90][90]*1e6 << "\t 90.95: " << m_covTotal_lh[90][95]*1e6 << endl;
    cout << "total ll - cols: " << m_covTotal_ll.GetNcols() << "\t rows: " << m_covTotal_ll.GetNrows() << "\t 90.90: " << m_covTotal_ll[90][90]*1e6 << "\t 90.95: " << m_covTotal_ll[90][95]*1e6 << endl;
  }
  // TString GetChi2Cor(TMatrixD m, vector<TString> vecS, VecD vecD, int nbins, vector<bool> skipBin)
  TString Chi2Cor_hh = GetChi2Cor(m_covTotal_hh, fits_hh, content_data_hh, 180, skipBins_hh, "hh");
  TString Chi2Cor_hl = GetChi2Cor(m_covTotal_hl, fits_hl, content_data_hl, 180, skipBins_hl, "hl");
  TString Chi2Cor_lh = GetChi2Cor(m_covTotal_lh, fits_lh, content_data_lh, 180, skipBins_lh, "lh");
  TString Chi2Cor_ll = GetChi2Cor(m_covTotal_ll, fits_ll, content_data_ll, 180, skipBins_ll, "ll");

  cout << "\nGet final correlated chi2 ..." << endl;
  TString Chi2Cor_final = Chi2Cor_hh+Chi2Cor_hl+Chi2Cor_lh+Chi2Cor_ll;
  cout << "\t ... calculate" << endl;
  TF2* chi2_JMS_cor = new TF2("JMS", Chi2Cor_final, -2, 2, -2, 2);

  cout << "Minimum " << chi2_JMS_cor->GetMinimum() << endl;
  cout << "Get Points from correlated chi2 (Minimum) ..." << endl;
  VecD nominal_JMS_cor = GetMinimumChi2(chi2_JMS_cor, "JMS");
  VecDD ellipse_JMS_cor = GetSigmaEllipse(chi2_JMS_cor, nominal_JMS_cor, "JMS", 1000, 0.1);
  cout << "Number of points for correlated 1\u03C3 ellipse: " << ellipse_JMS_cor.size() << endl;
  VecDD ellipse2_JMS_cor = GetSigmaEllipse(chi2_JMS_cor, nominal_JMS_cor, "JMS", 1000, 0.1, 6.18);
  cout << "Number of points for correlated 2\u03C3 ellipse: " << ellipse_JMS_cor.size() << endl;
  VecDD ellipse_JMS_cor_alt = GetSigmaEllipse_alt(chi2_JMS_cor, nominal_JMS_cor, "JMS", 2000, 0.1);
  cout << "Number of points for correlated 1\u03C3 ellipse_alt: " << ellipse_JMS_cor_alt.size() << endl;
  VecDD ellipse2_JMS_cor_alt = GetSigmaEllipse_alt(chi2_JMS_cor, nominal_JMS_cor, "JMS", 2000, 0.1, 6.18);
  cout << "Number of points for correlated 1\u03C3 ellipse_alt: " << ellipse2_JMS_cor_alt.size() << endl;

  VecDD AllExtrema_JMS_cor = GetExtreme(ellipse_JMS_cor, nominal_JMS_cor);
  VecD uu_cor = {AllExtrema_JMS_cor[0]};
  VecD dd_cor = {AllExtrema_JMS_cor[1]};
  VecD du_cor = {AllExtrema_JMS_cor[2]};
  VecD ud_cor = {AllExtrema_JMS_cor[3]};

  VecDD AllExtrema_JMS_cor_alt = GetExtreme(ellipse_JMS_cor_alt, nominal_JMS_cor);
  VecD uu_cor_alt = {AllExtrema_JMS_cor_alt[0]};
  VecD dd_cor_alt = {AllExtrema_JMS_cor_alt[1]};
  VecD du_cor_alt = {AllExtrema_JMS_cor_alt[2]};
  VecD ud_cor_alt = {AllExtrema_JMS_cor_alt[3]};

  // void PrintValues(VecD& nominal, VecD& uu, VecD& ud, VecD& du, VecD& dd, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms, TString name)
  PrintValues(nominal_JMS_cor, uu_cor, ud_cor, du_cor, dd_cor, "JMScorrelated.txt", zero, zero, zero, zero, vec0, vec0, zero);
  PrintValues(nominal_JMS_cor, uu_cor_alt, ud_cor_alt, du_cor_alt, dd_cor_alt, "JMScorrelated_alt.txt", zero, zero, zero, zero, vec0, vec0, zero);
  // VecDD ordered_ellipse_cor = upper_cor;
  // ordered_ellipse_cor.insert(ordered_ellipse_cor.end(), lower_cor.begin(), lower_cor.end());
  VecDD extrema_cor = {uu_cor, dd_cor, ud_cor, du_cor};
  VecDD extrema_cor_alt = {uu_cor_alt, dd_cor_alt, ud_cor_alt, du_cor_alt};
  Draw2DChi2(chi2_JMS_cor, ellipse_JMS_cor, extrema_cor, nominal_JMS_cor, "JMScor");
  Draw2DChi2(chi2_JMS_cor, ellipse_JMS_cor_alt, extrema_cor_alt, nominal_JMS_cor, "JMScor_alt");

  if(debug) cout << "\t ... Save hist in root file" << endl;
  TPolyMarker3D* zmin_point = new TPolyMarker3D(1);
  zmin_point->SetPoint(0, nominal_JMS_cor[0], nominal_JMS_cor[1], nominal_JMS_cor[2]);

  TPolyMarker3D* sigma_points = new TPolyMarker3D();
  if(ellipse_JMS_cor.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse_JMS_cor.size(); i++) sigma_points->SetPoint(i, ellipse_JMS_cor[i][0], ellipse_JMS_cor[i][1], ellipse_JMS_cor[i][2]);
  }
  TPolyMarker3D* sigma2_points = new TPolyMarker3D();
  if(ellipse2_JMS_cor.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse2_JMS_cor.size(); i++) sigma2_points->SetPoint(i, ellipse2_JMS_cor[i][0], ellipse2_JMS_cor[i][1], ellipse2_JMS_cor[i][2]);
  }

  TPolyMarker3D* sigma_points_alt = new TPolyMarker3D();
  if(ellipse_JMS_cor_alt.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse_JMS_cor_alt.size(); i++) sigma_points_alt->SetPoint(i, ellipse_JMS_cor_alt[i][0], ellipse_JMS_cor_alt[i][1], ellipse_JMS_cor_alt[i][2]);
  }
  TPolyMarker3D* sigma2_points_alt = new TPolyMarker3D();
  if(ellipse2_JMS_cor_alt.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse2_JMS_cor_alt.size(); i++) sigma2_points_alt->SetPoint(i, ellipse2_JMS_cor_alt[i][0], ellipse2_JMS_cor_alt[i][1], ellipse2_JMS_cor_alt[i][2]);
  }

  if(onlyLin && !addSyst){
    TString option = gSystem->AccessPathName("files/PaperPlots.root")?"recreate":"update";
    TFile *f_out = new TFile("files/PaperPlots.root", option);
    f_out->cd();
    if(f_out->GetDirectory("Functions")==0) gDirectory->mkdir("Functions");
    if(f_out->GetDirectory("Graphs")==0)    gDirectory->mkdir("Graphs");
    f_out->cd("Functions");
    chi2_JMS_cor->Write("JMS_Chi2", TObject::kOverwrite);
    f_out->cd("Graphs");
    zmin_point->Write("JMS_nominal", TObject::kOverwrite);
    double x,y,z;
    zmin_point->GetPoint(0, x,y,z);
    cout << x << "\t" << y << "\t" << z << endl;
    sigma_points->Write("JMS_ellipse", TObject::kOverwrite);
    sigma2_points->Write("JMS_ellipse_2sigma", TObject::kOverwrite);
    f_out->Close();
  }

  option = gSystem->AccessPathName(save_nfs+"/Ellipse.root")?"recreate":"update";
  TFile *f_out = new TFile(save_nfs+"/Ellipse.root", option);
  f_out->cd();
  if(f_out->GetDirectory("Functions")==0) gDirectory->mkdir("Functions");
  if(f_out->GetDirectory("Graphs")==0)    gDirectory->mkdir("Graphs");
  f_out->cd("Functions");
  chi2_JMS_cor->Write("JMS_Chi2", TObject::kWriteDelete);
  f_out->cd("Graphs");
  zmin_point->Write("JMS_nominal", TObject::kWriteDelete);
  double x,y,z;
  zmin_point->GetPoint(0, x,y,z);
  cout << x << "\t" << y << "\t" << z << endl;
  sigma_points->Write("JMS_ellipse", TObject::kWriteDelete);
  sigma2_points->Write("JMS_ellipse_2sigma", TObject::kWriteDelete);
  f_out->Close();

  cout << "====================================================" << endl;
  cout << "JMS factors uncorrelated" << endl;
  cout << "fJMS = " << setw(space) << nominal_JMS[0] << " + " << setw(space) << abs(nominal_JMS[0]-uu[0]) << " - " << setw(space) << abs(nominal_JMS[0]-dd[0]) << endl;
  cout << "fCOR = " << setw(space) << nominal_JMS[1] << " + " << setw(space) << abs(nominal_JMS[1]-uu[1]) << " - " << setw(space) << abs(nominal_JMS[1]-dd[1]) << endl;
  cout << "chi2 = " << setw(space) << nominal_JMS[2] << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "JMS factors correlated" << endl;
  cout << "fJMS = " << setw(space) << nominal_JMS_cor[0] << " + " << setw(space) << abs(nominal_JMS_cor[0]-uu_cor[0]) << " - " << setw(space) << abs(nominal_JMS_cor[0]-dd_cor[0]) << endl;
  cout << "fCOR = " << setw(space) << nominal_JMS_cor[1] << " + " << setw(space) << abs(nominal_JMS_cor[1]-uu_cor[1]) << " - " << setw(space) << abs(nominal_JMS_cor[1]-dd_cor[1]) << endl;
  cout << "chi2 = " << setw(space) << nominal_JMS_cor[2] << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "JMS factors correlated alt" << endl;
  cout << "fJMS = " << setw(space) << nominal_JMS_cor[0] << " + " << setw(space) << abs(nominal_JMS_cor[0]-uu_cor_alt[0]) << " - " << setw(space) << abs(nominal_JMS_cor[0]-dd_cor_alt[0]) << endl;
  cout << "fCOR = " << setw(space) << nominal_JMS_cor[1] << " + " << setw(space) << abs(nominal_JMS_cor[1]-uu_cor_alt[1]) << " - " << setw(space) << abs(nominal_JMS_cor[1]-dd_cor_alt[1]) << endl;
  cout << "chi2 = " << setw(space) << nominal_JMS_cor[2] << endl;
  cout << "====================================================" << endl <<endl;
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// Define functions

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
MapHH GetHistograms(TString process, TString h_name){
  VecTS collection_muon, collection_elec;

  if(debug) cout << "\t\t ... " << process << " - " << h_name << endl;

  if(process.EqualTo("data")){ collection_muon = data_muon; collection_elec = data_elec;}
  else if(process.EqualTo("ttbar")){ collection_muon = ttbar_muon; collection_elec = ttbar_elec;}
  else if(process.EqualTo("wjets")){ collection_muon = wjets_muon; collection_elec = wjets_elec;}
  else if(process.EqualTo("st")){ collection_muon = st_muon; collection_elec = st_elec;}
  else if(process.EqualTo("other")){ collection_muon = other_muon; collection_elec = other_elec;}
  else if(process.EqualTo("JECup")){ collection_muon = jec_up_muon; collection_elec = jec_up_elec;}
  else if(process.EqualTo("JECdown")){ collection_muon = jec_down_muon; collection_elec = jec_down_elec;}
  else if(process.EqualTo("CORup")){ collection_muon = cor_up_muon; collection_elec = cor_up_elec;}
  else if(process.EqualTo("CORdown")){ collection_muon = cor_down_muon; collection_elec = cor_down_elec;}
  else if(process.EqualTo("JERup")){ collection_muon = jer_up_muon; collection_elec = jer_up_elec;}
  else if(process.EqualTo("JERdown")){ collection_muon = jer_down_muon; collection_elec = jer_down_elec;}
  else if(process.EqualTo("hdampup")){ collection_muon = hdamp_up_muon; collection_elec = hdamp_up_elec;}
  else if(process.EqualTo("hdampdown")){ collection_muon = hdamp_down_muon; collection_elec = hdamp_down_elec;}
  else if(process.EqualTo("tuneup")){ collection_muon = tune_up_muon; collection_elec = tune_up_elec;}
  else if(process.EqualTo("tunedown")){ collection_muon = tune_down_muon; collection_elec = tune_down_elec;}
  else if(process.EqualTo("CRgm")){ collection_muon = cr_gm_muon; collection_elec = cr_gm_elec;}
  else if(process.EqualTo("CRqb")){ collection_muon = cr_qb_muon; collection_elec = cr_qb_elec;}
  else if(process.EqualTo("fsrup_sqrt2")){ collection_muon = fsr_up_sqrt2_muon; collection_elec = fsr_up_sqrt2_elec;}
  else if(process.EqualTo("fsrdown_sqrt2")){ collection_muon = fsr_down_sqrt2_muon; collection_elec = fsr_down_sqrt2_elec;}
  else if(process.EqualTo("isrup_sqrt2")){ collection_muon = isr_up_sqrt2_muon; collection_elec = isr_up_sqrt2_elec;}
  else if(process.EqualTo("isrdown_sqrt2")){ collection_muon = isr_down_sqrt2_muon; collection_elec = isr_down_sqrt2_elec;}
  else if(process.EqualTo("fsrup_2")){ collection_muon = fsr_up_2_muon; collection_elec = fsr_up_2_elec;}
  else if(process.EqualTo("fsrdown_2")){ collection_muon = fsr_down_2_muon; collection_elec = fsr_down_2_elec;}
  else if(process.EqualTo("isrup_2")){ collection_muon = isr_up_2_muon; collection_elec = isr_up_2_elec;}
  else if(process.EqualTo("isrdown_2")){ collection_muon = isr_down_2_muon; collection_elec = isr_down_2_elec;}
  else if(process.EqualTo("fsrup_4")){ collection_muon = fsr_up_4_muon; collection_elec = fsr_up_4_elec;}
  else if(process.EqualTo("fsrdown_4")){ collection_muon = fsr_down_4_muon; collection_elec = fsr_down_4_elec;}
  else if(process.EqualTo("isrup_4")){ collection_muon = isr_up_4_muon; collection_elec = isr_up_4_elec;}
  else if(process.EqualTo("isrdown_4")){ collection_muon = isr_down_4_muon; collection_elec = isr_down_4_elec;}
  else throw runtime_error("Check the process ´"+process+"´ to obtain the histograms");

  MapHH map;
  vector<TH1F*> h_muon = get_all_hists(collection_muon, h_name);
  vector<TH1F*> h_elec = get_all_hists(collection_elec, h_name);
  vector<TH1F*> h_combine = combine_channels(h_muon, h_elec);

  for(auto nhist:h_muon) nhist->SetTitle(h_name);
  for(auto nhist:h_elec) nhist->SetTitle(h_name);
  for(auto nhist:h_combine) nhist->SetTitle(h_name);

  map["muon"]["2016"] = h_muon[0];
  map["muon"]["2017"] = h_muon[1];
  map["muon"]["2018"] = h_muon[2];
  map["muon"]["combine"] = h_muon[3];

  map["elec"]["2016"] = h_elec[0];
  map["elec"]["2017"] = h_elec[1];
  map["elec"]["2018"] = h_elec[2];
  map["elec"]["combine"] = h_elec[3];

  map["combine"]["2016"] = h_combine[0];
  map["combine"]["2017"] = h_combine[1];
  map["combine"]["2018"] = h_combine[2];
  map["combine"]["combine"] = h_combine[3];

  return map;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
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
    // cout << bin << "\t" << sqrt(syserror2) << "\t" << sqrt(olderror2) << "\t" << sqrt(syserror2+olderror2) << endl;
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

void SubtractBackgroundsMap(MapHHH& m_data, MapHHH& m_st, MapHHH& m_wjets, MapHHH& m_other){
  // MapHHH map;
  for(auto hist: m_data){
    TString h = hist.first;
    for(auto channel: hist.second){
      TString c = channel.first;
      for(auto year: channel.second){
        TString y = year.first;
        vector<TH1F*> bkgs = {m_st[h][c][y], m_wjets[h][c][y], m_other[h][c][y]};
        vector<double> bkgrate;
        if(addition.Contains("noBKGsys")) bkgrate = {0.,0.,0.};
        else bkgrate = {0.23, 0.19, 1.00};
        TH1F* hist = SubtractBackgrounds(m_data[h][c][y], bkgs, bkgrate);
        m_data[h][c][y] = hist;
      }
    }
  }
}

void AddTwoMaps(MapHHH& map1, MapHHH& map2, int option){
  // MapHHH map_out = map1;
  for(auto hist: map1){
    TString h = hist.first;
    for(auto channel: map1[h]){
      TString c = channel.first;
      for(auto year: map1[h][c]){
        TString y = year.first;
        // if(y.Contains("combine")&&c.Contains("combine")) cout << h << " - " << map1[h][c][y]->GetBinContent(90) << " & " << map2[h][c][y]->GetBinContent(90);
        map1[h][c][y]->Add(map2[h][c][y], option);
        // if(y.Contains("combine")&&c.Contains("combine")) cout << " -> " << map1[h][c][y]->GetBinContent(90) << endl;;
      }
    }
  }
  // return map_out;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapM GetCovMatrixMap(MapH& map, const TString& save, const TString& process){
  MapM covs;
  for(auto hist: map){
    // auto start = high_resolution_clock::now(); // Calculation time - start

    TString h = hist.first; TString c = channel; TString y = year;
    TMatrixD mtemp = GetCovMatrix(map[h], debug);

    TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
    TString wbin = h; wbin.ReplaceAll(cut, "");
    DrawCov(htemp, save_nfs+"SYS/"+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    // bool onlyDiag = process.Contains("data")?false:true;
    // TMatrixD norm = NormCovMatrix(map[h], mtemp, false, onlyDiag);
    // TMatrixD* norm_pointer = new TMatrixD(norm);
    // TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
    // DrawCov(htemp_norm, save_nfs+"SYS/"+year+"/"+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
    //
    // auto stop = high_resolution_clock::now();  // Calculation time - stop
    // auto duration = duration_cast<seconds>(stop - start);
    // cout << process << " " << h << " took " << RED << duration.count() << "s" << RESET << endl;

    int dim = map[h]->GetNbinsX();
    covs[h].ResizeTo(dim,dim);
    covs[h] = mtemp;
  }
  return covs;
}

MapM GetCovMatrixNormMap(MapM& map, MapH& hists, const TString& save, const TString& process){
  MapM covs;
  for(auto hist: map){
    auto start = high_resolution_clock::now(); // Calculation time - start

    // TMatrixD mtemp = GetCovMatrix(map[h], debug);
    //
    // TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
    // DrawCov(htemp, save_nfs+"SYS/"+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    TString h = hist.first; TString c = channel; TString y = year;
    TString wbin = h; wbin.ReplaceAll(cut, "");
    TMatrixD mtemp = map[h];
    bool onlyDiag = process.Contains("data")?false:true;
    TMatrixD norm = NormCovMatrixAlt(hists[h], mtemp, false, onlyDiag);
    // TMatrixD* norm_pointer = new TMatrixD(norm);
    TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
    DrawCov(htemp_norm, save_nfs+"SYS/"+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    auto stop = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<seconds>(stop - start);
    if(!onlyDiag) cout << process << " " << h << " took " << RED << duration.count() << "s" << RESET << endl;

    int dim = hists[h]->GetNbinsX();
    covs[h].ResizeTo(dim,dim);
    covs[h] = norm;
  }
  return covs;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapH RebinAndNormalize(MapHHH& process, int& width){
  MapH map;
  for(auto hist: process){
    TString h = hist.first; TString c = channel; TString y = year;
    TH1F* h_rebin = (TH1F*) process[h][c][y]->Rebin(width);
    map[h] = Normalize(h_rebin); // already sets normalized error
  }
  return map;
}

MapH RebinMap(MapHHH& process, int& width){
  MapH map;
  for(auto hist: process){
    TString h = hist.first; TString c = channel; TString y = year;
    map[h] = (TH1F*) process[h][c][y]->Rebin(width);
  }
  return map;
}

MapH NormalizeMap(MapH& process){
  MapH map;
  for(auto hist: process){
    TString h = hist.first;
    map[h] = Normalize(((TH1F*) process[h])); // already sets normalized error
  }
  return map;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void printVector(VecD& vector, TString info) {
  cout << info << "; (Size: " << vector.size() << ") ";
  for(double value: vector) cout << value << ", ";
  cout << endl;
}

void printVector(VecI& vector, TString info) {
  cout << info << "; (Size: " << vector.size() << ") ";
  for(double value: vector) cout << value << ", ";
  cout << endl;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapSI2DGe Creat2DGraph(MapHH& map, MapVI& peak, VecTS& hists){
  if(debug) cout << "\t ... Inside Create2DGraph - " << channel << endl;
  VecD factor_x = {0.0,  1.0, -1.0,  0.0,  0.0};
  VecD factor_y = {0.0,  0.0,  0.0,  1.0, -1.0};
  VecD dummy = {0.0,  0.0,  0.0,  0.0,  0.0};
  VecD content, error;
  MapI2DGe storage;
  MapSI2DGe bins;

  for(TString hist: hists){ // defined in main()
    if(debug) cout << "new hist: " << hist << endl;
    for(unsigned int i=0; i<peak[hist].size(); i++)
    {
      // ----------------------------------------------
      // We only want to combine all years and
      // consider all channels (muon, elec and combine)

      int bin = peak[hist].at(i);
      TString name = to_string(bin);
      content.push_back(map["ttbar"][hist]->GetBinContent(bin));
      content.push_back(map["JECup"][hist]->GetBinContent(bin));
      content.push_back(map["JECdown"][hist]->GetBinContent(bin));
      content.push_back(map["CORup"][hist]->GetBinContent(bin));
      content.push_back(map["CORdown"][hist]->GetBinContent(bin));

      error.push_back(map["ttbar"][hist]->GetBinError(bin));
      error.push_back(map["JECup"][hist]->GetBinError(bin));
      error.push_back(map["JECdown"][hist]->GetBinError(bin));
      error.push_back(map["CORup"][hist]->GetBinError(bin));
      error.push_back(map["CORdown"][hist]->GetBinError(bin));

      TGraph2DErrors* one_bin = new TGraph2DErrors(5, &factor_x[0], &factor_y[0],&content[0],&dummy[0],&dummy[0],&error[0]);
      one_bin->SetName(name); // To avoid Warning: Replacing existing TGraph2D
      storage[bin] = one_bin;

      // void Plot2DGraph(TGraph2DErrors* bin_fit, int bin, TString title, VecD content, VecD error)
      Plot2DGraph(one_bin, bin, hist+channel+year+name, content, error);

      content.clear(); error.clear();
    }
    bins[hist] = storage;
  }
  return bins;
}

void Creat3DFit(MapHH& map, MapVI& peak, VecTS& hists){
  if(debug) cout << "\t ... Inside Create2DGraph - " << channel << endl;
  double factor_x[7] = {0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0};
  double factor_y[7] = {0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0};
  double factor_z[7] = {0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0};
  double dummy[7]    = {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
  double content[7], error[7];
  MapSI2DGe bins;

  for(TString hist: hists){ // defined in main()
    // if(debug) cout << "new hist: " << hist << endl;
    cout << "new hist: " << hist << endl;
    cout << "new hist size: " << peak[hist].size() << endl;
    cout << 1.E-2 << " " << 1.e-2 << endl;
    for(unsigned int i=0; i<peak[hist].size(); i++)
    {
      // ----------------------------------------------
      // We only want to combine all years and
      // consider all channels (muon, elec and combine)
      cout << "\t ... new bin: " << i << endl;

      cout << "\t ... content" << endl;
      int bin = peak[hist].at(i);
      TString name = to_string(bin);
      content[0] = (map["ttbar"][hist]->GetBinContent(bin));
      content[1] = (map["JECup"][hist]->GetBinContent(bin));
      content[2] = (map["JECdown"][hist]->GetBinContent(bin));
      content[3] = (map["CORup"][hist]->GetBinContent(bin));
      content[4] = (map["CORdown"][hist]->GetBinContent(bin));
      content[5] = (map["JERup"][hist]->GetBinContent(bin));
      content[6] = (map["JERdown"][hist]->GetBinContent(bin));

      cout << "\t ... error: "  << endl;
      error[0] = (map["ttbar"][hist]->GetBinError(bin));
      error[1] = (map["JECup"][hist]->GetBinError(bin));
      error[2] = (map["JECdown"][hist]->GetBinError(bin));
      error[3] = (map["CORup"][hist]->GetBinError(bin));
      error[4] = (map["CORdown"][hist]->GetBinError(bin));
      error[5] = (map["JERup"][hist]->GetBinError(bin));
      error[6] = (map["JERdown"][hist]->GetBinError(bin));

      cout << "\t ... BinData" << endl;
      // create a 3d binned data structure
      ROOT::Fit::BinData data(7,3);
      double xx[3];
      for(int i = 0; i < 7; ++i) {
        xx[0] = factor_x[i];
        xx[1] = factor_y[i];
        xx[2] = factor_z[i];
        // add the 3d-data coordinate, the predictor value (v[i])  and its errors
        data.Add(xx, content[i], error[i]);
      }

      cout << "\t ... TF3" << endl;
      TF3 * f3 = new TF3("f3","[0] + [1]*x + [2]*y + [3]*z",0,10,0,10,0,10);
      f3->SetParameters(2,2,2);
      ROOT::Fit::Fitter fitter;
      // wrapped the TF1 in a IParamMultiFunction interface for the Fitter class
      ROOT::Math::WrappedMultiTF1 wf(*f3,3);
      fitter.SetFunction(wf);

      bool ret = fitter.Fit(data);
      if (ret) {
        const ROOT::Fit::FitResult & res = fitter.Result();
        // print result (should be around 1)
        res.Print(std::cout);
        // copy all fit result info (values, chi2, etc..) in TF3
        f3->SetFitResult(res);
        // test fit p-value (chi2 probability)
        double prob = res.Prob();
        if (prob < 1.E-2)
        Error("exampleFit3D","Bad data fit - fit p-value is %f",prob);
        else
        std::cout << "Good fit : p-value  = " << prob << std::endl;

      }
      else
      Error("exampleFit3D","3D fit failed");

    }
    // return bins;
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

void Plot2DGraph(TGraph2DErrors* bin_fit, int bin, TString title, VecD content, VecD error){
  // -----------------------------------------------------------------------
  // Plot Points -----------------------------------------------------------
  double max_bin_content = GetMaxValue(content);
  double min_bin_content = GetMinValue(content);

  // Additional Points to plot the z axis properly -------------------------
  /*
  The additional points are necessary to scale the z axis properly. One can not zoom out of the z axis without
  cutting of the errors above the highest- and below the lowest point. The trick is to add two additional points
  and stretch the z axis. Afterward one can zoom in (and not out).
  */
  bin_fit->SetPoint(5, 0., 0., max_bin_content*1.8);
  bin_fit->SetPointError(5, 0., 0., 0.);
  bin_fit->SetPoint(6, 0., 0., min_bin_content*0.2);
  bin_fit->SetPointError(6, 0., 0., 0.);

  bin_fit->SetTitle("");
  // bin_fit->SetTitle(title);
  bin_fit->GetHistogram()->GetXaxis()->SetTitle("#it{f}^{JEC}");
  bin_fit->GetHistogram()->GetYaxis()->SetTitle("#it{f}^{XCone}");
  bin_fit->GetHistogram()->GetZaxis()->SetTitle("a.u.");
  bin_fit->GetHistogram()->GetXaxis()->SetTitleOffset(1.5);
  bin_fit->GetHistogram()->GetYaxis()->SetTitleOffset(2.0);
  bin_fit->GetHistogram()->GetZaxis()->SetTitleOffset(2.0);
  bin_fit->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
  bin_fit->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
  bin_fit->GetHistogram()->GetZaxis()->SetTitleSize(0.05);
  bin_fit->SetMarkerStyle(20);

  // void SetPolyLine(TPolyLine3D* line, double x1, double x2, double y1, double y2, double z1, double z2, int color)
  TPolyLine3D *xup = SetPolyLine(0, 0, 1, 1, min_bin_content*0.8, content[3]-error[3], kRed);
  TPolyLine3D *xdo = SetPolyLine(0, 0, -1, -1, min_bin_content*0.8, content[4]-error[4], kRed);
  TPolyLine3D *jup = SetPolyLine(1, 1, 0, 0, min_bin_content*0.8, content[1]-error[1], kRed);
  TPolyLine3D *jdo = SetPolyLine(-1, -1, 0, 0, min_bin_content*0.8, content[2]-error[2], kRed);
  TPolyLine3D *mid = SetPolyLine(0, 0, 0, 0, min_bin_content*0.8, content[0]-error[0], kRed);
  TPolyLine3D *xax = SetPolyLine(-1.05, 1.05, 0, 0, min_bin_content*0.8, min_bin_content*0.8, kBlack);
  TPolyLine3D *yax = SetPolyLine(0, 0, -1.05, 1.05, min_bin_content*0.8, min_bin_content*0.8, kBlack);

  TCanvas *B = new TCanvas(title, title, 800, 800); // name used to avoid Warning
  bin_fit->GetHistogram()->GetXaxis()->SetRangeUser(-2, 2);
  bin_fit->GetHistogram()->GetYaxis()->SetRangeUser(-2, 2);
  bin_fit->GetHistogram()->GetZaxis()->SetRangeUser(min_bin_content*0.8, max_bin_content*1.2);

  B->SetRightMargin(0.09);
  B->SetLeftMargin(0.20);
  bin_fit->Draw("err P");
  jup->Draw("same");
  jdo->Draw("same");
  mid->Draw("same");
  xup->Draw("same");
  xdo->Draw("same");
  yax->Draw("same");
  xax->Draw("same");

  TString sbin = to_string(bin);
  TString nbin;
  if(title.Contains("ll")) nbin = "Bin_"+sbin+"_ll.pdf";
  if(title.Contains("hl")) nbin = "Bin_"+sbin+"_hl.pdf";
  if(title.Contains("lh")) nbin = "Bin_"+sbin+"_lh.pdf";
  if(title.Contains("hh")) nbin = "Bin_"+sbin+"_hh.pdf";
  B->SaveAs(save+"/all_bins/"+nbin); // Double / (//) does not affect the dir
  B->Clear();

  // Remove Additional Points for fit --------------------------------------
  bin_fit->RemovePoint(6); // order is on purpose
  bin_fit->RemovePoint(5);
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TPolyLine3D* SetPolyLine(double x1, double x2, double y1, double y2, double z1, double z2, int color){
  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, x1, y1, z1);
  line->SetPoint(1, x2, y2, z2);
  line->SetLineStyle(2);
  line->SetLineColor(color);
  line->SetLineWidth(1);
  return line;
}


// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapF2 GetFits(MapSI2DGe& map, MapVI& peak, TString& hist){
  if(debug) cout << "\n\t ... Inside GetFits - " << channel << endl;

  MapF2 functions;
  vector<VecTS> name_fits = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
  VecII ndfs  = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
  VecII n_fit_para  = {{3}, {3, 3}, {4, 4}, {3}, {5}};
  int nlin = 0; int ntot = 0;
  TString replace = hist; replace.ReplaceAll("comparison_topjet_xcone_pass_rec/", "");

  cout << endl;
  cout << "Start fitting " << replace << " ... " << endl;
  for(unsigned int i=0; i<peak[hist].size(); i++)
  {
    ntot++;
    int bin = peak[hist].at(i);
    TString name = to_string(bin);
    TGraph2DErrors* bin_fit=map[hist][bin];

    TF2 *used_fit;
    TF2 *fit_lin  = new TF2(name, "[0] + [1]*x + [2]*y");
    TF2 *fit_xy2  = new TF2(name, "[0] + [1]*x + [2]*y*y");
    TF2 *fit_x2y  = new TF2(name, "[0] + [1]*x*x + [2]*y");
    TF2 *fit_xyy2 = new TF2(name, "[0] + [1]*x + [2]*y + [3]*y*y");
    TF2 *fit_xx2y = new TF2(name, "[0] + [1]*x + [2]*x*x + [3]*y");
    TF2 *fit_quad = new TF2(name, "[0] + [1]*x*x + [2]*y*y");
    TF2 *fit_poly = new TF2(name, "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y");
    vector<vector<TF2*>> fits = {{fit_lin}, {fit_xy2, fit_x2y}, {fit_xyy2, fit_xx2y}, {fit_quad}, {fit_poly}};

    for(unsigned int order=0; order<fits.size(); order++){
      if(onlyLin && order>0) continue;
      bool twoFits = fits[order].size()==2;
      double chi2a = -1; double chi2b = -1;
      double prob1 = -1; double prob2 = -1;

      bin_fit->Fit(fits[order][0], "Q");
      if(twoFits) bin_fit->Fit(fits[order][1], "Q");

      chi2a = fits[order][0]->GetChisquare();
      if(twoFits) chi2b = fits[order][1]->GetChisquare();

      prob1 = TMath::Prob(chi2a, ndfs[order][0]);
      if(twoFits) prob2 = TMath::Prob(chi2b, ndfs[order][1]);

      bool isProb2 = prob1<prob2;
      double prob = (isProb2)?prob2:prob1;
      double chi2 = (isProb2)?chi2b:chi2a;
      int index = (isProb2)?1:0;

      // -----------------------------------------------------------------------
      if(!onlyLin){
        if(prob>0.05 || ndfs[order][0]==0){ // prob > 0.05 is common in data science
          if(debug) cout << " -> Bin " << name << ": A " << GREEN << name_fits[order][index] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
          functions[name]=fits[order][index];
          if(order == 0) nlin++;
          // Information which can be extracted later.
          // This way, no further information then the function itself has to be provided
          // Necessary for Chi2 construction to decide which formula to use
          const char* info = name_fits[order][index];
          functions[name]->SetParName(0, info);
          break;
        }
      }
      else{
        functions[name]=fits[order][index];
        nlin++;
        const char* info = name_fits[order][index];
        functions[name]->SetParName(0, info);
        break; // extra caution
      }
    }
  }
  if(debug) cout << "\t ... " << nlin << " linear fits with " << ntot << " fits in total (" << dtos(100*nlin/(double)ntot, 2) << " %)" << endl;
  cout << "\t ... " << nlin << " linear fits with " << ntot << " fits in total (" << dtos(100*nlin/(double)ntot, 2) << " %)" << endl;
  return functions;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

VecTS GetFitFunctions(MapF2& fits, MapVI& peak, TString& hist){
  VecTS fit_functions;
  for(unsigned int i=0; i<peak[hist].size(); i++)
  {
    int bin = peak[hist].at(i);
    TString sbin = to_string(bin);
    TFormula *formula = fits[sbin]->GetFormula();
    int npar = formula->GetNpar();
    TString fname = formula->GetParName(0);

    // if(debug){
    //   cout << npar << "\t";
    //   for(unsigned int i=0; i<npar; i++)
    //   {
    //     cout << fits[sbin]->GetParameter(i) << "\t";
    //   }
    //   cout << endl;
    // }

    TString function;
    if(fname.EqualTo("linear")) function = "[0] + [1]*x + [2]*y";
    if(fname.EqualTo("mixed (xy2)")) function = "[0] + [1]*x + [2]*y*y";
    if(fname.EqualTo("mixed (x2y)")) function = "[0] + [1]*x*x + [2]*y";
    if(fname.EqualTo("mixed (xyy2)")) function = "[0] + [1]*x + [2]*y + [3]*y*y";
    if(fname.EqualTo("mixed (xx2y)")) function = "[0] + [1]*x + [2]*x*x + [3]*y";
    if(fname.EqualTo("quadratic")) function = "[0] + [1]*x*x + [2]*y*y";
    if(fname.EqualTo("polynomial of order 2")) function = "[0] + [1]*x + [2]*y + [3]*x*x + [4]*y*y";

    // if(debug) cout << function << endl;
    for(unsigned int i=0; i<npar; i++) function.ReplaceAll("["+to_string(i)+"]", dtos(fits[sbin]->GetParameter(i), precision));
    // if(debug) cout << function << endl;

    fit_functions.push_back(function);
  }
  return fit_functions;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

VecD StoreBinContentInVector(MapH &map, MapVI& peak, TString& hist){
  VecD content;
  for(unsigned int i=0; i<peak[hist].size(); i++){
    content.push_back(map[hist]->GetBinContent(peak[hist][i]));
  }
  return content;
}

vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist){
  vector<bool> skipBin;
  VecI v = peak[hist];
  for(unsigned int i=1; i<=map[hist][channel][year]->GetNbinsX(); i++){
    if(find(v.begin(), v.end(), i) != v.end()) skipBin.push_back(false);
    else skipBin.push_back(true);
  }
  return skipBin;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
TF2* GetChi2(MapF2& fits, MapHH& nhists, MapVI& peak, TString hist){
  int ipar = 0;
  TString fchi2 = "";
  for(unsigned int i=0; i<peak[hist].size(); i++)
  {
    int bin = peak[hist].at(i);
    TString sbin = to_string(bin);
    TFormula *formula = fits[sbin]->GetFormula();
    int npar = formula->GetNpar();
    TString fname = formula->GetParName(0);

    double data_content = nhists["data"][hist]->GetBinContent(bin);
    double data_error = nhists["data"][hist]->GetBinError(bin); // stat. error
    double tt_content = nhists["ttbar"][hist]->GetBinContent(bin); // fit error (estimated)
    double tt_error = nhists["ttbar"][hist]->GetBinError(bin); // fit error (estimated)
    double error2 = pow(data_error, 2) + pow(tt_error, 2);

    TString function = "("+dtos(data_content,50);
    if(fname.EqualTo("linear")) function += " - [0] - [1]*x - [2]*y";
    if(fname.EqualTo("mixed (xy2)")) function += " - [0] - [1]*x - [2]*y*y";
    if(fname.EqualTo("mixed (x2y)")) function += " - [0] - [1]*x*x - [2]*y";
    if(fname.EqualTo("mixed (xyy2)")) function += " - [0] - [1]*x - [2]*y - [3]*y*y";
    if(fname.EqualTo("mixed (xx2y)")) function += " - [0] - [1]*x - [2]*x*x - [3]*y";
    if(fname.EqualTo("quadratic")) function += " - [0] - [1]*x*x - [2]*y*y";
    if(fname.EqualTo("polynomial of order 2")) function += "[0] - [1]*x - [2]*y - [3]*x*x - [4]*y*y";
    function += ")^2/"+dtos(error2,50);
    for(int par=0; par<npar; par++)
    {
      TString before = "["+to_string(par)+"]";
      TString after = "("+dtos(formula->GetParameter(par),50)+")";
      function.ReplaceAll(before, after);
    }
    fchi2+=" + "+function;
  }
  // dummy parameter - SetParName only available if function has parameters!
  fchi2 += " + [0]"; // will be set to 0!
  TF2* chi2 = new TF2(hist, fchi2, -10, 10, -10, 10);

  chi2->SetParameter(0, 0);
  const char *info = fchi2.ReplaceAll(" + [0]", "");
  chi2->GetFormula()->SetParName(0, info); // to extract the formula as TString

  if(debug) cout << "Chi2(0,0) for " << hist << " is " << chi2->Eval(0, 0) << endl;
  return chi2;
}



TString GetChi2Cor(TMatrixD &m, vector<TString> &vecS, VecD &vecD, int nbins, vector<bool> &skipBin, TString wbin){
  TString term = "(di - (gi)) * (dj - (gj)) * Vij-1";
  TString function = "";
  TString temp_term = term;

  vector<bool> dummy;
  for(unsigned int i=0; i<nbins; i++){dummy.push_back(false);}
  TMatrixD m_trim = TrimMatrix(m, dummy); // Cut of diagonals with 0 for .Invert()
  vector<bool> skipBinTrim = TrimSkipBin(m, skipBin, true); // Keep same size as m_trim
  if(skipBinTrim.size()!=m_trim.GetNcols()) throw runtime_error("Trimmed Vector not equal to trimmed Matrix!");

  if(debug) printCovMatrix(m, "standard", 100000, 3);
  if(debug) printCovMatrix(m_trim, "trimmed", 100000, 3);

  TMatrixF minvert = TMatrixF(m_trim);
  minvert.SetTol(1.e-30);
  minvert.Invert();

  TMatrixF mId = TMatrixF(m_trim);
  mId *= minvert;

  if(debug) printCovMatrix(m_trim, "Nominal Matrix", 1.e8, 3);
  if(debug) printCovMatrix(minvert, "Inverted Matrix", 100000);
  if(debug) printCovMatrix(mId, "Identity Test", 1, 1);

  TMatrixD m_final = TrimMatrix(minvert, skipBinTrim);
  vector<bool> skipBinFinal = TrimSkipBin(minvert, skipBinTrim, false);
  TH2D* h_final = TMatrixDtoTH2D(m_final, m_final.GetNcols(), 0, 180, channel+year+wbin+"norm");
  DrawCov(h_final, save+"CovMatrix_"+wbin+"_"+channel+"_"+year, "#it{m}_{W}", 0.3);

  if(skipBinFinal.size()!=m_final.GetNcols()){
    cout << "Vector " << skipBinFinal.size() << " | Matrix " << m_final.GetNcols() << endl;
    throw runtime_error("Final trimmed Vector not equal to trimmed Matrix!");
  }

  if(debug) cout << term << endl;
  nbins=m_final.GetNcols();
  cout << "Final Check --- Matrix: " << nbins << "\t | skipBin: " << skipBinFinal.size() << "\t | data bin content: " << vecD.size() << "\t | terms: " << vecS.size() << endl;
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

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
VecD GetMinimumChi2(TF2* chi2, TString hist){
  double twoD_minX, twoD_minY;
  double twoD_minZ = chi2->GetMinimumXY(twoD_minX,twoD_minY);
  TString histID = hist; histID.ReplaceAll(cut, "");
  if(debug) cout << histID << " -  JEC: " << twoD_minX << "\t | XCone: " << twoD_minY << "\t (\u03A7^2 min: " << twoD_minZ << ")" << endl;
  return {twoD_minX, twoD_minY, twoD_minZ};
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
VecDD GetSigmaEllipse(TF2* chi2, VecD &minimum, TString hist, double runs, double acc, double at){
  if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";
  double z = minimum[2]; double y = minimum[1]; double x = minimum [0];
  auto start = high_resolution_clock::now(); // Calculation time - start
  VecDD points = FindXY(chi2, z+at, x-0.8, x+0.8, y-2.0, y+2.0, runs, acc);
  auto stop  = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  cout << "Numeric solution for "+hist.ReplaceAll(cut, "")+"s 1\u03C3 estimation took " << GREEN << duration.count() << "ms" << RESET << endl;
  if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;
  return points;
}

VecDD GetSigmaEllipse_alt(TF2* chi2, VecD &minimum, TString hist, double runs, double acc, double at){
  if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";
  double z = minimum[2]; double y = minimum[1]; double x = minimum [0];
  auto start = high_resolution_clock::now(); // Calculation time - start
  VecDD points = FindXY_alt(chi2, z+at, x-1.0, x+1.0, y-1.0, y+1.0, runs, acc);
  auto stop  = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  cout << "Numeric solution for "+hist.ReplaceAll(cut, "")+"s 1\u03C3 estimation took " << GREEN << duration.count() << "ms" << RESET << endl;
  if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;
  return points;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void Draw2DChi2(TF2* chi2, VecDD &points, VecDD &extrema, VecD &minimum, TString hist){

  bool isCor = hist.Contains("cor");

  TPolyMarker3D* zmin_point = new TPolyMarker3D(1);
  zmin_point->SetPoint(0, minimum[0], minimum[1], minimum[2]);
  zmin_point->SetMarkerColor(kBlack);
  zmin_point->SetMarkerStyle(kFullCircle);
  zmin_point->SetMarkerSize(0.4);

  TPolyMarker3D* sigma_points = new TPolyMarker3D();
  if(points.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
  }
  sigma_points->SetMarkerColor(kRed);
  sigma_points->SetMarkerStyle(kFullCircle);
  sigma_points->SetMarkerSize(0.1);

  TPolyMarker3D* extrema_points = new TPolyMarker3D();
  if(extrema.size()>0){
    extrema_points->SetPoint(0, extrema[0][0], extrema[0][1], extrema[0][2]);
    extrema_points->SetPoint(1, extrema[1][0], extrema[1][1], extrema[1][2]);
    extrema_points->SetPoint(2, extrema[2][0], extrema[2][1], extrema[2][2]);
    extrema_points->SetPoint(3, extrema[3][0], extrema[3][1], extrema[3][2]);
  }
  extrema_points->SetMarkerColor(kBlue);
  extrema_points->SetMarkerStyle(kFullCircle);
  extrema_points->SetMarkerSize(0.3);

  Int_t nb = 50;
  chi2->SetTitle("");
  chi2->SetFillStyle(1000);
  chi2->SetLineWidth(1);
  chi2->SetRange(-3, -3, 3, 3);
  chi2->SetContour(nb); // Contours

  const Int_t Number = 5;
  Double_t Red[Number]    = { 0.98, 0.00, 0.00, 0.00, 0.00};
  Double_t Green[Number]  = { 0.98, 0.84, 0.315, 0.00, 0.00};
  Double_t Blue[Number]   = { 0.98, 0.84, 0.56, 0.28, 0.00};
  Double_t Length_uncor[Number] = { 0.00, 0.19, 0.46, 0.73, 1.00};
  Double_t Length_cor[Number] = { 0.00, 0.05, 0.51, 0.86, 1.00};
  Double_t Length[Number] = { 0.00,
    isCor?0.20:0.19,
    isCor?0.56:0.46,
    isCor?0.86:0.73,
    1.00
  };
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  TCanvas *D = new TCanvas(hist,"D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  chi2->Draw("cont4z");
  chi2->GetHistogram()->GetXaxis()->SetTitle("#it{f}^{JEC}");
  chi2->GetHistogram()->GetYaxis()->SetTitle("#it{f}^{XCone}");
  chi2->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
  chi2->GetHistogram()->GetXaxis()->SetTitleOffset(1.0);
  chi2->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);
  chi2->GetHistogram()->GetZaxis()->SetTitleOffset(0.5);
  chi2->GetHistogram()->GetZaxis()->CenterTitle();
  chi2->SetMinimum(95);
  gPad->RedrawAxis();
  D->SetTheta(90);
  D->SetPhi(0);
  sigma_points->Draw("SAME P");
  zmin_point->Draw("SAME P");
  extrema_points->Draw("SAME P");
  CMSLabelOffset(0.11, 0.95, 0.105, -0.012);
  TString id = hist.ReplaceAll("comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_", "");
  D->SaveAs(save+"/chi2_"+id+".pdf");
  D->SaveAs(save_nfs+"/chi2_"+id+".pdf");
  D->Clear();
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
bool sortcolx( const VecD& v1, const VecD& v2 ) {
  return v1[0] < v2[0];
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
TEllipse* ApproxEllipse(const VecD& mid, const VecD& uu, const VecD& dd, const VecD& ud, const VecD& du){
  // Problems while drawing, but good alternative for drawing all points
  double sigma_uu = sqrt(pow(mid[0]-uu[0],2)+pow(mid[1]-uu[1],2));
  double sigma_dd = sqrt(pow(mid[0]-dd[0],2)+pow(mid[1]-dd[1],2));
  double sigma_ud = sqrt(pow(mid[0]-ud[0],2)+pow(mid[1]-ud[1],2));
  double sigma_du = sqrt(pow(mid[0]-du[0],2)+pow(mid[1]-du[1],2));

  double dist_min = (sigma_uu+sigma_dd)*0.5;
  double dist_max = (sigma_ud+sigma_du)*0.5;

  VecD sma = {uu[0]-mid[0], uu[1]-mid[1]}; // semi_major_axis
  double length = sqrt(sma[0]*sma[0]+sma[1]*sma[1]);
  double scalar = abs(sma[0]);
  double theta = acos(scalar/length)*(180/TMath::Pi());

  // First: Halbachse alonge xaxis, Second: Halbachse along yaxis
  TEllipse *ellipse = new TEllipse(mid[0], mid[1], dist_min, dist_max, 0, 360, theta);
  return ellipse;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
VecD ExtractJMSValues(const TF1* chi2function, const double& var){
  double minChi2   = chi2function->GetMinimum();
  double f_JMS     = chi2function->GetX(minChi2, -3, 3);
  double f_up      = chi2function->GetX(minChi2+var, f_JMS, 3);
  double f_down    = chi2function->GetX(minChi2+var, -3, f_JMS);

  double sigmaup   = f_up - f_JMS;
  double sigmadown = f_JMS - f_down;

  VecD values      = {f_JMS, sigmaup, sigmadown, f_up, f_down, minChi2};
  return values;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange){
  TCanvas* ccor = new TCanvas(xaxis, xaxis, 600, 600);
  chi2function->SetTitle("");
  chi2function->GetXaxis()->SetRangeUser(xRange[0], xRange[1]);
  chi2function->GetYaxis()->SetRangeUser(yRange[0], yRange[1]);
  chi2function->GetXaxis()->SetTitle("f^{JMS}_{"+xaxis+"}");
  chi2function->GetYaxis()->SetTitle("#chi^{2}");
  chi2function->GetYaxis()->SetTitleOffset(1.2);
  chi2function->Draw();
  ccor->SaveAs(save+"/chi2_projection_"+xaxis+".pdf");
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PoUinTXT(double& nom, double& up, double& down, TString variation){
  ofstream textfile;
  textfile.open(save+"/factor_projection_"+variation+".txt");
  textfile << right;
  textfile << setw(10) << "nom" << setw(10) << "up" << setw(10) << "down" << "\n";
  textfile << setw(10) <<  nom  << setw(10) <<  up  << setw(10) <<  down;
  textfile << "\n";
  textfile.close();
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
TF1* GetSplitFunction(VecDD& ellipse){
  if(debug) cout << "\t ... extract x&y values of points" << endl;
  VecD xJMS, yJMS;
  for(unsigned int i=0; i<ellipse.size(); i++){
    xJMS.push_back(ellipse[i][0]);
    yJMS.push_back(ellipse[i][1]);
  }
  if(debug) cout << "\t ... size xJMS and yJMS: " << xJMS.size() << "\t" << yJMS.size() << endl;
  if(debug) cout << "\t ... extract min and max values for x&y" << endl;
  // To construct function which splits ellipse in two regions.
  // Line goes through xmin to xmax.
  double xmin = *min_element(xJMS.begin(),xJMS.end());
  double xmax = *max_element(xJMS.begin(),xJMS.end());
  double ixmin = find(xJMS.begin(), xJMS.end(), xmin)-xJMS.begin();
  double ixmax = find(xJMS.begin(), xJMS.end(), xmax)-xJMS.begin();

  double ymin = *min_element(yJMS.begin(),yJMS.end());
  double ymax = *max_element(yJMS.begin(),yJMS.end()); ymax_global = ymax;
  double iymin = find(yJMS.begin(), yJMS.end(), ymin)-yJMS.begin();
  double iymax = find(yJMS.begin(), yJMS.end(), ymax)-yJMS.begin();

  VecD point_ymax = ellipse[iymax];
  VecD point_ymin = ellipse[iymin];
  VecD point_xmax = ellipse[ixmax];
  VecD point_xmin = ellipse[ixmin];

  if(debug){
    printVector(point_ymax, to_string((int) iymax));
    printVector(point_ymin, to_string((int) iymin));
    printVector(point_xmax, to_string((int) ixmax));
    printVector(point_xmin, to_string((int) ixmin));
  }

  double xmin_x = point_xmin[0]; double xmax_x = point_xmax[0];
  double xmin_y = point_xmin[1]; double xmax_y = point_xmax[1];

  if(debug){
    cout << "xmin (" << setw(space) << xmin_x << ", " << setw(space) << xmin_y << ")" << endl;
    cout << "xmax (" << setw(space) << xmax_x << ", " << setw(space) << xmax_y << ")" << endl;
  }

  if(debug) cout << "\t ... split ellipse in two regions" << endl;
  double m = (xmax_y-xmin_y)/(xmax_x-xmin_x);
  double bmax = xmax_y-m*xmax_x; double bmin = xmin_y-m*xmin_x;

  if(debug) cout << "Steigung: " << m << " | bmax: " << bmax << " | bmin: " << bmin << endl;
  TString nsplit = to_string(m)+"*x+"+to_string(bmax);
  TF1* fsplit = new TF1("fsplit", nsplit, 0, 2);
  return fsplit;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PrintValues(VecD& nominal, VecD& uu, VecD& ud, VecD& du, VecD& dd, TString name, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms){
  cout << endl;
  if(debug) cout << "\t ... print all information" << endl;
  cout << "====================================================" << endl;
  cout << "JMS factors ("+channel+", "+year+") --- "+name << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "nom (" << setw(space) << nominal[0] << ", " << setw(space) << nominal[1] << ") with z = " << nominal[2] << endl;
  cout << "uu  (" << setw(space) << uu[0] << ", " << setw(space) << uu[1] << ")" << endl;
  cout << "dd  (" << setw(space) << dd[0] << ", " << setw(space) << dd[1] << ")" << endl;
  cout << "du  (" << setw(space) << du[0] << ", " << setw(space) << du[1] << ")" << endl;
  cout << "ud  (" << setw(space) << ud[0] << ", " << setw(space) << ud[1] << ")" << endl;
  cout << "-----------------------------" << endl;
  cout << "fJMS = " << setw(space) << nominal[0] << " + " << setw(space) << abs(nominal[0]-uu[0]) << " - " << setw(space) << abs(nominal[0]-dd[0]) << endl;
  cout << "fCOR = " << setw(space) << nominal[1] << " + " << setw(space) << abs(nominal[1]-uu[1]) << " - " << setw(space) << abs(nominal[1]-dd[1]) << endl;
  cout << endl;
  if(!name.Contains("cor")){
    cout << "sigma_cor: " << sigma_cor_center << endl;
    cout << "\u03C1: " << -rho << endl;
    cout << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
    cout << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
    cout << "Error f_{JMS}: " << sigma_jms << endl;
  }
  cout << "====================================================" << endl <<endl;

  ofstream textfile;
  textfile.open(save+"/"+name);
  textfile << right;
  textfile << "nom (" << setw(space) << nominal[0] << ", " << setw(space) << nominal[1] << ") with z = " << nominal[2] << endl;
  textfile << "uu  (" << setw(space) << uu[0] << ", " << setw(space) << uu[1] << ")" << endl;
  textfile << "dd  (" << setw(space) << dd[0] << ", " << setw(space) << dd[1] << ")" << endl;
  textfile << "du  (" << setw(space) << du[0] << ", " << setw(space) << du[1] << ")" << endl;
  textfile << "ud  (" << setw(space) << ud[0] << ", " << setw(space) << ud[1] << ")" << endl;
  textfile << "-----------------------------" << endl;
  textfile << "fJMS = " << setw(space) << nominal[0] << " + " << setw(space) << abs(nominal[0]-uu[0]) << " - " << setw(space) << abs(nominal[0]-dd[0]) << endl;
  textfile << "fCOR = " << setw(space) << nominal[1] << " + " << setw(space) << abs(nominal[1]-uu[1]) << " - " << setw(space) << abs(nominal[1]-dd[1]) << endl;
  textfile << endl;
  if(!name.Contains("cor")){
    textfile << "sigma_cor: " << sigma_cor_center << endl;
    textfile << "\u03C1: " << -rho << endl;
    textfile << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
    textfile << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
    textfile << "Error f_{JMS}: " << sigma_jms << endl;
  }
  textfile.close();

  textfile.open(save_nfs+"/"+name);
  textfile << right;
  textfile << "nom (" << setw(space) << nominal[0] << ", " << setw(space) << nominal[1] << ") with z = " << nominal[2] << endl;
  textfile << "uu  (" << setw(space) << uu[0] << ", " << setw(space) << uu[1] << ")" << endl;
  textfile << "dd  (" << setw(space) << dd[0] << ", " << setw(space) << dd[1] << ")" << endl;
  textfile << "du  (" << setw(space) << du[0] << ", " << setw(space) << du[1] << ")" << endl;
  textfile << "ud  (" << setw(space) << ud[0] << ", " << setw(space) << ud[1] << ")" << endl;
  textfile << "-----------------------------" << endl;
  textfile << "fJMS = " << setw(space) << nominal[0] << " + " << setw(space) << abs(nominal[0]-uu[0]) << " - " << setw(space) << abs(nominal[0]-dd[0]) << endl;
  textfile << "fCOR = " << setw(space) << nominal[1] << " + " << setw(space) << abs(nominal[1]-uu[1]) << " - " << setw(space) << abs(nominal[1]-dd[1]) << endl;
  textfile << endl;
  if(!name.Contains("cor")){
    textfile << "sigma_cor: " << sigma_cor_center << endl;
    textfile << "\u03C1: " << -rho << endl;
    textfile << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
    textfile << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
    textfile << "Error f_{JMS}: " << sigma_jms << endl;
  }
  textfile.close();
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

VecDD GetExtreme(VecDD& ellipse_JMS, VecD& nominal_JMS){

  TF1* fsplit = GetSplitFunction(ellipse_JMS);

  VecDD upper, lower;
  for(unsigned int i=0;i<ellipse_JMS.size();i++){
    double xp = ellipse_JMS[i][0]; double yp = ellipse_JMS[i][1];
    double limit = fsplit->Eval(xp);
    if(yp<=limit) lower.push_back(ellipse_JMS[i]);
    else upper.push_back(ellipse_JMS[i]);
  }

  sort(upper.begin(), upper.end(), sortcolx); // Precaution
  sort(lower.begin(), lower.end(), sortcolx); // Precaution
  reverse(lower.begin(), lower.end()); // Prepare to attach both vectors

  if(debug) cout << "\t ... ellipse properties" << endl;
  VecD dist_upper, dist_lower;
  for(auto point: lower) dist_lower.push_back(sqrt(pow(nominal_JMS[0]-point[0], 2)+pow(nominal_JMS[1]-point[1], 2)));
  for(auto point: upper) dist_upper.push_back(sqrt(pow(nominal_JMS[0]-point[0], 2)+pow(nominal_JMS[1]-point[1], 2)));

  double lmin = *min_element(dist_lower.begin(),dist_lower.end());
  double lmax = *max_element(dist_lower.begin(),dist_lower.end());
  int ilmin = find(dist_lower.begin(), dist_lower.end(), lmin)-dist_lower.begin();
  int ilmax = find(dist_lower.begin(), dist_lower.end(), lmax)-dist_lower.begin();
  if(debug) cout << lmax << " ("<< ilmax << ")" << "\t" << lmin << " ("<< ilmin << ")" << endl;

  double umin = *min_element(dist_upper.begin(),dist_upper.end());
  double umax = *max_element(dist_upper.begin(),dist_upper.end());
  int iumin = find(dist_upper.begin(), dist_upper.end(), umin)-dist_upper.begin();
  int iumax = find(dist_upper.begin(), dist_upper.end(), umax)-dist_upper.begin();
  if(debug) cout << umax << " ("<< iumax << ")" << "\t" << umin << " ("<< iumin << ")" << endl;

  // =====================================================================================
  // === extract JMS factors - uncorrelated

  if(debug) cout << "\t ... extract JMS factors" << endl;
  if(debug) cout << "\t\t ... from ellipse" << endl;
  VecD uu = {upper[iumin][0], upper[iumin][1]};
  VecD dd = {lower[ilmin][0], lower[ilmin][1]};
  VecD du = {upper[iumax][0], upper[iumax][1]};
  VecD ud = {lower[ilmax][0], lower[ilmax][1]};

  return {uu,dd,du,ud};

}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

vector<vector<double>> FindXY_alt(TF2 *function, double zfix, double xmin, double xmax, double ymin, double ymax, int steps, double accuracy, bool count){
  vector<double> points;
  vector<vector<double>> all_points;
  double dx = (xmax-xmin)/steps;
  double dy = (ymax-ymin)/steps;
  double x, y, z;
  double xfirst=-1; double yfirst=-1;
  double xlast1=-1; double xlast2=-1;
  bool firstHit = false;
  bool RunY = true;
  bool HitInY = false;
  bool FoundX1 = false;

  auto start = high_resolution_clock::now(); // Calculation time - start
  for(int ypar=steps*0.05; ypar<steps; ypar++){
    // if(HitInY) cout << " ============================================================ " << ypar << endl;
    // if(HitInY)  cout << " ----- Last row ------ " << (int) xlast1 << "\t" << (int) xlast2 << "\t" << (int) xfirst << endl;
    // if(!HitInY&&firstHit) continue;
    // if(!firstHit) cout << "y = " << y << endl;
    HitInY=false;
    FoundX1=false;
    if(!RunY) continue; // Will be true if most upper point is hit
    if(count&&(ypar%1000==0)){
      auto stop = high_resolution_clock::now();  // Calculation time - stop
      auto duration = duration_cast<milliseconds>(stop - start);
      cout << "Passed y = " << ypar << " (" << duration.count()/1000 << "s)" << endl;
    }
    y=ymin+ypar*dy;
    for(int xpar=0; xpar<steps; xpar++){
      int distX1 = abs(xpar-xlast1); int distX2 = abs(xpar-xlast2);
      int lowerX1 = xlast1-xpar; int upperX2 = xpar-xlast1;
      if(ypar==yfirst) continue; // stop this line after first point is found
      if(firstHit && (distX1>5&&distX2>5)){continue;} // Next point only in certain intervall
      // if(firstHit && (xpar<(xlast1-10)||xpar>(xlast2+10))){continue;} // Next point only in certain intervall
      x=xmin+xpar*dx;
      z=function->Eval(x, y);
      if(abs(z-zfix)<accuracy){

        if(!firstHit){ // lowest point of ellipse - in theory
          firstHit = true;
          xfirst = xpar; yfirst = ypar;
          xlast1 = xpar; xlast2 = xpar;
          // cout << "First: " << xpar << "\t" << ypar << endl;
        }
        else{
          if(!FoundX1&&xpar<=xlast1){xlast1=xpar; FoundX1=true; /*cout << "below xlast1: " << xpar << "\t" << ypar << endl;*/} // left of xlast1 always first
          else if(xpar>=xlast2){xlast2=xpar; /*cout << "above xlast2: " << xpar << "\t" << ypar << endl;*/} // right of xlast2 always last
          else if(!FoundX1&&distX1<=distX2){xlast1=xpar; FoundX1=true; /*cout << "close xlast1: " << xpar << "\t" << ypar << endl;*/} // for top of ellipse
          else if(distX1>distX2){xlast2=xpar; /*cout << "close xlast2: " << xpar << "\t" << ypar << endl;*/} // for top of ellipse
          // else{ cout << "Check distances - " << xpar << "\t" << ypar << "\t" << distX1 << "\t" << distX2 << endl;}
          // else throw runtime_error("Check points in FindXY_alt");
        }

        HitInY=true;
        points.push_back(x);
        points.push_back(y);
        points.push_back(z);
        all_points.push_back(points);
        points={};
      }
    }
    // if(HitInY) cout << " ----- This row ----- " << (int) xlast1 << "\t" << (int) xlast2 << "\t" << (int) xfirst << endl;
  }

  return all_points;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

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
  TH2D* h_norm = TMatrixDtoTH2D(matrix, nbins, 0, 180, "cov_"+sysname+"_"+year+to_string(panel)+add);
  DrawCov(h_norm, save_nfs+"/SYS/"+year+"/cov_"+sysname+"_"+year+"_"+add, to_string(panel)+add);

  TCanvas *c = new TCanvas(sysname+add, "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("m_{W}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(-100, 100);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST ][");
  // c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  // if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/SYS/SYS_"+sysname+".pdf");
  c->SaveAs(save_nfs+"/SYS/"+year+"/SYS_"+sysname+"_"+add+".pdf");

  // c_sys->cd(panel);
  sys->Draw("HIST ][");

  delete c;
  return matrix;
}
