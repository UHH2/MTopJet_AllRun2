#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"
#include "../include/tdrstyle_all.h"

#include "TRandom2.h"
#include "TFormula.h"
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
MapF2 GetFits(MapSI2DGe& map, TString& hist);
MapHH GetHistograms(TString process, TString h_name);
MapH RebinAndNormalize(MapHHH& process, int& width);
MapH RebinMap(MapHHH& process, int& width);
MapH NormalizeMap(MapH& process);
MapH TrimMap(MapH& process, MapVI& range, TString name);
MapH TrimMatrixMap(MapM& Matrix);
MapM GetCovMatrixMap(MapH& map, const TString& process);
MapM GetCovMatrixNormMap(MapM& map, MapH& hists, const TString& process);
TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel=0, TString add = "");
MapM GetSYSCovMap(MapH& hists, MapH& vars, TString sysname, TString bin, int panel=0);
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TF1* GetSplitFunction(VecDD& ellipse);
TEllipse* ApproxEllipse(const VecD& mid, const VecD& uu, const VecD& dd, const VecD& ud, const VecD& du);
TPolyLine3D* SetPolyLine(double x1, double x2, double y1, double y2, double z1, double z2, int color);
VecD ExtractJMSValues(const TF1* chi2function, const double& var);
VecD GetMinimumChi2(TF2* chi2, TString hist);
VecDD GetSigmaEllipse(TF2* chi2, VecD& minimum, TString hist, double runs, double acc, double at=2.3);
VecDD GetSigmaEllipse_alt(TF2* chi2, VecD &minimum, TString hist, double runs, double acc, double at=2.3, double range = 1.0);
VecDD GetExtreme(VecDD& ellipseJMS, VecD& nominal_JMS);
vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist);
VecD StoreBinContentInVector(MapH& map, TString& hist);
MapSI2DGe Creat2DGraph(MapHH& map, VecTS& hists);
VecTS GetFitFunctions(MapF2& fits);

TString GetChi2Cor(TMatrixD &m, vector<TString> &vecS, VecD &vecD, int nbins, int nbins_hh, int nbins_hl, int nbins_lh, int nbins_ll, bool uncor);

vector<vector<double>> FindXY_alt(TF2 *function, double zfix, double xmin, double xmax, double ymin, double ymax, int steps = 1000, double accuracy = 0, bool count=false);

TCanvas *c_sys = new TCanvas("sys", "sys", 600, 600); // to Plot all systematics in on pdf

bool debug = false;
TString year, channel, save_nfs;
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

MapVI range = {
  {"comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hh", {70, 105}},
  {"comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hl", {75, 104}},
  {"comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_lh", {62, 98}},
  {"comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_ll", {63, 101}},
};

double ymax_global = 0; // TODO

inline TString wbin(TString h){return h.ReplaceAll(cut, "");}

int main(int argc, char* argv[]){

  printf("= ============================== =\n");
  printf("=   Usage: ./JEC_SYS <folder>    =\n");
  printf("= ============================== =\n\n");

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  // Input ---------------------------------------------------------------------
  year = "combine"; // argv[1]
  channel = "combine"; // argv[2]
  addition = argv[1];

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

  save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/"+year+"/"+channel+"/"+slimit+"/"+addition+"/";
  TString temp_path = save_nfs+"{SYS,all_bins}";
  CreateSavePath((string) temp_path);
  // CreateSavePath((string) temp_path);

  cout << "Start collecting Histograms ..." << endl;
  // ---------------------------------------------------------------------------
  // I load all histograms; This code was created after studies ended. Only one
  // option is available. If necessary include other options as well (no bins etc.)

  TString w_mass_hh = hists[0];
  TString w_mass_hl = hists[1];
  TString w_mass_lh = hists[2];
  TString w_mass_ll = hists[3];

  MapHHHH histograms; // DELETE
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
    throw runtime_error("<E> You are adding systematics to fit but not defined a correct FSR/ISR variation (sqrt2,2,4)");
  }

  if(debug) cout << "\t ... Get empty/peak bins" << endl;
  MapVI empty_bins, peak_bins;
  for(TString h: hists){
    TString c = channel; TString y = year;
    empty_bins[h] = bins_empty(m_data[h][c][y]);
    // peak_bins[h] = bins_upper_limit((TH1F*) m_data[h][c][y]->Rebin(bin_width), limit);
    vector<int> bin_pass_v;
    for(int bin=1; bin <= m_data[h][c][y]->GetNbinsX(); bin++){
      if((m_data[h][c][y]->GetBinContent(bin)>limit)) bin_pass_v.push_back(bin);
    }
    peak_bins[h] = bin_pass_v;
    if(debug){
      printVector(peak_bins[h], "PEAK:("+h+","+c+","+y+")");
      printVector(empty_bins[h], "EMPTY:("+h+","+c+","+y+")");
    }
  }

  // ===========================================================================
  // === Sub-study to check high uncertainty in data
  // === Result: large bin uncertatiny comes from data-sub
  // === Which sample causes this ? - Check all three

  MapH m_data_rebin_v0 = RebinMap(m_data, bin_width);
  MapH m_st_rebin_v0 = RebinMap(m_st, bin_width);
  MapH m_wjets_rebin_v0 = RebinMap(m_wjets, bin_width);
  MapH m_other_rebin_v0 = RebinMap(m_other, bin_width);

  // Study bin in Wmass
  MapM m_cov_data_v0 = GetCovMatrixMap(m_data_rebin_v0, "data_noSub");
  MapM m_cov_st_v0 = GetCovMatrixMap(m_st_rebin_v0, "st");
  MapM m_cov_wjets_v0 = GetCovMatrixMap(m_wjets_rebin_v0, "wjets");
  MapM m_cov_other_v0 = GetCovMatrixMap(m_other_rebin_v0, "other");

  // ===========================================================================

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

  // =================================================================================
  // === Only Peak

  if(debug) cout << "\t ... Trim hists" << endl;
  MapH m_data_peak = TrimMap(m_data_rebin, range, "m_data_trim");
  MapH m_ttbar_peak = TrimMap(m_ttbar_rebin, range, "m_ttbar_trim");
  MapH m_JECup_peak = TrimMap(m_JECup_rebin, range, "m_JECup_trim");
  MapH m_JECdown_peak = TrimMap(m_JECdown_rebin, range, "m_JECdown_trim");
  MapH m_CORup_peak = TrimMap(m_CORup_rebin, range, "m_CORup_trim");
  MapH m_CORdown_peak = TrimMap(m_CORdown_rebin, range, "m_CORdown_trim");
  MapH m_JERup_peak = TrimMap(m_JERup_rebin, range, "m_JERup_trim");
  MapH m_JERdown_peak = TrimMap(m_JERdown_rebin, range, "m_JERdown_trim");
  MapH m_tuneup_peak = TrimMap(m_tuneup_rebin, range, "m_tuneup_trim");
  MapH m_tunedown_peak = TrimMap(m_tunedown_rebin, range, "m_tunedown_trim");
  MapH m_hdampup_peak = TrimMap(m_hdampup_rebin, range, "m_hdampup_trim");
  MapH m_hdampdown_peak = TrimMap(m_hdampdown_rebin, range, "m_hdampdown_trim");
  MapH m_CRqb_peak = TrimMap(m_CRqb_rebin, range, "m_CRqb_trim");
  MapH m_CRgm_peak = TrimMap(m_CRgm_rebin, range, "m_CRgm_trim");
  MapH m_fsrup_peak = TrimMap(m_fsrup_rebin, range, "m_fsrup_trim");
  MapH m_fsrdown_peak = TrimMap(m_fsrdown_rebin, range, "m_fsrdown_trim");
  MapH m_isrup_peak = TrimMap(m_isrup_rebin, range, "m_isrup_trim");
  MapH m_isrdown_peak = TrimMap(m_isrdown_rebin, range, "m_isrdown_trim");

  for(auto h:m_data_peak){
    int nbins = h.second->GetNbinsX();
    int b_low = h.second->GetXaxis()->GetBinLowEdge(1);
    int b_up = h.second->GetXaxis()->GetBinUpEdge(nbins);
    int r_low = range[h.first][0];
    int r_up = range[h.first][1];
    if(b_low!=r_low || b_up!=r_up){
      printf("edge low %4i; edge up %4i; range low %4i; range up %4i\n",b_low,b_up,r_low,r_up);
      throw runtime_error("[ERROR] Bin Edges seem wrong");
    }
  }

  // =================================================================================
  // === Normalize

  MapH m_data_norm = NormalizeMap(m_data_peak);
  MapH m_ttbar_norm = NormalizeMap(m_ttbar_peak);
  MapH m_JECup_norm = NormalizeMap(m_JECup_peak);
  MapH m_JECdown_norm = NormalizeMap(m_JECdown_peak);
  MapH m_CORup_norm = NormalizeMap(m_CORup_peak);
  MapH m_CORdown_norm = NormalizeMap(m_CORdown_peak);
  MapH m_JERup_norm = NormalizeMap(m_JERup_peak);
  MapH m_JERdown_norm = NormalizeMap(m_JERdown_peak);
  MapH m_tuneup_norm = NormalizeMap(m_tuneup_peak);
  MapH m_tunedown_norm = NormalizeMap(m_tunedown_peak);
  MapH m_hdampup_norm = NormalizeMap(m_hdampup_peak);
  MapH m_hdampdown_norm = NormalizeMap(m_hdampdown_peak);
  MapH m_CRqb_norm = NormalizeMap(m_CRqb_peak);
  MapH m_CRgm_norm = NormalizeMap(m_CRgm_peak);
  MapH m_fsrup_norm = NormalizeMap(m_fsrup_peak);
  MapH m_fsrdown_norm = NormalizeMap(m_fsrdown_peak);
  MapH m_isrup_norm = NormalizeMap(m_isrup_peak);
  MapH m_isrdown_norm = NormalizeMap(m_isrdown_peak);

  // =================================================================================
  // === Store bin content of data in vector for chi2
  // === Define here for consitency checks & bin dimension definiton

  if(debug) cout << "\t ... Store Data bin Content" << endl;
  VecD content_data_hh = StoreBinContentInVector(m_data_norm, w_mass_hh);
  VecD content_data_hl = StoreBinContentInVector(m_data_norm, w_mass_hl);
  VecD content_data_lh = StoreBinContentInVector(m_data_norm, w_mass_lh);
  VecD content_data_ll = StoreBinContentInVector(m_data_norm, w_mass_ll);

  int dim_hh = content_data_hh.size();
  int dim_hl = content_data_hl.size();
  int dim_lh = content_data_lh.size();
  int dim_ll = content_data_ll.size();
  int dim_full = dim_hh+dim_hl+dim_lh+dim_ll;
  MapI dims = {{w_mass_hh, dim_hh},{w_mass_hl, dim_hl},{w_mass_lh, dim_lh},{w_mass_ll, dim_ll}};
  printf("dim hh: %4i hl: %4i lh: %4i ll: %4i\n",dim_hh,dim_hl,dim_lh,dim_ll);

  int n_bins_hh = m_data_norm[w_mass_hh]->GetNbinsX();
  int n_bins_hl = m_data_norm[w_mass_hl]->GetNbinsX();
  int n_bins_lh = m_data_norm[w_mass_lh]->GetNbinsX();
  int n_bins_ll = m_data_norm[w_mass_ll]->GetNbinsX();
  int n_bins_full = n_bins_hh+n_bins_hl+n_bins_lh+n_bins_ll;

  if(dim_hh != n_bins_hh){
    printf("#content %3i vs. #bins %3i in hh",dim_hh,m_data_norm[w_mass_hh]->GetNbinsX());
    throw runtime_error("[ERROR] Somethings wrong with dimension in hh");
  }
  if(dim_hl != n_bins_hl){
    printf("#content %3i vs. #bins %3i in hl",dim_hl,m_data_norm[w_mass_hl]->GetNbinsX());
    throw runtime_error("[ERROR] Somethings wrong with dimension in hl");
  }
  if(dim_lh != n_bins_lh){
    printf("#content %3i vs. #bins %3i in lh",dim_lh,m_data_norm[w_mass_lh]->GetNbinsX());
    throw runtime_error("[ERROR] Somethings wrong with dimension in lh");
  }
  if(dim_ll != n_bins_ll){
    printf("#content %3i vs. #bins %3i in ll",dim_ll,m_data_norm[w_mass_ll]->GetNbinsX());
    throw runtime_error("[ERROR] Somethings wrong with dimension in ll");
  }
  if(dim_full != n_bins_full){
    printf("#content %3i vs. #bins %3i in ll",dim_full,n_bins_full);
    throw runtime_error("[ERROR] Somethings wrong with full dimension");
  }

  // =================================================================================
  // === Get Fits

  MapHH m_all_norm;
  m_all_norm["data"] = m_data_norm;
  m_all_norm["ttbar"] = m_ttbar_norm;
  m_all_norm["JECup"] = m_JECup_norm;
  m_all_norm["JECdown"] = m_JECdown_norm;
  m_all_norm["CORup"] = m_CORup_norm;
  m_all_norm["CORdown"] = m_CORdown_norm;

  if(debug) cout << "Get single Chi2 ... " << endl;
  MapSI2DGe m_combine = Creat2DGraph(m_all_norm, hists);

  if(debug) cout << "\t ... Get Fit functions" << endl;
  MapF2 m_fit_combine_hh = GetFits(m_combine, w_mass_hh);
  MapF2 m_fit_combine_hl = GetFits(m_combine, w_mass_hl);
  MapF2 m_fit_combine_lh = GetFits(m_combine, w_mass_lh);
  MapF2 m_fit_combine_ll = GetFits(m_combine, w_mass_ll);

  space = cout_precision+4;
  cout << fixed << setprecision(cout_precision) << endl;
  cout << "Prepare results for JMS ..." << endl;

  if(debug) cout << "\t ... Construct FitFunctions" << endl;
  VecTS fits_hh = GetFitFunctions(m_fit_combine_hh);
  VecTS fits_hl = GetFitFunctions(m_fit_combine_hl);
  VecTS fits_lh = GetFitFunctions(m_fit_combine_lh);
  VecTS fits_ll = GetFitFunctions(m_fit_combine_ll);

  // =================================================================================
  // === Get covarianze matrix

  MapM m_cov_tt = GetCovMatrixMap(m_ttbar_peak, "tt");
  MapM m_cov_data = GetCovMatrixMap(m_data_peak, "data");

  MapM m_cov_tt_norm = GetCovMatrixNormMap(m_cov_tt, m_ttbar_peak, "tt_norm");
  MapM m_cov_data_norm = GetCovMatrixNormMap(m_cov_data, m_data_peak, "data_norm");

  MapH h_sys_tune, h_sys_hdamp, h_sys_CR, h_sys_fsr, h_sys_isr, h_sys_jer;
  MapM m_sys_tune, m_sys_hdamp, m_sys_CR, m_sys_fsr, m_sys_isr, m_sys_jer;
  for(auto h: hists){
    h_sys_tune[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_tune[h]->Add(m_tuneup_norm[h], -1);
    h_sys_hdamp[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_hdamp[h]->Add(m_hdampup_norm[h], -1);
    h_sys_CR[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_CR[h]->Add(m_CRqb_norm[h], -1);
    h_sys_fsr[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_fsr[h]->Add(m_fsrup_norm[h], -1);
    h_sys_isr[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_isr[h]->Add(m_isrup_norm[h], -1);
    h_sys_jer[h] = (TH1F*) m_ttbar_norm[h]->Clone(); h_sys_jer[h]->Add(m_JERup_norm[h], -1);

    // Get systematic cov matrix
    m_sys_tune[h].ResizeTo(dims[h], dims[h]); m_sys_tune[h] = GetSystCovMatrix(h_sys_tune[h], debug);
    m_sys_hdamp[h].ResizeTo(dims[h], dims[h]); m_sys_hdamp[h] = GetSystCovMatrix(h_sys_hdamp[h], debug);
    m_sys_CR[h].ResizeTo(dims[h], dims[h]); m_sys_CR[h] = GetSystCovMatrix(h_sys_CR[h], debug);
    m_sys_fsr[h].ResizeTo(dims[h], dims[h]); m_sys_fsr[h] = GetSystCovMatrix(h_sys_fsr[h], debug);
    m_sys_isr[h].ResizeTo(dims[h], dims[h]); m_sys_isr[h] = GetSystCovMatrix(h_sys_isr[h], debug);
    m_sys_jer[h].ResizeTo(dims[h], dims[h]); m_sys_jer[h] = GetSystCovMatrix(h_sys_jer[h], debug);
  }

  MapM m_covFit;
  for(auto h:hists){
    int dim = dims[h];
    cout << "\t ... Dimension of cov matricies is " << dim << endl;
    m_covFit[h].ResizeTo(dim, dim);
    for(int x=0; x<dim; x++){
      m_covFit[h][x][x] = m_cov_tt_norm[h][x][x];
    }
  }

  TMatrixD m_covFull_data = TMatrixD(dim_full, dim_full);
  TMatrixD m_covFull_fit = TMatrixD(dim_full, dim_full);

  int d_counter = 0;
  for(auto h:hists){
    int dim = dims[h];
    for(int k=0; k<dim; k++){
      for(int l=0; l<dim; l++){
        int r = k + d_counter;
        int c = l + d_counter;
        m_covFull_data[r][c] = m_cov_data_norm[h][k][l];
        m_covFull_fit[r][c] = m_covFit[h][k][l];
      }
    }
    d_counter += dim;
    cout << h << "  " << d_counter << endl;
  }
  if(d_counter != dim_full){
    printf("counter %3i vs. full dim %3i",d_counter,dim_full);
    throw runtime_error("[ERROR] Counter does not match full dim");
  }
  TH2D* temp_data = TMatrixDtoTH2D(m_covFull_data, m_covFull_data.GetNcols(), 0, dim_full, "datafullnormpeak"+channel+year);
  DrawCov(temp_data, save_nfs+"SYS/SYS_full_data_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  if(debug){
    cout << setw(7) << "stat"<< " Check if off diagonal is 0: | n0? " << setw(15) << m_covFull_data[20][20]*1e10;
    cout << " | n0? " << setw(15) << m_covFull_data[20][10]*1e10 << " | 0? " << setw(15) << m_covFull_data[100][20]*1e10;
    cout << " | 0? " << setw(15) << m_covFull_data[50][110]*1e10 << " | n0? " << setw(15) << m_covFull_data[40][40]*1e10;
    cout << " | n0? " << setw(15) << m_covFull_data[110][120]*1e10 << endl;
  }

  // === Correlated syst

  MapH h_full_sys;
  MapM m_covFull_sys;
  // h_sys_tune contains norm ttbar - norm systematic
  map<TString, MapH> h_sys = {{"tune", h_sys_tune}, {"hdamp", h_sys_hdamp}, {"CR", h_sys_CR}, {"isr", h_sys_isr}, {"fsr", h_sys_fsr}, {"jer", h_sys_jer}};
  for(auto s: h_sys){
    TH1F* hist = new TH1F(s.first, s.first, dim_full, 0, dim_full);
    d_counter = 0;
    for(auto h:hists){
      int dim = dims[h];
      for(int j=1; j<=dim; j++){
        // order hh - hl - lh - ll
        int bin = j + d_counter;
        hist->SetBinContent(bin, s.second[h]->GetBinContent(j));
      }
      d_counter += dim;
      cout << h << " " << s.first << " " << d_counter << endl;
    }
    if(d_counter != dim_full){
      printf("counter %3i vs. full dim %3i",d_counter,dim_full);
      throw runtime_error("[ERROR] Counter does not match full dim");
    }
    h_full_sys[s.first] = (TH1F*) hist->Clone();
    m_covFull_sys[s.first].ResizeTo(dim_full, dim_full);
    m_covFull_sys[s.first] = GetSystCovMatrix(h_full_sys[s.first], debug);
    if(debug) cout << s.first << " - Number bins " << hist->GetNbinsX() << " and " << h_full_sys[s.first]->GetNbinsX() << " and " << m_covFull_sys[s.first].GetNcols() << " and " << m_covFull_sys[s.first].GetNrows() << endl;
    TH2D* temp = TMatrixDtoTH2D(m_covFull_sys[s.first], m_covFull_sys[s.first].GetNcols(), 0, dim_full, s.first+"fullnormpeak"+channel+year);
    DrawCov(temp, save_nfs+"SYS/SYS_full_"+s.first+"_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  }

  // // === Unorrelated syst
  //
  // MapM m_covFull_sys_uncor;
  // MapMM m_sys = {{"tune", m_sys_tune}, {"hdamp", m_sys_hdamp}, {"CR", m_sys_CR}, {"isr", m_sys_isr}, {"fsr", m_sys_fsr}, {"jer", m_sys_jer}};
  // for(auto s: h_sys){
  //   int d_counter = 0;
  //   m_covFull_sys_uncor[s.first].ResizeTo(dim_full,dim_full);
  //   for(auto h:hists){
  //     int dim = dims[h];
  //     for(int k=0; k<dim; k++){
  //       for(int l=0; l<dim; l++){
  //         int r = k + d_counter;
  //         int c = l + d_counter;
  //         m_covFull_sys_uncor[s.first][r][c] = m_sys[s.first][h][k][l];
  //       }
  //     }
  //     d_counter += dim;
  //     cout << h << "  " << d_counter << endl;
  //   }
  // }
  // if(d_counter != dim_full){
  //   printf("counter %3i vs. full dim %3i",d_counter,dim_full);
  //   throw runtime_error("[ERROR] Counter does not match full dim");
  // }
  // TH2D* temp_data_uncor = TMatrixDtoTH2D(m_covFull_data, m_covFull_data.GetNcols(), 0, dim_full, "datafullnormpeak"+channel+year);
  // DrawCov(temp_data_uncor, save_nfs+"SYS/SYS_full_data_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  // if(debug){
  //   cout << setw(7) << "stat"<< " Check if off diagonal is 0: | n0? " << setw(15) << m_covFull_data[20][20]*1e10;
  //   cout << " | n0? " << setw(15) << m_covFull_data[20][10]*1e10 << " | 0? " << setw(15) << m_covFull_data[100][20]*1e10;
  //   cout << " | 0? " << setw(15) << m_covFull_data[50][110]*1e10 << " | n0? " << setw(15) << m_covFull_data[40][40]*1e10;
  //   cout << " | n0? " << setw(15) << m_covFull_data[110][120]*1e10 << endl;
  // }


  // === full matrix
  // Consider all systematics correlated in pt bins
  // Statistical and fit uncertainty uncorrelated

  TMatrixD m_covFull = TMatrixD(dim_full, dim_full);
  m_covFull += m_covFull_data + m_covFull_fit;
  for(auto m: m_covFull_sys) m_covFull += m.second;

  vector<double> content_data = content_data_hh; // Possible error if not in same order then matrix
  content_data.insert(content_data.end(), content_data_hl.begin(), content_data_hl.end());
  content_data.insert(content_data.end(), content_data_lh.begin(), content_data_lh.end());
  content_data.insert(content_data.end(), content_data_ll.begin(), content_data_ll.end());
  if(debug) cout << "total size content " << setw(10) << content_data.size() << " with hh " << setw(10) << content_data_hh.size() << " and hl " << setw(10) << content_data_hl.size() << " and lh " << setw(10) << content_data_lh.size() << " and ll " << setw(10) << content_data_ll.size() << endl;

  vector<TString> fits = fits_hh; // Possible error if not in same order then matrix
  fits.insert(fits.end(), fits_hl.begin(), fits_hl.end());
  fits.insert(fits.end(), fits_lh.begin(), fits_lh.end());
  fits.insert(fits.end(), fits_ll.begin(), fits_ll.end());
  if(debug) cout << "total size fits    " << setw(10) << fits.size() << " with hh " << setw(10) << fits_hh.size() << " and hl " << setw(10) << fits_hl.size() << " and lh " << setw(10) << fits_lh.size() << " and ll " << setw(10) << fits_ll.size() << endl;

  int c_hh = content_data_hh.size();
  int c_hl = content_data_hl.size();
  int c_lh = content_data_lh.size();
  int c_ll = content_data_ll.size();

  // ===========================================================================
  // === Get Chi2

  // === All syst Correlated

  TString Chi2Cor_full = GetChi2Cor(m_covFull, fits, content_data, dim_full, c_hh, c_hl, c_lh, c_ll, false);
  TF2* chi2_JMS_cor_full = new TF2("JMSfull", Chi2Cor_full, -2, 2, -2, 2);

  return 0;

  VecD nominal_JMS_cor_full = GetMinimumChi2(chi2_JMS_cor_full, "JMS");
  // VecDD ellipse_JMS_cor_full = GetSigmaEllipse_alt(chi2_JMS_cor_full, nominal_JMS_cor_full, "JMS bin correalted", 1500, 0.1);
  VecDD ellipse_JMS_cor_full = GetSigmaEllipse(chi2_JMS_cor_full, nominal_JMS_cor_full, "JMS bin correalted", 1000, 0.1);
  // VecDD ellipse2_JMS_cor_full = GetSigmaEllipse_alt(chi2_JMS_cor_full, nominal_JMS_cor_full, "JMS bin correalted", 3000, 0.1, 6.18, 1.5);
  VecDD ellipse2_JMS_cor_full = GetSigmaEllipse(chi2_JMS_cor_full, nominal_JMS_cor_full, "JMS bin correalted", 1000, 0.1, 6.18);
  cout << "Number of points for correlated 1\u03C3 ellipse_alt: " << ellipse_JMS_cor_full.size() << endl;
  cout << "Number of points for correlated 2\u03C3 ellipse_alt: " << ellipse2_JMS_cor_full.size() << endl;


  VecDD AllExtrema_JMS_cor_full = GetExtreme(ellipse_JMS_cor_full, nominal_JMS_cor_full);
  VecD uu_cor_full = {AllExtrema_JMS_cor_full[0]};
  VecD dd_cor_full = {AllExtrema_JMS_cor_full[1]};
  VecD du_cor_full = {AllExtrema_JMS_cor_full[2]};
  VecD ud_cor_full = {AllExtrema_JMS_cor_full[3]};

  TString full_chi2_jec = Chi2Cor_full;
  TString full_chi2_nan = Chi2Cor_full;
  TString full_chi2_cor = Chi2Cor_full;

  full_chi2_jec.ReplaceAll("y", "("+to_string(nominal_JMS_cor_full[1])+")");
  full_chi2_cor.ReplaceAll("x", to_string(nominal_JMS_cor_full[0]));
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

  VecD projection_cor_twosig = ExtractJMSValues(full_chi2_function_cor, 2.3);
  double minY = nominal_JMS_cor_full[1];
  double sigma_cor_center = projection_cor_twosig[1]+projection_cor_twosig[2];
  double sigma_cor_full = 2*(ymax_global-minY); // just to avoid comfusion
  double rho = sqrt(1-sigma_cor_center/sigma_cor_full);
  double zero = 0; VecD vec0 = {0};
  cout << setw(15) << minY << setw(15) << sigma_cor_center << setw(15) << sigma_cor_full << setw(15) << rho << endl;
  cout << setw(15) << projection_cor_twosig[1] << setw(15) << projection_cor_twosig[2] << setw(15) << ymax_global << endl;
  // PrintValues(nominal_JMS, uu, ud, du, dd, "JMS.txt", sigma_cor_center, rho, sigma_cor, sigma_jec, PoU_COR, PoU_JEC, sigma_jms);
  PrintValues(nominal_JMS_cor_full, uu_cor_full, ud_cor_full, du_cor_full, dd_cor_full, "JMScorrelated.txt", sigma_cor_center, rho, sigma_cor, sigma_jec, PoU_COR, PoU_JEC, sigma_jms);
  VecDD extrema_cor_full = {uu_cor_full, dd_cor_full, ud_cor_full, du_cor_full};
  Draw2DChi2(chi2_JMS_cor_full, ellipse_JMS_cor_full, extrema_cor_full, nominal_JMS_cor_full, "JMScor");

  if(debug) cout << "\t ... Save hist in root file" << endl;
  TPolyMarker3D* zmin_point_full = new TPolyMarker3D(1);
  zmin_point_full->SetPoint(0, nominal_JMS_cor_full[0], nominal_JMS_cor_full[1], nominal_JMS_cor_full[2]);

  TPolyMarker3D* sigma_points_full = new TPolyMarker3D();
  if(ellipse_JMS_cor_full.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse_JMS_cor_full.size(); i++) sigma_points_full->SetPoint(i, ellipse_JMS_cor_full[i][0], ellipse_JMS_cor_full[i][1], ellipse_JMS_cor_full[i][2]);
  }
  TPolyMarker3D* sigma2_points_full = new TPolyMarker3D();
  if(ellipse2_JMS_cor_full.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<ellipse2_JMS_cor_full.size(); i++) sigma2_points_full->SetPoint(i, ellipse2_JMS_cor_full[i][0], ellipse2_JMS_cor_full[i][1], ellipse2_JMS_cor_full[i][2]);
  }

  // ========================================================================================
  // === Write TObjects (chi2, minimum, ellipse) to root file for dedicated paper plot script

  TString option = gSystem->AccessPathName("files/PaperPlots_Peak.root")?"recreate":"update";
  TFile *f_out = new TFile("files/PaperPlots_Peak.root", option);
  f_out->cd();
  if(f_out->GetDirectory("Functions")==0) gDirectory->mkdir("Functions");
  if(f_out->GetDirectory("Graphs")==0)    gDirectory->mkdir("Graphs");
  f_out->cd("Functions");
  chi2_JMS_cor_full->Write("JMS_Chi2", TObject::kOverwrite);
  f_out->cd("Graphs");
  zmin_point_full->Write("JMS_nominal", TObject::kOverwrite);
  double x,y,z;
  zmin_point_full->GetPoint(0, x,y,z);
  cout << x << "\t" << y << "\t" << z << endl;
  sigma_points_full->Write("JMS_ellipse", TObject::kOverwrite);
  sigma2_points_full->Write("JMS_ellipse_2sigma", TObject::kOverwrite);
  f_out->Close();

  // ====================================================================
  // === Write extrame (minimum, ellipse) to root file ellipse comparison

  option = gSystem->AccessPathName(save_nfs+"/Ellipse.root")?"recreate":"update";
  f_out = new TFile(save_nfs+"/Ellipse.root", option);
  f_out->cd();
  if(f_out->GetDirectory("Functions")==0) gDirectory->mkdir("Functions");
  if(f_out->GetDirectory("Graphs")==0)    gDirectory->mkdir("Graphs");
  f_out->cd("Functions");
  chi2_JMS_cor_full->Write("JMS_Chi2", TObject::kWriteDelete);
  f_out->cd("Graphs");
  zmin_point_full->Write("JMS_nominal", TObject::kWriteDelete);
  // double x,y,z;
  zmin_point_full->GetPoint(0, x,y,z);
  cout << x << "\t" << y << "\t" << z << endl;
  sigma_points_full->Write("JMS_ellipse", TObject::kWriteDelete);
  sigma2_points_full->Write("JMS_ellipse_2sigma", TObject::kWriteDelete);
  f_out->Close();

  cout << "====================================================" << endl;
  cout << "===               'JMS Peak'                    ====" << endl;
  cout << "====================================================" << endl;
  cout << "JMS factors correlated" << endl;
  cout << "fJMS = " << setw(space) << nominal_JMS_cor_full[0] << " + " << setw(space) << abs(nominal_JMS_cor_full[0]-uu_cor_full[0]) << " - " << setw(space) << abs(nominal_JMS_cor_full[0]-dd_cor_full[0]) << endl;
  cout << "fCOR = " << setw(space) << nominal_JMS_cor_full[1] << " + " << setw(space) << abs(nominal_JMS_cor_full[1]-uu_cor_full[1]) << " - " << setw(space) << abs(nominal_JMS_cor_full[1]-dd_cor_full[1]) << endl;
  cout << "chi2 = " << setw(space) << nominal_JMS_cor_full[2] << endl;
  cout << "====================================================" << endl <<endl;

  return 0;


  // =================================================================================
  // === Continue normal

  // =================================================================================
  // === Example Code if histogram gets too big
  // === Run four times over nbins ->
  // === for 180 bins: 180*180*180*180 bins to loop over
  // === run ones and save in extra file

  // MapM m_cov_data = GetCovMatrixMap(m_data_rebin, save, "data");
  // MapM m_cov_norm_data; // skip long runtime
  // m_cov_norm_data[w_mass_hh].ResizeTo(mDim, mDim);
  // m_cov_norm_data[w_mass_hl].ResizeTo(mDim, mDim);
  // m_cov_norm_data[w_mass_lh].ResizeTo(mDim, mDim);
  // m_cov_norm_data[w_mass_ll].ResizeTo(mDim, mDim);
  //
  // option = gSystem->AccessPathName("files/CovMatrices_JMS.root")?"recreate":"update";
  // TFile *f = new TFile("files/CovMatrices_JMS.root", option);
  // f->cd();
  //
  // if(fast){
  //   cout << "Matricies for data " << RED << "from ROOT file" << RESET << endl;
  //   m_cov_norm_data[w_mass_hh] = TMatrixD(*((TMatrixD*) f->Get("covData_hh_norm_"+year+"_"+channel+add2)));
  //   m_cov_norm_data[w_mass_hl] = TMatrixD(*((TMatrixD*) f->Get("covData_hl_norm_"+year+"_"+channel+add2)));
  //   m_cov_norm_data[w_mass_lh] = TMatrixD(*((TMatrixD*) f->Get("covData_lh_norm_"+year+"_"+channel+add2)));
  //   m_cov_norm_data[w_mass_ll] = TMatrixD(*((TMatrixD*) f->Get("covData_ll_norm_"+year+"_"+channel+add2)));
  // }
  // else{
  //   m_cov_norm_data = GetCovMatrixNormMap(m_cov_data, m_data_rebin, save, "data");
  //   m_cov_norm_data[w_mass_hh].Write("covData_hh_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
  //   m_cov_norm_data[w_mass_hl].Write("covData_hl_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
  //   m_cov_norm_data[w_mass_lh].Write("covData_lh_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
  //   m_cov_norm_data[w_mass_ll].Write("covData_ll_norm_"+year+"_"+channel+add2, TObject::kOverwrite);
  // }
  //
  // TH2D* temp = new TH2D();
  // temp = TMatrixDtoTH2D(m_cov_norm_data[w_mass_hh], m_cov_norm_data[w_mass_hh].GetNcols(), 0, 180, "hhdatanorm"+channel+year); DrawCov(temp, save_nfs+"SYS/data_hh_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  // temp = TMatrixDtoTH2D(m_cov_norm_data[w_mass_hl], m_cov_norm_data[w_mass_hl].GetNcols(), 0, 180, "hldatanorm"+channel+year); DrawCov(temp, save_nfs+"SYS/data_hl_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  // temp = TMatrixDtoTH2D(m_cov_norm_data[w_mass_lh], m_cov_norm_data[w_mass_lh].GetNcols(), 0, 180, "lhdatanorm"+channel+year); DrawCov(temp, save_nfs+"SYS/data_lh_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);
  // temp = TMatrixDtoTH2D(m_cov_norm_data[w_mass_ll], m_cov_norm_data[w_mass_ll].GetNcols(), 0, 180, "lldatanorm"+channel+year); DrawCov(temp, save_nfs+"SYS/data_ll_norm_"+channel+"_"+year, "#it{m}_{W}", 0.3);

  // =================================================================================
  // === Get Skip bins - obsolet

  // vector<bool> skipBins_hh = GetSkipBins(m_data, peak_bins, w_mass_hh);
  // vector<bool> skipBins_hl = GetSkipBins(m_data, peak_bins, w_mass_hl);
  // vector<bool> skipBins_lh = GetSkipBins(m_data, peak_bins, w_mass_lh);
  // vector<bool> skipBins_ll = GetSkipBins(m_data, peak_bins, w_mass_ll);

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

MapM GetCovMatrixMap(MapH& map, const TString& process){
  MapM covs;
  for(auto hist: map){

    TString h = hist.first; TString c = channel; TString y = year;
    TMatrixD mtemp = GetCovMatrixDiag(map[h], debug);

    TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
    TString wbin = h; wbin.ReplaceAll(cut, "");
    DrawCov(htemp, save_nfs+"SYS/"+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    int dim = map[h]->GetNbinsX();
    covs[h].ResizeTo(dim,dim);
    covs[h] = mtemp;
  }
  return covs;
}

MapM GetCovMatrixNormMap(MapM& map, MapH& hists, const TString& process){
  MapM covs;
  for(auto hist: map){
    auto start = high_resolution_clock::now(); // Calculation time - start

    TString h = hist.first; TString c = channel; TString y = year;
    TString wbin = h; wbin.ReplaceAll(cut, "");
    TMatrixD mtemp = map[h];
    bool onlyDiag = process.Contains("data")?false:true;
    TMatrixD norm = NormCovMatrixAlt(hists[h], mtemp, false, onlyDiag);
    TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, norm.GetNcols(), h+c+y+process+"norm");
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

MapH TrimMap(MapH& process, MapVI& range, TString name){
  MapH map;
  for(auto h:hists){
    // TH1F* TrimHistogram(TH1F* hist_, double r_low, double r_high, TString name, bool debug){
    map[h] = TrimHistogram(process[h], range[h].at(0), range[h].at(1), h+"_trimmed_"+name, debug);
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

MapSI2DGe Creat2DGraph(MapHH& map, VecTS& hists){
  if(debug) cout << "\t ... Inside Create2DGraph - " << channel << endl;
  VecD factor_x = {0.0,  1.0, -1.0,  0.0,  0.0};
  VecD factor_y = {0.0,  0.0,  0.0,  1.0, -1.0};
  VecD dummy = {0.0,  0.0,  0.0,  0.0,  0.0};
  VecD content, error;
  MapSI2DGe bins;

  for(TString hist: hists){ // defined in main()
    if(debug) cout << "new hist: " << hist << endl;
    MapI2DGe storage;
    cout << "Create " << hist << " with " << map["data"][hist]->GetNbinsX() << " bins" << endl;
    for(unsigned int i=1; i<=map["data"][hist]->GetNbinsX(); i++)
    {
      cout << setw(3) << i;
      // ----------------------------------------------
      // We only want to combine all years and
      // consider all channels (muon, elec and combine)

      int bin = i;
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
    cout << endl;
    bins[hist] = storage;
    cout << "bin size " << bins[hist].size() << endl;
  }
  return bins;
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
  B->SaveAs(save_nfs+"/all_bins/"+nbin); // Double / (//) does not affect the dir
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

MapF2 GetFits(MapSI2DGe& map, TString& hist){
  if(debug) cout << "\n\t ... Inside GetFits - " << channel << endl;

  MapF2 functions;
  TString name_fit = "linear";
  int ndf  = 2; // ndf = n_points - n_parameters = 5 - n_para
  int n_fit_para  = 3;
  TString replace = hist; replace.ReplaceAll("comparison_topjet_xcone_pass_rec/", "");

  cout << "Start fitting " << replace << " ... " << "(size " << map[hist].size() << ")" << endl;
  cout << map[hist].size() << endl;
  for(unsigned int i=0; i<map[hist].size(); i++)
  {
    TString name = to_string(i+1);
    TGraph2DErrors* bin_fit=map[hist][i+1];

    TF2 *fit  = new TF2(name, "[0] + [1]*x + [2]*y");

    bin_fit->Fit(fit, "Q");
    double chi2 = fit->GetChisquare();
    double prob = TMath::Prob(chi2, ndf);

    // -----------------------------------------------------------------------
    functions[name]=fit;
    const char* info = name_fit;
    functions[name]->SetParName(0, info);
  }
  return functions;
}


// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

VecTS GetFitFunctions(MapF2& fits){
  VecTS fit_functions;
  for(unsigned int i=0; i<fits.size(); i++)
  {
    int bin = i+1;
    TString sbin = to_string(bin);
    TFormula *formula = fits[sbin]->GetFormula();
    int npar = formula->GetNpar();
    TString fname = formula->GetParName(0);

    TString function = "[0] + [1]*x + [2]*y";

    // cout << setw(3) << i << function << endl;
    for(unsigned int i=0; i<npar; i++) function.ReplaceAll("["+to_string(i)+"]", dtos(fits[sbin]->GetParameter(i), precision));
    // cout << setw(3) << "" << function << endl;

    fit_functions.push_back(function);
  }
  return fit_functions;
}


// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist){
  vector<bool> skipBin;
  VecI v = peak[hist];
  for(unsigned int i=1; i<=map[hist][channel][year]->GetNbinsX(); i++){
    if(find(v.begin(), v.end(), i) != v.end()) skipBin.push_back(false);
    else skipBin.push_back(true);
  }
  return skipBin;
}

VecD StoreBinContentInVector(MapH &map, TString& hist){
  VecD content;
  for(unsigned int i=1; i<=map[hist]->GetNbinsX(); i++){
    content.push_back(map[hist]->GetBinContent(i));
  }
  if(debug){
    for(unsigned int i = 0;i<content.size();i++) printf("(%3i, %2.5d)",i,content[i]);
    cout << endl;
  }
  return content;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

TString GetChi2Cor(TMatrixD &m, vector<TString> &vecS, VecD &vecD, int nbins, int nbins_hh, int nbins_hl, int nbins_lh, int nbins_ll, bool uncor){

  cout << "Size input vectors - fits " << setw(5) << vecS.size() << " | data " << setw(5) << vecS.size() << " | #bins " << setw(5) << nbins << endl;

  bool uncorrelated = uncor;
  int n = nbins_hh+nbins_hl+nbins_lh+nbins_ll;
  if(n!=vecD.size() && n!=0) throw runtime_error("WRONG BINNING!");

  TString term = "(di - (gi)) * (dj - (gj)) * Vij-1";
  TString i_term = "(di - (gi)) * (j_term)";
  TString j_term = "(dj - (gj)) * Vij-1";
  TString function = "";
  TString function_ij = "";
  TString temp_term = term;
  TString temp_term_i = term;
  TString temp_term_j = term;
  TString temp_term_ij = term;

  int skip_1 = rand()%nbins_hh;
  int skip_2 = rand()%nbins_hl+nbins_hh;
  int skip_3 = rand()%nbins_lh+(nbins_hh+nbins_hl);
  int skip_4 = rand()%nbins_ll+(nbins_hh+nbins_hl+nbins_lh);

  cout << "Skip four bins:" << setw(4) << skip_1 << setw(4) << skip_2 << setw(4) << skip_3 << setw(4) << skip_4 << endl;

  vector<int> bins = {};
  for(unsigned int i=0; i<nbins; i++){
    bool skip = false;
    if(i==skip_1||i==skip_2||i==skip_3||i==skip_4) skip = true;
    if(!skip) bins.push_back(i);
  }

  int N = bins.size();
  TMatrixD m_trim = TMatrixD(N,N);
  VecD data;
  vector<TString> fits;

  if(bins.size()!=m_trim.GetNcols()){
    cout << "Used bins " << N << " | Matrix " << m_trim.GetNcols() << endl;
    throw runtime_error("Used bins not equal to trimmed Matrix!");
  }
  for(int i=0; i<N; i++){
    data.push_back(vecD[bins[i]]);
    fits.push_back(vecS[bins[i]]);
    for(int j=0; j<N; j++){
      m_trim[i][j] = m[bins[i]][bins[j]];
    }
  }

  // printCovMatrix(m, "standard", 100000, 3);
  // printCovMatrix(m_trim, "trimmed", 100000, 3);

  TDecompLU lu(m);
  TMatrixD imat = TMatrixD(m);
  lu.Invert(imat);
  TMatrixD test(m, TMatrixD::kMult, imat);
  TMatrixD unit(TMatrixD::kUnit, m);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 1.e-8);
  if (!ok){
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<N; ++j)  printf("%7.2f ", m[i][j] );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("(%3d,%3d)   %9.8f",i,i,m[i][i]*1e8);
      if ( m[i][i] < 1.e-8 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // if(debug) printCovMatrix(imat, "Inverted Matrix", 100000);
  // if(debug) printCovMatrix(unit, "Identity Test", 1, 1);

  TH2D* h_final = TMatrixDtoTH2D(imat, N, 0, N, channel+year+"norm");
  DrawCov(h_final, save_nfs+"CovMatrix_inverted_"+channel+"_"+year, "#it{m}_{W}", 0.3);

  cout << "Before Check -- Matrix: " << m.GetNcols() << "\t | used bins: " << bins.size() << "\t | data bin content: " << vecD.size() << "\t | terms: " << vecS.size() << endl;
  cout << "Final Check --- Matrix: " << imat.GetNcols() << "\t | used bins: " << bins.size() << "\t | data bin content: " << data.size() << "\t | terms: " << fits.size() << endl;

  int b_hh = nbins_hh;
  int b_hl = b_hh+nbins_hl;
  int b_lh = b_hl+nbins_lh;
  int b_ll = b_lh+nbins_ll;
  cout << setw(5) << nbins_hh << setw(5) << nbins_hl << setw(5) << nbins_lh << setw(5) << nbins_ll << endl;
  cout << setw(5) << b_hh << setw(5) << b_hl << setw(5) << b_lh << setw(5) << b_ll << endl;

  for(unsigned int i = 0; i<N; i++){

    TString temp_comp_j = "";
    // double cold=0;
    // double cnew=0;
    for(unsigned int j = 0; j<N; j++){

      if(abs(imat[i][j])<=10e-30){cout << "check " << setw(3) << i << setw(3) << j << endl; continue;}

      temp_term_j = j_term;
      temp_term_j.ReplaceAll("dj", dtos(data[j], precision));
      temp_term_j.ReplaceAll("gj", fits[j]);
      temp_term_j.ReplaceAll("Vij-1", dtos(imat[i][j], precision));
      temp_comp_j += " + "+temp_term_j;

      // cout << "\t\t" << setw(3) << j << temp_term_j << endl;

    }

    temp_term_i = i_term;
    temp_term_i.ReplaceAll("di", dtos(data[i], precision));
    temp_term_i.ReplaceAll("gi", fits[i]);
    // cout << "\t" << setw(3) << i << temp_term_i << endl;
    temp_term_i.ReplaceAll("j_term", temp_comp_j);
    function_ij += " + "+temp_term_i;

  }
  return function_ij;
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
  VecDD points = FindXY(chi2, z+at, x-0.8, x+0.8, y-1.5, y+1.5, runs, acc);
  auto stop  = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  cout << "Numeric solution for "+hist.ReplaceAll(cut, "")+"s 1\u03C3 estimation took " << GREEN << duration.count() << "ms" << RESET << endl;
  if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;
  return points;
}

VecDD GetSigmaEllipse_alt(TF2* chi2, VecD &minimum, TString hist, double runs, double acc, double at, double range){
  if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";
  double z = minimum[2]; double y = minimum[1]; double x = minimum [0];
  auto start = high_resolution_clock::now(); // Calculation time - start
  VecDD points = FindXY_alt(chi2, z+at, x-range, x+range, y-range, y+range, runs, acc);
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
  D->SaveAs(save_nfs+"/chi2_"+id+".pdf");
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
  ccor->SaveAs(save_nfs+"/chi2_projection_"+xaxis+".pdf");
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void PoUinTXT(double& nom, double& up, double& down, TString variation){
  ofstream textfile;
  textfile.open(save_nfs+"/factor_projection_"+variation+".txt");
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
  // if(!name.Contains("cor")){
  cout << "sigma_cor: " << sigma_cor_center << endl;
  cout << "\u03C1: " << -rho << endl;
  cout << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
  cout << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
  cout << "Error f_{JMS}: " << sigma_jms << endl;
  // }
  cout << "====================================================" << endl <<endl;

  ofstream textfile;
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
  DrawCov(h_norm, save_nfs+"/SYS/cov_"+sysname+"_"+year+"_"+add, to_string(panel)+add);

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
  c->SaveAs(save_nfs+"/SYS/SYS_"+sysname+"_"+add+".pdf");

  // c_sys->cd(panel);
  sys->Draw("HIST ][");

  delete c;
  return matrix;
}
//
