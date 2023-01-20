#include "../include/CentralInclude.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>

using namespace std;

typedef map<int, TGraph2DErrors*> MapI2DGe;
typedef map<TString, MapI2DGe> MapSI2DGe;

typedef map<TString, TF2*> MapF2;
typedef map<TString, MapF2> MapFF2;
typedef map<int, TF2*> MapIF2;
typedef map<TString, MapIF2> MapSIF2;

// ---------------------------------------------------------------------------
// Declare functions

// MapHHH SetUpMap(VecTS hists);
void printVector(VecD vector, TString info);
void printVector(VecI vector, TString info);
void AddTwoMaps(MapHHH map1, MapHHH map2, int option);
void Plot2DGraph(TGraph2DErrors* bin_fit, int bin, TString title, VecD content, VecD error);
void DrawPoints(const VecDD& vec1, const VecDD& vec2);
void DrawTestEllipse(const VecDD& vec1);
void Draw2DChi2(TF2* chi2, VecDD points, VecD minimum, TString hist, TEllipse* ell);
void drawChi2Projection(TF1* chi2function, TString xaxis, VecD xRange, VecD yRange);
void PoUinTXT(double nom, double up, double down, TString variation);
bool sortcolx( const VecD& v1, const VecD& v2 );
void JMSinTXT(VecDD points, VecD PoU, double sig2cor, double sig1cor, double sig1jec, double sigJMS, TString variation);
MapF2 GetFits(MapSI2DGe map, MapVIII peak, TString hists);
MapHH GetHistograms(TString process, TString h_name);
MapHHH RebinAndNormalize(MapHHH process, int width);
MapSI2DGe Creat2DGraph(MapHHHH map, MapVIII peak, VecTS hists);
TF2* GetChi2(MapF2 fits, MapHHHH nhists, MapVIII peak, TString hist);
TEllipse* ApproxEllipse(const VecD& mid, const VecD& uu, const VecD& dd, const VecD& ud, const VecD& du);
TPolyLine3D* SetPolyLine(double x1, double x2, double y1, double y2, double z1, double z2, int color);
VecD ExtractJMSValues(const TF1* chi2function, const double& var);
VecD GetMinimumChi2(TF2* chi2, TString hist);
VecDD Get1sigmaEllipse(TF2* chi2, VecD minimum, TString hist, double runs, double acc);

bool debug = true;
TString save, year, channel;
TString schannel = "";
VecTS hists = {
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hl",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_lh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_ll"
};
TString cut = "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_";
int space;

int main(int argc, char* argv[]){

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  // Input ---------------------------------------------------------------------
  year = argv[1];
  channel = argv[2];
  if(channel.EqualTo("muon")||channel.EqualTo("elec")) schannel = "_"+channel;
  gErrorIgnoreLevel = kWarning;

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
  int limit = 100;
  // int limit = (channel.EqualTo("muon")||channel.EqualTo("elec"))?50:100;
  TString slimit = to_string(limit);


  // Create Directories --------------------------------------------------------
  if(debug) cout << "Create Directiories ..." << endl;

  save = get_save_path()+"/Plots";
  save = creat_folder_and_path(save, "JetCorrections");
  // save = creat_folder_and_path(save, "Dennis");
  save = creat_folder_and_path(save, "fit");
  save = creat_folder_and_path(save, year);
  save = creat_folder_and_path(save, channel);
  save = creat_folder_and_path(save, slimit);
  cout << save << endl;

  save = creat_folder_and_path(save, "BinWidth_"+to_string(bin_width));
  creat_folder(save, "all_bins");
  creat_folder(save, "projection");

  // =====================================================================================
  // === Get Histograms                                                                ===
  // =====================================================================================
  cout << "Start collecting Histograms ..." << endl;
  // ---------------------------------------------------------------------------
  // I load all histograms; This code was created after studies ended. Only one
  // option is available. If necessary include other options as well (no bins etc.)

  TString w_mass_hh = hists[0];
  TString w_mass_hl = hists[1];
  TString w_mass_lh = hists[2];
  TString w_mass_ll = hists[3];

  MapHHHH histograms;
  MapHHH m_data, m_ttbar, m_st, m_wjets, m_other, m_JECup, m_JECdown, m_CORup, m_CORdown;

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
  }

  if(debug) cout << "\t ... Get Background" << endl;
  MapHHH m_bkg = m_other;
  for(MapHHH map: {m_wjets, m_st}) AddTwoMaps(m_bkg, map, 1);

  AddTwoMaps(m_ttbar, m_bkg, 1);
  AddTwoMaps(m_JECup, m_bkg, 1);
  AddTwoMaps(m_JECdown, m_bkg, 1);
  AddTwoMaps(m_CORup, m_bkg, 1);
  AddTwoMaps(m_CORdown, m_bkg, 1);


  if(debug) cout << "\t ... Normalize and Rebin" << endl;
  MapHHH m_data_norm = RebinAndNormalize(m_data, bin_width);
  MapHHH m_ttbar_norm = RebinAndNormalize(m_ttbar, bin_width);
  MapHHH m_JECup_norm = RebinAndNormalize(m_JECup, bin_width);
  MapHHH m_JECdown_norm = RebinAndNormalize(m_JECdown, bin_width);
  MapHHH m_CORup_norm = RebinAndNormalize(m_CORup, bin_width);
  MapHHH m_CORdown_norm = RebinAndNormalize(m_CORdown, bin_width);

  // =====================================================================================
  // === Get considered bins (exclude empty ones in data and no peak bin)              ===
  // =====================================================================================

  if(debug) cout << "\t ... Get empty/peak bins" << endl;
  // TODO: adjust limits for non combined years
  // int limit = (channel.EqualTo("muon")||channel.EqualTo("elec"))?50:100;
  // int limit = 100;
  MapVIII empty_bins, peak_bins;
  for(TString h: hists){
    for(TString c: channels){
      for(TString y: years){
        if(!c.EqualTo(channel)||!y.EqualTo(year)) continue;
        empty_bins[h][c][y] = bins_empty(m_data_norm[h][c][y]);
        peak_bins[h][c][y] = bins_upper_limit((TH1F*) m_data[h][c][y]->Rebin(bin_width), limit);
        peak_bins[h][c][y] = bins_upper_limit((TH1F*) m_data[h][c][y]->Rebin(bin_width), limit);
        if(debug) printVector(peak_bins[h][c][y], "PEAK:("+h+","+c+","+y+")");
      }
    }
  }

  // =====================================================================================
  // === Creat TGraphErrors for each bin                                               ===
  // =====================================================================================
  if(debug) cout << "Create TGraph2DErrors ... " << endl;

  MapHHHH m_all_norm;
  m_all_norm["data"] = m_data_norm;
  m_all_norm["ttbar"] = m_ttbar_norm;
  m_all_norm["JECup"] = m_JECup_norm;
  m_all_norm["JECdown"] = m_JECdown_norm;
  m_all_norm["CORup"] = m_CORup_norm;
  m_all_norm["CORdown"] = m_CORdown_norm;

  // MapSI2DGe Creat2DGraph(MapHHHH map, MapVIII peak, VecTS hists, TString channel)
  MapSI2DGe m_combine = Creat2DGraph(m_all_norm, peak_bins, hists);

  // =====================================================================================
  // === Creat Chi2 functions for each bin                                              ===
  // =====================================================================================

  if(debug) cout << "Get single Chi2 ... " << endl;
  // MapSIF2 GetFits(MapSI2DGe map, MapVIII peak, VecTS hists){
  if(debug) cout << "\t ... Get Fit functions" << endl;
  MapF2 m_fit_combine_hh = GetFits(m_combine, peak_bins, w_mass_hh);
  MapF2 m_fit_combine_hl = GetFits(m_combine, peak_bins, w_mass_hl);
  MapF2 m_fit_combine_lh = GetFits(m_combine, peak_bins, w_mass_lh);
  MapF2 m_fit_combine_ll = GetFits(m_combine, peak_bins, w_mass_ll);

  if(debug) cout << "\t ... Construct Chi2" << endl;
  MapFF2 m_fit_combine;
  m_fit_combine[w_mass_hh]=m_fit_combine_hh; m_fit_combine[w_mass_ll]=m_fit_combine_ll;
  m_fit_combine[w_mass_hl]=m_fit_combine_hl; m_fit_combine[w_mass_lh]=m_fit_combine_lh;

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
  VecDD ellipse_hh = Get1sigmaEllipse(chi2_hh, minimum_hh, w_mass_hh, 2000, 0.1);
  VecDD ellipse_hl = Get1sigmaEllipse(chi2_hl, minimum_hl, w_mass_hl, 2000, 0.1);
  VecDD ellipse_lh = Get1sigmaEllipse(chi2_lh, minimum_lh, w_mass_lh, 2000, 0.1);
  VecDD ellipse_ll = Get1sigmaEllipse(chi2_ll, minimum_ll, w_mass_ll, 2000, 0.1);

  TEllipse *ell_hh = new TEllipse();
  TEllipse *ell_hl = new TEllipse();
  TEllipse *ell_lh = new TEllipse();
  TEllipse *ell_ll = new TEllipse();

  if(debug) cout << "\t ... Draw Chi2" << endl;
  gStyle->SetPalette(kDeepSea); // kDeepSea kGreyScale kRust
  TColor::InvertPalette();

  Draw2DChi2(chi2_hh, ellipse_hh, minimum_hh, w_mass_hh, ell_hh);
  Draw2DChi2(chi2_hl, ellipse_hl, minimum_hl, w_mass_hl, ell_hl);
  Draw2DChi2(chi2_lh, ellipse_lh, minimum_lh, w_mass_lh, ell_lh);
  Draw2DChi2(chi2_ll, ellipse_ll, minimum_ll, w_mass_ll, ell_ll);

  // =====================================================================================
  // === Creat chi2 function for JMS                                                   ===
  // =====================================================================================

  // precision ---------------------------------------------------------------------
  int precision = 6;
  space = precision+4;
  cout << fixed << setprecision(precision);

  cout << endl;
  cout << "-------------------" << endl;
  cout << "-       JMS       -" << endl;
  cout << "-------------------" << endl;
  cout << "Prepare results for JMS ..." << endl;

  TString function_hh = chi2_hh->GetFormula()->GetParName(0);
  TString function_hl = chi2_hl->GetFormula()->GetParName(0);
  TString function_lh = chi2_lh->GetFormula()->GetParName(0);
  TString function_ll = chi2_ll->GetFormula()->GetParName(0);
  TString function_JMS = function_hh+function_hl+function_lh+function_ll;
  TF2* chi2_JMS = new TF2("JMS", function_JMS, -3, 3, -3, 3);

  VecD nominal_JMS = GetMinimumChi2(chi2_JMS, "JMS");
  VecDD ellipse_JMS = Get1sigmaEllipse(chi2_JMS, nominal_JMS, "JMS", 4000, 0.01);

  cout << "Extract factors from 1\u03C3 ellipse ..." << endl;
  if(debug) cout << "\t ... extract x&y values of points" << endl;
  VecD xJMS, yJMS;
  for(unsigned int i=0; i<ellipse_JMS.size(); i++){
    xJMS.push_back(ellipse_JMS[i][0]);
    yJMS.push_back(ellipse_JMS[i][1]);
  }

  if(debug) cout << "\t ... extract min and max values for x&y" << endl;
  // To construct function which splits ellipse in two regions.
  // Line goes through xmin to xmax.
  double xmin = *min_element(xJMS.begin(),xJMS.end());
  double xmax = *max_element(xJMS.begin(),xJMS.end());
  double ixmin = find(xJMS.begin(), xJMS.end(), xmin)-xJMS.begin();
  double ixmax = find(xJMS.begin(), xJMS.end(), xmax)-xJMS.begin();

  double ymin = *min_element(yJMS.begin(),yJMS.end());
  double ymax = *max_element(yJMS.begin(),yJMS.end());
  double iymin = find(yJMS.begin(), yJMS.end(), ymin)-yJMS.begin();
  double iymax = find(yJMS.begin(), yJMS.end(), ymax)-yJMS.begin();

  VecD point_ymax = ellipse_JMS[iymax];
  VecD point_ymin = ellipse_JMS[iymin];
  VecD point_xmax = ellipse_JMS[ixmax];
  VecD point_xmin = ellipse_JMS[ixmin];

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

  VecDD upper, lower;
  for(unsigned int i=0;i<ellipse_JMS.size();i++){
    double xp = ellipse_JMS[i][0]; double yp = ellipse_JMS[i][1];
    double limit = fsplit->Eval(xp);
    if(yp<=limit) lower.push_back(ellipse_JMS[i]);
    else upper.push_back(ellipse_JMS[i]);
  }
  // DrawPoints(upper, lower);

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
  // cout << lmax << " ("<< ilmax << ")" << "\t" << lmin << " ("<< ilmin << ")" << endl;

  double umin = *min_element(dist_upper.begin(),dist_upper.end());
  double umax = *max_element(dist_upper.begin(),dist_upper.end());
  int iumin = find(dist_upper.begin(), dist_upper.end(), umin)-dist_upper.begin();
  int iumax = find(dist_upper.begin(), dist_upper.end(), umax)-dist_upper.begin();
  // cout << umax << " ("<< iumax << ")" << "\t" << umin << " ("<< iumin << ")" << endl;


  // =====================================================================================
  // === extract JMS factors                                                           ===
  // =====================================================================================

  if(debug) cout << "\t ... extract JMS factors" << endl;
  if(debug) cout << "\t\t ... from ellipse" << endl;
  VecD uu = {upper[iumin][0], upper[iumin][1]};
  VecD dd = {lower[ilmin][0], lower[ilmin][1]};
  VecD du = {upper[iumax][0], upper[iumax][1]};
  VecD ud = {lower[ilmax][0], lower[ilmax][1]};

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
  double sigma_cor_full = ymax; // just to avoid comfusion
  double rho = sqrt(1-sigma_cor_center/sigma_cor_full);

  cout << endl;
  if(debug) cout << "\t ... print all information" << endl;
  cout << "================================================" << endl;
  cout << "JMS factors ("+channel+", "+year+"):" << endl;
  cout << "----------------------------" << endl;
  cout << "nom (" << setw(space) << nominal_JMS[0] << ", " << setw(space) << nominal_JMS[1] << ")" << endl;
  cout << "uu  (" << setw(space) << uu[0] << ", " << setw(space) << uu[1] << ")" << endl;
  cout << "dd  (" << setw(space) << dd[0] << ", " << setw(space) << dd[1] << ")" << endl;
  cout << "du  (" << setw(space) << du[0] << ", " << setw(space) << du[1] << ")" << endl;
  cout << "ud  (" << setw(space) << ud[0] << ", " << setw(space) << ud[1] << ")" << endl;
  cout << endl;
  cout << "sigma_cor: " << sigma_cor_center << endl;
  cout << "\u03C1: " << -rho << endl;
  cout << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
  cout << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
  cout << "Error f_{JMS}: " << sigma_jms << endl;
  cout << "================================================" << endl;

  ofstream textfile;
  textfile.open(save+"/JMS.txt");
  textfile << right;
  textfile << "nom (" << setw(space) << nominal_JMS[0] << ", " << setw(space) << nominal_JMS[1] << ")" << endl;
  textfile << "uu  (" << setw(space) << uu[0] << ", " << setw(space) << uu[1] << ")" << endl;
  textfile << "dd  (" << setw(space) << dd[0] << ", " << setw(space) << dd[1] << ")" << endl;
  textfile << "du  (" << setw(space) << du[0] << ", " << setw(space) << du[1] << ")" << endl;
  textfile << "ud  (" << setw(space) << ud[0] << ", " << setw(space) << ud[1] << ")" << endl;
  textfile << endl;
  textfile << "sigma_cor: " << sigma_cor_center << endl;
  textfile << "\u03C1: " << -rho << endl;
  textfile << "PoU COR: " << pow(sigma_cor/PoU_COR[0],2) << endl;
  textfile << "PoU JEC: " << pow(sigma_jec/PoU_JEC[0],2) << endl;
  textfile << "Error f_{JMS}: " << sigma_jms << endl;
  textfile.close();

  VecDD ordered_ellipse = upper;
  ordered_ellipse.insert(ordered_ellipse.end(), lower.begin(), lower.end());
  // DrawTestEllipse(ordered_ellipse);
  TEllipse* ellipse = ApproxEllipse(nominal_JMS, upper[iumin], lower[ilmin], upper[iumax], lower[ilmax]);
  Draw2DChi2(chi2_JMS, ordered_ellipse, nominal_JMS, "JMS", ellipse);

  double te = 0.324328574805972405743285992347054;
  cout << te << "\t" << to_string(te) << "\t" << dtos(te, 10) << endl;
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
  if(process.EqualTo("data")){ collection_muon = data_muon; collection_elec = data_elec;}
  else if(process.EqualTo("ttbar")){ collection_muon = ttbar_muon; collection_elec = ttbar_elec;}
  else if(process.EqualTo("wjets")){ collection_muon = wjets_muon; collection_elec = wjets_elec;}
  else if(process.EqualTo("st")){ collection_muon = st_muon; collection_elec = st_elec;}
  else if(process.EqualTo("other")){ collection_muon = other_muon; collection_elec = other_elec;}
  else if(process.EqualTo("JECup")){ collection_muon = jec_up_muon; collection_elec = jec_up_elec;}
  else if(process.EqualTo("JECdown")){ collection_muon = jec_down_muon; collection_elec = jec_down_elec;}
  else if(process.EqualTo("CORup")){ collection_muon = cor_up_muon; collection_elec = cor_up_elec;}
  else if(process.EqualTo("CORdown")){ collection_muon = cor_down_muon; collection_elec = cor_down_elec;}
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
void AddTwoMaps(MapHHH map1, MapHHH map2, int option){
  // MapHHH map;
  for(auto hist: map1){
    TString h = hist.first;
    for(auto channel: hist.second){
      TString c = channel.first;
      for(auto year: channel.second){
        TString y = year.first;
        // map[h][c][y] = AddHists(year.second, bkg[h][c][y], 1);
        year.second->Add(map2[h][c][y], option);
      }
    }
  }
  // return map;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
MapHHH RebinAndNormalize(MapHHH process, int width){
  MapHHH map;
  for(auto hist: process){
    TString h = hist.first;
    for(auto channel: hist.second){
      TString c = channel.first;
      for(auto year: channel.second){
        TString y = year.first;
        TH1F* h_rebin = (TH1F*) year.second->Rebin(width);
        map[h][c][y] = Normalize(h_rebin); // already sets normalized error
      }
    }
  }
  return map;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void printVector(VecD vector, TString info) {
  cout << info << "; (Size: " << vector.size() << ") ";
  for(double value: vector) cout << value << ", ";
  cout << endl;
}

void printVector(VecI vector, TString info) {
  cout << info << "; (Size: " << vector.size() << ") ";
  for(double value: vector) cout << value << ", ";
  cout << endl;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
MapSI2DGe Creat2DGraph(MapHHHH map, MapVIII peak, VecTS hists){
  if(debug) cout << "\t ... Inside Create2DGraph - " << channel << endl;
  VecD factor_x = {0.0,  1.0, -1.0,  0.0,  0.0};
  VecD factor_y = {0.0,  0.0,  0.0,  1.0, -1.0};
  VecD dummy = {0.0,  0.0,  0.0,  0.0,  0.0};
  VecD content, error;
  MapI2DGe storage;
  MapSI2DGe bins;

  for(TString hist: hists){ // defined in main()
    if(debug) cout << "new hist: " << hist << endl;
    for(unsigned int i=0; i<peak[hist][channel][year].size(); i++)
    {
      // ----------------------------------------------
      // We only want to combine all years and
      // consider all channels (muon, elec and combine)

      int bin = peak[hist][channel][year].at(i);
      TString name = to_string(bin);
      content.push_back(map["ttbar"][hist][channel][year]->GetBinContent(bin));
      content.push_back(map["JECup"][hist][channel][year]->GetBinContent(bin));
      content.push_back(map["JECdown"][hist][channel][year]->GetBinContent(bin));
      content.push_back(map["CORup"][hist][channel][year]->GetBinContent(bin));
      content.push_back(map["CORdown"][hist][channel][year]->GetBinContent(bin));

      error.push_back(map["ttbar"][hist][channel][year]->GetBinError(bin));
      error.push_back(map["JECup"][hist][channel][year]->GetBinError(bin));
      error.push_back(map["JECdown"][hist][channel][year]->GetBinError(bin));
      error.push_back(map["CORup"][hist][channel][year]->GetBinError(bin));
      error.push_back(map["CORdown"][hist][channel][year]->GetBinError(bin));

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

  bin_fit->SetTitle(title);
  bin_fit->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
  bin_fit->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
  bin_fit->GetHistogram()->GetZaxis()->SetTitle("a.u.");
  bin_fit->GetHistogram()->GetXaxis()->SetTitleOffset(2.5);
  bin_fit->GetHistogram()->GetYaxis()->SetTitleOffset(2.5);
  bin_fit->GetHistogram()->GetZaxis()->SetTitleOffset(2.6);
  bin_fit->GetHistogram()->GetXaxis()->SetTitleSize(0.026);
  bin_fit->GetHistogram()->GetYaxis()->SetTitleSize(0.026);
  bin_fit->GetHistogram()->GetZaxis()->SetTitleSize(0.026);
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
  B->SetLeftMargin(0.15);
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
MapF2 GetFits(MapSI2DGe map, MapVIII peak, TString hist){
  if(debug) cout << "\n\t ... Inside GetFits - " << channel << endl;

  MapF2 functions;
  vector<VecTS> name_fits = {{"linear"}, {"mixed (xy2)", "mixed (x2y)"}, {"mixed (xyy2)", "mixed (xx2y)"}, {"quadratic"}, {"polynomial of order 2"}};
  VecII ndfs  = {{2}, {2, 2}, {1, 1}, {2}, {0}}; // ndf = n_points - n_parameters = 5 - n_para
  VecII n_fit_para  = {{3}, {3, 3}, {4, 4}, {3}, {5}};

  TString replace = hist; replace.ReplaceAll("comparison_topjet_xcone_pass_rec/", "");
  cout << "Start fitting " << replace << " ... " << endl;
  for(unsigned int i=0; i<peak[hist][channel][year].size(); i++)
  {
    int bin = peak[hist][channel][year].at(i);
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
      if(prob>0.05 || ndfs[order][0]==0){ // prob > 0.05 is common in data science
        if(debug) cout << " -> Bin " << name << ": A " << GREEN << name_fits[order][index] << RESET << " fit fullfills the criterion (" << prob << ")" << endl;
        functions[name]=fits[order][index];

        // Information which can be extracted later.
        // This way, no further information then the function itself has to be provided
        // Necessary for Chi2 construction to decide which formula to use
        const char* info = name_fits[order][index];
        functions[name]->SetParName(0, info);
        break;
      }
    }
  }
  return functions;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
TF2* GetChi2(MapF2 fits, MapHHHH nhists, MapVIII peak, TString hist){
  int ipar = 0;
  TString fchi2 = "";
  for(unsigned int i=0; i<peak[hist][channel][year].size(); i++)
  {
    int bin = peak[hist][channel][year].at(i);
    TString sbin = to_string(bin);
    TFormula *formula = fits[sbin]->GetFormula();
    int npar = formula->GetNpar();
    TString fname = formula->GetParName(0);

    double data_content = nhists["data"][hist][channel][year]->GetBinContent(bin);
    double data_error = nhists["data"][hist][channel][year]->GetBinError(bin); // stat. error
    double tt_content = nhists["ttbar"][hist][channel][year]->GetBinContent(bin); // fit error (estimated)
    double tt_error = nhists["ttbar"][hist][channel][year]->GetBinError(bin); // fit error (estimated)
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
VecDD Get1sigmaEllipse(TF2* chi2, VecD minimum, TString hist, double runs, double acc){
  if(debug) cout << "Chi2 - Get 1\u03C3 Ellipse\n";
  double z = minimum[2]; double y = minimum[1]; double x = minimum [0];
  auto start = high_resolution_clock::now(); // Calculation time - start
  VecDD points = FindXY(chi2, z+2.3, x-1.8, x+1.8, y-1.8, y+1.8, runs, acc);
  auto stop  = high_resolution_clock::now();  // Calculation time - stop
  auto duration = duration_cast<milliseconds>(stop - start);
  // cout << "Numeric solution for "+hist.ReplaceAll(cut, "")+"s 1\u03C3 estimation took " << GREEN << duration.count()/1000 << "s" << RESET << endl;
  cout << "Numeric solution for "+hist.ReplaceAll(cut, "")+"s 1\u03C3 estimation took " << GREEN << duration.count() << "ms" << RESET << endl;
  if(debug) cout << "Number Points at 1\u03C3: " << points.size() << endl;
  return points;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void Draw2DChi2(TF2* chi2, VecDD points, VecD minimum, TString hist, TEllipse* ell){

  TGraph* test = new TGraph();
  test->SetPoint(0, 0.3, 0.3);
  test->SetMarkerColor(kBlack);

  TGraph2D *zmin_point = new TGraph2D();
  zmin_point->SetPoint(0, minimum[0], minimum[1], minimum[2]);
  zmin_point->SetMarkerColor(kBlack);
  zmin_point->SetMarkerStyle(kFullCircle);
  zmin_point->SetMarkerSize(0.5);

  TGraph2D *sigma_points = new TGraph2D();
  if(points.size()>0){ // precaution, ellipse for bins is not important
    for(unsigned int i=0; i<points.size(); i++) sigma_points->SetPoint(i, points[i][0], points[i][1], points[i][2]);
  }
  sigma_points->SetMarkerColor(kRed);
  sigma_points->SetMarkerStyle(kFullCircle);
  sigma_points->SetMarkerSize(0.1);

  // VecD x, y, z;
  // for(auto vec:points){
  //   x.push_back(vec[0]);
  //   y.push_back(vec[1]);
  //   z.push_back(vec[2]);
  // }
  // x.push_back(points[0][0]); y.push_back(points[0][1]); z.push_back(points[0][2]);
  // TPolyLine3D* line = new TPolyLine3D(points.size()+1, &x[0], &y[0], &z[0]);
  // line->SetLineColor(kRed);
  // line->SetLineWidth(1);

  //TEllipse* ApproxEllipse(const VecD& mid, const VecD& uu, const VecD& dd, const VecD& ud, const VecD& du)
  ell->SetLineColor(kGreen+2);
  ell->SetFillStyle(0);

  chi2->SetTitle("");
  chi2->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
  chi2->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
  chi2->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
  chi2->GetHistogram()->GetZaxis()->CenterTitle();
  chi2->GetHistogram()->GetZaxis()->SetTitleOffset(0.5);
  chi2->SetFillStyle(1000);
  chi2->SetLineWidth(1);
  chi2->SetRange(-3, -3, 3, 3);
  chi2->SetContour(50);        // Contours

  TCanvas *D = new TCanvas(hist,"D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  // test->Draw("P");
  chi2->Draw("cont4z");
  D->SetTheta(90);
  D->SetPhi(0);
  zmin_point->Draw("SAME P");
  sigma_points->Draw("SAME P");
  // line->Draw("SAME LINE");
  // ell->Draw("SAME");
  TString id = hist.ReplaceAll("comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_", "");
  D->SaveAs(save+"/chi2_"+id+".pdf");
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
void PoUinTXT(double nom, double up, double down, TString variation){
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
void DrawPoints(const VecDD& vec1, const VecDD& vec2){

  TGraph2D *range = new TGraph2D();
  range->SetPoint(0, 0, -0.6, vec1[0][2]);
  range->SetPoint(1, 0, 1.1, vec1[0][2]);
  range->SetPoint(2, 1.0, -0.6, vec1[0][2]);
  range->SetPoint(3, 1.0, 1.1, vec1[0][2]);
  range->SetMarkerColor(kGray);
  range->SetMarkerStyle(kFullCircle);
  range->SetMarkerSize(0.1);

  TGraph2D *points1 = new TGraph2D();
  for(unsigned int i=0; i<vec1.size(); i++) points1->SetPoint(i, vec1[i][0], vec1[i][1], vec1[i][2]);
  points1->SetMarkerColor(kRed);
  points1->SetMarkerStyle(kFullCircle);
  points1->SetMarkerSize(0.1);

  TGraph2D *points2 = new TGraph2D();
  for(unsigned int i=0; i<vec2.size(); i++) points2->SetPoint(i, vec2[i][0], vec2[i][1], vec2[i][2]);
  points2->SetMarkerColor(kBlue);
  points2->SetMarkerStyle(kFullCircle);
  points2->SetMarkerSize(0.1);

  points1->SetTitle("");
  // points1->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
  // points1->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
  points1->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
  points1->GetHistogram()->GetZaxis()->CenterTitle();
  points1->GetHistogram()->GetZaxis()->SetTitleOffset(0.5);
  points1->SetFillStyle(1000);
  points1->SetLineWidth(1);

  TCanvas *D = new TCanvas("testpoints","D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  D->SetTheta(90);
  D->SetPhi(0);
  range->Draw("P");
  points1->Draw("SAME P");
  points2->Draw("SAME P");
  D->SaveAs(save+"/test_points_sort.pdf");
  D->Clear();
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void DrawTestEllipse(const VecDD& vec1){

  TGraph2D *range = new TGraph2D();
  range->SetPoint(0, 0, -0.6, vec1[0][2]);
  range->SetPoint(1, 0, 1.1, vec1[0][2]);
  range->SetPoint(2, 1.0, -0.6, vec1[0][2]);
  range->SetPoint(3, 1.0, 1.1, vec1[0][2]);
  range->SetMarkerColor(kGray);
  range->SetMarkerStyle(kFullCircle);
  range->SetMarkerSize(0.1);

  TGraph2D *points1 = new TGraph2D();
  for(unsigned int i=0; i<vec1.size(); i++) points1->SetPoint(i, vec1[i][0], vec1[i][1], vec1[i][2]);
  points1->SetMarkerColor(kRed);
  points1->SetMarkerStyle(kFullCircle);
  points1->SetMarkerSize(0.1);

  points1->SetTitle("");
  // points1->GetHistogram()->GetXaxis()->SetTitle("JEC factor");
  // points1->GetHistogram()->GetYaxis()->SetTitle("Additional XCone correction factor");
  points1->GetHistogram()->GetZaxis()->SetTitle("#chi^{2}");
  points1->GetHistogram()->GetZaxis()->CenterTitle();
  points1->GetHistogram()->GetZaxis()->SetTitleOffset(0.5);
  points1->SetFillStyle(1000);
  points1->SetLineWidth(1);

  VecD x, y, z;
  for(auto vec:vec1){
    x.push_back(vec[0]);
    y.push_back(vec[1]);
    z.push_back(vec[2]);
  }
  x.push_back(vec1[0][0]); y.push_back(vec1[0][1]); z.push_back(vec1[0][2]);
  TPolyLine3D* line = new TPolyLine3D(vec1.size()+1, &x[0], &y[0], &z[0]);
  line->SetLineColor(kBlue);

  TCanvas *D = new TCanvas("testpoints","D", 600, 600);
  D->SetRightMargin(0.12);
  D->SetLogz();
  D->SetTheta(90);
  D->SetPhi(0);
  range->Draw("P");
  // points1->Draw("SAME LINE");
  points1->Draw("SAME P");
  // line->Draw("SAME");
  D->SaveAs(save+"/test_ellipse_line.pdf");
  D->Clear();
}
