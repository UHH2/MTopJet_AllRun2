#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"

#include <fstream>
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include <TMath.h>
#include <TString.h>

// #include <io.h>
using namespace std;

void draw2(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color, vector<TString> names, TString name);
void draw(TF2* dummy, vector<TPolyMarker3D*> min, vector<TPolyMarker3D*> graph, vector<int> color, vector<int> style,  vector<TString> names, TString name);
void CompareAreas(vector<TEllipse*> ellipse, vector<TString> names);
MapVD ExtractFromFile(TString path, TString c);
TPolyMarker3D* ExtractFromROOT(TString name, TString year, TString channel, TString sys, TString bins = "100");
TF2* ExtractFunctionFromROOT(TString year, TString channel, TString sys);
VecD ExtractPoints(TString line, const TString jms);
TEllipse* build_ellipse(MapVD points, TString name="");
TGraph* graph_multiple_points(MapVD points);
inline void set_ellipse(TEllipse *ellipse, int color);
inline void set_graph(TGraph *graph, int color);
inline void set_graph(TPolyMarker3D *graph, int color, double size);
inline double angle_to_xaxis(const VecD vec);
inline double distance_points(const VecD vec1, const VecD vec2);
inline double length_vector(const VecD vec);
inline VecD subtract_two_vectors(const VecD vec1, const VecD vec2);

bool debug = false;
TString save;

vector<double> lumis = {138., 36., 41., 60.};

int main(int argc, char* argv[]){

  gErrorIgnoreLevel = kWarning;
  TString save_path = get_save_path();

  TString dir = save_path+"/JetCorrections/fit/combine/";
  save = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/Ellipses/";
  VecTS channels = {"combine", "muon", "elec"};
  TString width = "/100/BinWidth_1/";
  TString file_txt = "JMS.txt";

  MapVDD JMSfactor;

  cout << "\nExtract Points ... " << endl;
  // JMSfactor["combine_uncor"]        = ExtractFromFile(save_path+"/JetCorrections/fit/combine/combine/100/BinWidth_1/JMS.txt", "combine");
  // JMSfactor["combine_cor"]          = ExtractFromFile(save_path+"/JetCorrections/fit/combine/combine/100/BinWidth_1/JMScorrelated.txt", "combine");
  // JMSfactor["combine_cor_syst"]     = ExtractFromFile(save_path+"/JetCorrections/fit/combine/combine/100/BinWidth_1/syst/JMScorrelated.txt", "combine");
  // JMSfactor["combine_lin_cor"]      = ExtractFromFile(save_path+"/JetCorrections/fit/combine/combine/100/BinWidth_1/lin/JMScorrelated.txt", "combine");
  // JMSfactor["combine_cor_noBKGsys"] = ExtractFromFile(save_path+"/JetCorrections/fit/combine/combine/100/BinWidth_1/noBKGsys/JMScorrelated.txt", "elec");
  //
  // JMSfactor["combine_add"] = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/combine/100/BinWidth_1/JMS.txt", "combine");
  // JMSfactor["muon_add"]    = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/muon/100/BinWidth_1/JMS.txt", "muon");
  // JMSfactor["elec_add"]    = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/elec/100/BinWidth_1/JMS.txt", "elec");
  //
  // JMSfactor["combine_V20"] = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/combine/100/BinWidth_1/JMS.txt", "elec");
  // JMSfactor["muon_V20"]    = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/muon/100/BinWidth_1/JMS.txt", "elec");
  // JMSfactor["elec_V20"]    = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/elec/100/BinWidth_1/JMS.txt", "elec");
  //
  // JMSfactor["2016_combine"]     = ExtractFromFile(save_path+"/JetCorrections/fit/2016/combine/40/BinWidth_1/JMS.txt", "combine");
  // JMSfactor["2017_combine"]     = ExtractFromFile(save_path+"/JetCorrections/fit/2017/combine/40/BinWidth_1/JMS.txt", "combine");
  // JMSfactor["2018_combine"]     = ExtractFromFile(save_path+"/JetCorrections/fit/2018/combine/40/BinWidth_1/JMS.txt", "combine");
  // JMSfactor["2016_combine_cor"] = ExtractFromFile(save_path+"/JetCorrections/fit/2016/combine/40/BinWidth_1/lin/JMScorrelated.txt", "combine");
  // JMSfactor["2017_combine_cor"] = ExtractFromFile(save_path+"/JetCorrections/fit/2017/combine/40/BinWidth_1/lin/JMScorrelated.txt", "combine");
  // JMSfactor["2018_combine_cor"] = ExtractFromFile(save_path+"/JetCorrections/fit/2018/combine/40/BinWidth_1/lin/JMScorrelated.txt", "combine");

  // for(TString c: channels){
  //   if(debug) cout << "\t ... channel "+c << endl;
  //   TString file = dir+c+width+file_txt;
  //   JMSfactor[c] = ExtractFromFile(file, c);
  // }


  map<TString, TF2*> JMS_chi2;
  // JMS_chi2["combine_cor"] = ExtractFunctionFromROOT("combine", "combine", "");
  // JMS_chi2["combine_cor_syst_2"] = ExtractFunctionFromROOT("combine", "combine", "syst_2");

  map<TString, TPolyMarker3D*> JMS_min;
  JMS_min["combine_cor"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "");
  JMS_min["combine_cor_syst"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst");
  JMS_min["combine_cor_syst_noJER"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst_noJER");
  JMS_min["combine_cor_fitall"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "_fitall");
  JMS_min["combine_cor_syst_peak"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst_peak");
  JMS_min["combine_cor_syst_peak_m4"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst_peak_m4");
  JMS_min["combine_cor_syst_peak_uncor"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst_peak_uncor");
  JMS_min["combine_cor_stat_peak"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "stat_peak");
  JMS_min["epjc"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "syst_peak_m4");
  JMS_min["arxiv"] = ExtractFromROOT("Graphs/JMS_nominal", "combine", "combine", "_fitall");

  map<TString, TPolyMarker3D*> JMS_ellipse;
  // === For different methods
  JMS_ellipse["combine_cor"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "");
  JMS_ellipse["combine_cor_syst"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst");
  JMS_ellipse["combine_cor_syst_noJER"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst_noJER");
  JMS_ellipse["combine_cor_fitall"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "_fitall");
  JMS_ellipse["combine_cor_syst_peak"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst_peak");
  JMS_ellipse["combine_cor_syst_peak_m4"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst_peak_m4");
  JMS_ellipse["combine_cor_syst_peak_uncor"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst_peak_uncor");
  JMS_ellipse["combine_cor_stat_peak"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "stat_peak");
  JMS_ellipse["epjc"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "syst_peak_m4");
  JMS_ellipse["arxiv"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "_fitall");

  // === For single bins
  JMS_ellipse["jms"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "quick_test_2_1");
  JMS_ellipse["all"] = ExtractFromROOT("Graphs/JMS_ellipse_all", "combine", "combine", "quick_test_2_1");
  JMS_ellipse["hh"] = ExtractFromROOT("Graphs/JMS_ellipse_hh", "combine", "combine", "quick_test_2_1");
  JMS_ellipse["hl"] = ExtractFromROOT("Graphs/JMS_ellipse_hl", "combine", "combine", "quick_test_2_1");
  JMS_ellipse["lh"] = ExtractFromROOT("Graphs/JMS_ellipse_lh", "combine", "combine", "quick_test_2_1");
  JMS_ellipse["ll"] = ExtractFromROOT("Graphs/JMS_ellipse_ll", "combine", "combine", "quick_test_2_1");

  // === For bin width = 3
  JMS_ellipse["bin_width=3"] = ExtractFromROOT("Graphs/JMS_ellipse", "combine", "combine", "tmp", "3");


  // Start drawing ---------------------------------------------------------------------------------
  cout << "Start plotting ... " << endl;
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  TF2* dummy = new TF2("dummy","x*y", 1.5, 1.5);
  // draw(dummy, {JMS_min["combine_cor"], JMS_min["combine_cor_syst"], JMS_min["combine_cor_syst_noJER"]}, {JMS_ellipse["combine_cor"], JMS_ellipse["combine_cor_syst"], JMS_ellipse["combine_cor_syst_noJER"]}, {kRed+2, kBlue+2, kGreen+2}, {kOpenTriangleUp, kOpenTriangleDown, kOpenCircle}, {"No syst", "syst", "syst no JER"}, "syst_comparison_cor");
  // draw(dummy, {JMS_min["combine_cor"], JMS_min["combine_cor_syst"]}, {JMS_ellipse["combine_cor"], JMS_ellipse["combine_cor_syst"]}, {kRed+2, kBlue+2}, {kOpenTriangleUp, kOpenTriangleDown}, {"No syst", "syst"}, "syst_comparison_cor");
  // draw(dummy, {JMS_min["combine_cor_syst_peak"], JMS_min["combine_cor_fitall"]}, {JMS_ellipse["combine_cor_syst_peak"], JMS_ellipse["combine_cor_syst"]}, {kRed+2, kBlue+2}, {kOpenTriangleUp, kOpenTriangleDown}, {"Only Peak", "Peak"}, "syst_comparison_peak");
  draw(dummy, {JMS_min["combine_cor_syst_peak"], JMS_min["combine_cor_syst_peak_m4"]}, {JMS_ellipse["combine_cor_syst_peak"], JMS_ellipse["combine_cor_syst_peak_m4"]}, {kRed+2, kBlue+2}, {kOpenTriangleUp, kOpenTriangleDown}, {"Only Peak", "Peak"}, "syst_comparison_peak");
  // draw(dummy, {JMS_min["combine_cor_syst_peak"], JMS_min["combine_cor_syst_peak_uncor"], JMS_min["combine_cor_stat_peak"]}, {JMS_ellipse["combine_cor_syst_peak"], JMS_ellipse["combine_cor_syst_peak_uncor"], JMS_ellipse["combine_cor_stat_peak"]}, {kRed+2, kBlue+2, kGreen+2}, {kOpenTriangleUp, kOpenTriangleDown, kOpenCircle}, {"Correlated", "Uncorrelated", "Only stat"}, "syst_comparison_peak_uncerts");
  // draw(dummy, {JMS_min["combine_cor"], JMS_min["combine_cor_fitall"]}, {JMS_ellipse["combine_cor"], JMS_ellipse["combine_cor_fitall"]}, {kRed+2, kBlue+2}, {kOpenTriangleUp, kOpenTriangleDown}, {"only tt", "largest"}, "fit_comparison");

  // void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> leg, TString name)
  // draw({extrema_combine, extrema_combine_cor}, {ellipse_combine, ellipse_combine_cor}, {kBlue+2, kRed+2}, {"uncorrelated", "correlated"}, "correlation_comparison");
  // draw({extrema_combine_cor, extrema_combine_cor_noBKGsys}, {ellipse_combine_cor, ellipse_combine_cor_noBKGsys}, {kBlue+2, kRed+2}, {"with bkg rate", "no bkg rate"}, "BKGsys_cor_comparison");
  // draw({extrema_combine_cor, extrema_2016_comb, extrema_2017_comb, extrema_2018_comb}, {ellipse_combine_cor, ellipse_2016_comb, ellipse_2017_comb, ellipse_2018_comb}, {kBlack, kBlue+2, kRed+2, kGreen+2}, {"Combine", "2016", "2017", "2018"}, "Year_combine_comparison");
  // draw({extrema_combine_lin_cor, extrema_2016_comb_cor, extrema_2017_comb_cor, extrema_2018_comb_cor}, {ellipse_combine_lin_cor, ellipse_2016_comb_cor, ellipse_2017_comb_cor, ellipse_2018_comb_cor}, {kBlack, kBlue+2, kRed+2, kGreen+2}, {"Combine", "2016", "2017", "2018"}, "Year_comparison_cor");
  // draw({extrema_combine_cor, extrema_combine_cor_syst}, {ellipse_combine_cor, ellipse_combine_cor_syst}, {kBlack, kBlue+2}, {"No syst", "syst"}, "syst_comparison_cor");
  // draw({extrema_combine_cor, extrema_combine_lin_cor}, {ellipse_combine_cor, ellipse_combine_lin_cor}, {kBlue+2, kRed+2}, {"with quadratic", "only linear"}, "linear_comparison");
  // draw({extrema_combine, extrema_combine_cor}, {ellipse_combine, ellipse_combine_cor}, {kBlue+2, kRed+2}, {"uncorrelated", "correlated"}, "correlation_comparison");
  // draw({extrema_combine, extrema_muon, extrema_elec}, {ellipse_combine, ellipse_muon, ellipse_elec}, {kBlue+2, kRed+2, kGreen+2}, {"combine", "muon", "elec"}, "channel_comparison");
  // draw({extrema_combine, extrema_add_comb}, {ellipse_combine, ellipse_add_comb}, {kBlue+2, kRed+2}, {"data-bkg", "tt+bkg"}, "bkg_comparison");
  // draw({extrema_combine, extrema_V20_comb}, {ellipse_combine, ellipse_V20_comb}, {kBlue+2, kRed+2}, {"JEC V11", "JEC V20"}, "JECversion_comparison");
  // draw({extrema_add_comb, extrema_add_muon, extrema_add_elec}, {ellipse_add_comb, ellipse_add_muon, ellipse_add_elec}, {kBlue+2, kRed+2, kGreen+2}, {"combine", "muon", "elec"}, "AddBKGChannelComparison");

  // CompareAreas({ellipse_combine_lin_cor, ellipse_2016_comb_cor, ellipse_2017_comb_cor, ellipse_2018_comb_cor}, {"Combine", "2016", "2017", "2018"});

}

// ======================================================================================================
// ===  Compare Areas                                                                                 ===
// ======================================================================================================

void CompareAreas(vector<TEllipse*> ellipse, vector<TString> names){
  if(ellipse.size()!=names.size()) throw runtime_error("ERROR in CompareAreas; Not the same size.");
  vector<double> areas;
  for(unsigned int i=0; i<ellipse.size();i++){
    double area = ellipse[i]->GetR1()*ellipse[i]->GetR2()*TMath::Pi();
    areas.push_back(area);
    cout << "Area for "+names[i]+" is\tA = " << area << endl;
  }
  cout << names[0]+"/"+names[1] << "\t area " << 1/(areas[0]/areas[1]) << "\t lumi " << lumis[0]/lumis[1] << endl;
  cout << names[0]+"/"+names[2] << "\t area " << 1/(areas[0]/areas[2]) << "\t lumi " << lumis[0]/lumis[2] << endl;
  cout << names[0]+"/"+names[3] << "\t area " << 1/(areas[0]/areas[3]) << "\t lumi " << lumis[0]/lumis[3] << endl;
  cout << names[1]+"/"+names[2] << "\t area " << 1/(areas[1]/areas[2]) << "\t lumi " << lumis[1]/lumis[2] << endl;
  cout << names[2]+"/"+names[3] << "\t area " << 1/(areas[2]/areas[3]) << "\t lumi " << lumis[2]/lumis[3] << endl;
  cout << names[1]+"/"+names[3] << "\t area " << 1/(areas[1]/areas[3]) << "\t lumi " << lumis[1]/lumis[3] << endl;

}


// ======================================================================================================
// ===  Extract points                                                                                ===
// ======================================================================================================

MapVD ExtractFromFile(TString path, TString c){
  if(debug) cout << "\t\t ... get file " << path << endl;
  ifstream file(path);
  TString error = path+" Does NOT exists!";
  if(!file.good()) throw runtime_error(error);

  // Go over every line ----------------------------------------------------------------------------
  if(debug) cout << "\t\t ... loop over lines" << endl;
  MapVD JMSfactor;
  string line;
  while (std::getline(file, line))
  {
    TString sline = (TString) line;
    bool nom = line.find("nom") != string::npos;
    bool uu = line.find("uu") != string::npos;
    bool dd = line.find("dd") != string::npos;
    bool ud = line.find("ud") != string::npos;
    bool du = line.find("du") != string::npos;

    if(!(nom||uu||dd||ud||du)) continue;
    // line will be converted to TString
    if(nom) JMSfactor["nom"] = ExtractPoints(line, "nom");
    if(uu) JMSfactor["uu"] = ExtractPoints(line, "uu");
    if(dd) JMSfactor["dd"] = ExtractPoints(line, "dd");
    if(du) JMSfactor["du"] = ExtractPoints(line, "du");
    if(ud) JMSfactor["ud"] = ExtractPoints(line, "ud");
  }

  if(debug){
    cout << fixed << right << "\t" << setw(7) << "factor" << setw(10) << "x" << setw(10) << "y" << " - " << c << endl;
    cout << "\t" << "------------------------------" << endl;
    cout << "\t" << setw(7) << "nom" << setw(10) << JMSfactor["nom"][0] << setw(10) << JMSfactor["nom"][1] <<  endl;
    cout << "\t" << setw(7) << "uu" << setw(10) << JMSfactor["uu"][0] << setw(10) << JMSfactor["uu"][1] <<  endl;
    cout << "\t" << setw(7) << "dd" << setw(10) << JMSfactor["dd"][0] << setw(10) << JMSfactor["dd"][1] <<  endl;
    cout << "\t" << setw(7) << "ud" << setw(10) << JMSfactor["ud"][0] << setw(10) << JMSfactor["ud"][1] <<  endl;
    cout << "\t" << setw(7) << "du" << setw(10) << JMSfactor["du"][0] << setw(10) << JMSfactor["du"][1] <<  endl;
    cout << endl;
  }

  return JMSfactor;
}

VecD ExtractPoints(TString line, const TString jms){
  line.ReplaceAll(jms, ""); // jms_(_x,_y)
  line.ReplaceAll(" ", ""); // _(_x,_y)
  line.ReplaceAll("(", ""); // (x,y)
  line.ReplaceAll(")", ""); // x,y)
  line.ReplaceAll(",", " "); // x,y -> x y
  stringstream number((string) line);
  float x; number >> x;
  float y; number >> y;
  return {x, y};
}

TPolyMarker3D* ExtractFromROOT(TString name, TString year, TString channel, TString sys, TString bins){
  TString dir = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/";
  dir += year+"/"+channel+"/"+bins+"/";
  if(!sys.EqualTo("")) dir += sys+"/";
  TFile *file = new TFile(dir+"Ellipse.root");
  // TPolyMarker3D* ellipse = (TPolyMarker3D*) file->Get("Graphs/JMS_ellipse");
  TPolyMarker3D* graph;
  if(name.Contains("ellipse")) graph = (TPolyMarker3D*) file->Get(name);
  else if(name.Contains("nominal")) graph = (TPolyMarker3D*) file->Get(name);
  return graph;
}

TF2* ExtractFunctionFromROOT(TString year, TString channel, TString sys){
  TString dir = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/";
  dir += year+"/"+channel+"/100/";
  if(sys.EqualTo("")) dir += sys+"/";
  TFile *file = new TFile(dir+"Ellipse.root");
  TF2* func = (TF2*) file->Get("Functions/JMS_Chi2");
  return func;
}

// ======================================================================================================
// ===  Smaller calculations                                                                          ===
// ======================================================================================================


// ------------------------------------------------------------------------------------------------
inline double distance_points(const VecD vec1, const VecD vec2){
  if(vec1.size()!=2&&vec2.size()!=2) throw runtime_error("distance_points: Point-Vector must look like {x,y}");
  return sqrt(pow(vec1[0]-vec2[0],2)+pow(vec1[1]-vec2[1],2));
}

// ------------------------------------------------------------------------------------------------
inline VecD subtract_two_vectors(const VecD vec1, const VecD vec2){
  if(vec1.size()!=2&&vec2.size()!=2) throw runtime_error("subtract_two_vectors: Point-Vector must look like {x,y}");
  return {vec1[0]-vec2[0], vec1[1]-vec2[1]};
}

// ------------------------------------------------------------------------------------------------
inline double length_vector(const VecD vec){
  if(vec.size()!=2) throw runtime_error("length_vector: Point-Vector must look like {x,y}");
  return {sqrt(vec[0]*vec[0]+vec[1]*vec[1])};
}

// ------------------------------------------------------------------------------------------------
inline double angle_to_xaxis(const VecD vec){
  // look at vector xaxis - (1,0)
  if(vec.size()!=2) throw runtime_error("angle_to_xaxis: Point-Vector must look like {x,y}");
  double length = length_vector(vec);
  double scalar = abs(vec[0]);
  return {acos(scalar/length)*(180/TMath::Pi())};
}

// ======================================================================================================
// ===  Build Graphs and ellipse                                                                      ===
// ======================================================================================================

// ------------------------------------------------------------------------------------------------
TGraph* graph_multiple_points(MapVD points){
  TGraph *graph = new TGraph();
  VecD nom = points["nom"];
  graph->SetPoint(0, points["nom"][0], points["nom"][1]);
  graph->SetPoint(1, points["uu"][0], points["uu"][1]);
  graph->SetPoint(2, points["dd"][0], points["dd"][1]);
  graph->SetPoint(3, points["ud"][0], points["ud"][1]);
  graph->SetPoint(4, points["du"][0], points["du"][1]);
  return graph;
}

// ------------------------------------------------------------------------------------------------
TEllipse* build_ellipse(MapVD points, TString name){

  double distance_min_1 = distance_points(points["nom"], points["uu"]);
  double distance_min_2 = distance_points(points["nom"], points["dd"]);
  double distance_max_1 = distance_points(points["nom"], points["ud"]);
  double distance_max_2 = distance_points(points["nom"], points["du"]);

  double dist_min = (distance_min_1+distance_min_2)*0.5;
  double dist_max = (distance_max_1+distance_max_2)*0.5;

  VecD angle_vec = subtract_two_vectors(points["uu"], points["nom"]);
  double theta = angle_to_xaxis(angle_vec);

  // First: Halbachse alonge xaxis, Second: Halbachse along yaxis
  TEllipse *ellipse = new TEllipse(points["nom"][0], points["nom"][1], dist_min, dist_max, 0, 360, theta);
  return ellipse;
}

// ======================================================================================================
// ===  Draw Settings                                                                                 ===
// ======================================================================================================

// ------------------------------------------------------------------------------------------------
inline void set_graph(TGraph *graph, int color){
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(0.5);
  graph->GetXaxis()->SetLimits(-3, 3);
  graph->GetYaxis()->SetRangeUser(-3, 3);
  graph->GetHistogram()->GetXaxis()->SetTitle("JEC");
  graph->GetHistogram()->GetYaxis()->SetTitle("XCone");
}

inline void set_graph(TPolyMarker3D *graph, int color, int style, double size){
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(style);
  graph->SetMarkerSize(size);
}
// ------------------------------------------------------------------------------------------------
inline void set_ellipse(TEllipse *ellipse, int color){
  ellipse->SetLineColor(color);
  ellipse->SetFillStyle(0);
}

// ------------------------------------------------------------------------------------------------
void draw2(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> names, TString name){
  TString save_path = get_save_path();

  for(unsigned int i=0; i<color.size(); i++){
    set_graph(graph[i], color[i]);
    set_ellipse(ellipse[i], color[i]);
  }

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  graph[0]->GetXaxis()->SetTitle("JEC factor");
  graph[0]->GetYaxis()->SetTitle("Additional XCone correction factor");

  TCanvas *A = new TCanvas(name, "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  graph[0]->Draw("AP");
  ellipse[0]->Draw("same");
  for(unsigned int i=1; i<color.size(); i++){
    graph[i]->Draw("same P");
    ellipse[i]->Draw("same");
  }
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetTextSize(0.03);
  for(unsigned int i=0; i<color.size(); i++) leg->AddEntry(ellipse[i], names[i],"l");
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path+"/JetCorrections/Ellipse/"+name+".pdf");
}


// ------------------------------------------------------------------------------------------------
void draw(TF2* dummy, vector<TPolyMarker3D*> min, vector<TPolyMarker3D*> graph, vector<int> color, vector<int> style,  vector<TString> names, TString name){
  cout << "\t ... " << name << endl;

  TString save_path = get_save_path();

  // TH2F* dummy = new TH2F("dummy","dummy", 100, -1.5, 1.5, 100, -1.5, 1.5);

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  dummy->GetXaxis()->SetTitle("JEC factor");
  dummy->GetYaxis()->SetTitle("Additional XCone correction factor");

  SetupCanvas(false);
  m_can->cd();
  m_can->SetLeftMargin(0.15);
  m_can->SetBottomMargin(0.15);
  m_can->SetRightMargin(0.05);

  Int_t nb = 50;
  dummy->SetContour(nb); // Contours
  dummy->SetTitle("");
  dummy->SetFillStyle(0);
  dummy->SetLineWidth(4);
  dummy->SetRange(-1.7, -1.7, 1.7, 1.7);
  // if(name.Contains("AllBins")) dummy->SetRange(-4, -4, 4, 4);

  const Int_t Number = 2;
  // Double_t Red[Number]    = { 0.99, 0.49, 0.10, 0.00};
  // Double_t Green[Number]  = { 0.99, 0.80, 0.40, 0.00};
  // Double_t Blue[Number]   = { 0.99, 0.88, 0.66, 0.33};
  // Double_t Length[Number] = { 0.00, 0.33, 0.66, 1.00};
  Double_t Red[Number]    = { 1.00, 1.00};
  Double_t Green[Number]  = { 1.00, 1.00};
  Double_t Blue[Number]   = { 1.00, 1.00};
  Double_t Length[Number] = { 0.00, 1.00};

  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  dummy->Draw("cont4z");
  Cosmetics(dummy->GetHistogram(), "#it{f}^{ JEC}", 0, 0, 0, 0, false);
  dummy->GetHistogram()->SetFillColorAlpha(kWhite, 0);
  dummy->GetHistogram()->GetYaxis()->SetTitle("#it{f}^{ XCone}");
  dummy->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  dummy->GetHistogram()->GetZaxis()->SetTitleOffset(1.1);
  dummy->GetHistogram()->GetZaxis()->SetTitleSize(0.06);
  dummy->GetHistogram()->GetZaxis()->SetLabelSize(0.04);
  dummy->GetHistogram()->GetZaxis()->SetTickLength(0.02);
  dummy->GetHistogram()->GetZaxis()->SetLabelOffset(0.011);
  dummy->GetHistogram()->GetZaxis()->SetTitle("#chi^{ 2}");
  dummy->GetHistogram()->GetZaxis()->CenterTitle();

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)dummy->GetHistogram()->GetListOfFunctions()->FindObject("palette");

  // the following lines move the palette. Choose the values you need for the position.
  palette->SetX1NDC(1.5);
  palette->SetX2NDC(1.51);
  palette->SetY1NDC(1.5);
  palette->SetY2NDC(1.51);

  gPad->Modified();
  gPad->Update();

  // dummy->SetMinimum(190);
  dummy->SetMinimum(0);
  gPad->RedrawAxis();
  m_can->SetTheta(90);
  m_can->SetPhi(0);

  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetTextSize(0.03);
  for(unsigned int i=0; i<color.size(); i++){
    TH1F* filler = new TH1F();
    filler->SetMarkerColor(color[i]);
    filler->SetMarkerStyle(kFullCircle);
    leg->AddEntry(filler, names[i],"p");
    if(graph.size()>0){
      set_graph(graph[i], color[i], style[i],0.05);
      graph[i]->Draw("same P");
    }
    if(min.size()>0){
      set_graph(min[i], color[i], style[i],0.4);
      min[i]->Draw("same P");
    }
  }
  DrawLumi(138., false, true, false, "JMS");
  leg->Draw();
  // gPad->RedrawAxis();
  m_can->SaveAs(save+name+".pdf");
}
