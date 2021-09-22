#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

#include <fstream>
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include <TMath.h>
#include <TString.h>

// #include <io.h>
using namespace std;

inline void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> names, TString name);
inline void set_ellipse(TEllipse *ellipse, int color);
inline void set_graph(TGraph *graph, int color);
MapVD ExtractFromFile(TString path, TString c);
inline double angle_to_xaxis(const VecD vec);
inline double distance_points(const VecD vec1, const VecD vec2);
inline double length_vector(const VecD vec);
inline VecD ExtractPoints(TString line, const TString jms);
inline VecD subtract_two_vectors(const VecD vec1, const VecD vec2);
inline TEllipse* build_ellipse(MapVD points);
inline TGraph* graph_multiple_points(MapVD points);

bool debug = false;

int main(int argc, char* argv[]){

  gErrorIgnoreLevel = kWarning;
  TString save_path = get_save_path();

  TString dir = save_path+"/JetCorrections/fit/combine/";
  VecTS channels = {"combine", "muon", "elec"};
  TString width = "/100/BinWidth_1/";
  TString file_txt = "JMS.txt";

  MapVDD JMSfactor;

  cout << "\nExtract Points ... " << endl;
  JMSfactor["combine_add"] = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/combine/100/BinWidth_1/JMS.txt", "combine");
  JMSfactor["muon_add"] = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/muon/100/BinWidth_1/JMS.txt", "muon");
  JMSfactor["elec_add"] = ExtractFromFile(save_path+"/JetCorrections/fit_t+b/combine/elec/100/BinWidth_1/JMS.txt", "elec");
  JMSfactor["combine_V20"] = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/combine/100/BinWidth_1/JMS.txt", "elec");
  JMSfactor["muon_V20"] = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/muon/100/BinWidth_1/JMS.txt", "elec");
  JMSfactor["elec_V20"] = ExtractFromFile(save_path+"/JetCorrections/fit_V20/combine/elec/100/BinWidth_1/JMS.txt", "elec");

  for(TString c: channels){
    if(debug) cout << "\t ... channel "+c << endl;
    TString file = dir+c+width+file_txt;
    JMSfactor[c] = ExtractFromFile(file, c);
  }

  // Fill Graph and Ellipse ------------------------------------------------------------------------
  cout << "Create Graphs ... " << endl;
  TGraph* extrema_combine = graph_multiple_points(JMSfactor["combine"]);
  TGraph* extrema_muon = graph_multiple_points(JMSfactor["muon"]);
  TGraph* extrema_elec = graph_multiple_points(JMSfactor["elec"]);
  TGraph* extrema_add_comb = graph_multiple_points(JMSfactor["combine_add"]);
  TGraph* extrema_add_muon = graph_multiple_points(JMSfactor["muon_add"]);
  TGraph* extrema_add_elec = graph_multiple_points(JMSfactor["elec_add"]);
  TGraph* extrema_V20_comb = graph_multiple_points(JMSfactor["combine_V20"]);
  TGraph* extrema_V20_muon = graph_multiple_points(JMSfactor["muon_V20"]);
  TGraph* extrema_V20_elec = graph_multiple_points(JMSfactor["elec_V20"]);

  cout << "Create Ellipses ... " << endl;
  TEllipse* ellipse_combine = build_ellipse(JMSfactor["combine"]);
  TEllipse* ellipse_muon = build_ellipse(JMSfactor["muon"]);
  TEllipse* ellipse_elec = build_ellipse(JMSfactor["elec"]);
  TEllipse* ellipse_add_comb = build_ellipse(JMSfactor["combine_add"]);
  TEllipse* ellipse_add_muon = build_ellipse(JMSfactor["muon_add"]);
  TEllipse* ellipse_add_elec = build_ellipse(JMSfactor["elec_add"]);
  TEllipse* ellipse_V20_comb = build_ellipse(JMSfactor["combine_V20"]);
  TEllipse* ellipse_V20_muon = build_ellipse(JMSfactor["muon_V20"]);
  TEllipse* ellipse_V20_elec = build_ellipse(JMSfactor["elec_V20"]);

  // Start drawing ---------------------------------------------------------------------------------
  cout << "Start plotting ... " << endl;
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  // void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> leg, TString name)
  draw({extrema_combine, extrema_muon, extrema_elec}, {ellipse_combine, ellipse_muon, ellipse_elec}, {kBlue+2, kRed+2, kGreen+2}, {"combine", "muon", "elec"}, "channel_comparison");
  draw({extrema_combine, extrema_add_comb}, {ellipse_combine, ellipse_add_comb}, {kBlue+2, kRed+2}, {"data-bkg", "tt+bkg"}, "bkg_comparison");
  draw({extrema_combine, extrema_V20_comb}, {ellipse_combine, ellipse_V20_comb}, {kBlue+2, kRed+2}, {"JEC V11", "JEC V20"}, "JECversion_comparison");
  draw({extrema_add_comb, extrema_add_muon, extrema_add_elec}, {ellipse_add_comb, ellipse_add_muon, ellipse_add_elec}, {kBlue+2, kRed+2, kGreen+2}, {"combine", "muon", "elec"}, "AddBKGChannelComparison");
}


// ======================================================================================================
// ===  Extract points                                                                                ===
// ======================================================================================================

MapVD ExtractFromFile(TString path, TString c){
  if(debug) cout << "\t\t ... get file" << endl;
  cout << path << endl;
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

inline VecD ExtractPoints(TString line, const TString jms){
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
inline TGraph* graph_multiple_points(MapVD points){
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
inline TEllipse* build_ellipse(MapVD points){

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

// ------------------------------------------------------------------------------------------------
inline void set_ellipse(TEllipse *ellipse, int color){
  ellipse->SetLineColor(color);
  ellipse->SetFillStyle(0);
}

// ------------------------------------------------------------------------------------------------
inline void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> names, TString name){
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
