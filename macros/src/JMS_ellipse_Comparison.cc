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

/*
. ██████ ██    ██ ████████
.██      ██    ██    ██
.██      ██    ██    ██
.██      ██    ██    ██
. ██████  ██████     ██
*/

inline vector<double> cut_points(TString line, const TString jms, const TString seperater){
  line.ReplaceAll(jms+": ( " , "");
  line.ReplaceAll(")", "");
  line.ReplaceAll(" ", "");
  int index_sep = line.Index(seperater);
  TString x_jms = line;
  TString y_jms = line;
  x_jms.Replace(index_sep, line.Length(), "");
  y_jms.Replace(0, index_sep+1, "");
  return {x_jms.Atof(), y_jms.Atof()};
}

inline TString cut_points(TString line, TString coord){
  line.ReplaceAll(" ", "");
  line.ReplaceAll(":", "");
  line.ReplaceAll(coord, "");
  return line;
}

/*
.██████   █████  ███████ ██  ██████ ███████
.██   ██ ██   ██ ██      ██ ██      ██
.██████  ███████ ███████ ██ ██      ███████
.██   ██ ██   ██      ██ ██ ██           ██
.██████  ██   ██ ███████ ██  ██████ ███████
*/

// ------------------------------------------------------------------------------------------------
inline double distance_points(const vector<double> vec1, const vector<double> vec2){
  if(vec1.size()!=2&&vec2.size()!=2) throw runtime_error("distance_points: Point-Vector must look like {x,y}");
  return sqrt(pow(vec1[0]-vec2[0],2)+pow(vec1[1]-vec2[1],2));
}

// ------------------------------------------------------------------------------------------------
inline vector<double> subtract_two_vectors(const vector<double> vec1, const vector<double> vec2){
  if(vec1.size()!=2&&vec2.size()!=2) throw runtime_error("subtract_two_vectors: Point-Vector must look like {x,y}");
  return {vec1[0]-vec2[0], vec1[1]-vec2[1]};
}

// ------------------------------------------------------------------------------------------------
inline double length_vector(const vector<double> vec){
  if(vec.size()!=2) throw runtime_error("length_vector: Point-Vector must look like {x,y}");
  return {sqrt(vec[0]*vec[0]+vec[1]*vec[1])};
}

// ------------------------------------------------------------------------------------------------
inline double angle_to_xaxis(const vector<double> vec){
  // look at vector xaxis - (1,0)
  if(vec.size()!=2) throw runtime_error("angle_to_xaxis: Point-Vector must look like {x,y}");
  double length = length_vector(vec);
  double scalar = abs(vec[0]);
  return {acos(scalar/length)*(180/TMath::Pi())};
}

/*
.██████  ██    ██ ██ ██      ██████
.██   ██ ██    ██ ██ ██      ██   ██
.██████  ██    ██ ██ ██      ██   ██
.██   ██ ██    ██ ██ ██      ██   ██
.██████   ██████  ██ ███████ ██████
*/

// ------------------------------------------------------------------------------------------------
inline TGraph* graph_one_point(const vector<double> point){
  if(point.size()!=2) throw runtime_error("graph_one_point: Point-Vector must look like {x,y}");
  TGraph *graph = new TGraph();
  graph->SetPoint(0, point[0], point[1]);
  return graph;
}

// ------------------------------------------------------------------------------------------------
inline TGraph* graph_multiple_points(const vector<vector<double>> points, int number){
  if(points.size()<number) throw runtime_error("graph_multiple_points: Vector with Points too small");
  TGraph *graph = new TGraph();
  for(int i=0; i<number; i++){
    if(points[i].size()!=2) throw runtime_error("graph_multiple_points: Point-Vector must look like {x,y} - Point: "+to_string(i));
    graph->SetPoint(i, points[i][0], points[i][1]);
  }
  return graph;
}

// ------------------------------------------------------------------------------------------------
inline TEllipse* build_ellipse(const vector<vector<double>> points){
  // 1st entry of vector<vector> - mid point;
  // 2nd entry of vector<vector> - uu  point;
  // 3rd entry of vector<vector> - dd  point;
  // 4th entry of vector<vector> - ud  point;
  // 5th entry of vector<vector> - du  point;

  double distance_min_1 = distance_points(points[0], points[1]);
  double distance_min_2 = distance_points(points[0], points[2]);
  double distance_max_1 = distance_points(points[0], points[3]);
  double distance_max_2 = distance_points(points[0], points[4]);

  double dist_min = (distance_min_1+distance_min_2)*0.5;
  double dist_max = (distance_max_1+distance_max_2)*0.5;

  vector<double> angle_vec = subtract_two_vectors(points[1], points[0]);
  double theta = angle_to_xaxis(angle_vec);

  // First: Halbachse alonge xaxis, Second: Halbachse along yaxis
  TEllipse *ellipse = new TEllipse(points[0][0], points[0][1], dist_min, dist_max, 0, 360, theta);
  return ellipse;
}

/*
.███████ ███████ ████████
.██      ██         ██
.███████ █████      ██
.     ██ ██         ██
.███████ ███████    ██
*/

// ------------------------------------------------------------------------------------------------
inline void set_graph(TGraph *graph, int color){
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(0.5);
}

// ------------------------------------------------------------------------------------------------
inline void set_main_graph(TGraph *graph, int color){
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(0.5);
  graph->GetXaxis()->SetLimits(-2, 2);
  graph->GetYaxis()->SetRangeUser(-2, 2);
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
  for(unsigned int i=0; i<color.size(); i++){
    if(i==0) set_main_graph(graph[i], color[i]);
    else     set_graph(graph[i], color[i]);
    set_ellipse(ellipse[i], color[i]);
  }

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

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
  leg->SetTextSize(0.02);
  for(unsigned int i=0; i<color.size(); i++) leg->AddEntry(ellipse[i], names[i],"l");
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/ellipse/"+name+".pdf");
}

/*
.███    ███  █████  ██ ███    ██
.████  ████ ██   ██ ██ ████   ██
.██ ████ ██ ███████ ██ ██ ██  ██
.██  ██  ██ ██   ██ ██ ██  ██ ██
.██      ██ ██   ██ ██ ██   ████
*/

int main(int argc, char* argv[]){

  bool debug = true;
  bool show_points = true;
  if(debug) show_points=true;

  // /*
  // ███    ███ ████████  ██████  ██████
  // ████  ████    ██    ██    ██ ██   ██
  // ██ ████ ██    ██    ██    ██ ██████
  // ██  ██  ██    ██    ██    ██ ██
  // ██      ██    ██     ██████  ██
  // */
  //
  // // #################################################################################################
  // // Get Path ########################################################################################
  //
  // TString dir = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/";
  // vector<TString> uncerts = {"only_data_uncert/", "mc_data_uncert/", "mc_factor_data_uncert/"};
  // TString ptbins = "pt_bins/";
  // vector<TString> masses = {"combined/", "1695/combined/", "1755/combined/"};
  // TString subdirs = "btag/rebin45/";
  // vector<TString> fits = {"masspeak/", "linear/masspeak/"};
  // TString file_txt = "jec_factor_all.txt";
  //
  // // #################################################################################################
  // // Get Points ######################################################################################
  //
  // vector<vector<double>> nominal_data, nominal_mcdata, nominal_data_lin, nominal_mcdata_lin;
  // vector<vector<double>> min_data, min_mcdata, min_data_lin, min_mcdata_lin, min_factor, min_factor_lin; // 1695
  // vector<vector<double>> max_data, max_mcdata, max_data_lin, max_mcdata_lin, max_factor, max_factor_lin; // 1755
  //
  // // Order: data, data_lin, mcdata, mcdata_lin, factor, factor_lin
  // vector<TGraph*>   graphs_mid_nominal, graphs_mid_min, graphs_mid_max;
  // vector<TGraph*>   graphs_all_nominal, graphs_all_min, graphs_all_max;
  // vector<TEllipse*> ellipse_nominal,    ellipse_min,    ellipse_max;
  //
  // vector<TString> points_name = {"min", "uu", "dd", "ud", "du"};
  //
  // for(TString mass: masses){
  //   for(TString unc: uncerts){
  //     for(TString fit: fits){
  //
  //       // Get File --------------------------------------------------------------------------------------
  //       TString file_points = dir+unc+ptbins+mass+subdirs+fit+file_txt;
  //       ifstream file(file_points);
  //       if(!file.good()) continue;
  //
  //       if(show_points){
  //         cout << "\n#########################################################################################\n";
  //         cout << file_points << endl;
  //       }
  //
  //       // Get iterarotes --------------------------------------------------------------------------------
  //       // uncerts = 0:only_data_uncert   1:mc_data_uncert   2:mc_factor_data_uncert
  //       // masses  = 0:combined           1:1695             2:1755
  //       // fits    = 0:linear             1:all
  //       std::vector<TString>::iterator it_unc_  = std::find(uncerts.begin(), uncerts.end(), unc);
  //       std::vector<TString>::iterator it_mass_ = std::find(masses.begin(), masses.end(), mass);
  //       std::vector<TString>::iterator it_fit_  = std::find(fits.begin(), fits.end(), fit);
  //
  //       int it_unc  = std::distance(uncerts.begin(), it_unc_);
  //       int it_mass = std::distance(masses.begin(), it_mass_);
  //       int it_fit  = std::distance(fits.begin(), it_fit_);
  //       if(show_points) cout << "unc: " << it_unc << " | mass: " << it_mass << " | fit: " << it_fit << endl;
  //       // if(it_unc==0&&it_mass==0&&it_fit==0) continue;
  //
  //       // Go over every line ----------------------------------------------------------------------------
  //       int i = 0;
  //       string line;
  //       vector<double> points_min, points_uu, points_dd, points_ud, points_du;
  //       while (std::getline(file, line))
  //       {
  //         i++;
  //         bool skip=false;
  //         TString line_new = (TString) line;
  //         vector<TString> points;
  //
  //         if(i==1||i==2||i==3||i==6) continue; // Bins, blanck, zmin, blanck
  //         if(i==4)   line_new    = cut_points(line_new, "x"); // xmin
  //         if(i==5)   line_new    = cut_points(line_new, "y"); // ymin
  //         if(i==7)   points_uu   = cut_points(line_new, "uu", ",");
  //         if(i==8)   points_dd   = cut_points(line_new, "dd", ",");
  //         if(i==9)   points_ud   = cut_points(line_new, "ud", ",");
  //         if(i==10)  points_du   = cut_points(line_new, "du", ",");
  //
  //         if(i==4) points_min.push_back(Round(line_new.Atof(), 3));
  //         if(i==5) points_min.push_back(Round(line_new.Atof(), 3));
  //       }
  //
  //       cout << endl;
  //       if(show_points){
  //         cout << setw(11) << "x" << setw(8) << "y" << endl;
  //         cout << setw(20) << "--------------------";
  //         cout << "\n" << setw(3) << points_name[0] << ": " ;
  //         for(unsigned int i=0; i<points_min.size(); i++) cout << setw(6) << points_min[i] << "  ";
  //         cout << "\n" << setw(3) << points_name[1] << ": " ;
  //         for(unsigned int i=0; i<points_uu.size(); i++) cout << setw(6) << points_uu[i] << "  ";
  //         cout << "\n" << setw(3) << points_name[2] << ": " ;
  //         for(unsigned int i=0; i<points_dd.size(); i++) cout << setw(6) << points_dd[i] << "  ";
  //         cout << "\n" << setw(3) << points_name[3] << ": " ;
  //         for(unsigned int i=0; i<points_ud.size(); i++) cout << setw(6) << points_ud[i] << "  ";
  //         cout << "\n" << setw(3) << points_name[4] << ": " ;
  //         for(unsigned int i=0; i<points_du.size(); i++) cout << setw(6) << points_du[i] << "  ";
  //         cout << endl;
  //       }
  //       // Fill Graph and Ellipse ------------------------------------------------------------------------
  //       vector<vector<double>> all_points = {points_min, points_uu, points_dd, points_du, points_ud};
  //       TGraph*   graph_one = graph_one_point(all_points[0]);
  //       TGraph*   graph_all = graph_multiple_points(all_points, 5);
  //       TEllipse* ellipse   = build_ellipse(all_points);
  //
  //       if(it_mass==0){
  //         graphs_mid_nominal.push_back(graph_one);
  //         graphs_all_nominal.push_back(graph_all);
  //         ellipse_nominal.push_back(ellipse);
  //       } else if(it_mass==1){
  //         graphs_mid_min.push_back(graph_one);
  //         graphs_all_min.push_back(graph_all);
  //         ellipse_min.push_back(ellipse);
  //       } else if(it_mass==2){
  //         graphs_mid_max.push_back(graph_one);
  //         graphs_all_max.push_back(graph_all);
  //         ellipse_max.push_back(ellipse);
  //       }
  //     }
  //   }
  // }
  //
  // if(debug){
  //   cout << endl;
  //   cout << graphs_mid_nominal.size() << "   " << graphs_all_nominal.size() << "   " << ellipse_nominal.size() << endl;
  //   cout << graphs_mid_min.size()     << "   " << graphs_all_min.size()     << "   " << ellipse_min.size() << endl;
  //   cout << graphs_mid_max.size()     << "   " << graphs_all_max.size()     << "   " << ellipse_max.size() << endl;
  // }
  // cout << endl;
  //
  // // #################################################################################################
  // // #################################################################################################
  // // Order: data, data_lin, mcdata, mcdata_lin, factor, factor_lin
  // //           0         1       2           3       4           5
  // gStyle->SetPadTickY(1);
  // gStyle->SetPadTickX(1);
  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetLegendBorderSize(0);
  //
  // // void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> leg, TString name)
  // vector<TGraph*>   graphs;
  // vector<TEllipse*> ellipse;
  // vector<int>       colors;
  // vector<TString>   names;

  // #################################################################################################
  // #################################################################################################
  // // Nominal, data, fits -----------------------------------------------------------------------------
  // if(debug) cout << "Combined - data fits\n";
  // graphs  = {graphs_all_nominal[0], graphs_all_nominal[1]};
  // ellipse = {ellipse_nominal[0], ellipse_nominal[1]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_nom_data_fits");

  // // Nominal, mcdata, fits -----------------------------------------------------------------------------
  // if(debug) cout << "Combined - mcdata fits\n";
  // graphs  = {graphs_all_nominal[2], graphs_all_nominal[3]};
  // ellipse = {ellipse_nominal[2], ellipse_nominal[3]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_nom_mcdata_fits");

  // Nominal, uncerts, all ----------------------------------------------------------------------------
  // if(debug) cout << "Combined - uncerts\n";
  // graphs  = {graphs_all_nominal[0], graphs_all_nominal[2]};
  // ellipse = {ellipse_nominal[0], ellipse_nominal[2]};
  // colors  = {kBlue, kRed};
  // names   = {"DATA", "MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_nom_uncerts");

  // // Nominal, uncerts, lin ----------------------------------------------------------------------------
  // if(debug) cout << "Combined - uncerts lin\n";
  // graphs  = {graphs_all_nominal[1], graphs_all_nominal[3]};
  // ellipse = {ellipse_nominal[1], ellipse_nominal[3]};
  // colors  = {kBlue, kRed};
  // names   = {"DATA", "MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_nom_uncerts_lin");

  // #################################################################################################
  // #################################################################################################
  if(debug) cout << "1695\n";

  // // 1695, data, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_min[0], graphs_all_min[1]};
  // ellipse = {ellipse_min[0], ellipse_min[1]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1695_data_fits");

  // // 1695, mcdata, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_min[2], graphs_all_min[3]};
  // ellipse = {ellipse_min[2], ellipse_min[3]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1695_mcdata_fits");
  //
  // // 1695, factor, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_min[4], graphs_all_min[5]};
  // ellipse = {ellipse_min[4], ellipse_min[5]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1695_factor_fits");
  //
  // // 1695, uncerts, lin ----------------------------------------------------------------------------
  // graphs  = {graphs_all_min[1], graphs_all_min[3], graphs_all_min[5]};
  // ellipse = {ellipse_min[1], ellipse_min[3], ellipse_min[5]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"DATA", "MC+DATA", "factor*MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_1695_uncerts_lin");

  // // 1695, uncerts, all ----------------------------------------------------------------------------
  // graphs  = {graphs_all_min[0], graphs_all_min[2], graphs_all_min[4]};
  // ellipse = {ellipse_min[0], ellipse_min[2], ellipse_min[4]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"DATA", "MC+DATA", "factor*MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_1695_uncerts");

  // #################################################################################################
  // #################################################################################################
  if(debug) cout << "1755\n";

  // // 1755, data, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_max[0], graphs_all_max[1]};
  // ellipse = {ellipse_max[0], ellipse_max[1]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1755_data_fits");
  //
  // // 1755, mcdata, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_max[2], graphs_all_max[3]};
  // ellipse = {ellipse_max[2], ellipse_max[3]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1755_mcdata_fits");
  //
  // // 1755, factor, fits -----------------------------------------------------------------------------
  // graphs  = {graphs_all_max[4], graphs_all_max[5]};
  // ellipse = {ellipse_max[4], ellipse_max[5]};
  // colors  = {kBlue, kRed};
  // names   = {"all fits", "only lin fits"};
  // draw(graphs, ellipse, colors, names, "compare_1755_factor_fits");
  //
  // // 1755, uncerts, lin ----------------------------------------------------------------------------
  // graphs  = {graphs_all_max[1], graphs_all_max[3], graphs_all_max[5]};
  // ellipse = {ellipse_max[1], ellipse_max[3], ellipse_max[5]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"DATA", "MC+DATA", "factor*MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_1755_uncerts_lin");

  // 1755, uncerts, all ----------------------------------------------------------------------------
  // graphs  = {graphs_all_max[0], graphs_all_max[2], graphs_all_max[4]};
  // ellipse = {ellipse_max[0], ellipse_max[2], ellipse_max[4]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"DATA", "MC+DATA", "factor*MC+DATA"};
  // draw(graphs, ellipse, colors, names, "compare_1755_uncerts");

  // #################################################################################################
  // #################################################################################################
  // // nom, 1695, 1755, data, all ----------------------------------------------------------------------
  // if(debug) cout << "All\n";
  // graphs  = {graphs_all_nominal[0], graphs_all_min[0], graphs_all_max[0]};
  // ellipse = {ellipse_nominal[0], ellipse_min[0], ellipse_max[0]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_data");

  // // nom, 1695, 1755, data, lin ----------------------------------------------------------------------
  // graphs  = {graphs_all_nominal[1], graphs_all_min[1], graphs_all_max[1]};
  // ellipse = {ellipse_nominal[1], ellipse_min[1], ellipse_max[1]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_data_lin");

  // nom, 1695, 1755, mc+data, all ----------------------------------------------------------------------
  // graphs  = {graphs_all_nominal[0], graphs_all_min[0], graphs_all_max[0]};
  // ellipse = {ellipse_nominal[0], ellipse_min[0], ellipse_max[0]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_mcdata_45");

  // // nom, 1695, 1755, mc+data, lin ----------------------------------------------------------------------
  // graphs  = {graphs_all_nominal[3], graphs_all_min[3], graphs_all_max[3]};
  // ellipse = {ellipse_nominal[3], ellipse_min[3], ellipse_max[3]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_mcdata_lin");

  // // nom, 1695, 1755, factor, all ----------------------------------------------------------------------
  // graphs  = {graphs_all_nominal[2], graphs_all_min[4], graphs_all_max[4]};
  // ellipse = {ellipse_nominal[2], ellipse_min[4], ellipse_max[4]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_factor");

  // // nom, 1695, 1755, factor, lin ----------------------------------------------------------------------
  // graphs  = {graphs_all_nominal[3], graphs_all_min[5], graphs_all_max[5]};
  // ellipse = {ellipse_nominal[3], ellipse_min[5], ellipse_max[5]};
  // colors  = {kBlue, kRed, kGreen};
  // names   = {"m_{t}=172.5 GeV", "m_{t}=169.5 GeV", "m_{t}=175.5 GeV"};
  // draw(graphs, ellipse, colors, names, "compare_mass_factor_lin");

  /*
  ██████  ██████  ███    ███ ██████  ██
  ██      ██    ██ ████  ████ ██   ██ ██
  ██      ██    ██ ██ ████ ██ ██████  ██
  ██      ██    ██ ██  ██  ██ ██   ██ ██
  ██████  ██████  ██      ██ ██████  ██
  */

  // #################################################################################################
  // Get Path ########################################################################################

  TString dir          = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/mc_data_uncert/pt_bins/combined/btag/rebin180/";
  TString file_txt     = "jec_factor_all.txt";
  vector<TString> fits = {"masspeak/", "mtop_combi/masspeak/"};

  // #################################################################################################
  // Get Points ######################################################################################

  // Order: data, data_lin, mcdata, mcdata_lin, factor, factor_lin
  vector<TGraph*>   graphs_all_single, graphs_all_added;
  vector<TEllipse*> ellipse_single   , ellipse_added   ;

  vector<TString> points_name = {"min", "uu", "dd", "ud", "du"};

  for(TString fit: fits){

    // Get File --------------------------------------------------------------------------------------
    TString file_points = dir+fit+file_txt;
    ifstream file(file_points);
    if(!file.good()) continue;

    if(show_points){
      cout << "\n#########################################################################################\n";
      cout << file_points << endl;
    }

    // Get iterarotes --------------------------------------------------------------------------------
    // fits    = 0: nominal             1: all samples

    std::vector<TString>::iterator it_fit_  = std::find(fits.begin(), fits.end(), fit);
    int it_fit  = std::distance(fits.begin(), it_fit_);
    if(show_points) cout << "fit: " << it_fit << endl;
    // if(it_unc==0&&it_mass==0&&it_fit==0) continue;

    // Go over every line ----------------------------------------------------------------------------
    int i = 0;
    string line;
    vector<double> points_min, points_uu, points_dd, points_ud, points_du;
    while (std::getline(file, line))
    {
      i++;
      bool skip=false;
      TString line_new = (TString) line;

      if(i==1||i==2||i==3||i==6) continue; // Bins, blanck, zmin, blanck
      if(i==4)   line_new    = cut_points(line_new, "x"); // xmin
      if(i==5)   line_new    = cut_points(line_new, "y"); // ymin
      if(i==7)   points_uu   = cut_points(line_new, "uu", ",");
      if(i==8)   points_dd   = cut_points(line_new, "dd", ",");
      if(i==9)   points_ud   = cut_points(line_new, "ud", ",");
      if(i==10)  points_du   = cut_points(line_new, "du", ",");

      if(i==4) points_min.push_back(Round(line_new.Atof(), 3));
      if(i==5) points_min.push_back(Round(line_new.Atof(), 3));
    }

    cout << endl;
    if(show_points){
      cout << setw(11) << "x" << setw(8) << "y" << endl;
      cout << setw(20) << "--------------------";
      cout << "\n" << setw(3) << points_name[0] << ": " ;
      for(unsigned int i=0; i<points_min.size(); i++) cout << setw(6) << points_min[i] << "  ";
      cout << "\n" << setw(3) << points_name[1] << ": " ;
      for(unsigned int i=0; i<points_uu.size(); i++) cout << setw(6) << points_uu[i] << "  ";
      cout << "\n" << setw(3) << points_name[2] << ": " ;
      for(unsigned int i=0; i<points_dd.size(); i++) cout << setw(6) << points_dd[i] << "  ";
      cout << "\n" << setw(3) << points_name[3] << ": " ;
      for(unsigned int i=0; i<points_ud.size(); i++) cout << setw(6) << points_ud[i] << "  ";
      cout << "\n" << setw(3) << points_name[4] << ": " ;
      for(unsigned int i=0; i<points_du.size(); i++) cout << setw(6) << points_du[i] << "  ";
      cout << endl;
    }
    // Fill Graph and Ellipse ------------------------------------------------------------------------
    vector<vector<double>> all_points = {points_min, points_uu, points_dd, points_du, points_ud};
    TGraph*   graph_all = graph_multiple_points(all_points, 5);
    TEllipse* ellipse   = build_ellipse(all_points);

    if(it_fit==0){
      graphs_all_single.push_back(graph_all);
      ellipse_single.push_back(ellipse);
    } else if (it_fit==1){
      graphs_all_added.push_back(graph_all);
      ellipse_added.push_back(ellipse);
    }
  }


  if(debug){
    cout << endl;
    cout << graphs_all_single.size() << "   " << ellipse_single.size() << endl;
    cout << graphs_all_added.size()  << "   " << ellipse_added.size()  << endl;
  }
  cout << endl;

  // #################################################################################################
  // #################################################################################################
  // Order: data, data_lin, mcdata, mcdata_lin, factor, factor_lin
  //           0         1       2           3       4           5
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  // void draw(vector<TGraph*> graph, vector<TEllipse*> ellipse, vector<int> color,  vector<TString> leg, TString name)
  vector<TGraph*>   graphs;
  vector<TEllipse*> ellipse;
  vector<int>       colors;
  vector<TString>   names;

  // compare single + added --------------------------------------------------------------------------
  graphs  = {graphs_all_single[0], graphs_all_added[0]};
  ellipse = {ellipse_single[0], ellipse_added[0]};
  colors  = {kBlue, kRed};
  names   = {"m_{t}=172.5 GeV", "m_{t}=added"};
  draw(graphs, ellipse, colors, names, "compare_single_added");
}
