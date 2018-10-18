#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TVectorD.h>
#include <TF1.h>
#include <vector>
#include <TLatex.h>

class plotter{

 public:
  plotter(TString);
  void draw_matrix(TH2* hist, TString file_name, bool zlog);
  void draw_output(TH1* output, TH1D* truth, bool norm, TString file_name);
  void draw_output_smear(std::vector<TH1*> output, TH1D* truth, TString file_name);
  void draw_output_stat(TH1* output_, TH1* stat_, TH1D* truth_, bool norm, TString file_name);
  void draw_output_mass(TH1* output, std::vector<TH1D*> mtop_templates, std::vector<bool> show, bool norm, TString file_name);
  void draw_output_pseudo(TH1* output, TH1D* pseudotruth, TH1D* mctruth, bool norm, TString file_name);
  void draw_lcurve(TGraph *lcurve, double x1, double y1, TString file_name);
  void draw_projection(TH1D* proj, TH1D* compare, TString file_name );
  void draw_1D_hist(TH1D* hist, TString file_name);
  void draw_delta(TH1* hist, TString file_name);
  void draw_rec(TH1D* data, TH1D* sig, TH1D* bgr, TString file_name);
  void draw_purity(TH1D* numerator, TH1D* denominator, TString file_name);
  TH1D* get_difference(TH1D* hist1, TH1D* hist2);
  void draw_chi2(TF1 * fit, std::vector<double> masses, std::vector<double> chi2, double mass, double uncert, TString file_name);
  void draw_delta_comparison( TH1* total_, std::vector<TH1*> MODEL_DELTA, std::vector<TString> UncertNames, TString category, TString file_name);
  void draw_bias(TH1* output_, TH1D* truth_, TH1* bias_, TString file_name);
  void draw_smearFit(TH1D* variation, TF1* fit_, TString file_name);

 private:
  TH1* add_error_bar(TH1* hist, std::vector<double> errors);
  TString directory;

};
