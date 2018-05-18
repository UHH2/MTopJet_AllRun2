#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
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
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldDensity.h"
#include <TString.h>

class unfolding{

 public:
  unfolding(TH1D* input, std::vector<TH1D*> backgrounds, std::vector<TString>bgr_name, TH1D* signal, TH2* migration_matrix, std::vector<TH2*>sys_matrix, std::vector<TString>sys_name, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan);
  TH1* get_output(bool);
  TH2* get_prob_matrix();
  TH2* get_cor_matrix();
  TH2* get_cov_matrix();
  std::vector<TH1*> get_delta();
  std::vector<TH2*> get_sys_matrix();
  TGraph* get_lcurve();
  double get_best_point(TString xy);
  double get_tau();


 private:
  TH2* CreateCovMatrixFromDelta(TH1* delta);
  TH1 *output;
  TH1 *output_all;
  TH2 *CorM;
  TH2 *CovM;
  TH2 *ProbM;
  TGraph *lcurve = 0;
  double tau;
  double lcurveX, lcurveY;
  std::vector<TH2*>SysCov;
  std::vector<TH1*> SysDelta;
};
