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
  unfolding(TH1D* input, std::vector<TH1D*> backgrounds,  std::vector<TString> bgr_name, TH1D* signal, TH2* migration_matrix, std::vector< std::vector<TH2*> > sys_matrix, std::vector< std::vector<TString> > sys_name, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan, double tau_);
  TH1* get_output(bool);
  TH2* get_prob_matrix();
  TH2* get_cor_matrix();
  TH1* GetBiasDistribution();
  TSpline* GetRhoTau();
  std::vector< std::vector<TH1*> > get_sys_delta();
  std::vector< std::vector<TH2*> > GetSysCov();
  // TH2* GetTotalCov();

  std::vector<TH1*> get_bgr_delta();
  TGraph* get_lcurve();
  double get_best_point(TString xy);
  double get_tau();

  TH2* GetInputStatCov();
  TH2* GetMatrixStatCov();
  std::vector<TH2*> GetBgrStatCov();
  std::vector<TH2*> GetBgrScaleCov();


 private:
  TH2* CreateCovMatrixFromDelta(TH1* delta);
  TH1 *output;
  TH1 *output_all;
  TH2 *CorM;
  TH2 *ProbM;
  TH1 *bias;
  TSpline *rhoLogTau=0;
  TGraph *lcurve = 0;
  double tau;
  double lcurveX, lcurveY;
  std::vector< std::vector<TH1*> > SysDelta;
  std::vector< std::vector<TH2*> > SysCov;

  std::vector<TH2*> BgrCov;
  std::vector<TH1*> BgrDelta;

  std::vector<TH2*> CovBgrStat;
  std::vector<TH2*> CovBgrScale;
  TH2* CovInputStat;
  TH2* CovMatrixStat;
  // TH2* CovTotal;

};
