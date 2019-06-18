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
#include <vector>
#include "TROOT.h"
#include <TVectorD.h>
#include <TDecompLU.h>
#include <iostream>
#include <iomanip>
#include <TF1.h>

class chi2fit{
 public:
  chi2fit(TH1D* data, TH2D* cov_matrix, std::vector<TH1D*> masspoints, std::vector<double> masses, double lower, double upper, bool width);

  void CalculateChi2();

  TF1* GetChi2Fit();
  std::vector<double> GetChi2Values();
  double GetMass();
  double GetMin();
  double GetUncertainty();

 private:
  void SetMasses(std::vector<double> masses);
  double ComputeChi2(TH1D* MC, std::vector<int> bins);
  void SetBins(std::vector<TH1D*> masspoints);
  void SetRange(double lower, double upper);
  void SetMasspoints(std::vector<TH1D*> masspoints);
  void SetData(TH1D* data);
  void SetCovMatrix(TH2D* cov_matrix);


  std::vector<TH1D*> masspoints_;
  TH1D* data_;
  int lower_, upper_;
  bool width_;
  TH2D* cov_;
  std::vector< std::vector<int> > bins_;
  TVectorD chi2_;
  TVectorD masses_;
  TF1 * fit_;

};
