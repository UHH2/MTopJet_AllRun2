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

using namespace std;

class chi2sys{
 public:
  chi2sys(TH1* data_, TH2* cov_, TH1D* truth_);
  TH2D* GetCov(){return DELTAcov;};
  TH1D* GetDelta(){return DELTA;};
  double GetFactor(){return factor;};

 private:
  void CreateDeltaCov(double delta);
  double GetChi2NDFValue(TH2D* inputcov);
  TH1D* data;
  TH2D* oldcov;
  TH2D* DELTAcov;
  TH1D* DELTA;
  TH1D* truth;
  int Nbins;
  double factor;
};
