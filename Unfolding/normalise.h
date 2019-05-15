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

class Normalise{
 public:
  Normalise(TH1* hist, TH2* oldMatrix, double lower, double upper, bool width);
  Normalise(TH1* hist, double lower, double upper, bool width);
  ~Normalise();
  TH1D* GetHist();
  TH2D* GetMatrix();

 private:
  void SetHist(TH1* hist);
  void SetMatrix(TH2* matrix);
  void SetRange(double lower, double upper);
  void Norm();
  TH2D* MatrixErrorPropagation(TH2*old_cov);
  TH2D* CalcMatrix();
  void UpdateHist();

  TH1D* hist_;
  TH2D* matrix_;
  int lower_, upper_;
  bool width_;
};
