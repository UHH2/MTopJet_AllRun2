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

class ModelSys{
 public:
  ModelSys(TH1* unfolded, TH1D* truth);
  TH2* GetCov();

 private:
  void CalculateCov(TH1D* hist);
  TH1D* GetDiff(TH1* unfolded, TH1D* truth);
  TH2* cov;

};
