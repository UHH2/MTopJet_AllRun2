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

/* TString directory; */
/* TString input_file; */
/* TString output_file; */
/* TString binning_xml; */

class unfolding{

 public:
  unfolding(TH1D* input, TH1D* background, TH1D* signal, TH2* migration_matrix, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan);
  TH1* get_output(bool);
  TH2* get_prob_matrix();
  TH2* get_cor_matrix();
  TH2* get_cov_matrix();
  TGraph* get_lcurve();
  double get_best_point(TString xy);
  double get_tau();


 private:
  TH1 *output;
  TH1 *output_all;
  TH2 *CorM;
  TH2 *CovM;
  TH2 *ProbM;
  TGraph *lcurve = 0;
  double tau;
  double lcurveX, lcurveY;
};
