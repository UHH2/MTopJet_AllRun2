#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;
using namespace uhh2;



class CorrectionFactor_BestFit: public uhh2::AnalysisModule{
public:

  explicit CorrectionFactor_BestFit(uhh2::Context &,  const std::string & , Year year_);
  virtual bool process(uhh2::Event & event) override;
  double get_mass_BestFit(uhh2::Event &, vector<TopJet>, vector<double>);

private:

  uhh2::Event::Handle<std::vector<TopJet>>h_oldjets;
  Year year;
  TString str_year;
  vector<double> points;
  bool isMC;
  std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;

  // XCone ---------------------------------------------------------------------
  void get_function(TString);
  void get_additionalSYS();
  double get_factor_XCone(double, double, double);

  //double arr[6][12]; // values from table
  //double par[12][3]; // values with function parameters in 12 eta bins
  double eta_bins[13] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191};
  std::vector<TF1*> CentralCorrectionFunctions;
  std::vector<TGraph*> UpDownCorrectionGraphs;
  TGraph *AdditionalSys;

  // JEC -----------------------------------------------------------------------

  JetCorrectionUncertainty* jec_uncertainty;
  int direction =0; // -1 = down, +1 = up, 0 = nominal
  std::unique_ptr<FactorizedJetCorrector> corrector;

};
