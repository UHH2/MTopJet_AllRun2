#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"

#include "UHH2/MTopJet/include/GenericSubJetCorrector_flavor.h"
#include "UHH2/MTopJet/include/Vector_utils.h"

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

class FactorizedJetCorrector;

class CorrectionFactor_JMS: public uhh2::AnalysisModule{
public:

  explicit CorrectionFactor_JMS(uhh2::Context &,  const std::string &, const std::string &, Year year_);
  virtual bool process(uhh2::Event & event) override;
  double get_mass_BestFit(uhh2::Event &, const vector<double>&);

private:

  uhh2::Event::Handle<std::vector<TopJet>>    h_oldjets;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_genjets;
  TString str_year;
  Year year;
  vector<double> points;
  bool isMC;
  std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;

  // XCone ---------------------------------------------------------------------
  void get_function(TString &);
  void get_additionalSYS();
  double get_factor_XCone(double, double, double);

  //double arr[6][12]; // values from table
  //double par[12][3]; // values with function parameters in 12 eta bins
  double eta_bins[13] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191};
  std::vector<TF1*> CentralCorrectionFunctions;
  std::vector<TGraph*> UpDownCorrectionGraphs;
  TGraph *AdditionalSys;

  // JEC -----------------------------------------------------------------------
  vector<double> get_factor_JEC(FactorizedJetCorrector & corrector, vector<Jet> subjets, const Event & event, JetCorrectionUncertainty* jec_unc, double point_J);
  JetCorrectionUncertainty* jec_uncertainty;
  int direction =0; // -1 = down, +1 = up, 0 = nominal
  std::unique_ptr<FactorizedJetCorrector> corrector;
  std::unique_ptr<FactorizedJetCorrector> corrector_MC_2016, corrector_MC_2017, corrector_MC_2018;
  JetCorrectionUncertainty *uncertainty_MC_2016, *uncertainty_MC_2017, *uncertainty_MC_2018;

};
