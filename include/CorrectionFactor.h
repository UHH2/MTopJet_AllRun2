#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/AnalysisModule.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
using namespace uhh2;

Particle GetLepton(uhh2::Event & event);

class CorrectionFactor: public uhh2::AnalysisModule{
public:

  explicit CorrectionFactor(uhh2::Context &,  const std::string &);
  virtual bool process(uhh2::Event & ) override; 
  /* void get_values_table(); */
  void get_values();
  /* double get_factor_table(double, double); */
  double get_factor(double, double);

private:

  uhh2::Event::Handle<std::vector<TopJet>>h_oldjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_newjets;
  double arr[6][12]; // values from table
  double par[12][3]; // values with function parameters in 12 eta bins
  double pt_bins[7] = {0, 80, 130, 180, 250, 350, 500};
  double eta_bins[13] = {-4, -1.5, -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 4};
};
