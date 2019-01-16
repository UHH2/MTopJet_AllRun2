#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/AnalysisModule.h"

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


class NonClosureUncertainty: public uhh2::AnalysisModule{
public:

  explicit NonClosureUncertainty(uhh2::Context &);
  virtual bool process(uhh2::Event & ) override;

private:
  void get_function();
  double get_factor(double);
  Particle GetLepton(uhh2::Event & event);
  uhh2::Event::Handle<std::vector<TopJet>>h_oldjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_newjets;
  TGraph *NonClosureFunction;


};
