#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/AnalysisModule.h"
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;
using namespace uhh2;

class StoreBJet: public uhh2::AnalysisModule{
public:
  StoreBJet(uhh2::Context & ctx, const std::string & collectionname);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_bjets;
};
