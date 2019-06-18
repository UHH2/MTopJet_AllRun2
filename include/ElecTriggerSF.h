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
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;
using namespace uhh2;


class ElecTriggerSF: public uhh2::AnalysisModule{
public:

  explicit ElecTriggerSF(uhh2::Context & ctx, std::string, TString pe);
  virtual bool process(uhh2::Event & ) override;

private:
  TH1F *h_sf_lo;
  TH1F *h_sf_me;
  TH1F *h_sf_hi;
  bool isMC;
  TString pteta;

};
