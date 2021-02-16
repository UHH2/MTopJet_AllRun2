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
#include <TGraphAsymmErrors.h>
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

// -----------------------------------------------------------------------------
// new method from Alex F

class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString syst_direction_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, syst_direction;
  Year year;
  std::vector<float> pt_bins;
  std::unique_ptr<TGraphAsymmErrors> g_sf_pt1, g_sf_pt2;
  uhh2::Event::Handle<float> h_ele_weight, h_ele_weight_up, h_ele_weight_down;

};
