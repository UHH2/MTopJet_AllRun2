#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/TTbarGen.h"

#include <iostream>
#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
using namespace std;

// -----------------------------------------------------------------------------
class WeightHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  WeightHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

protected:
  TH1F *h_weight_w, *h_weight_1;
};

// -----------------------------------------------------------------------------
class CountingEventHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  CountingEventHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

protected:
  TH1F *events_w, *events_nw;
};

// -----------------------------------------------------------------------------
class WeightRangeHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  WeightRangeHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

protected:
  TH1F *h_weight_big_w, *h_weight_big_1;
  TH1F *h_weight_medium_w, *h_weight_medium_1;
  TH1F *h_weight_small_w, *h_weight_small_1;

  TH1F *h_weight_negativ_big_w, *h_weight_negativ_big_1;
  TH1F *h_weight_negativ_medium_w, *h_weight_negativ_medium_1;
  TH1F *h_weight_negativ_small_w, *h_weight_negativ_small_1;

  TH1F *events, *events_weight;

};

// -----------------------------------------------------------------------------
class GenWeightRangeHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  GenWeightRangeHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *h_weight_big_w, *h_weight_big_1;
  TH1F *h_weight_medium_w, *h_weight_medium_1;
  TH1F *h_weight_small_w, *h_weight_small_1;

  TH1F *h_weight_negativ_big_w, *h_weight_negativ_big_1;
  TH1F *h_weight_negativ_medium_w, *h_weight_negativ_medium_1;
  TH1F *h_weight_negativ_small_w, *h_weight_negativ_small_1;

  TH1F *events, *events_weight;

};
