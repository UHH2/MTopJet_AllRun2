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

/*
This collection contains HistClasses which are to small for seperated files.
They are used for ControlHists like event counting or HistClasses which are only
creating one hists for one variable. If the Classes are getting to big, they
will be shifted into an individual file.
*/

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

// -----------------------------------------------------------------------------
class JetMassScaleHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  JetMassScaleHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;
  virtual void fill_mass(const uhh2::Event & ev, const double mass, const double wmass);

private:
  double mass, wmass;

  TH1F *h_mass, *h_wmass;
  TH1F *events, *events_weight;
};

// -----------------------------------------------------------------------------
class MissingPtHist: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  MissingPtHist(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

private:

  TH1F *MPT;
};

// -----------------------------------------------------------------------------
class PositionBTagHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  PositionBTagHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

private:

  TH1F *Position1stBtag, *Position2ndBtag;
  TH1F *PositionBoth;

  uhh2::Event::Handle<std::vector<TopJet>>h_hadjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_lepjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;
};
