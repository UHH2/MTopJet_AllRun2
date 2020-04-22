#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

class RecoGenHists_xcone_topjet: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_xcone_topjet(uhh2::Context & ctx, const std::string & dirname, bool isTTbar, double masscut);

  virtual void fill(const uhh2::Event & ev) override;

private:
  bool isTTbar_;
  double masscut_;

protected:

  TH1F *HadJetPt, *HadJetMass, *HadJetEta, *HadJetPhi;
  TH1F *XCone_HadJetEta, *XCone_HadJetPhi;
  TH1F *distance_hadx_hadak;
  TH1F *HadJetTau1, *HadJetTau2, *HadJetTau3, *HadJetTau32, *HadJetTau23;

  TH1F *number_ak8_jets, *number_xcone_jets, *number_ak8_had_jets;

  TH1F *Number_matched_all, *Number_matched_top, *Number_matched_q1, *Number_matched_q2, *Number_matched_bottom;
  TH1F *HadJetMass_fullymerged, *HadJetMass_semimerged, *HadJetMass_notmerged;
  TH1F *HadJetTau32_fullymerged, *HadJetTau32_semimerged, *HadJetTau32_notmerged;

  uhh2::Event::Handle<std::vector<TopJet>>h_recjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_xcone;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

};
