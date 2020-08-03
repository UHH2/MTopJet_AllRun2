#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/JetIds.h"

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

  TH1F *h_number_events;

  // Tau Hists
  TH1F *h_HadJetPt, *h_HadJetMass, *h_HadJetEta, *h_HadJetPhi;
  TH1F *h_XCone_HadJetEta, *h_XCone_HadJetPhi;
  TH1F *h_distance_hadx_hadak;
  TH1F *h_HadJetTau1, *h_HadJetTau2, *h_HadJetTau3, *h_HadJetTau32, *h_HadJetTau23;

  TH1F *h_number_ak8_jets, *h_number_xcone_had_jet, *h_number_ak8_had_jets, *h_number_ak8_not_matched;

  // Gen Hists
  TH1F *h_Number_matched_all, *h_Number_matched_top, *h_Number_matched_q1, *h_Number_matched_q2, *h_Number_matched_bottom;
  TH1F *h_HadJetMass_fullymerged, *h_HadJetMass_semimerged, *h_HadJetMass_notmerged;
  TH1F *h_HadJetTau32_fullymerged, *h_HadJetTau32_semimerged, *h_HadJetTau32_notmerged;

  // Wjet Hists
  TH1F *h_number_ak4_jets, *h_number_ak4_matched_xcone_had, *h_number_equal_distance;
  TH1F *h_ak4_btag, *h_ak4_btag_high, *h_ak4_btag_high_btagcut;
  TH1F *h_wjet_pt_match, *h_wmass_match;
  TH1F *h_wmass_match_ptbin_low, *h_wmass_match_ptbin_midlow, *h_wmass_match_ptbin_midhigh, *h_wmass_match_ptbin_high;
  TH1F *h_distance_ak4_xcone_subjet_btag, *h_difference_distance_small_mid, *h_difference_distance_small_high, *h_difference_distance_mid_high;
  TH1F *h_btag_difference_high_mid, *h_btag_percentage_high_mid, *h_btag_difference_high_mid_sel_btag_min, *h_btag_percentage_high_mid_sel_btag_min;

  TH1F *h_events_no_ak4, *h_events_ak4;
  TH1F *h_wjet_pt_btagcut, *h_wmass_btagcut;
  TH1F *h_btagcut_missed, *h_btagcut_passed;

  TH1F *h_wjet_pt_one_btag, *h_wmass_one_btag;
  TH1F *h_wjet_pt_btagcut_one_btag, *h_wmass_btagcut_one_btag;
  TH1F *h_pass_one_btag, *h_pass_btagcut_one_btag;

  TH1F *h_wjet_pt_one_subjet, *h_wmass_one_subjet;
  TH1F *h_wjet_pt_btagcut_one_subjet, *h_wmass_btagcut_one_subjet;
  TH1F *h_pass_one_subjet, *h_pass_btagcut_one_subjet;

  TH1F *h_wjet_pt_one_btag_subjet, *h_wmass_one_btag_subjet;
  TH1F *h_wjet_pt_btagcut_one_btag_subjet, *h_wmass_btagcut_one_btag_subjet;
  TH1F *h_pass_one_btag_subjet, *h_pass_btagcut_one_btag_subjet;

  TH1F *h_number_no_close_subjets, *h_number_one_close_subjets, *h_number_multiple_close_subjets;
  TH1F *h_wmass_combined_methods, *h_wjet_pt_combined_methods, *h_wjet_combination;
  TH1F *h_wmass_min_mass, *h_wjet_pt_min_mass, *h_wmass_compare, *h_wjet_pt_compare;

  uhh2::Event::Handle<std::vector<TopJet>> h_recjets, h_hadjet, h_lepjet, h_fatjets;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

};
