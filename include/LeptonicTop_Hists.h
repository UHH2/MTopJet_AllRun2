#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <UHH2/common/include/MCWeight.h>
#include "UHH2/core/include/LorentzVector.h"
#include <UHH2/common/include/TTbarReconstruction.h>
#include "UHH2/MTopJet/include/Vector_utils.h"

#include <math.h>
#include <vector>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;
using namespace std;

class LeptonicTop_Hists: public Hists {
public:

  LeptonicTop_Hists(Context & ctx, const string & dirname);
  virtual void fill(const Event & ev) override;

protected:

  TH1F *Full_LepJetMass, *Lead_LepJetMass;
  TH1F *muon_pt, *muon_perp_full_momentum, *muon_perp_lead_momentum, *muon_para_full_momentum, *muon_para_lead_momentum;
  TH1F *muon_boost_perp_full_momentum, *muon_boost_perp_lead_momentum, *muon_boost_para_full_momentum, *muon_boost_para_lead_momentum;
  TH1F *top_boost_full_pt, *top_boost_lead_pt, *muon_boost_full_pt, *muon_boost_lead_pt;
  TH1F *leading_subjet_position;

  uhh2::Event::Handle<std::vector<TopJet>> h_fatjets;

};
