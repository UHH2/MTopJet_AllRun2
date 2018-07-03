#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;



class RecoHists_puppi: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  RecoHists_puppi(uhh2::Context & ctx, const std::string & dirname);
  virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *HadJetMass, *HadJetMass_rebin, *HadJetMass_B, *LepJetMass, *HadMassLepMass;
  TH1F *HadJetEta, *HadJetPhi, *LepJetEta, *LepJetPhi;
  TH1F *HadJetPT, *LepJetPT;

  uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;

};
