#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <UHH2/common/include/MCWeight.h>
#include "UHH2/core/include/LorentzVector.h"
#include <UHH2/common/include/TTbarReconstruction.h>


#include <math.h>
#include <vector>

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
  TH1F *MuonPt, *MuonPerpendicularFullPt, *MuonPerpendicularLeadPt, *MuonParallelFullPt, *MuonParallelLeadPt;
  TH1F *MuonBoostFullPt, *MuonBoostLeadPt, *TopBoostFullPt, *TopBoostLeadPt;
  TH1F *leading_subjet_position;

  uhh2::Event::Handle<std::vector<TopJet>> h_fatjets;

};

// Lorentz to TLorentzVector
TLorentzVector lorentz_to_tlorentz(const LorentzVector v4);
