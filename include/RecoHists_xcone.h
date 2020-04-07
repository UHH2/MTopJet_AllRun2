#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <UHH2/common/include/MCWeight.h>


#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;



class RecoHists_xcone: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  RecoHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type);
  virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *HadJetMass, *HadJetMass_rebin, *HadJetMass_B, *LepJetMass, *HadMassLepMass, *SoftdropMass_had, *SoftdropMass_lep, *SoftdropMass_Sel;
  TH1F *HadJetMass_ptbin_1, *HadJetMass_ptbin_2, *HadJetMass_ptbin_3, *HadJetMass_ptbin_4;
  TH1F *HadJetEta, *HadJetPhi, *LepJetEta, *LepJetPhi;
  TH1F *HadJetPT, *LepJetPT, *FatJetPT_had, *FatJetPT_lep;
  TH1F *DeltaRDiff;
  TH1F *number_hadjet, *number_lepjet;
  TH1F *FatJetPTDiff_had, *FatJetMassDiff_had, *FatJetPTDiff_lep, *FatJetMassDiff_lep;
  TH1F *RhoA, *RhoA_fat, *RhoA_diff, *E_diff;
  TH2F *Mass_Vertices;
  TH1F *JER_factor;
  TH1F *DeltaR_btagT_nextjet, *DeltaR_btagM_nextjet, *DeltaR_btagT_xcone;
  TH1F *pt_probe1, *pt_probe2;
  TH1F *csvmax;
  TH1F *number_smalldR_all,*number_smalldR_pass,*number_largedR_all,*number_largedR_pass;

  uhh2::Event::Handle<std::vector<TopJet>>h_hadjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_lepjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;

};
