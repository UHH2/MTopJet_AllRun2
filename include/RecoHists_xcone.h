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



class RecoHists_xcone: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *HadJetMass, *HadJetMass_rebin, *LepJetMass, *HadMassLepMass, *SoftdropMass_had, *SoftdropMass_lep, *SoftdropMass_Sel;
    TH1F *HadJetEta, *HadJetPhi, *LepJetEta, *LepJetPhi;
    TH1F *HadJetPT, *LepJetPT, *FatJetPT_had, *FatJetPT_lep;
    TH1F *DeltaRDiff;
    TH1F *number_hadjet, *number_lepjet;
    TH1F *FatJetPTDiff_had, *FatJetMassDiff_had, *FatJetPTDiff_lep, *FatJetMassDiff_lep;
    TH1F *RhoA, *RhoA_fat, *RhoA_diff, *E_diff;
    TH2F *Mass_Vertices;

    uhh2::Event::Handle<std::vector<Jet>>h_hadjets;
    uhh2::Event::Handle<std::vector<Jet>>h_lepjets;
    uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;

};


