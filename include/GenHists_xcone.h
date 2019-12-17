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



class GenHists_xcone: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_xcone(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *Mass_HadJet23, *Mass_HadJet23_rebin, *Mass_LepJet23, *MassDiff23;
    TH1F *PT_HadJet23, *PT_LepJet23;
    TH1F *number_HadJet23, *number_LepJet23;

    TH1F *Mass_HadJet33, *Mass_HadJet33_B, *Mass_HadJet33_C, *Mass_HadJet33_rebin, *Mass_HadJet33_unfold, *Mass_LepJet33, *MassDiff33;
    TH1F *PT_HadJet33, *PT_LepJet33;
    TH1F *number_HadJet33, *number_LepJet33;

    TH1F *DeltaR_23_33, *DeltaMass_23_33, *DeltaPT_23_33;
    TH1F *DeltaR_jet2_23_33, *DeltaMass_jet2_23_33, *DeltaPT_jet2_23_33;

    TH1F *SoftDropMass;

    TH1F *min_mass_Wjet, *min_mass_Wjet_zoom, *min_mass_Wjet_ptlow, *min_mass_Wjet_pthigh, *min_mass_Wjet_ptmix, *all_mass_Wjet;
    TH1F *pt_Wjet, *pt_Wjet_i, *pt_Wjet_j, *pt_Wjet_diff;
    TH2F *pt_Wjets;
    TH1F* Wjet_combination;

    TH1F* pt_subjets_had, *pt_subjets_had1, *pt_subjets_had2, *pt_subjets_had3;

    uhh2::Event::Handle<std::vector<GenTopJet>>h_hadjets23;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_lepjets23;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_hadjets33;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_lepjets33;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_fatjets33;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_softdrop;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
};
