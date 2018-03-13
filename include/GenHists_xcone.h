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


    uhh2::Event::Handle<std::vector<Particle>>h_hadjets23;
    uhh2::Event::Handle<std::vector<Particle>>h_lepjets23;
    uhh2::Event::Handle<std::vector<Particle>>h_hadjets33;
    uhh2::Event::Handle<std::vector<Particle>>h_lepjets33;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_softdrop;
};
