#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/TTbarGen.h"


class ClusteringHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  ClusteringHists(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH2F *FatJet1, *FatJet2, *DecayProducts, *Lepton, *NotClustered, *AllPF, *SubJet1_1, *SubJet1_2, *SubJet1_3, *SubJet2_1, *SubJet2_2,  *NotClustered_sub1,  *NotClustered_sub2;

    uhh2::Event::Handle<std::vector<int>>h_list_fat;
    uhh2::Event::Handle<std::vector<int>>h_list_sub1;
    uhh2::Event::Handle<std::vector<int>>h_list_sub2;
    uhh2::Event::Handle<TTbarGen>h_ttbargen;
};
