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
    TH2F *FatJet1, *FatJet2, *DecayProducts, *DecayProducts_lep, *DecayProducts_had, *Lepton, *NotClustered, *AllPF, *SubJet1_1, *SubJet1_2, *SubJet1_3, *SubJet2_1, *SubJet2_2,  *NotClustered_sub1,  *NotClustered_sub2, *Top, *Top_lep, *Top_had, *Neutrino, *Bottom_had, *Bottom_lep, *LepJet, *HadJet;
    TH1F *PT, *Mass, *List_incjet, *List_subjet1, *List_subjet2, *Number_incjets, *LepJetBool, *Eta;

    uhh2::Event::Handle<std::vector<int>>h_list_fat;
    uhh2::Event::Handle<std::vector<int>>h_list_sub1;
    uhh2::Event::Handle<std::vector<int>>h_list_sub2;
    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<Jet>>h_jets;
    uhh2::Event::Handle<std::vector<Jet>>h_incjets;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_fat1;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_fat2;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub1_1;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub1_2;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub1_3;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub2_1;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub2_2;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_fat0;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub1_0;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_sub2_0;
    uhh2::Event::Handle<std::vector<Jet>>h_pf_all;

};
