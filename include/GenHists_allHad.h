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



class GenHists_allHad: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_allHad(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *number_top, *top_pt, *top_Dphi;
    TH1F *bjet_index, *bjet_index_nomatch_sub, *bjet_index_nomatch_fat, *b_in_leading_sub;
    TH1F *index_top_antitop;
    TH1F *deltaR_b_bjet;

    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_fatjets_gen;
};
