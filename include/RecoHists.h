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


class RecoHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *RecoJetNumber, *NLepton;
    TH1F *RecoJet1Mass, *RecoJet2Mass, *Mass1Mass2;
    TH1F *RecoJet1PT, *RecoJet2PT, *RecoJet1Jet2PT, *RecoJet3PT, *LeptonPT, *RecoJetPT;
    TH1F *RecoJet2Eta;

    uhh2::Event::Handle<std::vector<Jet>>h_jets;

};


