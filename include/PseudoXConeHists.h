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

class PseudoXConeHists: public uhh2::Hists {

  public:
    PseudoXConeHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

  protected:

    TH1F *pt_diff_fat, *eta_diff_fat, *dR_fat;
    TH1F *pt_diff_sub, *eta_diff_sub, *dR_sub;

    uhh2::Event::Handle<std::vector<TopJet>>h_xcone;
    uhh2::Event::Handle<std::vector<TopJet>>h_pseudoxcone;

};
