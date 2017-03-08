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



class GenHists_particles: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_particles(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *W_pt, *top_pt, *elec_pt, *muon_pt;

};
