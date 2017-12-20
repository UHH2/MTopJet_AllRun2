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



class GenHists_particles: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_particles(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *number_top, *hadtop_pt, *leptop_pt, *hadtop_mass, *leptop_mass,  *hadtop_phi, *leptop_phi, *lepton_pt,*deltaR_hadtop_b, *deltaR_leptop_b, *deltaR_lep_b, *deltaR_lep_neu, *deltaR_hadtop_leptop, *deltaPhi_hadtop_leptop;
    TH1F *deltaR_hadtop_jet1, *deltaPT_hadtop_jet1;
    TH1F *deltaR_hadtop_genjet1, *deltaPT_hadtop_genjet1;

    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<Jet>>h_hadjets;
    uhh2::Event::Handle<std::vector<Particle>>h_hadjets_gen;


};
