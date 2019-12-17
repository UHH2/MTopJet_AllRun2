#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/core/include/GenParticle.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace uhh2;
using namespace std;


class GenHists_process: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_process(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *number_each_particle, *number_muon, *number_top, *number_bottom, *number_W;
    TH1F *number_each_particle_cut, *number_muon_cut, *number_top_cut, *number_bottom_cut, *number_W_cut;
    TH1F *mother_muon_1, *mother_muon_1_cut, *mother_muon_1_outside, *daughter_muon, *daughter_muon_cut, *mother_muon_2;
    TH1F *number_genparticles;
    TH1F *only_one_lepton, *only_one_lepton_cut;
    TH1F *all_status, *status_muon, *status_muon_cut;
    TH1F *deltaR_muon1_muon2, *deltaR_muon1_muon3, *deltaR_muon1_muon4, *deltaR_muon1_muon5, *deltaR_muon1_muon6, *deltaR_muon1_muon7 ;

    TH2F *index_status;

};
