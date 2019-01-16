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
    TH1F *deltaR_hadtop_jet1, *deltaPT_hadtop_jet1, *deltaR_bot_jet1;
    TH1F *deltaR_hadtop_genjet1, *deltaPT_hadtop_genjet1;
    TH1F *bjet_index, *bjet_index_nomatch, *MinDeltaR_bot_subjet, *InsideFat_MinDeltaR_bot_subjet, *b_inside_fatjet;

    TH1F *closest_particle_0,*closest_particle_50, *closest_particle_100, *closest_particle_200, *closest_particle_300, *closest_particle_400, *closest_particle_500;
    TH2F *dR_pt_b_subjet, *dR_pt_ud_subjet, *dR_pt_g_subjet;

    TH1F *fraction_ud_A, *fraction_s_A, *fraction_c_A, *fraction_b_A, *fraction_g_A, *fraction_nomatch_A;
    TH1F *fraction_ud_B, *fraction_s_B, *fraction_c_B, *fraction_b_B, *fraction_g_B;
    TH1F *fraction_ud_C, *fraction_s_C, *fraction_c_C, *fraction_b_C, *fraction_g_C;
    TH1F *N_matches;

    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<TopJet>>h_hadjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_hadjets_gen;


};
