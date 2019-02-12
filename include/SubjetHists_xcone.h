#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <UHH2/common/include/JetIds.h>
#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

class SubjetHists_xcone: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  SubjetHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type);

    virtual void fill(const uhh2::Event & ev) override;

protected:


    TH1F *pt_had_subjets, *pt_had_subjets_fine, *pt_lep_subjets, *pt_had_subjet1_fine, *pt_had_subjet2_fine, *pt_had_subjet3_fine;
    TH1F *pt_had_subjet1, *pt_had_subjet2, *pt_had_subjet3;
    TH1F *eta_had_subjets, *eta_lep_subjets, *eta_abs_had_subjets;
    TH1F *area_had_subjets, *area_had1_subjet, *area_had2_subjet, *area_had3_subjet;
    TH1F *area_lep_subjets, *area_lep1_subjet, *area_lep2_subjet, *area_lep3_subjet;
    TH1F *area_all_subjets, *pt_all_subjets, *eta_all_subjets;
    TH1F *area_ak4, *pt_ak4, *eta_ak4;
    TH1F *match_to_subjet;

    TH1F *pt_had_combine, *pt_lep_combine;
    TH1F *mass_had_combine, *mass_had_combine_cut, *mass_lep_combine;
    TH1F *min_mass_Wjet, *min_mass_Wjet_zoom, *min_mass_Wjet_ptlow, *min_mass_Wjet_pthigh, *min_mass_Wjet_ptmix, *all_mass_Wjet;
    TH1F *pt_Wjet, *pt_Wjet_i, *pt_Wjet_j, *pt_Wjet_diff;
    TH2F *pt_Wjets;
    TH1F* Wjet_combination;

    TH1F *JEC_all_subjets, *JEC_L1_all_subjets, *JEC_L2L3_all_subjets;
    TH1F *JEC_had_subjets, *JEC_L1_had_subjets, *JEC_L2L3_had_subjets;
    TH1F *JEC_had_bjet, *JEC_L1_had_bjet, *JEC_L2L3_had_bjet;
    TH1F *JEC_had_lightjets, *JEC_L1_had_lightjets, *JEC_L2L3_had_lightjets;
    TH1F *JEC_lep_subjets, *JEC_L1_lep_subjets, *JEC_L2L3_lep_subjets;
    TH1F *JEC_ak4, *JEC_L1_ak4, *JEC_L2L3_ak4;

    TH1F *JEC_L1_had_subjets_ETA_00to05;
    TH1F *JEC_L1_had_subjets_ETA_05to10;
    TH1F *JEC_L1_had_subjets_ETA_10to24;
    TH1F *JEC_L1_had_subjets_ETA_24to28;
    TH1F *JEC_L1_had_subjets_ETA_28to40;


    TH1F *pt_check;
    TH1F *area_final;
    TH2F *area_pt_final;

    TH2F *AreaVsPT, *AreaVsPT_ak4;

    uhh2::Event::Handle<std::vector<TopJet>>h_recfatjets;

};
