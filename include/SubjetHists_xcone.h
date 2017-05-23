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
 
class SubjetHists_xcone: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  SubjetHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:


    TH1F *pt_had_subjets, *pt_had_subjets_fine, *pt_lep_subjets;
    TH1F *eta_had_subjets, *eta_lep_subjets;
    TH1F *area_had_subjets, *area_had1_subjet, *area_had2_subjet, *area_had3_subjet;
    TH1F *area_lep_subjets, *area_lep1_subjet, *area_lep2_subjet, *area_lep3_subjet;
    TH1F *area_all_subjets, *pt_all_subjets, *eta_all_subjets;
    TH1F *area_ak4, *pt_ak4, *eta_ak4;

    TH1F *pt_had_combine, *pt_lep_combine;
    TH1F *mass_had_combine, *mass_had_combine_cut, *mass_lep_combine;
    TH1F *min_mass_Wjet, *min_mass_Wjet_zoom, *all_mass_Wjet;

    TH1F *JEC_all_subjets, *JEC_L1_all_subjets, *JEC_L2L3_all_subjets;
    TH1F *JEC_ak4, *JEC_L1_ak4, *JEC_L2L3_ak4;

    TH1F *pt_check;
    TH1F *area_final;
    TH2F *area_pt_final;

    uhh2::Event::Handle<std::vector<TopJet>>h_recfatjets;

};


