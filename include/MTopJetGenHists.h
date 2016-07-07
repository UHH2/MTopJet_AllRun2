#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/TTbarGen.h"
 

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class MTopJetGenHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    MTopJetGenHists(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *GenEvNumber, *GenEvNumberStable, *GenPT_old, *GenNumber_old, *GenEta_old, *GenJet1Mass_old, *GenPT_ak06, *GenNumber_ak06, *GenEta_ak06, *GenJet1Mass_ak06, *GenJet1Mass_ak07, *GenJet1Mass_ak08, *GenJet1Mass_ak06_unmatched, *GenJet1Mass_ak07_unmatched, *GenJet1Mass_ak08_unmatched, *GenJet1Mass_ak06_matched, *GenJet1Mass_ak07_matched, *GenJet1Mass_ak08_matched;
    TH1F *RecoEvNumber, *RecoPT_old, *RecoNumber_old, *RecoEta_old, *RecoJet1Mass_old, *RecoPT_ak06, *RecoNumber_ak06, *RecoEta_ak06, *RecoJet1Mass_ak06;

    TH1F *TopHadPT, *Number_Elec, *Number_Muon;

    TH1F *deltaR_bot_jet1, *deltaR_q1_jet1, *deltaR_q2_jet1, *deltaR_bot_lep_jet1, *deltaR_lep1_jet1, *deltaR_lep2_jet1;

    uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;

};


