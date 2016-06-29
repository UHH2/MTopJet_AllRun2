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

    TH1F *GenEvNumber, *GenPT_old, *GenNumber_old, *GenEta_old, *GenJet1Mass_old, *GenPT_ak08, *GenNumber_ak08, *GenEta_ak08, *GenJet1Mass_ak08, *GenPT_ak10, *GenNumber_ak10, *GenEta_ak10, *GenJet1Mass_ak10, *GenPT_ak12, *GenNumber_ak12, *GenEta_ak12, *GenJet1Mass_ak12, *GenPT_ak14, *GenNumber_ak14, *GenEta_ak14, *GenJet1Mass_ak14;

    TH1F *RecoEvNumber, *RecoPT_old, *RecoNumber_old, *RecoEta_old, *RecoJet1Mass_old, *RecoPT_ak08, *RecoNumber_ak08, *RecoEta_ak08, *RecoJet1Mass_ak08;
    uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;

};


