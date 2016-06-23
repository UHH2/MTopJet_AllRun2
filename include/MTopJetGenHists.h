#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/MTopJet/include/GenJetProps.h"
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
    TH1F *GenPartNumber, *GenNumber1_new, *GenNumber2_new, *GenNumber1_old, *GenNumber2_old, *GenPT1_new, *GenPT2_new, *GenPT1_old, *GenPT2_old, *GenJetMass_new, *GenJetMass_old;

    uhh2::Event::Handle<double> h_ht;
    uhh2::Event::Handle<double> h_st;
    uhh2::Event::Handle<FlavorParticle> h_primlep; 


};


