#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/MTopJet/include/JetCluster.h"


/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class GenHists_xconeN5: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_xconeN5(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *Combined_Mass1, *Combined_PT1, *Combined_Mass2, *Combined_PT2;
    TH1F *deltaR_lep_jet1, *deltaR_lep_jet2, *deltaR_lep_jet3, *deltaR_lep_jet4, *deltaR_lep_combinedjet1, *deltaR_lep_combinedjet2, *deltaR_combinedjet1_combinedjet2;


    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<Jet>>h_jets;
};


