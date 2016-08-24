#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
/* #include "UHH2/core/include/PFParticle.h" */
/* #include "UHH2/common/include/TTbarGen.h" */
/* #include "UHH2/MTopJet/include/JetCluster.h" */
 

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class MTopJetRecoHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  MTopJetRecoHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *RecoJetNumber, *NLepton;
    TH1F *RecoJet1Mass, *Mass1Mass2;
    TH1F *RecoJet1PT, *RecoJet2PT, *RecoJet1Jet2PT, *RecoJet3PT, *LeptonPT, *RecoJetPT;
    TH1F *RecoJet2Eta;

    uhh2::Event::Handle<std::vector<Jet>>h_jets;

};


