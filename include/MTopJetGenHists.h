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
class MTopJetGenHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  MTopJetGenHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *GenJet1Mass;
    TH1F *GenJet1PT, *GenJet2PT, *LeptonPT;
    TH1F *TopHadPT, *TopLepPT;
    TH1F *deltaR_bot_jet1, *deltaR_bot_lep_jet1,*deltaR_q1_jet1, *deltaR_q2_jet1, *deltaR_lep1_jet1, *deltaR_lep2_jet1,  *deltaR_lep1_jet2, *deltaR_lep2_jet2, *deltaR_tophad_jet1, *deltaR_toplep_jet1;

    uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;
    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<Jet>>h_jets;

};


