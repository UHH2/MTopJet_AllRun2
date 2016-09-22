#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/TTbarGen.h"
/* #include "UHH2/MTopJet/include/JetCluster.h" */


/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class GenHists_topjet: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  GenHists_topjet(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *GenJetNumber;
    TH1F *GenJet1Mass, *GenJet2Mass, *GenJet2MassBool, *GenJet2LepMassBool, *GenJet2LepMass, *GenJet2Mass_scale, *GenJet2LepMass_scale, *Mass1Mass2;
    TH2F *GenJet1MassJet2LepMass;
    TH1F *GenJet1PT, *GenJet2PT, *GenJet1Jet2PT, *GenJet3PT, *LeptonPT, *GenJetPT;
    TH1F *GenJet2Eta;
    TH1F *TopHadPT, *TopLepPT;
    TH1F *deltaR_bot_jet1, *deltaR_botlep_jet1,*deltaR_q1_jet1, *deltaR_q2_jet1, *deltaR_lep1_jet1, *deltaR_lep2_jet1, *deltaR_tophad_jet1, *deltaR_toplep_jet1;
    TH1F *deltaR_bot_jet2, *deltaR_botlep_jet2,*deltaR_q1_jet2, *deltaR_q2_jet2, *deltaR_lep1_jet2, *deltaR_lep2_jet2, *deltaR_tophad_jet2, *deltaR_toplep_jet2;
    /* TH1F *deltaR_bot_jet3, *deltaR_botlep_jet3,*deltaR_q1_jet3, *deltaR_q2_jet3, *deltaR_lep1_jet3, *deltaR_lep2_jet3, *deltaR_tophad_jet3, *deltaR_toplep_jet3; */
    /* TH1F *deltaR_bot_jet4, *deltaR_botlep_jet4,*deltaR_q1_jet4, *deltaR_q2_jet4, *deltaR_lep1_jet4, *deltaR_lep2_jet4, *deltaR_tophad_jet4, *deltaR_toplep_jet4; */
    /* TH1F *deltaR_bot_jet5, *deltaR_botlep_jet5,*deltaR_q1_jet5, *deltaR_q2_jet5, *deltaR_lep1_jet5, *deltaR_lep2_jet5, *deltaR_tophad_jet5, *deltaR_toplep_jet5; */
    /* TH1F *deltaR_bot_jet6, *deltaR_botlep_jet6,*deltaR_q1_jet6, *deltaR_q2_jet6, *deltaR_lep1_jet6, *deltaR_lep2_jet6, *deltaR_tophad_jet6, *deltaR_toplep_jet6; */
    
    uhh2::Event::Handle<TTbarGen>h_ttbargen;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_gentopjets;
};


