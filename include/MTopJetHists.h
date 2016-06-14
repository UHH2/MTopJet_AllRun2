#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

 

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class MTopJetHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    MTopJetHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *N_PrimVertices, *N_TrueInteractions, *Weights, *WeightsLogBins, *MET, *HT, *HTLep, *ST, *BTAG_L, *BTAG_M, *BTAG_T, *deltaR_lep_topjet1, *deltaR_lep_topjet2, *deltaR_lep_jet1, *deltaR_lep_jet2, *deltaPhi_lep_topjet1, *deltaPhi_lep_topjet2, *TopNumber, *TopPT1, *TopPT2, *TopJetMass;
    TH2F *TopJetMass1_TopJetMass2;

    uhh2::Event::Handle<double> h_ht;
    uhh2::Event::Handle<double> h_st;
    uhh2::Event::Handle<FlavorParticle> h_primlep; 

};


