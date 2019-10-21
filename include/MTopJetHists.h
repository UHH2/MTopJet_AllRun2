#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/JetIds.h"



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
    TH1F *N_PrimVertices, *N_TrueInteractions, *Weights, *WeightsLogBins, *MET, *HT, *HTLep, *ST, *deltaR_lep_Bjet, *deltaR_lep_topjet1, *deltaR_lep_topjet2, *deltaR_lep_jet1, *deltaR_lep_jet2, *deltaPhi_lep_topjet1, *deltaPhi_lep_topjet2, *TopNumber, *TopPT1, *TopPT2, *TopJetMass;
    TH1F *BTAG_L_CSV, *BTAG_M_CSV, *BTAG_T_CSV, *BTAG_L_DJ, *BTAG_M_DJ, *BTAG_T_DJ, *BTAG_VALUE_DJ;
    TH2F *TopJetMass1_TopJetMass2;

    uhh2::Event::Handle<double> h_ht;
    uhh2::Event::Handle<double> h_st;
    uhh2::Event::Handle<FlavorParticle> h_primlep;

};
