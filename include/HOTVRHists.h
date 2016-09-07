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
class HOTVRHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname, double rho);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:


    TH1F  *Reff_HOTVR, *Reff2_HOTVR;
    TH2F *Mass_Reff_HOTVR, *Mass2_Reff2_HOTVR;
    uhh2::Event::Handle<std::vector<Jet>>h_jets;
    double rho_;

};


