#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;
 
class RecoGenHists_topjet: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_topjet(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *MassReso, *PtReso;
    TH1F *DeltaR_Rec1_Gen1, *DeltaR_Rec2_Gen2, *GenJetMass;
    TH1F *pt_rec1_gen1_beforeMatching, *pt_rec2_gen2_beforeMatching, *pt_rec1_gen1_afterMatching, *pt_rec2_gen2_afterMatching;

    uhh2::Event::Handle<std::vector<TopJet>>h_recjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_genjets;

};


