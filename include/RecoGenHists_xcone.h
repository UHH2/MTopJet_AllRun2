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
 
class RecoGenHists_xcone: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *MassReso, *PtReso;
    TH1F *DeltaR_Rec1_Gen1, *DeltaR_Rec2_Gen2, *GenJetMass;
    TH1F *pt_rec1_gen1_beforeMatching, *pt_rec2_gen2_beforeMatching, *pt_rec1_gen1_afterMatching, *pt_rec2_gen2_afterMatching;
    TH2F *RecGenMass, *RecGenPT;
    TH1F *MassReso_mass1, *MassReso_mass2, *MassReso_mass3, *MassReso_mass4, *MassReso_mass5, *MassReso_mass6;
    TH1F *PtReso_mass1, *PtReso_mass2, *PtReso_mass3, *PtReso_mass4, *PtReso_mass5, *PtReso_mass6;
    TH1F *MassReso_pt1, *MassReso_pt2, *MassReso_pt3, *MassReso_pt4, *MassReso_pt5, *MassReso_pt6;
    TH1F *PtReso_pt1, *PtReso_pt2, *PtReso_pt3, *PtReso_pt4, *PtReso_pt5, *PtReso_pt6;

    uhh2::Event::Handle<std::vector<Jet>>h_recjets_had;
    uhh2::Event::Handle<std::vector<Jet>>h_recjets_lep;
    uhh2::Event::Handle<std::vector<Particle>>h_genjets_had;
    uhh2::Event::Handle<std::vector<Particle>>h_genjets_lep;
 
};


