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
 
class RecoGenHists_subjets: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_subjets(uhh2::Context & ctx, const std::string & dirname,  const std::string & type);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *MassReso, *PtReso, *WMassReso;
    TH1F *MassReso_iso, *PtReso_iso;
    TH1F *PtReso_1, *PtReso_2, *PtReso_3, *PtReso_4, *PtReso_5, *PtReso_6;
    TH1F *PtReso_iso_1, *PtReso_iso_2, *PtReso_iso_3, *PtReso_iso_4, *PtReso_iso_5, *PtReso_iso_6;
    TH1F *area_all, *area_iso;
    TH1F *min_mass_Wjet_rec, *min_mass_Wjet_gen;

    uhh2::Event::Handle<std::vector<TopJet>>h_recjets;
    uhh2::Event::Handle<std::vector<Jet>>h_hadjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_genjets;

};


