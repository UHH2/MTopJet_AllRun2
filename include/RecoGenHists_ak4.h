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
 
class RecoGenHists_ak4: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_ak4(uhh2::Context & ctx, const std::string & dirname, bool);
    
    virtual void fill(const uhh2::Event & ev) override;

protected:

    TH1F *MassReso, *PtReso, *WMassReso;
    TH1F *PtReso_1, *PtReso_2, *PtReso_3, *PtReso_4, *PtReso_5, *PtReso_6;
    TH1F *area_all;
    TH1F *Energy_Check;
    bool use_JEC;
 
};


