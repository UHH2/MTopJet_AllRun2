#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;
using namespace std;

class RecoGenHists_allHad: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  RecoGenHists_allHad(uhh2::Context & ctx, const std::string & dirname,  const std::string & type, const std::string & subjet_selection, double ptreccut_, double ptgencut_);
  virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *MassReso, *PtReso, *WMassReso;
  TH1F *PtReso_1, *PtReso_2, *PtReso_3, *PtReso_4, *PtReso_5, *PtReso_6, *PtReso_7, *PtReso_8, *PtReso_9, *PtReso_10;
  TH1F *PtReso_area1, *PtReso_area2, *PtReso_area3, *PtReso_area4, *PtReso_area5, *PtReso_area6, *PtReso_area7, *PtReso_area8;
  TH1F *PtReso_rec1, *PtReso_rec2, *PtReso_rec3, *PtReso_rec4, *PtReso_rec5, *PtReso_rec6, *PtReso_rec7, *PtReso_rec8, *PtReso_rec9, *PtReso_rec10;
  TH1F *PtRec_1, *PtRec_2, *PtRec_3, *PtRec_4, *PtRec_5, *PtRec_6, *PtRec_7, *PtRec_8, *PtRec_9, *PtRec_10;
  TH1F *area_all, *area_iso;
  TH1F *min_mass_Wjet_rec, *min_mass_Wjet_gen;

  uhh2::Event::Handle<std::vector<TopJet>>h_recjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_genjets;
  uhh2::Event::Handle<TTbarGen>h_ttbargen;

  double ptgencut, ptreccut;
  bool only_one_jet;
  int index_onlyjet;
};
