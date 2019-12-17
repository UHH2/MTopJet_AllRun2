#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/GenTopJet.h"
#include "UHH2/core/include/GenJet.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"

#include <vector>
#include <iostream>
#include <math.h>

using namespace std;
using namespace uhh2;

class RemoveLepton: public uhh2::AnalysisModule{

 public:
   RemoveLepton(uhh2::Context & ctx, const std::string & name_jet);
   virtual bool process(uhh2::Event & ) override;

 private:
   uhh2::Event::Handle<std::vector<TopJet>>h_topjets;
};


class RemoveLeptonGen: public uhh2::AnalysisModule{

 public:
   RemoveLeptonGen(uhh2::Context & ctx, const std::string & name_jet);
   virtual bool process(uhh2::Event & ) override;

 private:
   uhh2::Event::Handle<std::vector<GenTopJet>>h_topjets;
   uhh2::Event::Handle<TTbarGen> h_ttbargen;
};
