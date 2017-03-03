#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/GenJetWithParts.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"

#include <vector>
#include <iostream> 
#include <math.h>

using namespace std;
using namespace uhh2;


class CombineXCone{

 public:
  Particle GetLepton(uhh2::Event &);
  Jet AddSubjets(vector<Jet> subjets, double ptmin);
  GenParticle GetLepton_gen(uhh2::Event &);
  Particle AddSubjets_gen(vector<Particle> subjets, double ptmin);
};


class CombineXCone33: public uhh2::AnalysisModule{
public:

  explicit CombineXCone33(uhh2::Context &);
  virtual bool process(uhh2::Event & ) override; 
    
private:

  uhh2::Event::Handle<std::vector<Jet>>h_xcone33hadjets;
  uhh2::Event::Handle<std::vector<Jet>>h_xcone33lepjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;

};

class CombineXCone33_gen: public uhh2::AnalysisModule{
public:

  explicit CombineXCone33_gen(uhh2::Context &);
  virtual bool process(uhh2::Event & ) override; 
    
private:

  uhh2::Event::Handle<std::vector<Particle>>h_GENxcone33hadjets;
  uhh2::Event::Handle<std::vector<Particle>>h_GENxcone33lepjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENfatjets;

};
