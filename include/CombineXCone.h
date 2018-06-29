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

  bool FindLepton(uhh2::Event &);
  bool FindLepton_gen(uhh2::Event &);
  Particle GetLepton(uhh2::Event &);
  TopJet CreateTopJetFromSubjets(vector<Jet> subjets, double ptmin);
  GenParticle GetLepton_gen(uhh2::Event &);
  GenTopJet CreateTopJetFromSubjets_gen(vector<Particle> subjets, double ptmin);
};


class CombineXCone33: public uhh2::AnalysisModule{
public:

  explicit CombineXCone33(uhh2::Context &,  const std::string &, const std::string & , const std::string &);
  virtual bool process(uhh2::Event & ) override;

private:

  uhh2::Event::Handle<std::vector<TopJet>>h_xcone33hadjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_xcone33lepjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_fatjets;

};

class CombineXCone33_gen: public uhh2::AnalysisModule{
public:

  explicit CombineXCone33_gen(uhh2::Context &);
  virtual bool process(uhh2::Event & ) override;

private:

  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENxcone33hadjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENxcone33lepjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENfatjets;

};


class CombineXCone23_gen: public uhh2::AnalysisModule{
public:

  explicit CombineXCone23_gen(uhh2::Context &);
  virtual bool process(uhh2::Event & ) override;

private:

  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENxcone23hadjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_GENxcone23lepjets;
  uhh2::Event::Handle<std::vector<GenTopJet>>h_GEN23fatjets;

};


class CopyJets: public uhh2::AnalysisModule{
public:

  explicit CopyJets(uhh2::Context &, const std::string &, const std::string &);
  virtual bool process(uhh2::Event & ) override;

private:

  uhh2::Event::Handle<std::vector<TopJet>>h_new;
  uhh2::Event::Handle<std::vector<TopJet>>h_old;

};
