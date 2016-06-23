#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/MTopJet/include/GenJetProps.h"
#include "UHH2/common/include/TTbarGen.h"
#include <fastjet/PseudoJet.hh> 
#include "UHH2/common/include/Utils.h"

class JetCluster{

 private:

  std::vector<fastjet::PseudoJet> _hadrons;
  bool IsStableHadron(GenParticle* p);
  fastjet::PseudoJet convert_particle(GenParticle* genparticle);

 public:

  enum E_algorithm { 
    e_ca, 
    e_akt,
    e_kt,
  };

  std::vector<fastjet::PseudoJet> get_genjets(std::vector<GenParticle>* genparts, enum JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);
  std::vector<fastjet::PseudoJet> particle_in;
  fastjet::ClusterSequence* clust_seq;
  std::vector<fastjet::PseudoJet> new_jets;
  fastjet::JetDefinition *jetdef;
};
