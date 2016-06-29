#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include <fastjet/PseudoJet.hh> 
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/GenJetWithParts.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/common/include/FJet.h"

class JetCluster{

 private:

  std::vector<fastjet::PseudoJet> _hadrons;
  bool IsStable(GenParticle* p);
  fastjet::PseudoJet convert_particle(GenParticle* genparticle);
  fastjet::PseudoJet convert_recoparticle(PFParticle* pfparticle);

 public:

  enum E_algorithm { 
    e_ca, 
    e_akt,
    e_kt,
  };

  std::vector<fastjet::PseudoJet> get_genjets(std::vector<GenParticle>* genparts, enum JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);
  std::vector<fastjet::PseudoJet> particle_in;
  std::vector<fastjet::PseudoJet> get_recojets(std::vector<PFParticle>* pfparts, enum JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);
  std::vector<fastjet::PseudoJet> particle_in_reco;
  fastjet::ClusterSequence* clust_seq;
  fastjet::ClusterSequence* clust_seq_reco;
  std::vector<fastjet::PseudoJet> new_jets;
  std::vector<fastjet::PseudoJet> new_recojets;
  fastjet::JetDefinition *jetdef;
  fastjet::JetDefinition *jetdef_reco;
};
