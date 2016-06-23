#include "UHH2/MTopJet/include/JetCluster.h"
#include <math.h>
#include <vector>
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/MTopJet/include/MTopJetGenHists.h"
//#include "UHH2/MTopJet/include/GenJetProps.h"
#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>


#include <iostream> 


bool JetCluster::IsStableHadron(GenParticle* p)
{
	int st = p->status();
	if (st==1) return true;
	else return false;
}

fastjet::PseudoJet JetCluster::convert_particle(GenParticle* genparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(genparticle->pt(),genparticle->eta(),genparticle->phi(),genparticle->energy());
  fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return gen_particle;

}

std::vector<fastjet::PseudoJet> JetCluster::get_genjets(std::vector<GenParticle>* genparts, enum  JetCluster::E_algorithm algorithm, double jet_radius, double ptmin){

  for (unsigned int i=0; i<(genparts->size()); ++i){
      GenParticle* part = &(genparts->at(i));
      if (IsStableHadron(part)){
      	  particle_in.push_back(convert_particle(part));
      	  continue;
      	}
  }

  if(algorithm==e_ca) jetdef= new fastjet::JetDefinition(fastjet::cambridge_algorithm,jet_radius);
  if(algorithm==e_akt)jetdef= new fastjet::JetDefinition(fastjet::antikt_algorithm,jet_radius);

  clust_seq=new fastjet::ClusterSequence(particle_in, *jetdef);

  new_jets = sorted_by_pt(clust_seq->inclusive_jets(ptmin));

  delete jetdef;
  return new_jets;

}

