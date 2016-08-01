#include "UHH2/MTopJet/include/JetCluster.h"
#include <math.h>
#include <vector>
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/MTopJet/include/MTopJetGenHists.h"
#include "UHH2/common/include/TTbarGen.h"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/VariableRPlugin.hh"
#include "UHH2/core/include/Utils.h"
#include <math.h>
#include <vector>
#include "UHH2/core/include/PFParticle.h"


#include <iostream> 


using namespace std;
using namespace fastjet;
using namespace contrib;

// Check if Gen particle is stable
bool JetCluster::IsStable(GenParticle* p)
{
	int st = p->status();
	if (st==1) return true;
	else return false;
}

// Check if Gen particle is Electron/Muon
bool JetCluster::IsLepton(GenParticle* p)
{
  int id = abs(p->pdgId());
  if (id==11 || id==13) return true;
  else return false;
}

// convert Gen Particles to PseudoJet
fastjet::PseudoJet JetCluster::convert_particle(GenParticle* genparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(genparticle->pt(),genparticle->eta(),genparticle->phi(),genparticle->energy());
  fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return gen_particle;

}

// Substract Lepton
std::vector<fastjet::PseudoJet> JetCluster::substract_lepton(GenParticle* genparticle, std::vector<fastjet::PseudoJet> fjets, double jet_radius){
  std::vector<fastjet::PseudoJet> new_fjets;
  fastjet::PseudoJet new_fjet;
  fastjet::PseudoJet fjet;
  fastjet::PseudoJet lepton = convert_particle(genparticle);
  double px, py, pz, E;
  double dR;
  int n1=0, n2=0;
  for(unsigned int i=0; i<fjets.size(); ++i){
    fjet = fjets[i];
    dR = sqrt((fjet.pseudorapidity()-lepton.pseudorapidity ())*(fjet.pseudorapidity()-lepton.pseudorapidity ()) + (fjet.phi()-lepton.phi())*(fjet.phi()-lepton.phi()));
    if(dR < jet_radius){
      // new_fjet = fjet - lepton;
      px = fjet.px()-lepton.px();
      py = fjet.py()-lepton.py();
      pz = fjet.pz()-lepton.pz();
      E = fjet.E()-lepton.E();
      ++n1;
    }
    else{
      px = fjet.px();
      py = fjet.py();
      pz = fjet.pz();
      E = fjet.E();
      ++n2;
    }
    new_fjet.reset (px, py, pz, E);
    new_fjets.push_back(new_fjet);
  }
  //cout<<"jets in event with lepton: "<<n1<<endl;
  //cout<<"jets in event without lepton: "<<n2<<endl;

  return new_fjets;
}

// Get clustered Gen Jet from Gen Particles 
std::vector<fastjet::PseudoJet> JetCluster::get_genjets(std::vector<GenParticle>* genparts, enum  JetCluster::E_algorithm algorithm, double jet_radius, double ptmin){

  fastjet::PseudoJet lepton;
  for (unsigned int i=0; i<(genparts->size()); ++i){
      GenParticle* part = &(genparts->at(i));
      if(IsStable(part)){
	  particle_in.push_back(convert_particle(part));
	  if(IsLepton(part)){
	      lepton = convert_particle(part);
	  }
	}
  }

  if(algorithm==e_ca) jetdef= new fastjet::JetDefinition(fastjet::cambridge_algorithm,jet_radius);
  if(algorithm==e_akt)jetdef= new fastjet::JetDefinition(fastjet::antikt_algorithm,jet_radius);

  clust_seq=new fastjet::ClusterSequence(particle_in, *jetdef);
  new_jets = sorted_by_pt(clust_seq->inclusive_jets(ptmin));

  // // run Jet-Lepton Cleaner
  // for (unsigned int i=0; i<(genparts->size()); ++i){
  //     GenParticle* part = &(genparts->at(i));
  //     if(IsLepton(part)){
 	
  // 	new_jets_cleaned = sorted_by_pt(substract_lepton(part, new_jets, jet_radius));
  // 	// new_jets_cleaned = sorted_by_pt(new_jets);
  //     }
  //     else{
  // 	new_jets_cleaned = sorted_by_pt(new_jets);

  //     }
  // }
  // // sort Jets by pT
  // new_jets_cleaned = sorted_by_pt(new_jets_cleaned);

  delete jetdef;
  return new_jets; //change to new_jets_cleaned if cleaning is used here (does not work now)
}

// // Get clustered HOTVR Jets from Gen Particles
// std::vector<fastjet::PseudoJet> JetCluster::get_hotvr_jets(std::vector<GenParticle>* genparts,enum  JetCluster::E_algorithm algorithm, double rho, double min_r, double max_r, double mu, double theta, double pt_cut){ 

//    for (unsigned int i=0; i<(genparts->size()); ++i){
//       GenParticle* part = &(genparts->at(i));
// 	if(IsStable(part)){
// 	  particle_in2.push_back(convert_particle(part));
// 	}
//   }

//   HOTVR::ClusterType clustertype;
//   if(algorithm==e_akt) clustertype=HOTVR::ClusterType::AKTLIKE;
//   if(algorithm==e_ca) clustertype=HOTVR::ClusterType::CALIKE;
//   if(algorithm==e_kt) clustertype=HOTVR::ClusterType::KTLIKE;
//   HOTVR plugin_hotvr(mu, theta, min_r, max_r, rho, pt_cut, clustertype);  //call HOTVR algorithm
//   fastjet::JetDefinition jet_def(&plugin_hotvr);
//   clust_seq2=new fastjet::ClusterSequence(particle_in2, jet_def);
    
    
//   std::vector<fastjet::PseudoJet> hotvr_jets,rejected_jets,soft_jets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets
	
//   //get vector from the plugin
//   hotvr_jets=plugin_hotvr.get_jets();
//   rejected_jets=plugin_hotvr.get_rejected_cluster();
//   soft_jets=plugin_hotvr.get_soft_cluster();
//   return hotvr_jets;
// }

// convert PF Particles to PseudoJet
fastjet::PseudoJet JetCluster::convert_recoparticle(PFParticle* pfparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(pfparticle->pt(),pfparticle->eta(),pfparticle->phi(),pfparticle->energy());
  fastjet::PseudoJet pf_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return pf_particle;
}

// Get clustered Reco Jet from PF Particles 
std::vector<fastjet::PseudoJet> JetCluster::get_recojets(std::vector<PFParticle>* pfparts, enum  JetCluster::E_algorithm algorithm, double jet_radius, double ptmin){

  for (unsigned int i=0; i<(pfparts->size()); ++i){
      PFParticle* part = &(pfparts->at(i));
      particle_in_reco.push_back(convert_recoparticle(part));
  }

  if(algorithm==e_ca) jetdef_reco= new fastjet::JetDefinition(fastjet::cambridge_algorithm,jet_radius);
  if(algorithm==e_akt)jetdef_reco= new fastjet::JetDefinition(fastjet::antikt_algorithm,jet_radius);

  clust_seq_reco=new fastjet::ClusterSequence(particle_in_reco, *jetdef_reco);

  new_recojets = sorted_by_pt(clust_seq_reco->inclusive_jets(ptmin));

  // new_recojets = particle_in_reco;

  delete jetdef;
  return new_recojets;

}

// convert Pseudo Jet to Jet
std::vector<Jet> JetCluster::convert_pseudojet_to_jet(std::vector<fastjet::PseudoJet> fjet){
  vector<Jet> new_jets;
  Jet new_jet;
  // TLorentzVector jet;
  for(unsigned i=0; i<fjet.size(); ++i){
    new_jet.set_pt(fjet[i].pt());
    new_jet.set_eta(fjet[i].eta());
    new_jet.set_phi(fjet[i].phi());
    new_jet.set_energy(fjet[i].E());
    new_jets.push_back(new_jet);
  }
  return new_jets;
}


JetProducer::JetProducer(uhh2::Context & ctx, const std::string & name, float ptmin, float jet_radius):
  h_newgenjets(ctx.declare_event_output<std::vector<Jet>>(name)),
  ptmin_(ptmin),
  jet_radius_(jet_radius) {}

bool JetProducer::process(uhh2::Event & event){
  
  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc = new JetCluster();
  std::vector<fastjet::PseudoJet> gen;
  gen = jetc->get_genjets(genparts, JetCluster::e_akt, jet_radius_, ptmin_); 
  // cout<<gen.size()<<endl;

  event.set(h_newgenjets, jetc->convert_pseudojet_to_jet(gen));

  delete jetc;
  return true;
}

RecoJetProducer::RecoJetProducer(uhh2::Context & ctx, const std::string & name, float ptmin, float jet_radius):
  h_newrecojets(ctx.declare_event_output<std::vector<Jet>>(name)),
  h_pfpart(ctx.get_handle<vector<PFParticle>>("PFParticles")),
  ptmin_(ptmin),
  jet_radius_(jet_radius) {}

bool RecoJetProducer::process(uhh2::Event & event){

  std::vector<PFParticle> pfparts = event.get(h_pfpart);
  JetCluster* jetc_reco = new JetCluster();
  std::vector<fastjet::PseudoJet> reco;
  reco = jetc_reco->get_recojets(&pfparts, JetCluster::e_akt, jet_radius_, ptmin_); 
  event.set(h_newrecojets, jetc_reco->convert_pseudojet_to_jet(reco));

  delete jetc_reco;
  return true;
}
