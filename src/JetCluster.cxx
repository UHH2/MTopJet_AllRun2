#include "UHH2/MTopJet/include/JetCluster.h"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/XConePlugin.hh"
// #include "fastjet/contrib/ClusteringVetoPlugin.hh"

using namespace std;
using namespace fastjet;
using namespace contrib;



// --------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------ Jet Clustering ----------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------

// ------------------------------------ Get clustered Reco Jet from PF Particles ---------------------------------------------------------------------------------
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

// ------------------------------------ Get clustered Gen Jet from Gen Particles ---------------------------------------------------------------------------------
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

   delete jetdef;
  return new_jets; 
}


// ------------------------------------ Get clustered HOTVR Jets from Gen Particles ---------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> JetCluster::get_hotvr_jets(std::vector<GenParticle>* genparts,enum  JetCluster::E_algorithm algorithm, double rho, double min_r, double max_r, double mu, double theta, double pt_cut){ 

   for (unsigned int i=0; i<(genparts->size()); ++i){
      GenParticle* part = &(genparts->at(i));
	if(IsStable(part)){
	  particle_in2.push_back(convert_particle(part));
	}
  }

  HOTVR::ClusterType clustertype;
  if(algorithm==e_akt) clustertype=HOTVR::ClusterType::AKTLIKE;
  if(algorithm==e_ca) clustertype=HOTVR::ClusterType::CALIKE;
  if(algorithm==e_kt) clustertype=HOTVR::ClusterType::KTLIKE;
  HOTVR plugin_hotvr(mu, theta, min_r, max_r, rho, pt_cut, clustertype);  //call HOTVR algorithm
  fastjet::JetDefinition jet_def_hotvr(&plugin_hotvr);
  clust_seq_hotvr=new fastjet::ClusterSequence(particle_in2, jet_def_hotvr);

  std::vector<fastjet::PseudoJet> hotvr_jets,rejected_jets,soft_jets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets

  //get vector from the plugin
  hotvr_jets=plugin_hotvr.get_jets();
  rejected_jets=plugin_hotvr.get_rejected_cluster();
  soft_jets=plugin_hotvr.get_soft_cluster();
  return hotvr_jets;
}


// ------------------------------------ Get clustered HOTVR Jets from PFParticles ---------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> JetCluster::get_hotvr_recojets(std::vector<PFParticle>* pfparts,enum  JetCluster::E_algorithm algorithm, double rho, double min_r, double max_r, double mu, double theta, double pt_cut){ 

 for (unsigned int i=0; i<(pfparts->size()); ++i){
      PFParticle* part = &(pfparts->at(i));
      particle_in_reco2.push_back(convert_recoparticle(part));
  }

  HOTVR::ClusterType clustertype;
  if(algorithm==e_akt) clustertype=HOTVR::ClusterType::AKTLIKE;
  if(algorithm==e_ca) clustertype=HOTVR::ClusterType::CALIKE;
  if(algorithm==e_kt) clustertype=HOTVR::ClusterType::KTLIKE;
  HOTVR plugin_hotvr(mu, theta, min_r, max_r, rho, pt_cut, clustertype);  //call HOTVR algorithm
  fastjet::JetDefinition jet_def_hotvr(&plugin_hotvr);
  clust_seq_hotvr=new fastjet::ClusterSequence(particle_in_reco2, jet_def_hotvr);

  std::vector<fastjet::PseudoJet> hotvr_recojets,rejected_recojets,soft_recojets ; //vector of hotvr_jets, jets that were rejcted durning the clustering procedure and soft jets
	
  //get vector from the plugin
  hotvr_recojets=plugin_hotvr.get_jets();
  rejected_recojets=plugin_hotvr.get_rejected_cluster();
  soft_recojets=plugin_hotvr.get_soft_cluster();
  return hotvr_recojets;
}

// ------------------------------------ Get clustered XCone Jets from Gen Particles ---------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> JetCluster::get_xcone_jets(std::vector<GenParticle>* genparts, int N, double R0, double beta, double ptmin){ 

  std::vector<fastjet::PseudoJet> xcone_jets; 
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle* part = &(genparts->at(i));
    if(IsStable(part)){
      particle_in2.push_back(convert_particle(part));
    }
  }
  XConePlugin plugin_xcone(N, R0, beta);
  fastjet::JetDefinition jet_def_xcone(&plugin_xcone);
  clust_seq_xcone=new fastjet::ClusterSequence(particle_in2, jet_def_xcone);
  xcone_jets = sorted_by_pt(clust_seq_xcone->inclusive_jets(ptmin));



  return xcone_jets;
}
// ------------------------------------ Get clustered XCone Jets from PFParticles ---------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> JetCluster::get_xcone_recojets(std::vector<PFParticle>* pfparts, int N, double R0, double beta, double ptmin){ 
 
 std::vector<fastjet::PseudoJet> xcone_recojets; 
 for (unsigned int i=0; i<(pfparts->size()); ++i){
      PFParticle* part = &(pfparts->at(i));
      particle_in_reco2.push_back(convert_recoparticle(part));
  }
  XConePlugin plugin_xcone(N, R0, beta);
  fastjet::JetDefinition jet_def_xcone(&plugin_xcone);
  clust_seq_xcone=new fastjet::ClusterSequence(particle_in_reco2, jet_def_xcone);
  xcone_recojets = sorted_by_pt(clust_seq_xcone->inclusive_jets(ptmin));

  return xcone_recojets;
  
}
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------



// --------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------ additional functions ----------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------

// ------------------------------------ Check if Gen particle is stable ---------------------------------------------------------------------------------
bool JetCluster::IsStable(GenParticle* p)
{
	int st = p->status();
	if (st==1) return true;
	else return false;
}

// ------------------------------------ Check if Gen particle is Electron/Muon ---------------------------------------------------------------------------------
bool JetCluster::IsLepton(GenParticle* p)
{
  int id = abs(p->pdgId());
  if (id==11 || id==13) return true;
  else return false;
}

// ------------------------------------ convert Gen Particles to PseudoJet ---------------------------------------------------------------------------------
fastjet::PseudoJet JetCluster::convert_particle(GenParticle* genparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(genparticle->pt(),genparticle->eta(),genparticle->phi(),genparticle->energy());
  fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return gen_particle;

}

// ------------------------------------ convert PF Particles to PseudoJet ---------------------------------------------------------------------------------
fastjet::PseudoJet JetCluster::convert_recoparticle(PFParticle* pfparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(pfparticle->pt(),pfparticle->eta(),pfparticle->phi(),pfparticle->energy());
  fastjet::PseudoJet pf_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return pf_particle;
}

// ------------------------------------ convert Pseudo Jet to Jet ---------------------------------------------------------------------------------
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

// ------------------------------------ Substract Lepton ---------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------


// --------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------ Jet Producer ------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------


// ------------------------------------ Gen Jets ----------------------------------------------------------------------------------
JetProducer::JetProducer(uhh2::Context & ctx, const std::string & name, float ptmin, float jet_radius):
  h_newgenjets(ctx.declare_event_output<std::vector<Jet>>(name)),
  ptmin_(ptmin),
  jet_radius_(jet_radius) {}

bool JetProducer::process(uhh2::Event & event){
  
  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc = new JetCluster();
  std::vector<fastjet::PseudoJet> gen;
  gen = jetc->get_genjets(genparts, JetCluster::e_akt, jet_radius_, ptmin_); 

  event.set(h_newgenjets, jetc->convert_pseudojet_to_jet(gen));

  delete jetc;
  return true;
}

// ------------------------------------ Reco Jets ---------------------------------------------------------------------------------
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
// ------------------------------------ HOTVR Reco Jets -----------------------------------------------------------------------------
RecoHOTVRJetProducer::RecoHOTVRJetProducer(uhh2::Context & ctx, const std::string & name, double rho):
  h_newrecohotvrjets(ctx.declare_event_output<std::vector<Jet>>(name)),
  h_pfpart(ctx.get_handle<vector<PFParticle>>("PFParticles")),
  rho_(rho) {}

bool RecoHOTVRJetProducer::process(uhh2::Event & event){
  // standard values for HOTVR:
  double mu(30.),     // massjump threshold
    theta(0.7),       // massjump parameter
    max_r(1.5),       // maximum allowed distance R
    min_r(0.1),       // minimum allowed distance R
    // rho(600),         // cone shrinking parameter
    pt_cut(30.);      // minimum pT of subjets

  std::vector<PFParticle> pfparts = event.get(h_pfpart);
  JetCluster* jetc_reco = new JetCluster();
  std::vector<fastjet::PseudoJet> reco;

  reco = jetc_reco->get_hotvr_recojets(&pfparts, JetCluster::e_ca, rho_, min_r, max_r, mu, theta, pt_cut);

  event.set(h_newrecohotvrjets, jetc_reco->convert_pseudojet_to_jet(reco));

  delete jetc_reco;
  return true;
}


// ------------------------------------ HOTVR Gen Jets -----------------------------------------------------------------------------
GenHOTVRJetProducer::GenHOTVRJetProducer(uhh2::Context & ctx, const std::string & name, double rho):
  h_newgenhotvrjets(ctx.declare_event_output<std::vector<Jet>>(name)),
  rho_(rho) {}

bool GenHOTVRJetProducer::process(uhh2::Event & event){
  // standard values for HOTVR:
  double mu(30.),     // massjump threshold
    theta(0.7),       // massjump parameter
    max_r(1.5),       // maximum allowed distance R
    min_r(0.1),       // minimum allowed distance R
    // rho(600),         // cone shrinking parameter
    pt_cut(30.);      // minimum pT of subjets

  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc = new JetCluster();
  std::vector<fastjet::PseudoJet> gen;
  gen = jetc->get_hotvr_jets(genparts, JetCluster::e_ca, rho_, min_r, max_r, mu, theta, pt_cut);

  event.set(h_newgenhotvrjets, jetc->convert_pseudojet_to_jet(gen));

  delete jetc;
  return true;
}

// ------------------------------------ XCone Reco Jets ----------------------------------------------------------------------------
RecoXCONEJetProducer::RecoXCONEJetProducer(uhh2::Context & ctx, const std::string & name, int N, double R0, double beta, double ptmin):
  h_newrecoxconejets(ctx.declare_event_output<std::vector<Jet>>(name)),
  h_pfpart(ctx.get_handle<vector<PFParticle>>("PFParticles")),
  N_(N),
  R0_(R0),
  beta_(beta),
  ptmin_(ptmin) {}

bool RecoXCONEJetProducer::process(uhh2::Event & event){


  std::vector<PFParticle> pfparts = event.get(h_pfpart);
  JetCluster* jetc_reco = new JetCluster();
  std::vector<fastjet::PseudoJet> reco;
  reco = jetc_reco->get_xcone_recojets(&pfparts, N_, R0_, beta_, ptmin_);

  event.set(h_newrecoxconejets, jetc_reco->convert_pseudojet_to_jet(reco));

  delete jetc_reco;
  return true;
}

// ------------------------------------ XCone Gen Jets -----------------------------------------------------------------------------
GenXCONEJetProducer::GenXCONEJetProducer(uhh2::Context & ctx, const std::string & name, int N, double R0, double beta, double ptmin):
  h_newgenxconejets(ctx.declare_event_output<std::vector<Jet>>(name)),
  N_(N),
  R0_(R0),
  beta_(beta),
  ptmin_(ptmin) {}

bool GenXCONEJetProducer::process(uhh2::Event & event){


  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc = new JetCluster();
  std::vector<fastjet::PseudoJet> gen;
  gen = jetc->get_xcone_jets(genparts, N_, R0_, beta_, ptmin_);

  event.set(h_newgenxconejets, jetc->convert_pseudojet_to_jet(gen));

  delete jetc;
  return true;
}

// ------------------------------------ merging XCone Gen Jets -------------------------------------------------------------------
MergeXConeGen::MergeXConeGen(uhh2::Context & ctx, const std::string & in_jet, const std::string & out_jet):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  h_injets(ctx.get_handle<std::vector<Jet>>(in_jet)),
  h_newjets(ctx.declare_event_output<std::vector<Jet>>(out_jet)) {}

bool MergeXConeGen::process(uhh2::Event & event){

  // get Lepton from TTbar
  const auto & ttbargen = event.get(h_ttbargen);
  GenParticle lep1, lep2, lepton;
  if(ttbargen.IsTopHadronicDecay()){
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
  }

  // get distance to jets
  std::vector<Jet> jets = event.get(h_injets);
  TLorentzVector jet1_v4, jet2_v4, jet3_v4, jet4_v4, combinedjet1_v4, combinedjet2_v4;
  Jet jet1;
  jet1 = jets.at(0);
  float dR1, dR2, dR3, dR4;
  if(jets.size() > 3){
    jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
    jet2_v4.SetPxPyPzE(jets.at(1).v4().Px(), jets.at(1).v4().Py(), jets.at(1).v4().Pz(), jets.at(1).v4().E()); 
    jet3_v4.SetPxPyPzE(jets.at(2).v4().Px(), jets.at(2).v4().Py(), jets.at(2).v4().Pz(), jets.at(2).v4().E()); 
    jet4_v4.SetPxPyPzE(jets.at(3).v4().Px(), jets.at(3).v4().Py(), jets.at(3).v4().Pz(), jets.at(3).v4().E()); 
  
    // claculate distance to Lepton
    dR1 = uhh2::deltaR(jet1, lepton);
    dR2 = uhh2::deltaR(jets.at(1), lepton);
    dR3 = uhh2::deltaR(jets.at(2), lepton);
    dR4 = uhh2::deltaR(jets.at(3), lepton);
    if(dR1 < dR2 && dR1 < dR3 && dR1 < dR4){
      combinedjet1_v4 = jet2_v4 + jet3_v4 + jet4_v4;
      combinedjet2_v4 = jet1_v4;
    }
    if(dR2 < dR1 && dR2 < dR3 && dR2 < dR4){
      combinedjet1_v4 = jet1_v4 + jet3_v4 + jet4_v4;
      combinedjet2_v4 = jet2_v4;
    }    
    if(dR3 < dR1 && dR3 < dR2 && dR3 < dR4){
      combinedjet1_v4 = jet1_v4 + jet2_v4 + jet4_v4;
      combinedjet2_v4 = jet3_v4;
    }
    if(dR4 < dR1 && dR4 < dR2 && dR4 < dR3){
      combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4;
      combinedjet2_v4 = jet4_v4;
    }
  }

  std::vector<Jet> new_jets;
  Jet new_jet1, new_jet2;
  // TLorentzVector jet;
  new_jet1.set_pt(combinedjet1_v4.Pt());
  new_jet1.set_eta(combinedjet1_v4.Eta());
  new_jet1.set_phi(combinedjet1_v4.Phi());
  new_jet1.set_energy(combinedjet1_v4.E());
  new_jet2.set_pt(combinedjet2_v4.Pt());
  new_jet2.set_eta(combinedjet2_v4.Eta());
  new_jet2.set_phi(combinedjet2_v4.Phi());
  new_jet2.set_energy(combinedjet2_v4.E());

  new_jets.push_back(new_jet1);
  new_jets.push_back(new_jet2);
  event.set(h_newjets, new_jets);

  return true;
}

// ------------------------------------ merging XCone Reco Jets -------------------------------------------------------------------
MergeXConeReco::MergeXConeReco(uhh2::Context & ctx, const std::string & in_jet, const std::string & out_jet):
  h_injets(ctx.get_handle<std::vector<Jet>>(in_jet)),
  h_newjets(ctx.declare_event_output<std::vector<Jet>>(out_jet)) {}

bool MergeXConeReco::process(uhh2::Event & event){

  // get Reco Lepton
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }

  // get distance to jets
  std::vector<Jet> jets = event.get(h_injets);
  TLorentzVector jet1_v4, jet2_v4, jet3_v4, jet4_v4, combinedjet1_v4, combinedjet2_v4;
  Jet jet1;
  jet1 = jets.at(0);
  float dR1, dR2, dR3, dR4;
  if(jets.size() > 3){
    jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
    jet2_v4.SetPxPyPzE(jets.at(1).v4().Px(), jets.at(1).v4().Py(), jets.at(1).v4().Pz(), jets.at(1).v4().E()); 
    jet3_v4.SetPxPyPzE(jets.at(2).v4().Px(), jets.at(2).v4().Py(), jets.at(2).v4().Pz(), jets.at(2).v4().E()); 
    jet4_v4.SetPxPyPzE(jets.at(3).v4().Px(), jets.at(3).v4().Py(), jets.at(3).v4().Pz(), jets.at(3).v4().E()); 
  
    // claculate distance to Lepton
    dR1 = uhh2::deltaR(jet1, lepton);
    dR2 = uhh2::deltaR(jets.at(1), lepton);
    dR3 = uhh2::deltaR(jets.at(2), lepton);
    dR4 = uhh2::deltaR(jets.at(3), lepton);
    if(dR1 < dR2 && dR1 < dR3 && dR1 < dR4){
      combinedjet1_v4 = jet2_v4 + jet3_v4 + jet4_v4;
      combinedjet2_v4 = jet1_v4;
    }
    if(dR2 < dR1 && dR2 < dR3 && dR2 < dR4){
      combinedjet1_v4 = jet1_v4 + jet3_v4 + jet4_v4;
      combinedjet2_v4 = jet2_v4;
    }    
    if(dR3 < dR1 && dR3 < dR2 && dR3 < dR4){
      combinedjet1_v4 = jet1_v4 + jet2_v4 + jet4_v4;
      combinedjet2_v4 = jet3_v4;
    }
    if(dR4 < dR1 && dR4 < dR2 && dR4 < dR3){
      combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4;
      combinedjet2_v4 = jet4_v4;
    }
  }

  std::vector<Jet> new_jets;
  Jet new_jet1, new_jet2;
  // TLorentzVector jet;
  new_jet1.set_pt(combinedjet1_v4.Pt());
  new_jet1.set_eta(combinedjet1_v4.Eta());
  new_jet1.set_phi(combinedjet1_v4.Phi());
  new_jet1.set_energy(combinedjet1_v4.E());
  new_jet2.set_pt(combinedjet2_v4.Pt());
  new_jet2.set_eta(combinedjet2_v4.Eta());
  new_jet2.set_phi(combinedjet2_v4.Phi());
  new_jet2.set_energy(combinedjet2_v4.E());

  new_jets.push_back(new_jet1);
  new_jets.push_back(new_jet2);
  event.set(h_newjets, new_jets);

  return true;
}

// ------------------------------------ merging XCone Gen Jets N=6  -------------------------------------------------------------------
MergeXConeN6Gen::MergeXConeN6Gen(uhh2::Context & ctx, const std::string & in_jet, const std::string & out_jet):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  h_injets(ctx.get_handle<std::vector<Jet>>(in_jet)),
  h_newjets(ctx.declare_event_output<std::vector<Jet>>(out_jet)) {}

bool MergeXConeN6Gen::process(uhh2::Event & event){

  // get Lepton from TTbar
  const auto & ttbargen = event.get(h_ttbargen);
  GenParticle lep1, lep2, lepton;
  if(ttbargen.IsTopHadronicDecay()){
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
  }

  // get distance to jets
  std::vector<Jet> jets = event.get(h_injets);
  TLorentzVector jet1_v4, jet2_v4, jet3_v4, jet4_v4, jet5_v4, combinedjet1_v4, combinedjet2_v4, combinedjet3_v4;
  Jet jet1;
  jet1 = jets.at(0);
  float dR1, dR2, dR3, dR4, dR5;
  if(jets.size() > 4){
    jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
    jet2_v4.SetPxPyPzE(jets.at(1).v4().Px(), jets.at(1).v4().Py(), jets.at(1).v4().Pz(), jets.at(1).v4().E()); 
    jet3_v4.SetPxPyPzE(jets.at(2).v4().Px(), jets.at(2).v4().Py(), jets.at(2).v4().Pz(), jets.at(2).v4().E()); 
    jet4_v4.SetPxPyPzE(jets.at(3).v4().Px(), jets.at(3).v4().Py(), jets.at(3).v4().Pz(), jets.at(3).v4().E()); 
    jet5_v4.SetPxPyPzE(jets.at(4).v4().Px(), jets.at(4).v4().Py(), jets.at(4).v4().Pz(), jets.at(4).v4().E()); 

    // claculate distance to Lepton
    dR1 = uhh2::deltaR(jet1, lepton);
    dR2 = uhh2::deltaR(jets.at(1), lepton);
    dR3 = uhh2::deltaR(jets.at(2), lepton);
    dR4 = uhh2::deltaR(jets.at(3), lepton);
    dR5 = uhh2::deltaR(jets.at(4), lepton);

    // calculate distance between jets;
    float dR12, dR13, dR14, dR15, dR23, dR24, dR25, dR34, dR35, dR45;
    dR12 = uhh2::deltaR(jets.at(0), jets.at(1));
    dR13 = uhh2::deltaR(jets.at(0), jets.at(2));
    dR14 = uhh2::deltaR(jets.at(0), jets.at(3));
    dR15 = uhh2::deltaR(jets.at(0), jets.at(4));
    dR23 = uhh2::deltaR(jets.at(1), jets.at(2));
    dR24 = uhh2::deltaR(jets.at(1), jets.at(3));
    dR25 = uhh2::deltaR(jets.at(1), jets.at(4));
    dR34 = uhh2::deltaR(jets.at(2), jets.at(3));
    dR35 = uhh2::deltaR(jets.at(2), jets.at(4));
    dR45 = uhh2::deltaR(jets.at(3), jets.at(4));


    // first get nearest jet to lep -> combined jet 2 / if one jet has big distance to all other 3 jets (3 dRs of this jet > other dRs) -> combine other 3 jets / else combine 4 jets
    if(dR1 < dR2 && dR1 < dR3 && dR1 < dR4 && dR1 < dR5){
      combinedjet2_v4 = jet1_v4;
      if((dR23 > dR34 && dR23 > dR35 && dR23 > dR45) && (dR24 > dR34 && dR24 > dR35 && dR24 > dR45) && (dR25 > dR34 && dR25 > dR35 && dR25 > dR45)){
	combinedjet1_v4 = jet3_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet2_v4;
      }
      else if((dR23 > dR24 && dR23 > dR25 && dR23 > dR45) && (dR34 > dR24 && dR34 > dR25 && dR34 > dR45) && (dR35 > dR24 && dR35 > dR25 && dR35 > dR45)){
	combinedjet1_v4 = jet2_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet3_v4;
      }	
      else if((dR24 > dR23 && dR24 > dR25 && dR24 > dR35) && (dR34 > dR23 && dR34 > dR25 && dR34 > dR35) && (dR45 > dR23 && dR45 > dR25 && dR45 > dR35)){
	combinedjet1_v4 = jet2_v4 + jet3_v4 + jet5_v4;
	combinedjet3_v4 = jet4_v4;
      }
      else if((dR25 > dR23 && dR25 > dR24 && dR25 > dR34) && (dR35 > dR23 && dR35 > dR24 && dR35 > dR34) && (dR45 > dR23 && dR45 > dR24 && dR45 > dR34)){
	combinedjet1_v4 = jet2_v4 + jet3_v4 + jet4_v4;
	combinedjet3_v4 = jet5_v4;
      }
      else{
	combinedjet1_v4 = jet2_v4 + jet3_v4 + jet4_v4 + jet5_v4;
      }
    }

    if(dR2 < dR1 && dR2 < dR3 && dR2 < dR4 && dR2 < dR5){
      combinedjet2_v4 = jet2_v4;
      if((dR13 > dR34 && dR13 > dR35 && dR13 > dR45) && (dR14 > dR34 && dR14 > dR35 && dR14 > dR45) && (dR15 > dR34 && dR15 > dR35 && dR15 > dR45)){
	combinedjet1_v4 = jet3_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet1_v4;
      }
      else if((dR13 > dR14 && dR13 > dR15 && dR13 > dR45) && (dR34 > dR14 && dR34 > dR15 && dR34 > dR45) && (dR35 > dR14 && dR35 > dR15 && dR35 > dR45)){
	combinedjet1_v4 = jet1_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet3_v4;
      }	
      else if((dR14 > dR13 && dR14 > dR15 && dR14 > dR35) && (dR34 > dR13 && dR34 > dR15 && dR34 > dR35) && (dR45 > dR13 && dR45 > dR15 && dR45 > dR35)){
	combinedjet1_v4 = jet1_v4 + jet3_v4 + jet5_v4;
	combinedjet3_v4 = jet4_v4;
      }
      else if((dR15 > dR13 && dR15 > dR14 && dR15 > dR34) && (dR35 > dR13 && dR35 > dR14 && dR35 > dR34) && (dR45 > dR13 && dR45 > dR14 && dR45 > dR34)){
	combinedjet1_v4 = jet1_v4 + jet3_v4 + jet4_v4;
	combinedjet3_v4 = jet5_v4;
      }
      else{
	combinedjet1_v4 = jet1_v4 + jet3_v4 + jet4_v4 + jet5_v4;
      }
    }    
    if(dR3 < dR1 && dR3 < dR2 && dR3 < dR4 && dR3 < dR5){
      combinedjet2_v4 = jet3_v4;
      if((dR12 > dR24 && dR12 > dR25 && dR12 > dR45) && (dR14 > dR24 && dR14 > dR25 && dR14 > dR45) && (dR15 > dR24 && dR15 > dR25 && dR15 > dR45)){
	combinedjet1_v4 = jet2_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet1_v4;
      }
      else if((dR12 > dR14 && dR12 > dR15 && dR12 > dR45) && (dR24 > dR14 && dR24 > dR15 && dR24 > dR45) && (dR25 > dR14 && dR25 > dR15 && dR25 > dR45)){
	combinedjet1_v4 = jet1_v4 + jet4_v4 + jet5_v4;
	combinedjet3_v4 = jet2_v4;
      }	
      else if((dR14 > dR12 && dR14 > dR15 && dR14 > dR25) && (dR24 > dR12 && dR24 > dR15 && dR24 > dR25) && (dR45 > dR12 && dR45 > dR15 && dR45 > dR25)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet5_v4;
	combinedjet3_v4 = jet4_v4;
      }
      else if((dR15 > dR12 && dR15 > dR14 && dR15 > dR24) && (dR25 > dR12 && dR25 > dR14 && dR25 > dR24) && (dR45 > dR12 && dR45 > dR14 && dR45 > dR24)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet4_v4;
	combinedjet3_v4 = jet5_v4;
      }
      else{
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet4_v4 + jet5_v4;
      }
    }
    if(dR4 < dR1 && dR4 < dR2 && dR4 < dR3 && dR4 < dR5){
      combinedjet2_v4 = jet4_v4;
      if((dR12 > dR23 && dR12 > dR25 && dR12 > dR35) && (dR13 > dR23 && dR13 > dR25 && dR13 > dR35) && (dR15 > dR23 && dR15 > dR25 && dR15 > dR35)){
	combinedjet1_v4 = jet2_v4 + jet3_v4 + jet5_v4;
	combinedjet3_v4 = jet1_v4;
      }
      else if((dR12 > dR13 && dR12 > dR15 && dR12 > dR35) && (dR23 > dR13 && dR23 > dR15 && dR23 > dR35) && (dR25 > dR13 && dR25 > dR15 && dR25 > dR35)){
	combinedjet1_v4 = jet1_v4 + jet3_v4 + jet5_v4;
	combinedjet3_v4 = jet2_v4;
      }	
      else if((dR13 > dR12 && dR13 > dR15 && dR13 > dR25) && (dR23 > dR12 && dR23 > dR15 && dR23 > dR25) && (dR35 > dR12 && dR35 > dR15 && dR35 > dR25)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet5_v4;
	combinedjet3_v4 = jet3_v4;
      }
      else if((dR15 > dR12 && dR15 > dR13 && dR15 > dR23) && (dR25 > dR12 && dR25 > dR13 && dR25 > dR23) && (dR35 > dR12 && dR35 > dR13 && dR35 > dR23)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4;
	combinedjet3_v4 = jet5_v4;
      }
      else{
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4 + jet5_v4;
      }
    }
    if(dR5 < dR1 && dR5 < dR2 && dR5 < dR3 && dR5 < dR4){
      combinedjet2_v4 = jet5_v4;
      if((dR12 > dR23 && dR12 > dR24 && dR12 > dR34) && (dR13 > dR23 && dR13 > dR24 && dR13 > dR34) && (dR14 > dR23 && dR14 > dR24 && dR14 > dR34)){
	combinedjet1_v4 = jet2_v4 + jet3_v4 + jet4_v4;
	combinedjet3_v4 = jet1_v4;
      }
      else if((dR12 > dR13 && dR12 > dR14 && dR12 > dR34) && (dR23 > dR13 && dR23 > dR14 && dR23 > dR34) && (dR24 > dR13 && dR24 > dR14 && dR24 > dR34)){
	combinedjet1_v4 = jet1_v4 + jet3_v4 + jet4_v4;
	combinedjet3_v4 = jet2_v4;
      }	
      else if((dR13 > dR12 && dR13 > dR14 && dR13 > dR24) && (dR23 > dR12 && dR23 > dR14 && dR23 > dR24) && (dR34 > dR12 && dR34 > dR14 && dR34 > dR24)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet4_v4;
	combinedjet3_v4 = jet3_v4;
      }
      else if((dR14 > dR12 && dR14 > dR13 && dR14 > dR23) && (dR24 > dR12 && dR24 > dR13 && dR24 > dR23) && (dR34 > dR12 && dR34 > dR13 && dR34 > dR23)){
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4;
	combinedjet3_v4 = jet4_v4;
      }
     else{
	combinedjet1_v4 = jet1_v4 + jet2_v4 + jet3_v4 + jet4_v4;
      }
    }
  }

  std::vector<Jet> new_jets;
  Jet new_jet1, new_jet2, new_jet3;
  // TLorentzVector jet;
  new_jet1.set_pt(combinedjet1_v4.Pt());
  new_jet1.set_eta(combinedjet1_v4.Eta());
  new_jet1.set_phi(combinedjet1_v4.Phi());
  new_jet1.set_energy(combinedjet1_v4.E());
  new_jet2.set_pt(combinedjet2_v4.Pt());
  new_jet2.set_eta(combinedjet2_v4.Eta());
  new_jet2.set_phi(combinedjet2_v4.Phi());
  new_jet2.set_energy(combinedjet2_v4.E());
  new_jet3.set_pt(combinedjet3_v4.Pt());
  new_jet3.set_eta(combinedjet3_v4.Eta());
  new_jet3.set_phi(combinedjet3_v4.Phi());
  new_jet3.set_energy(combinedjet3_v4.E());
  new_jets.push_back(new_jet1);
  new_jets.push_back(new_jet2);
  new_jets.push_back(new_jet3);
  event.set(h_newjets, new_jets);

  return true;
}
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
