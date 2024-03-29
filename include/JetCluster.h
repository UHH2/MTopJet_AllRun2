#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include "fastjet/AreaDefinition.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "UHH2/core/include/PFParticle.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/GenJet.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/AnalysisModule.h"

#include <vector>
#include <iostream>
#include <math.h>
#include <random>

using namespace std;
using namespace fastjet;

class JetCluster{

 private:

  std::vector<fastjet::PseudoJet> _hadrons;
  bool IsLepton(GenParticle* p);
  bool IsNeutrino(GenParticle* p);
  fastjet::PseudoJet convert_recoparticle(PFParticle* pfparticle);

 public:
  bool IsStable(GenParticle* p);
  fastjet::PseudoJet convert_particle(GenParticle* genparticle);

  enum E_algorithm {
    e_ca,
    e_akt,
    e_kt,
  };

  std::vector<fastjet::PseudoJet> get_genjets(std::vector<GenParticle>* genparts, enum JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);
  std::vector<fastjet::PseudoJet> get_recojets(std::vector<PFParticle>* pfparts, enum JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);
  std::vector<fastjet::PseudoJet> get_hotvr_jets(std::vector<GenParticle>* genparts, enum  JetCluster::E_algorithm algorithm, double rho, double min_r, double max_r, double mu, double theta, double pt_cut);
  std::vector<fastjet::PseudoJet> get_hotvr_recojets(std::vector<PFParticle>* pfparts, enum  JetCluster::E_algorithm algorithm, double rho, double min_r, double max_r, double mu, double theta, double pt_cut);
  std::vector<fastjet::PseudoJet> get_xcone_jets(std::vector<GenParticle>* genparts, int N, double R0, double beta, double ptmin);
  std::vector<fastjet::PseudoJet> get_xcone23_jets(uhh2::Event & event, uhh2::Event::Handle<vector<int>> h_list, uhh2::Event::Handle<vector<Jet>> h_particle_fatjet1, uhh2::Event::Handle<vector<Jet>> h_particle_fatjet2, uhh2::Event::Handle<vector<Jet>> h_particle_subjet1_1, uhh2::Event::Handle<vector<Jet>> h_particle_subjet1_2, uhh2::Event::Handle<vector<Jet>> h_particle_subjet1_3, uhh2::Event::Handle<vector<Jet>> h_particle_subjet2_1, uhh2::Event::Handle<vector<Jet>> h_particle_subjet2_2, uhh2::Event::Handle<vector<Jet>> h_particle_fatjet_0, uhh2::Event::Handle<vector<Jet>> h_particle_subjet1_0, uhh2::Event::Handle<vector<Jet>> h_particle_subjet2_0, uhh2::Event::Handle<vector<Jet>> h_particle_all, std::vector<GenParticle>* genparts, double ptmin, double ptmin_sub, int choose_jet);

  std::vector<fastjet::PseudoJet> get_xcone33_genjets(uhh2::Event & event, std::vector<GenParticle>* genparts, uhh2::Event::Handle<vector<Jet>>, uhh2::Event::Handle<vector<Jet>>, uhh2::Event::Handle<vector<Jet>>);
  std::vector<fastjet::PseudoJet> get_xcone33_recojets(uhh2::Event & event, std::vector<PFParticle>* pfparts, uhh2::Event::Handle<vector<Jet>>, uhh2::Event::Handle<vector<Jet>>, uhh2::Event::Handle<vector<Jet>>);

  std::vector<fastjet::PseudoJet> get_xcone_recojets(std::vector<PFParticle>* pfparts, int N, double R0, double beta, double ptmin);


  std::vector<Jet> convert_pseudojet_to_jet(std::vector<fastjet::PseudoJet> fjet);
  fastjet::PseudoJet substract_lepton(GenParticle*, fastjet::PseudoJet);
  std::vector<fastjet::PseudoJet> add_ghosts(std::vector<fastjet::PseudoJet> gen_in);

  void write_genjets(uhh2::Event & event, enum  JetCluster::E_algorithm algorithm, double jet_radius, double ptmin);

  std::vector<fastjet::PseudoJet> particle_in;
  std::vector<fastjet::PseudoJet> particle_in_subjet1;
  std::vector<fastjet::PseudoJet> particle_in_subjet2;
  std::vector<fastjet::PseudoJet> particle_in_subjet0;
  std::vector<fastjet::PseudoJet> particle_in2;
  std::vector<fastjet::PseudoJet> particle_in_noGhost;
  std::vector<fastjet::PseudoJet> particle_in_reco2;
  std::vector<fastjet::PseudoJet> particle_in_reco;
  fastjet::ClusterSequence* clust_seq;
  fastjet::ClusterSequence* clust_seq_hotvr;
  fastjet::ClusterSequence* clust_seq_reco;
  fastjet::ClusterSequence* clust_seq_xcone;
  fastjet::ClusterSequence* clust_seq_sub1;
  fastjet::ClusterSequence* clust_seq_sub2;
  std::vector<fastjet::PseudoJet> new_jets;
  std::vector<fastjet::PseudoJet> new_jets_cleaned;
  std::vector<fastjet::PseudoJet> new_recojets;
  fastjet::JetDefinition *jetdef;
  fastjet::JetDefinition *jetdef_reco;

  //double _rho, _mu, _theta, _rmin, _rmax, _ptmin,_radius,_pt_cut;
  string _clustering_algorithmus;

};


class JetProducer: public uhh2::AnalysisModule{
public:

  explicit JetProducer(uhh2::Context & , const std::string & , float, float);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_newgenjets;
  float ptmin_;
  float jet_radius_;
};


class RecoJetProducer: public uhh2::AnalysisModule{
public:

  explicit RecoJetProducer(uhh2::Context&, const std::string &, float, float);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_newrecojets;
  uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;
  float ptmin_;
  float jet_radius_;
};


class GenXCONEJetProducer: public uhh2::AnalysisModule{
public:

  explicit GenXCONEJetProducer(uhh2::Context&, const std::string &, int, double, double, double);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_newgenxconejets;
  int N_;
  double R0_;
  double beta_;
  double ptmin_;
};


class RecoXCONEJetProducer: public uhh2::AnalysisModule{
public:

  explicit RecoXCONEJetProducer(uhh2::Context&, const std::string &, int, double, double, double);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_newrecoxconejets;
  uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;
  int N_;
  double R0_;
  double beta_;
  double ptmin_;
};

class MergeXConeGen: public uhh2::AnalysisModule{
public:

  explicit MergeXConeGen(uhh2::Context&, const std::string &, const std::string &);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<TTbarGen>h_ttbargen;
  uhh2::Event::Handle<std::vector<Jet>>h_injets;
  uhh2::Event::Handle<std::vector<Jet>>h_newjets;
};

class MergeXConeReco: public uhh2::AnalysisModule{
public:

  explicit MergeXConeReco(uhh2::Context&, const std::string &, const std::string &);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_injets;
  uhh2::Event::Handle<std::vector<Jet>>h_newjets;
};

class MergeXConeN6Gen: public uhh2::AnalysisModule{
public:

  explicit MergeXConeN6Gen(uhh2::Context&, const std::string &, const std::string &);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<TTbarGen>h_ttbargen;
  uhh2::Event::Handle<std::vector<Jet>>h_injets;
  uhh2::Event::Handle<std::vector<Jet>>h_newjets;
};

class GenXCONE23JetProducer: public uhh2::AnalysisModule{
public:

  explicit GenXCONE23JetProducer(uhh2::Context&, const std::string &, const std::string &, double, double, double);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_xcone23seedjets;
  uhh2::Event::Handle<std::vector<Jet>>h_xcone23fatjets;
  uhh2::Event::Handle<std::vector<Jet>>h_xcone23subjets;
  uhh2::Event::Handle<std::vector<int>>h_particle_fatjet;
  uhh2::Event::Handle<std::vector<int>>h_particle_subjets1;
  uhh2::Event::Handle<std::vector<int>>h_particle_subjets2;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_fatjet1;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_fatjet2;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet1_1;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet1_2;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet1_3;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet2_1;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet2_2;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_fatjet0;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet1_0;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_subjet2_0;
  uhh2::Event::Handle<std::vector<Jet>>h_particle_all;



  double ptmin_;
  double ptmin_sub1_;
  double ptmin_sub2_;
};

class GenXCONE33JetProducer: public uhh2::AnalysisModule{

public:

  explicit GenXCONE33JetProducer(uhh2::Context&);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33fatjets;
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33subjets_1;
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33subjets_2;
};

class RecoXCONE33JetProducer: public uhh2::AnalysisModule{

public:

  explicit RecoXCONE33JetProducer(uhh2::Context&);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33fatjets;
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33subjets_1;
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33subjets_2;
  uhh2::Event::Handle<std::vector<PFParticle>> h_pfpart;
};


class XCone33Merge_gen: public uhh2::AnalysisModule{

public:

  explicit XCone33Merge_gen(uhh2::Context&);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33jets;
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33subjets_1;
  uhh2::Event::Handle<std::vector<Jet>>h_gen_xcone33subjets_2;
  uhh2::Event::Handle<TTbarGen>h_ttbargen;
};

class XCone33Merge_reco: public uhh2::AnalysisModule{

public:

  explicit XCone33Merge_reco(uhh2::Context&);
  virtual bool process(uhh2::Event & ) override;

private:
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33jets;
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33subjets_1;
  uhh2::Event::Handle<std::vector<Jet>>h_reco_xcone33subjets_2;
};
