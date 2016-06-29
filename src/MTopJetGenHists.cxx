#include "UHH2/MTopJet/include/MTopJetGenHists.h"
#include "UHH2/MTopJet/include/JetCluster.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/PFParticle.h"

#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

MTopJetGenHists::MTopJetGenHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  
  RecoPT_old = book<TH1F>("p_{T} old ak8 Jets", "p_{T} [GeV/c]", 20, 0, 1000);
  RecoNumber_old = book<TH1F>("Number old ak8 Jets ", "number", 15, 0, 15);
  RecoEta_old = book<TH1F>("#eta old ak8 Jets", "#eta", 20, -3, 3);
  RecoJet1Mass_old = book<TH1F>("Jet Mass old ak8 Jet", "M_{jet}", 15, 0, 300);

  RecoEvNumber = book<TH1F>("Number Reco Parts per Event", "number", 30, 0, 300);

  RecoPT_ak08 = book<TH1F>("p_{T} AK8 Jets", "p_{T} [GeV/c]", 20, 0, 1000);
  RecoNumber_ak08 = book<TH1F>("Number AK8 Jets", "number", 15, 0, 15);
  RecoEta_ak08 = book<TH1F>("#eta AK8 Jets", "#eta", 20, -3, 3);
  RecoJet1Mass_ak08 = book<TH1F>("Jet Mass AK8", "M_{jet}", 15, 0, 300);

  GenEvNumber = book<TH1F>("Number Gen Parts per Event", "number", 30, 0, 300);

  GenPT_ak08 = book<TH1F>("p_{T} AK8 GenJets", "p_{T} [GeV/c]", 20, 0, 1000);
  GenNumber_ak08 = book<TH1F>("Number AK8 GenJets", "number", 15, 0, 15);
  GenEta_ak08 = book<TH1F>("#eta AK8 GenJets", "#eta", 20, -3, 3);
  GenJet1Mass_ak08 = book<TH1F>("GenJet Mass AK8", "M_{jet}", 15, 0, 300);

  GenPT_ak10 = book<TH1F>("p_{T} AK10 GenJets", "p_{T} [GeV/c]", 20, 0, 1000);
  GenNumber_ak10 = book<TH1F>("Number AK10 GenJets", "number", 15, 0, 15);
  GenEta_ak10 = book<TH1F>("#eta AK10 GenJets", "#eta", 20, -3, 3);
  GenJet1Mass_ak10 = book<TH1F>("GenJet Mass AK10", "M_{jet}", 15, 0, 300);

  GenPT_ak12 = book<TH1F>("p_{T} AK12 GenJets", "p_{T} [GeV/c]", 20, 0, 1000);
  GenNumber_ak12 = book<TH1F>("Number AK12 GenJets", "number", 15, 0, 15);
  GenEta_ak12 = book<TH1F>("#eta AK12 GenJets", "#eta", 20, -3, 3);
  GenJet1Mass_ak12 = book<TH1F>("GenJet Mass AK12", "M_{jet}", 15, 0, 300);

  GenPT_ak14 = book<TH1F>("p_{T} AK14 GenJets", "p_{T} [GeV/c]", 20, 0, 1000);
  GenNumber_ak14 = book<TH1F>("Number AK14 GenJets", "number", 15, 0, 15);
  GenEta_ak14 = book<TH1F>("#eta AK14 GenJets", "#eta", 20, -3, 3);
  GenJet1Mass_ak14 = book<TH1F>("GenJet Mass AK14", "M_{jet}", 15, 0, 300);

   h_pfpart=ctx.get_handle<vector<PFParticle>>("PFParticles");
}



void MTopJetGenHists::fill(const Event & event){

  // Cluster Gen Jets
  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc=new JetCluster();
  std::vector<fastjet::PseudoJet> gen_ak08, gen_ak10, gen_ak12, gen_ak14;
  gen_ak08 = jetc->get_genjets(genparts, JetCluster::e_akt, 0.8, 0);
  gen_ak10 = jetc->get_genjets(genparts, JetCluster::e_akt, 1.0, 0);
  gen_ak12 = jetc->get_genjets(genparts, JetCluster::e_akt, 1.2, 0);
  gen_ak14 = jetc->get_genjets(genparts, JetCluster::e_akt, 1.4, 0);
  ////

  // Cluster Reco Jets
  std::vector<PFParticle> pfparts = event.get(h_pfpart);
  JetCluster* jetc_reco=new JetCluster();
  std::vector<fastjet::PseudoJet> reco_ak08 = jetc_reco->get_recojets(&pfparts, JetCluster::e_akt, 0.8, 200);
  ////

  // get weight
  double weight = event.weight;
  ////

  //Reco parts Number per Event
  float reconumber = pfparts.size();
  RecoEvNumber->Fill(reconumber, weight);
  ////

  //// Number of Reco Jets
  RecoNumber_old->Fill(event.topjets->size(), weight);
  RecoNumber_ak08->Fill(reco_ak08.size(), weight);
  ////

  // Reco PT & Eta
  if((event.topjets->size())>0){
    for(unsigned int i=0; i<event.topjets->size();i++){
    const Particle* recjet = &event.topjets->at(i);
    RecoPT_old->Fill(recjet->v4().pt(),weight);
    RecoEta_old->Fill(recjet->v4().eta(),weight);
    }
  } 
  if((reco_ak08.size())>0){
    for(unsigned int i=0; i<reco_ak08.size();i++){
    fastjet::PseudoJet recjet = reco_ak08[i];
    RecoPT_ak08->Fill(recjet.pt(),weight);
    RecoEta_ak08->Fill(recjet.eta(),weight);
   }
  } 
 
  //Mass Reco Jet
  if((event.topjets->size())>0){
    const Particle* recjet = &event.topjets->at(0);
    RecoJet1Mass_old->Fill(recjet->v4().M(),weight);
  } 
  if((reco_ak08.size())>0){
    fastjet::PseudoJet recjet = reco_ak08[0];
    TLorentzVector jet;
    jet.SetPtEtaPhiE(recjet.pt(),recjet.eta(),recjet.phi(),recjet.E());
    RecoJet1Mass_ak08->Fill(jet.M(),weight);
  } 
  ////

  //Gen parts Number per Event
  float gennumber = pfparts.size();
  GenEvNumber->Fill(gennumber, weight);
  ////

  //// Number of Gen Jets
  GenNumber_ak08->Fill(gen_ak08.size(), weight);
  GenNumber_ak10->Fill(gen_ak10.size(), weight);
  GenNumber_ak12->Fill(gen_ak12.size(), weight);
  GenNumber_ak14->Fill(gen_ak14.size(), weight);
  ////

  // Gen PT & Eta
  if((gen_ak08.size())>0){
    for(unsigned int i=0; i<gen_ak08.size();i++){
    fastjet::PseudoJet genjet = gen_ak08[i];
    GenPT_ak08->Fill(genjet.pt(),weight);
    GenEta_ak08->Fill(genjet.eta(),weight);
   }
  } 
  if((gen_ak10.size())>0){
    for(unsigned int i=0; i<gen_ak10.size();i++){
    fastjet::PseudoJet genjet = gen_ak10[i];
    GenPT_ak10->Fill(genjet.pt(),weight);
    GenEta_ak10->Fill(genjet.eta(),weight);
    }
  } 
  if((gen_ak12.size())>0){
    for(unsigned int i=0; i<gen_ak12.size();i++){
    fastjet::PseudoJet genjet = gen_ak12[i];
    GenPT_ak12->Fill(genjet.pt(),weight);
    GenEta_ak12->Fill(genjet.eta(),weight);
   }
  } 
  if((gen_ak14.size())>0){
    for(unsigned int i=0; i<gen_ak14.size();i++){
    fastjet::PseudoJet genjet = gen_ak14[i];
    GenPT_ak14->Fill(genjet.pt(),weight);
    GenEta_ak14->Fill(genjet.eta(),weight);
    }
  } 


  //Mass Gen Jet
  if((gen_ak08.size())>0){
    fastjet::PseudoJet genjet = gen_ak08[0];
    TLorentzVector jet;
    jet.SetPtEtaPhiE(genjet.pt(),genjet.eta(),genjet.phi(),genjet.E());
    GenJet1Mass_ak08->Fill(jet.M(),weight);
  } 
  if((gen_ak10.size())>0){
    fastjet::PseudoJet genjet = gen_ak10[0];
    TLorentzVector jet;
    jet.SetPtEtaPhiE(genjet.pt(),genjet.eta(),genjet.phi(),genjet.E());
    GenJet1Mass_ak10->Fill(jet.M(),weight);
  } 
  if((gen_ak12.size())>0){
    fastjet::PseudoJet genjet = gen_ak12[0];
    TLorentzVector jet;
    jet.SetPtEtaPhiE(genjet.pt(),genjet.eta(),genjet.phi(),genjet.E());
    GenJet1Mass_ak12->Fill(jet.M(),weight);
  } 
  if((gen_ak14.size())>0){
    fastjet::PseudoJet genjet = gen_ak14[0];
    TLorentzVector jet;
    jet.SetPtEtaPhiE(genjet.pt(),genjet.eta(),genjet.phi(),genjet.E());
    GenJet1Mass_ak14->Fill(jet.M(),weight);
  } 
  ////


}


