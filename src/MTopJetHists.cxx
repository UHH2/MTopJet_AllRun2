#include "UHH2/MTopJet/include/MTopJetHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

MTopJetHists::MTopJetHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  // Event Hists
  N_PrimVertices = book<TH1F>("N_PrimVertices", "number of primary vertices", 56, -0.5, 55.5);
  //N_TrueInteractions = book<TH1F>("N_TrueInteractions", "number of true interactions", 50, 0, 50);
  //Weights = book<TH1F>("Weights", "weights", 100,0,2);
  MET = book<TH1F>("MET", "missing E_{T}", 200,0,1000);
  HT = book<TH1F>("HT", "H_{T} Jets", 100, 0, 3500);
  HTLep = book<TH1F>("HTLep", "H_{T} Lep", 100, 0, 1000);
  ST = book<TH1F>("ST", "S_{T}", 100, 0, 5000);

  // b-tag
  BTAG_L = book<TH1F>("BTAG_L", "N b-tags loose", 10, 0, 10);
  BTAG_M = book<TH1F>("BTAG_M", "N b-tags medium", 10, 0, 10);
  BTAG_T = book<TH1F>("BTAG_T", "N b-tags tight", 10, 0, 10);

  deltaR_lep_Bjet = book<TH1F>("deltaR_lep_Bjet", "#Delta R(lep, b-tagged Jet)", 80, 0, 4.0);


  // delta R (lep, jet)
  deltaR_lep_jet1 = book<TH1F>("deltaR_lep_jet1", "#Delta R(lep,1st Jet)", 80, 0, 4.0);
  deltaR_lep_jet2 = book<TH1F>("deltaR_lep_jet2", "#Delta R(lep,2nd Jet)", 80, 0, 4.0);
 
  // delta R (lep, topjet2)
  deltaR_lep_topjet1 = book<TH1F>("deltaR_lep_topjet1", "#Delta R(lep,1st Top Jet)", 80, 0, 4.0);
  deltaR_lep_topjet2 = book<TH1F>("deltaR_lep_topjet2", "#Delta R(lep,2nd Top Jet)", 80, 0, 4.0);

  // delta Phi (lep, TopJet)
  deltaPhi_lep_topjet1 = book<TH1F>("deltaPhi_lep_topjet1", "#Delta Phi(lep,1st Top Jet)", 80, -7.0, 7.0);
  deltaPhi_lep_topjet2 = book<TH1F>("deltaPhi_lep_topjet2", "#Delta Phi(lep,2nd Top Jet)", 80, -7.0, 7.0);

  // TopJet number, pt, Mass
  TopNumber = book<TH1F>("Number Top Jets", "number", 10, 0, 10);
  TopPT1 = book<TH1F>("1st Top Jet p_{T}", "p_{T}^{topjet1} [GeV/c]", 20, 0, 1000);
  TopPT2 = book<TH1F>("2nd Top Jet p_{T}", "p_{T}^{topjet2} [GeV/c]", 20, 0, 1000);
  TopJetMass = book<TH1F>("1st Top Jet Mass","M^{topjet1} [GeV/c^{2}]", 15, 0, 300);

  // TopJetMass1 vs TopJetMass2
  TopJetMass1_TopJetMass2 = book<TH2F>("TopJetMass1_TopJetMass2", "x=M_Top1 y=M_Top2", 30, 0, 300., 30, 0, 300.);
}

void MTopJetHists::fill(const Event & event){

  // get weight
  double weight = event.weight;
  N_PrimVertices->Fill(event.pvs->size(), event.weight);

  //HT
  auto met = event.met->pt();
  hist("MET")->Fill(met, weight);
  double st = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  double ht_elec = 0.0;
  double ht_muon = 0.0;

  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
  }

  for(const auto & electron : *event.electrons){
    ht_elec += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_muon += muon.pt();
  }
  
  ht_lep = ht_elec + ht_muon + met;
  st = ht_elec + ht_muon + ht_jets + met;
  hist("HT")->Fill(ht_jets,weight);
  hist("HTLep")->Fill(ht_lep,weight);
  hist("ST")->Fill(st, weight);

  // b-tag

  std::vector<float> jets_CSV;
  jets_CSV.reserve(event.jets->size());
  for(const auto& j : *event.jets) jets_CSV.push_back(j.btag_combinedSecondaryVertex());
  std::sort(jets_CSV.begin(), jets_CSV.end(), [](const float s1, const float s2){return s1 > s2;});

  int jetN__CSVL(0), jetN__CSVM(0), jetN__CSVT(0);
  for(unsigned int i=0; i<jets_CSV.size(); ++i){

    const float& csv = jets_CSV.at(i);

    if(csv > 0.460) ++jetN__CSVL;
    if(csv > 0.800) ++jetN__CSVM;
    if(csv > 0.935) ++jetN__CSVT;

  }

  hist("BTAG_L")->Fill(jetN__CSVL, weight);
  hist("BTAG_M")->Fill(jetN__CSVM, weight);
  hist("BTAG_T")->Fill(jetN__CSVT, weight);

  for(const auto & jet : *event.jets){
    if(jet.btag_combinedSecondaryVertex() > 0.935){
      if((event.muons->size())!= 0){
	deltaR_lep_Bjet->Fill(deltaR(event.muons->at(0), jet),weight);
      }
      else deltaR_lep_Bjet->Fill(deltaR(event.electrons->at(0), jet),weight);
    }
  }
  //

  // delta R1 Jet
  if((event.jets->size())>0){
  double dR1;
  if((event.muons->size())!= 0){
    dR1 = deltaR(event.muons->at(0), event.jets->at(0));
  }
  else dR1 = deltaR(event.electrons->at(0), event.jets->at(0));

  deltaR_lep_jet1->Fill(dR1, weight);
  }
  //

  // delta R2  Jet
  if((event.jets->size())>1){
  double dR2;
  if((event.muons->size())!= 0){
    dR2 = deltaR(event.muons->at(0), event.jets->at(1));
  }
  else dR2 = deltaR(event.electrons->at(0), event.jets->at(1));

  deltaR_lep_jet2->Fill(dR2, weight);
  }
  //

  // delta R2 Top Jet
  if((event.topjets->size())>1){
  double dR2;
  if((event.muons->size())!= 0){
    dR2 = deltaR(event.muons->at(0), event.topjets->at(1));
  }
  else dR2 = deltaR(event.electrons->at(0), event.topjets->at(1));

  deltaR_lep_topjet2->Fill(dR2, weight);
  }
  //

  // delta R1 Top Jet
  if((event.topjets->size())>0){
  double dR1;
  if((event.muons->size())!= 0){
    dR1 = deltaR(event.muons->at(0), event.topjets->at(0));
  }
  else dR1 = deltaR(event.electrons->at(0), event.topjets->at(0));

  deltaR_lep_topjet1->Fill(dR1, weight);
  }
  //

  // delta Phi Top Jet
  if((event.topjets->size())>0){
    double dphi;
    const Particle* TopJet = &event.topjets->at(0);
    if((event.muons->size())!= 0){
      const Particle* Lep = &event.muons->at(0);
      dphi = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    else {
      const Particle* Lep = &event.electrons->at(0);
      dphi = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
      deltaPhi_lep_topjet1->Fill(dphi, weight);
  }
  //

  // delta Phi Top Jet2
  if((event.topjets->size())>1){
  double dphi2;
 const Particle* TopJet = &event.topjets->at(1);
    if((event.muons->size())!= 0){
      const Particle* Lep = &event.muons->at(0);
      dphi2 = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    else {
      const Particle* Lep = &event.electrons->at(0);
      dphi2 = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
      deltaPhi_lep_topjet2->Fill(dphi2, weight);
  }
  //

  //Top Number
  float number;
  number = event.topjets->size();
  TopNumber->Fill(number, weight);
  //

  //Top pT
  if((event.topjets->size())>0){
  assert(event.jets);
  float toppt1;
  const Particle* Top1 = &event.topjets->at(0);
  toppt1 = Top1->v4().pt();
  TopPT1->Fill(toppt1, weight);
  }
  if((event.topjets->size())>1){
  assert(event.jets);
  float toppt2;
  const Particle* Top2 = &event.topjets->at(1);
  toppt2 = Top2->v4().pt();
  TopPT2->Fill(toppt2, weight);
  }
  //

  //TopJetMass
  if((event.topjets->size())>0){
  assert(event.jets);
  float mass1;
  const Particle* Top1 = &event.topjets->at(0);
  mass1 = Top1->v4().M();
  TopJetMass->Fill(mass1, weight);
  }
  //

  //TopJetMass M1 vs M2
  if((event.topjets->size())>1){
  assert(event.jets);
  float mass1, mass2;
  const Particle* Top1 = &event.topjets->at(0);
  const Particle* Top2 = &event.topjets->at(1);
  mass1 = Top1->v4().M();
  mass2 = Top2->v4().M();
  TopJetMass1_TopJetMass2->Fill(mass1, mass2, weight);
  }
}
