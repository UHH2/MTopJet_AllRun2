#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/common/include/Utils.h>
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
  MET   = book<TH1F>("MET", "missing p_{T} [GeV]", 200,0,1000);
  HT    = book<TH1F>("HT", "H_{T} Jets [GeV]", 100, 0, 3500);
  HTLep = book<TH1F>("HTLep", "H_{T} Lep [GeV]", 100, 0, 1000);
  ST    = book<TH1F>("ST", "S_{T} [GeV]", 100, 0, 5000);

  // b-tag for CSVv2
  BTAG_L_CSV = book<TH1F>("BTAG_L_CSV", "N b-tags loose", 10, 0, 10);
  BTAG_M_CSV = book<TH1F>("BTAG_M_CSV", "N b-tags medium", 10, 0, 10);
  BTAG_T_CSV = book<TH1F>("BTAG_T_CSV", "N b-tags tight", 10, 0, 10);

  // b-tag for Deepjet
  BTAG_L_DJ = book<TH1F>("BTAG_L_DJ", "N b-tags loose", 10, 0, 10);
  BTAG_M_DJ = book<TH1F>("BTAG_M_DJ", "N b-tags medium", 10, 0, 10);
  BTAG_T_DJ = book<TH1F>("BTAG_T_DJ", "N b-tags tight", 10, 0, 10);

  BTAG_VALUE_DJ = book<TH1F>("BTAG_VALUE_DJ", "Probability", 20, 0, 1);

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
  TopNumber  = book<TH1F>("Number Top Jets", "number", 10, 0, 10);
  TopPT1     = book<TH1F>("1st Top Jet p_{T}", "p_{T}^{topjet1} [GeV/c]", 20, 0, 1000);
  TopPT2     = book<TH1F>("2nd Top Jet p_{T}", "p_{T}^{topjet2} [GeV/c]", 20, 0, 1000);
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

  hist("BTAG_L_CSV")->Fill(jetN__CSVL, weight);
  hist("BTAG_M_CSV")->Fill(jetN__CSVM, weight);
  hist("BTAG_T_CSV")->Fill(jetN__CSVT, weight);



  int jetN__DJL(0), jetN__DJM(0), jetN__DJT(0);
  for(const auto& j : *event.jets){
    if(DeepJetBTag(DeepJetBTag::WP_LOOSE)(j, event)) ++jetN__DJL;
    if(DeepJetBTag(DeepJetBTag::WP_MEDIUM)(j, event)) ++jetN__DJM;
    if(DeepJetBTag(DeepJetBTag::WP_TIGHT)(j, event)) ++jetN__DJT;
  }

  hist("BTAG_L_DJ")->Fill(jetN__DJL, weight);
  hist("BTAG_M_DJ")->Fill(jetN__DJM, weight);
  hist("BTAG_T_DJ")->Fill(jetN__DJT, weight);

  for(const auto& j : *event.jets) hist("BTAG_VALUE_DJ")->Fill(j.btag_DeepJet(), weight);

  for(const auto & jet : *event.jets){
    if(jet.btag_combinedSecondaryVertex() > 0.935){
      if((event.muons->size())!= 0){
        deltaR_lep_Bjet->Fill(deltaR(event.muons->at(0), jet),weight);
      }
      else if((event.electrons->size())!= 0) deltaR_lep_Bjet->Fill(deltaR(event.electrons->at(0), jet),weight);
    }
  }
  //

  // delta R1 Jet
  if((event.jets->size())>0){
    double dR1 = -1;
    if((event.muons->size())!= 0){
      dR1 = deltaR(event.muons->at(0), event.jets->at(0));
    }
    else if((event.electrons->size())!= 0) dR1 = deltaR(event.electrons->at(0), event.jets->at(0));

    deltaR_lep_jet1->Fill(dR1, weight);
  }
  //

  // delta R2  Jet
  if((event.jets->size())>1){
    double dR2 = -1;
    if((event.muons->size())!= 0){
      dR2 = deltaR(event.muons->at(0), event.jets->at(1));
    }
    else if((event.electrons->size())!= 0) dR2 = deltaR(event.electrons->at(0), event.jets->at(1));

    deltaR_lep_jet2->Fill(dR2, weight);
  }
  //

  // delta R2 Top Jet
  if((event.topjets->size())>1){
    double dR2 = -1;
    if((event.muons->size())!= 0){
      dR2 = deltaR(event.muons->at(0), event.topjets->at(1));
    }
    else if((event.electrons->size())!= 0) dR2 = deltaR(event.electrons->at(0), event.topjets->at(1));

    deltaR_lep_topjet2->Fill(dR2, weight);
  }
  //

  // delta R1 Top Jet
  if((event.topjets->size())>0){
    double dR1 = -1;
    if((event.muons->size())!= 0){
      dR1 = deltaR(event.muons->at(0), event.topjets->at(0));
    }
    else if((event.electrons->size())!= 0) dR1 = deltaR(event.electrons->at(0), event.topjets->at(0));

    deltaR_lep_topjet1->Fill(dR1, weight);
  }
  //

  // delta Phi Top Jet
  if((event.topjets->size())>0){
    double dphi = -1;
    const Particle* TopJet = &event.topjets->at(0);
    if((event.muons->size())!= 0){
      const Particle* Lep = &event.muons->at(0);
      dphi = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    else if((event.electrons->size())!= 0) {
      const Particle* Lep = &event.electrons->at(0);
      dphi = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    deltaPhi_lep_topjet1->Fill(dphi, weight);
  }
  //

  // delta Phi Top Jet2
  if((event.topjets->size())>1){
    double dphi2 = -1;
    const Particle* TopJet = &event.topjets->at(1);
    if((event.muons->size())!= 0){
      const Particle* Lep = &event.muons->at(0);
      dphi2 = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    else if((event.electrons->size())!= 0) {
      const Particle* Lep = &event.electrons->at(0);
      dphi2 = (TopJet->v4().Phi()) - (Lep->v4().Phi());
    }
    deltaPhi_lep_topjet2->Fill(dphi2, weight);
  }
}
