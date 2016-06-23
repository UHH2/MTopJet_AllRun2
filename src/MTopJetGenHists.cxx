#include "UHH2/MTopJet/include/MTopJetGenHists.h"
#include "UHH2/MTopJet/include/JetCluster.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/MTopJet/include/GenJetProps.h"
#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

MTopJetGenHists::MTopJetGenHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  // GenJet number, pt, Mass
  GenPartNumber = book<TH1F>("Number Gen Parts", "number", 3, 0, 3);
  GenNumber1_new = book<TH1F>("new Number Gen Jets 1", "number", 10, 0, 10);
  GenNumber2_new = book<TH1F>("new Number Gen Jets 2", "number", 10, 0, 10);
  GenNumber1_old = book<TH1F>("old Number Gen Jets 1", "number", 10, 0, 10);
  GenNumber2_old = book<TH1F>("old Number Gen Jets 2", "number", 10, 0, 10);  
  GenPT1_new = book<TH1F>("new ak4 Gen Jet p_{T}", "p_{T} [GeV/c]", 20, 0, 1000);
  GenPT2_new = book<TH1F>("new ak8 Gen Jet p_{T}", "p_{T} [GeV/c]", 20, 0, 1000);
  GenJetMass_new = book<TH1F>("new 1st Gen Jet Mass","M^{topjet1} [GeV/c^{2}]", 15, 0, 300);
  GenPT1_old = book<TH1F>("old ak4 Gen Jet p_{T}", "p_{T} [GeV/c]", 20, 0, 1000);
  GenPT2_old = book<TH1F>("old ak8 Gen Jet p_{T}", "p_{T} [GeV/c]", 20, 0, 1000);
  GenJetMass_old = book<TH1F>("old 1st Gen Jet Mass","M^{topjet1} [GeV/c^{2}]", 15, 0, 300);   


}




void MTopJetGenHists::fill(const Event & event){
  // Cluster Jets
  std::vector<GenParticle>* genparts = event.genparticles;
  JetCluster* jetc=new JetCluster();
  std::vector<fastjet::PseudoJet> fatjet1 = jetc->get_genjets(genparts, JetCluster::e_akt, 0.4, 10);
  std::vector<fastjet::PseudoJet> fatjet2 = jetc->get_genjets(genparts, JetCluster::e_akt, 0.8, 10);

  // get weight
  double weight = event.weight;

  //Gen parts Number
  float gennumber=0;
  gennumber = genparts->size();
  GenPartNumber->Fill(1, gennumber);

  //Gen Number
  float number1;
  number1 = fatjet1.size();
  GenNumber1_new->Fill(number1,weight);

  float number2;
  number2 = fatjet2.size();
  GenNumber2_new->Fill(number2,weight);

  float number1_;
  number1_ = event.genjets->size();
  GenNumber1_old->Fill(number1_,weight);

  // float number2_;
  // number2_ = event.gentopjets->size();
  // GenNumber2_old->Fill(number2_,weight);

  //

  //Gen pT
  if((fatjet1.size())>0){
    // assert(event.jets);
    for(unsigned int i=0; i<fatjet1.size();i++){
    fastjet::PseudoJet Gen = fatjet1[i];
    float toppt1 = Gen.pt();
    GenPT1_new->Fill(toppt1,weight);
    }
  }
  if((fatjet2.size())>0){
    // assert(event.jets);
    for(unsigned int i=0; i<fatjet2.size();i++){
    fastjet::PseudoJet Gen = fatjet2[i];
    float toppt2 = Gen.pt();
    GenPT2_new->Fill(toppt2,weight);
    }
  }

  if((event.genjets->size())>0){
    for(unsigned int i=0; i<event.genjets->size();i++){
    const Particle* gen = &event.genjets->at(i);
    float pt = gen->v4().pt();
    GenPT1_old->Fill(pt,weight);
    }
  }

  // if((event.gentopjets->size())>0){
  //   for(unsigned int i=0; i<event.gentopjets->size();i++){
  //   const Particle* gen = &event.gentopjets->at(i);
  //   float pt2 = gen->v4().pt();
  //   GenPT2_old->Fill(pt2,weight);
  //   }
  // }

  

  // GenJetMass
  // if((fatjet.size())>0){
  //   // assert(event.jets);
  //   float mass1 = 0;
  //   fastjet::PseudoJet Gen1 = fatjet[0];
  //   mass1 = Gen1.v4();
  //   GenJetMass->Fill(mass1, weight);
  // }
  
}


