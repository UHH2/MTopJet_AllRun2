#include "UHH2/MTopJet/include/ClusteringHists.h"
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

ClusteringHists::ClusteringHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  FatJet1 = book<TH2F>("JetDisplay_pf0", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  FatJet2 = book<TH2F>("JetDisplay_pf1", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  SubJet1_1 = book<TH2F>("JetDisplay_sub1_1", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  SubJet1_2 = book<TH2F>("JetDisplay_sub1_2", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  SubJet1_3 = book<TH2F>("JetDisplay_sub1_3", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  SubJet2_1 = book<TH2F>("JetDisplay_sub2_1", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  SubJet2_2 = book<TH2F>("JetDisplay_sub2_2", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  NotClustered = book<TH2F>("JetDisplay_radiation", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  NotClustered_sub1 = book<TH2F>("JetDisplay_sub1_radiation", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  NotClustered_sub2 = book<TH2F>("JetDisplay_sub2_radiation", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  AllPF = book<TH2F>("JetDisplay_pf_all", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  DecayProducts = book<TH2F>("JetDisplay_decay", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);
  Lepton = book<TH2F>("JetDisplay_lepton", "x=#eta y=#phi", 50, -M_PI , M_PI, 50, -M_PI, M_PI);


  // handle for particle-jet-lists
  h_list_fat=ctx.get_handle<std::vector<int>>("particle_fatjet_list");
  h_list_sub1=ctx.get_handle<std::vector<int>>("particle_subjets1_list");
  h_list_sub2=ctx.get_handle<std::vector<int>>("particle_subjets2_list");
  // handle for TTbarGen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
}



void ClusteringHists::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get needed objects  ---------------------------------
  //---------------------------------------------------------------------------------------
  const auto & ttbargen = event.get(h_ttbargen);
  vector<int> list_fat = event.get(h_list_fat);
  vector<int> list_sub1 = event.get(h_list_sub1);
  vector<int> list_sub2 = event.get(h_list_sub2);

  // get List with Ghost Particles
  JetCluster* jetc = new JetCluster();
  std::vector<GenParticle>* genparts = event.genparticles;
  std::vector<fastjet::PseudoJet> particle_in, particle_in_noGhost, fat1, fat2, notclustered, sub1_1, sub1_2, sub1_3, notclustered_sub1, sub2_1, sub2_2, notclustered_sub2; 
  particle_in_noGhost.clear();
  particle_in.clear();
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle* part = &(genparts->at(i));
    if(jetc->IsStable(part)){
      particle_in_noGhost.push_back(jetc->convert_particle(part));
    }
  }
  particle_in = jetc->add_ghosts(particle_in_noGhost);

  // make particle list for each jet:
  // fat jets
  for(unsigned i=0; i<particle_in.size(); i++){
    if(list_fat[i] == 0){
      fat1.push_back(particle_in[i]);
    }
    else if(list_fat[i] == 1){
      fat2.push_back(particle_in[i]);
    }
    else if(list_fat[i] == -1){
      notclustered.push_back(particle_in[i]);
    }
  }
  //sub jets
  for(unsigned i=0; i<fat1.size(); i++){
    if(list_sub1[i] == 0){
      sub1_1.push_back(fat1[i]);
    }
    else if(list_sub1[i] == 1){
      sub1_2.push_back(fat1[i]);
    }
    else if(list_sub1[i] == 2){
      sub1_3.push_back(fat1[i]);
    }
    else if(list_sub1[i] == -1){
      notclustered_sub1.push_back(fat1[i]);
    }
  }
  for(unsigned i=0; i<fat2.size(); i++){
    if(list_sub2[i] == 0){
      sub2_1.push_back(fat2[i]);
    }
    else if(list_sub2[i] == 1){
      sub2_2.push_back(fat2[i]);
    }
    else if(list_sub2[i] == -1){
      notclustered_sub2.push_back(fat2[i]);
    }
  }

  //---------------------------------------------------------------------------------------
  //--------------------------------- get decay products  ---------------------------------
  //---------------------------------------------------------------------------------------
  // bool matched = false;
  // get stable particles from ttbar decay and sort them into leptonic and hadronic
  GenParticle bot, q1, q2, bot_lep, lep1, lep2, lepton; //leptons already defined above
  if(ttbargen.IsTopHadronicDecay()){
    bot = ttbargen.bTop();
    q1 = ttbargen.Wdecay1();
    q2 = ttbargen.Wdecay2();
    bot_lep = ttbargen.bAntitop();
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    bot = ttbargen.bAntitop();
    q1 = ttbargen.WMinusdecay1();
    q2 = ttbargen.WMinusdecay2();
    bot_lep = ttbargen.bTop();
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  else if(!(ttbargen.IsTopHadronicDecay() || ttbargen.IsAntiTopHadronicDecay())){
    throw runtime_error("no hadronic Top found");
  }
  //check which lep is neutrino and which is elec/muon
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
  }

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
  

  // get weight
  // double weight = event.weight;
  ////
  for(unsigned int j=0; j<fat1.size(); j++){
    FatJet1->Fill(fat1.at(j).eta(), fat1.at(j).phi_std(), fat1.at(j).E());
  }
  for(unsigned int j=0; j<fat2.size(); j++){
    FatJet2->Fill(fat2.at(j).eta(), fat2.at(j).phi_std(), fat2.at(j).E());
  }  
  for(unsigned int j=0; j<notclustered.size(); j++){
    NotClustered->Fill(notclustered.at(j).eta(), notclustered.at(j).phi_std(), notclustered.at(j).E());
  } 
  for(unsigned int j=0; j<particle_in.size(); j++){
    AllPF->Fill(particle_in.at(j).eta(), particle_in.at(j).phi_std(), particle_in.at(j).E());
  }  
  for(unsigned int j=0; j<sub1_1.size(); j++){
    SubJet1_1->Fill(sub1_1.at(j).eta(), sub1_1.at(j).phi_std(), sub1_1.at(j).E());
  }  
  for(unsigned int j=0; j<sub1_2.size(); j++){
    SubJet1_2->Fill(sub1_2.at(j).eta(), sub1_2.at(j).phi_std(), sub1_2.at(j).E());
  } 
  for(unsigned int j=0; j<sub1_3.size(); j++){
    SubJet1_3->Fill(sub1_3.at(j).eta(), sub1_3.at(j).phi_std(), sub1_3.at(j).E());
  } 
  for(unsigned int j=0; j<sub2_1.size(); j++){
    SubJet2_1->Fill(sub2_1.at(j).eta(), sub2_1.at(j).phi_std(), sub2_1.at(j).E());
  }  
  for(unsigned int j=0; j<sub2_2.size(); j++){
    SubJet2_2->Fill(sub2_2.at(j).eta(), sub2_2.at(j).phi_std(), sub2_2.at(j).E());
  } 
  for(unsigned int j=0; j<notclustered_sub1.size(); j++){
    NotClustered_sub1->Fill(notclustered_sub1.at(j).eta(), notclustered_sub1.at(j).phi_std(), notclustered_sub1.at(j).E());
  } 
  for(unsigned int j=0; j<notclustered_sub2.size(); j++){
    NotClustered_sub2->Fill(notclustered_sub2.at(j).eta(), notclustered_sub2.at(j).phi_std(), notclustered_sub2.at(j).E());
  } 
  DecayProducts->Fill(bot.eta(), bot.phi(), bot.energy());
  DecayProducts->Fill(bot_lep.eta(), bot_lep.phi(), bot_lep.energy());
  DecayProducts->Fill(q1.eta(), q1.phi(), q1.energy());
  DecayProducts->Fill(q2.eta(), q2.phi(), q2.energy());
  DecayProducts->Fill(bot.eta(), bot.phi(), bot.energy());
  DecayProducts->Fill(lep1.eta(), lep1.phi(), lep1.energy());
  DecayProducts->Fill(lep2.eta(), lep2.phi(), lep2.energy());

  Lepton->Fill(lepton.eta(), lepton.phi(), lepton.energy());
}
