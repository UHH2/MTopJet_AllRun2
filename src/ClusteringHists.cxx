#include "UHH2/MTopJet/include/ClusteringHists.h"


ClusteringHists::ClusteringHists(uhh2::Context & ctx, const std::string & dirname, double phishift_): Hists(ctx, dirname){
  // book all histograms here

  Top = book<TH2F>("JetDisplay_top", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Top_lep = book<TH2F>("JetDisplay_toplep", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Top_had = book<TH2F>("JetDisplay_tophad", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  FatJet1 = book<TH2F>("JetDisplay_pf0", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  FatJet2 = book<TH2F>("JetDisplay_pf1", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  SubJet1_1 = book<TH2F>("JetDisplay_sub1_1", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  SubJet1_2 = book<TH2F>("JetDisplay_sub1_2", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  SubJet1_3 = book<TH2F>("JetDisplay_sub1_3", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  SubJet2_1 = book<TH2F>("JetDisplay_sub2_1", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  SubJet2_2 = book<TH2F>("JetDisplay_sub2_2", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  NotClustered = book<TH2F>("JetDisplay_radiation", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  NotClustered_sub1 = book<TH2F>("JetDisplay_sub1_radiation", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  NotClustered_sub2 = book<TH2F>("JetDisplay_sub2_radiation", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  AllPF = book<TH2F>("JetDisplay_pf_all", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  DecayProducts = book<TH2F>("JetDisplay_decay", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  DecayProducts_lep = book<TH2F>("JetDisplay_decaylep", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  DecayProducts_had = book<TH2F>("JetDisplay_decayhad", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Bottom_lep = book<TH2F>("JetDisplay_bottomlep", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Bottom_had = book<TH2F>("JetDisplay_bottomhad", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Lepton = book<TH2F>("JetDisplay_lepton", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  Neutrino = book<TH2F>("JetDisplay_neutrino", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  LepJet = book<TH2F>("LepJet", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  HadJet = book<TH2F>("HadJet", "x=#eta y=#phi", 200, -M_PI, M_PI, 200, -M_PI, M_PI);
  PT = book<TH1F>("JetpT", "p_T(jet)", 2, 0.5, 2.5);
  Mass = book<TH1F>("Jetmass", "M(jet)", 2, 0.5, 2.5);
  Number_incjets = book<TH1F>("Number_incjets", "number of inc jets", 4, -0.5, 3.5);
  List_incjet = book<TH1F>("List_incjet1", "entry in particle list", 5, -1.5, 3.5);
  List_subjet1 = book<TH1F>("List_subjet1", "entry in particle list", 5, -1.5, 3.5);
  List_subjet2 = book<TH1F>("List_subjet2", "entry in particle list", 5, -1.5, 3.5);
  LepJetBool = book<TH1F>("LepJetBool", "lepton is in jet x", 2, 0.5, 2.5);
  Eta = book<TH1F>("Eta", "genpartidcle eta", 160, -4, 4);

  // handle for particle-jet-lists
  h_list_fat=ctx.get_handle<std::vector<int>>("particle_fatjet_list");
  h_list_sub1=ctx.get_handle<std::vector<int>>("particle_subjets1_list");
  h_list_sub2=ctx.get_handle<std::vector<int>>("particle_subjets2_list");
  // handle for TTbarGen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
  // handle for clustered jets
  h_jets=ctx.get_handle<std::vector<Jet>>("xcone23_gen_fatjets");
  h_subjets=ctx.get_handle<std::vector<Jet>>("xcone23_gen_subjets");
  h_incjets=ctx.get_handle<std::vector<Jet>>("xcone23_gen_seedjet");
  h_pf_fat1=ctx.get_handle<std::vector<Jet>>("particle_fatjet1");
  h_pf_fat2=ctx.get_handle<std::vector<Jet>>("particle_fatjet2");
  h_pf_sub1_1=ctx.get_handle<std::vector<Jet>>("particle_subjet1_1");
  h_pf_sub1_2=ctx.get_handle<std::vector<Jet>>("particle_subjet1_2");
  h_pf_sub1_3=ctx.get_handle<std::vector<Jet>>("particle_subjet1_3");
  h_pf_sub2_1=ctx.get_handle<std::vector<Jet>>("particle_subjet2_1");
  h_pf_sub2_2=ctx.get_handle<std::vector<Jet>>("particle_subjet2_2");
  // h_pf_fat0=ctx.get_handle<std::vector<Jet>>("particle_fatjet0");
  h_pf_sub1_0=ctx.get_handle<std::vector<Jet>>("particle_subjet1_0");
  h_pf_sub2_0=ctx.get_handle<std::vector<Jet>>("particle_subjet2_0");
  h_pf_all=ctx.get_handle<std::vector<Jet>>("particle_all");

  phishift = phishift_;

}



void ClusteringHists::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get needed objects  ---------------------------------
  //---------------------------------------------------------------------------------------
  const auto & ttbargen = event.get(h_ttbargen);
  vector<int> list_fat = event.get(h_list_fat);
  vector<int> list_sub1 = event.get(h_list_sub1);
  vector<int> list_sub2 = event.get(h_list_sub2);
  vector<Jet> jets = event.get(h_jets);
  vector<Jet> incjets = event.get(h_incjets);
  vector<Jet> subjets = event.get(h_subjets);
  vector<Jet> pf_fat1 = event.get(h_pf_fat1);
  vector<Jet> pf_fat2 = event.get(h_pf_fat2);
  vector<Jet> pf_sub1_1 = event.get(h_pf_sub1_1);
  vector<Jet> pf_sub1_2 = event.get(h_pf_sub1_2);
  vector<Jet> pf_sub1_3 = event.get(h_pf_sub1_3);
  vector<Jet> pf_sub2_1 = event.get(h_pf_sub2_1);
  vector<Jet> pf_sub2_2 = event.get(h_pf_sub2_2);
  // vector<Jet> pf_fat0 = event.get(h_pf_fat0);
  vector<Jet> pf_sub1_0 = event.get(h_pf_sub1_0);
  vector<Jet> pf_sub2_0 = event.get(h_pf_sub2_0);
  vector<Jet> pf_all = event.get(h_pf_all);

  //---------------------------------------------------------------------------------------
  //--------------------------------- get decay products  ---------------------------------
  //---------------------------------------------------------------------------------------
  // bool matched = false;
  // get stable particles from ttbar decay and sort them into leptonic and hadronic
  GenParticle tophad, toplep, bot, q1, q2, bot_lep, lep1, lep2, lepton, neutrino; //leptons already defined above
  if(ttbargen.IsTopHadronicDecay()){
    tophad = ttbargen.Top();
    toplep = ttbargen.Antitop();
    bot = ttbargen.bTop();
    q1 = ttbargen.Wdecay1();
    q2 = ttbargen.Wdecay2();
    bot_lep = ttbargen.bAntitop();
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    tophad = ttbargen.Antitop();
    toplep = ttbargen.Top();
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
    neutrino = lep2;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
    neutrino = lep1;
  }

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

  Number_incjets->Fill(incjets.size(), 1);

  // get List with Ghost Particles
  std::vector<GenParticle>* genparts = event.genparticles;
  std::vector<GenParticle> genparts2;
  std::vector<fastjet::PseudoJet> particle_in, particle_in_noGhost, fat1, fat2, notclustered, sub1_1, sub1_2, sub1_3, notclustered_sub1, sub2_1, sub2_2, notclustered_sub2;
  particle_in_noGhost.clear();
  particle_in.clear();
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle* part = &(genparts->at(i));
    if(part->status() == 1){
      genparts2.push_back(genparts->at(i));
      particle_in_noGhost.push_back(convert_particle(part));
    }
  }
  particle_in = add_ghosts(particle_in_noGhost);

  // make particle list for each jet (+ look for lepton):
  // fat jets
  bool lep_in_jet1 = false;
  bool lep_in_jet2 = false;
  for(unsigned i=0; i<particle_in.size(); i++){
    List_incjet->Fill(list_fat[i], 1);
    if(list_fat[i] == 0){
      fat1.push_back(particle_in[i]);
      if(i < genparts2.size() && (abs(genparts2.at(i).pdgId()) == 11 || abs(genparts2.at(i).pdgId()) == 13)){
	lep_in_jet1 = true;
      }
    }
    else if(list_fat[i] == 1){
      fat2.push_back(particle_in[i]);
      if(i < genparts2.size() && (abs(genparts2.at(i).pdgId()) == 11 || abs(genparts2.at(i).pdgId()) == 13)){
	lep_in_jet2 = true;
      }
    }
    else if(list_fat[i] == -1){
      notclustered.push_back(particle_in[i]);
    }
  }
  // if(lep_in_jet1) cout<<"Lepton is in Jet 1"<<endl;
  // if(lep_in_jet2) cout<<"Lepton is in Jet 2"<<endl;

  // if lepton is not clustered into jet, find nearest jet
  // double dR1, dR2;
  // if(!(lep_in_jet1 || lep_in_jet2)){
  //   cout<<"Lepton is not in particle list"<<endl;
  //   dR1 = uhh2::deltaR(lepton, incjets.at(0));
  //   dR2 = uhh2::deltaR(lepton, incjets.at(1));
  //   if(dR1 < dR2) lep_in_jet1 = true;
  //   if(dR1 > dR2) lep_in_jet2 = true;
  // }

  if(lep_in_jet2){
    //sub jets
    for(unsigned i=0; i<fat1.size(); ++i){
      List_subjet1->Fill(list_sub1[i], 1);
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
    for(unsigned i=0; i<fat2.size(); ++i){
      List_subjet2->Fill(list_sub2[i], 1);
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
  }
  if(lep_in_jet1){
    //sub jets
    for(unsigned i=0; i<fat2.size(); ++i){
      List_subjet1->Fill(list_sub1[i], 1);
      if(list_sub1[i] == 0){
	sub1_1.push_back(fat2[i]);
      }
      else if(list_sub1[i] == 1){
	sub1_2.push_back(fat2[i]);
      }
      else if(list_sub1[i] == 2){
	sub1_3.push_back(fat2[i]);
      }
      else if(list_sub1[i] == -1){
	notclustered_sub1.push_back(fat2[i]);
      }
    }
    for(unsigned i=0; i<fat1.size(); ++i){
      List_subjet2->Fill(list_sub2[i], 1);
      if(list_sub2[i] == 0){
	sub2_1.push_back(fat1[i]);
      }
      else if(list_sub2[i] == 1){
	sub2_2.push_back(fat1[i]);
      }
      else if(list_sub2[i] == -1){
	notclustered_sub2.push_back(fat1[i]);
      }
    }
  }

  bool use_int_list = false;
  // get weight
  // double weight = event.weight;
  ////
  if(use_int_list){
    if(lep_in_jet2){
      for(unsigned int j=0; j<fat1.size(); ++j){
	FatJet1->Fill(fat1.at(j).eta(), ShiftPhi(fat1.at(j).phi_std()), fat1.at(j).E());
      }
      for(unsigned int j=0; j<fat2.size(); ++j){
	FatJet2->Fill(fat2.at(j).eta(), ShiftPhi(fat2.at(j).phi_std()), fat2.at(j).E());
      }
    }
    if(lep_in_jet1){
      for(unsigned int j=0; j<fat2.size(); ++j){
	FatJet1->Fill(fat2.at(j).eta(), ShiftPhi(fat2.at(j).phi_std()), fat2.at(j).E());
      }
      for(unsigned int j=0; j<fat1.size(); ++j){
	FatJet2->Fill(fat1.at(j).eta(), ShiftPhi(fat1.at(j).phi_std()), fat1.at(j).E());
      }
    }
    for(unsigned int j=0; j<notclustered.size(); ++j){
      NotClustered->Fill(notclustered.at(j).eta(), ShiftPhi(notclustered.at(j).phi_std()), notclustered.at(j).E());
    }
    for(unsigned int j=0; j<particle_in.size(); ++j){
      AllPF->Fill(particle_in.at(j).eta(), ShiftPhi(particle_in.at(j).phi_std()), particle_in.at(j).E());
    }
    for(unsigned int j=0; j<sub1_1.size(); ++j){
      SubJet1_1->Fill(sub1_1.at(j).eta(), ShiftPhi(sub1_1.at(j).phi_std()), sub1_1.at(j).E());
    }
    for(unsigned int j=0; j<sub1_2.size(); ++j){
      SubJet1_2->Fill(sub1_2.at(j).eta(), ShiftPhi(sub1_2.at(j).phi_std()), sub1_2.at(j).E());
    }
    for(unsigned int j=0; j<sub1_3.size(); ++j){
      SubJet1_3->Fill(sub1_3.at(j).eta(), ShiftPhi(sub1_3.at(j).phi_std()), sub1_3.at(j).E());
    }
    for(unsigned int j=0; j<sub2_1.size(); ++j){
      SubJet2_1->Fill(sub2_1.at(j).eta(), ShiftPhi(sub2_1.at(j).phi_std()), sub2_1.at(j).E());
    }
    for(unsigned int j=0; j<sub2_2.size(); ++j){
      SubJet2_2->Fill(sub2_2.at(j).eta(), ShiftPhi(sub2_2.at(j).phi_std()), sub2_2.at(j).E());
    }
    for(unsigned int j=0; j<notclustered_sub1.size(); ++j){
      NotClustered_sub1->Fill(notclustered_sub1.at(j).eta(), ShiftPhi(notclustered_sub1.at(j).phi_std()), notclustered_sub1.at(j).E());
    }
    for(unsigned int j=0; j<notclustered_sub2.size(); ++j){
      NotClustered_sub2->Fill(notclustered_sub2.at(j).eta(), ShiftPhi(notclustered_sub2.at(j).phi_std()), notclustered_sub2.at(j).E());
    }
  }
  if(!use_int_list){
    double phi;
    for(unsigned int j=0; j<pf_fat1.size(); ++j){
      if( pf_fat1.at(j).phi() > M_PI) phi = pf_fat1.at(j).phi() - 2*M_PI;
      else phi = pf_fat1.at(j).phi();
      FatJet1->Fill(pf_fat1.at(j).eta(), ShiftPhi(phi), pf_fat1.at(j).energy());
    }
    for(unsigned int j=0; j<pf_fat2.size(); ++j){
      if( pf_fat2.at(j).phi() > M_PI) phi = pf_fat2.at(j).phi() - 2*M_PI;
      else phi = pf_fat2.at(j).phi();
      FatJet2->Fill(pf_fat2.at(j).eta(), ShiftPhi(phi), pf_fat2.at(j).energy());
    }
    // for(unsigned int j=0; j<notclustered.size(); ++j){
    //   NotClustered->Fill(notclustered.at(j).eta(), notclustered.at(j).phi(), notclustered.at(j).energy());
    // }
    for(unsigned int j=0; j<pf_all.size(); ++j){
      if( pf_all.at(j).phi() > M_PI) phi = pf_all.at(j).phi() - 2*M_PI;
      else phi = pf_all.at(j).phi();
      AllPF->Fill(pf_all.at(j).eta(), ShiftPhi(phi), pf_all.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub1_1.size(); ++j){
      if( pf_sub1_1.at(j).phi() > M_PI) phi = pf_sub1_1.at(j).phi() - 2*M_PI;
      else phi = pf_sub1_1.at(j).phi();
      SubJet1_1->Fill(pf_sub1_1.at(j).eta(), ShiftPhi(phi), pf_sub1_1.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub1_2.size(); ++j){
      if( pf_sub1_2.at(j).phi() > M_PI) phi = pf_sub1_2.at(j).phi() - 2*M_PI;
      else phi = pf_sub1_2.at(j).phi();
      SubJet1_2->Fill(pf_sub1_2.at(j).eta(), ShiftPhi(phi), pf_sub1_2.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub1_3.size(); ++j){
      if( pf_sub1_3.at(j).phi() > M_PI) phi = pf_sub1_3.at(j).phi() - 2*M_PI;
      else phi = pf_sub1_3.at(j).phi();
      SubJet1_3->Fill(pf_sub1_3.at(j).eta(), ShiftPhi(phi), pf_sub1_3.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub2_1.size(); ++j){
      if( pf_sub2_1.at(j).phi() > M_PI) phi = pf_sub2_1.at(j).phi() - 2*M_PI;
      else phi = pf_sub2_1.at(j).phi();
      SubJet2_1->Fill(pf_sub2_1.at(j).eta(), ShiftPhi(phi), pf_sub2_1.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub2_2.size(); ++j){
      if( pf_sub2_2.at(j).phi() > M_PI) phi = pf_sub2_2.at(j).phi() - 2*M_PI;
      else phi = pf_sub2_2.at(j).phi();
      SubJet2_2->Fill(pf_sub2_2.at(j).eta(), ShiftPhi(phi), pf_sub2_2.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub1_0.size(); ++j){
      if( pf_sub1_0.at(j).phi() > M_PI) phi = pf_sub1_0.at(j).phi() - 2*M_PI;
      else phi = pf_sub1_0.at(j).phi();
      NotClustered_sub1->Fill(pf_sub1_0.at(j).eta(), ShiftPhi(phi), pf_sub1_0.at(j).energy());
    }
    for(unsigned int j=0; j<pf_sub2_0.size(); ++j){
      if( pf_sub2_0.at(j).phi() > M_PI) phi = pf_sub2_0.at(j).phi() - 2*M_PI;
      else phi = pf_sub2_0.at(j).phi();
      NotClustered_sub2->Fill(pf_sub2_0.at(j).eta(), ShiftPhi(phi), pf_sub2_0.at(j).energy());
    }
  }
  DecayProducts->Fill(bot.eta(), ShiftPhi(bot.phi()), bot.energy());
  DecayProducts->Fill(bot_lep.eta(), ShiftPhi(bot_lep.phi()), bot_lep.energy());
  DecayProducts->Fill(q1.eta(), ShiftPhi(q1.phi()), q1.energy());
  DecayProducts->Fill(q2.eta(), ShiftPhi(q2.phi()), q2.energy());
  DecayProducts->Fill(lep1.eta(), ShiftPhi(lep1.phi()), lep1.energy());
  DecayProducts->Fill(lep2.eta(), ShiftPhi(lep2.phi()), lep2.energy());
  DecayProducts_had->Fill(bot.eta(), ShiftPhi(bot.phi()), bot.energy());
  DecayProducts_had->Fill(q1.eta(), ShiftPhi(q1.phi()), q1.energy());
  DecayProducts_had->Fill(q2.eta(), ShiftPhi(q2.phi()), q2.energy());
  //DecayProducts_lep->Fill(lep1.eta(), ShiftPhi(lep1.phi()), lep1.energy());
  //DecayProducts_lep->Fill(lep2.eta(), ShiftPhi(lep2.phi()), lep2.energy());
  DecayProducts_lep->Fill(bot_lep.eta(), ShiftPhi(bot_lep.phi()), bot_lep.energy());

  Bottom_had->Fill(bot.eta(), ShiftPhi(bot.phi()), bot.energy());
  Bottom_lep->Fill(bot_lep.eta(), ShiftPhi(bot_lep.phi()), bot_lep.energy());

  Lepton->Fill(lepton.eta(), ShiftPhi(lepton.phi()), lepton.energy());
  Neutrino->Fill(neutrino.eta(), ShiftPhi(neutrino.phi()), neutrino.energy());

  Top->Fill(toplep.eta(), ShiftPhi(toplep.phi()), toplep.energy());
  Top->Fill(tophad.eta(), ShiftPhi(tophad.phi()), tophad.energy());
  Top_lep->Fill(toplep.eta(), ShiftPhi(toplep.phi()), toplep.energy());
  Top_had->Fill(tophad.eta(), ShiftPhi(tophad.phi()), tophad.energy());

  for(unsigned int i = 0; i< genparts2.size(); ++i){
    Eta->Fill(genparts2.at(i).eta(), 1);
  }

  double phi1, phi2;
  if(incjets.at(0).phi() > M_PI) phi1 = M_PI - incjets.at(0).phi();
  if(incjets.at(0).phi() <= M_PI) phi1 = incjets.at(0).phi();
  if(incjets.at(1).phi() > M_PI) phi2 = M_PI - incjets.at(1).phi();
  if(incjets.at(1).phi() <= M_PI) phi2 = incjets.at(1).phi();

  cout << "  top had pt = " << tophad.pt() << endl;
  cout << "  top lep pt = " << toplep.pt() << endl;
  cout << "  ----------------------------" << endl;

  if(lep_in_jet2){
    HadJet->Fill(incjets.at(0).eta(), ShiftPhi(phi1), incjets.at(0).energy());
    LepJet->Fill(incjets.at(1).eta(), ShiftPhi(phi2), incjets.at(1).energy());
    cout << "  had jet pt   = " << incjets.at(0).pt() << endl;
    cout << "  had jet mass = " << incjets.at(0).v4().M() << endl;
    cout << "  lep jet pt   = " << incjets.at(1).pt() << endl;
    cout << "  lep jet mass = " << incjets.at(1).v4().M() << endl;
  }

  if(lep_in_jet1){
    HadJet->Fill(incjets.at(1).eta(), ShiftPhi(phi2), incjets.at(1).energy());
    LepJet->Fill(incjets.at(0).eta(), ShiftPhi(phi1), incjets.at(0).energy());
    cout << "  had jet pt   = " << incjets.at(1).pt() << endl;
    cout << "  had jet mass = " << incjets.at(1).v4().M() << endl;
    cout << "  lep jet pt   = " << incjets.at(0).pt() << endl;
    cout << "  lep jet mass = " << incjets.at(0).v4().M() << endl;
  }

  for(auto sub: subjets){
    cout << "  ----------------------------" << endl;
    cout << "  subjet pt    = " << sub.pt() << endl;
    cout << "  subjet eta   = " << sub.eta() << endl;
    cout << "  subjet phi   = " << ShiftPhi(sub.phi()) << endl;
  }
  cout << "  ----------------------------" << endl;
  cout << "  lepton pt  = " << lepton.pt() << endl;
  cout << "  lepton eta = " << lepton.eta() << endl;
  cout << "----------------------------------------------------" << endl;


  PT->Fill(1, jets.at(0).pt());
  PT->Fill(2, jets.at(1).pt());

  if(lep_in_jet1) LepJetBool->Fill(1,1);
  if(lep_in_jet2) LepJetBool->Fill(2,1);


  Mass->Fill(1, jets.at(0).v4().M());
  Mass->Fill(2, jets.at(1).v4().M());

  // cout<<"(px, py, pz, E) = ("<< jets.at(0).v4().Px() <<", "<< jets.at(0).v4().Py() <<", "<< jets.at(0).v4().Pz() <<", "<< jets.at(0).v4().E() <<")"<<endl;
  // cout<<"(pt, eta, phi, E) = ("<< jets.at(0).v4().Pt() <<", "<< jets.at(0).v4().Eta() <<", "<< jets.at(0).v4().Phi() <<", "<< jets.at(0).v4().E() <<")"<<endl;
  // cout<<"Mass = "<< jets.at(0).v4().M() <<endl;
  // cout<<"Mass^2 = "<< jets.at(0).v4().M2() <<endl;

}

fastjet::PseudoJet ClusteringHists::convert_particle(GenParticle* genparticle){
  TLorentzVector particle;
  particle.SetPtEtaPhiE(genparticle->pt(),genparticle->eta(),genparticle->phi(),genparticle->energy());
  fastjet::PseudoJet gen_particle(particle.Px(),particle.Py(),particle.Pz(),particle.E());
  return gen_particle;
}

std::vector<fastjet::PseudoJet> ClusteringHists::add_ghosts(std::vector<fastjet::PseudoJet> gen_in){
  //fastjet::PseudoJet ghost;
  double pt, eta, phi, E, p;
  TLorentzVector ghost_v4;
  for(unsigned int i=0; i < 200; ++i){
    for(unsigned int i_=0; i_ < 200; ++i_ ){
      phi = -M_PI + (i+0.5)*(2*M_PI/200);
      eta = -M_PI + (i_+0.5)*(2*M_PI/200);
      // create random energy
      std::mt19937 rng( std::random_device{}() );
      std::uniform_real_distribution<> dist(0.1, 1);
      double randnr = dist(rng);
      E = randnr*0.0001;
      p = sqrt(E*E);
      pt = p/cosh(eta);
      ghost_v4.SetPtEtaPhiE(pt, eta, phi, E);
      fastjet::PseudoJet ghost( ghost_v4.Px(), ghost_v4.Py(), ghost_v4.Pz(), ghost_v4.E() );
      gen_in.push_back(ghost);
    }
  }
  return gen_in;
}

double ClusteringHists::ShiftPhi(double phi){
  double shift = phishift;
  double newphi = phi + shift;
  if(newphi >  M_PI) newphi = newphi - 2*M_PI;
  if(newphi < -M_PI) newphi = newphi + 2*M_PI;
  return newphi;
}
