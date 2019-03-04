#include "UHH2/MTopJet/include/GenHists_particles.h"


GenHists_particles::GenHists_particles(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // book all histograms here
  number_top = book<TH1F>("number_top", "number of tops", 10, 0, 10);
  hadtop_pt = book<TH1F>("hadtop_pt", "p_{T}", 100, 0, 1000);
  leptop_pt = book<TH1F>("leptop_pt", "p_{T}", 100, 0, 1000);
  hadtop_mass = book<TH1F>("hadtop_mass", "Mass", 150, 100, 250);
  leptop_mass = book<TH1F>("leptop_mass", "Mass", 150, 100, 250);
  hadtop_phi = book<TH1F>("hadtop_phi", "#Phi", 30, -6, 6);
  leptop_phi = book<TH1F>("leptop_phi", "#Phi", 30, -6, 6);
  lepton_pt = book<TH1F>("lepton_pt", "p_{T}", 100, 0, 1000);
  deltaR_hadtop_b = book<TH1F>("deltaR_hadtop_b", "#Delta R(had. top, b)", 30, 0, 6);
  deltaR_leptop_b = book<TH1F>("deltaR_leptop_b", "#Delta R(lep. top, b)", 30, 0, 6);
  deltaR_lep_b = book<TH1F>("deltaR_lep_b", "#Delta R(lep, b)", 30, 0, 6);
  deltaR_lep_neu = book<TH1F>("deltaR_lep_neu", "#Delta R(lepton, neutrino)", 30, 0, 6);
  deltaR_hadtop_leptop = book<TH1F>("deltaR_hadtop_leptop", "#Delta R(had. top, lep. top)", 30, 0, 6);
  deltaPhi_hadtop_leptop = book<TH1F>("deltaPhi_hadtop_leptop", "#Delta #Phi(had. top, lep. top)", 60, 0, 6);
  deltaR_bot_jet1 = book<TH1F>("deltaR_bot_jet1", "#Delta R(bot, jet1)", 60, 0, 6);
  deltaR_hadtop_jet1 = book<TH1F>("deltaR_hadtop_jet1", "#Delta R(had. top, jet1)", 60, 0, 6);
  deltaPT_hadtop_jet1 = book<TH1F>("deltaPT_hadtop_jet1", "p_{T,jet} - p_{T,top}", 50, -100, 100);
  deltaR_hadtop_genjet1 = book<TH1F>("deltaR_hadtop_genjet1", "#Delta R(had. top, gen jet)", 60, 0, 6);
  deltaPT_hadtop_genjet1 = book<TH1F>("deltaPT_hadtop_genjet1", "p_{T,gen jet} - p_{T,top}", 50, -100, 100);
  MinDeltaR_bot_subjet = book<TH1F>("MinDeltaR_bot_subjet", "min #Delta R(subjet, b)", 30, 0, 6);
  InsideFat_MinDeltaR_bot_subjet = book<TH1F>("InsideFat_MinDeltaR_bot_subjet", "min #Delta R(subjet, b) only if b in fatjet", 30, 0, 6);
  b_inside_fatjet = book<TH1F>("b_inside_fatjet", "b inside fatjet", 2, -0.5, 1.5);
  bjet_index = book<TH1F>("bjet_index", "jet index including b-quark", 4, -1.5, 2.5);
  bjet_index_nomatch = book<TH1F>("bjet_index_nomatch", "no match found", 1, -0.5, 0.5);

  pthadtop_ptleptop = book<TH2F>("pthadtop_ptleptop", "x=had. top p_{T} y=lep top p_{T}", 100, 0, 1000, 100, 0, 1000);


  dR_pt_b_subjet = book<TH2F>("dR_pt_b_subjet", "x=#Delta R(b, subjet) y=p_{T} subjet", 60, 0, 6, 60, 0, 600);
  dR_pt_ud_subjet = book<TH2F>("dR_pt_ud_subjet", "x=#Delta R(ud, subjet) y=p_{T} subjet", 60, 0, 6, 60, 0, 600);
  dR_pt_g_subjet = book<TH2F>("dR_pt_g_subjet", "x=#Delta R(g, subjet) y=p_{T} subjet", 60, 0, 6, 60, 0, 600);

  closest_particle_0 = book<TH1F>("closest_particle_0", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_50 = book<TH1F>("closest_particle_50", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_100 = book<TH1F>("closest_particle_100", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_200 = book<TH1F>("closest_particle_200", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_300 = book<TH1F>("closest_particle_300", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_400 = book<TH1F>("closest_particle_400", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);
  closest_particle_500 = book<TH1F>("closest_particle_500", "closest particle 0:none, 1:gluon, 2:light, 3:b", 4, -0.5, 3.5);

  fraction_ud_A = book<TH1F>("fraction_ud_A", "subjet p_T", 120, 0, 600);
  fraction_s_A = book<TH1F>("fraction_s_A", "subjet p_T", 120, 0, 600);
  fraction_c_A = book<TH1F>("fraction_c_A", "subjet p_T", 120, 0, 600);
  fraction_b_A = book<TH1F>("fraction_b_A", "subjet p_T", 120, 0, 600);
  fraction_g_A = book<TH1F>("fraction_g_A", "subjet p_T", 120, 0, 600);
  fraction_nomatch_A = book<TH1F>("fraction_nomatch_A", "subjet p_T", 120, 0, 600);

  fraction_ud_B = book<TH1F>("fraction_ud_B", "subjet p_T", 120, 0, 600);
  fraction_s_B = book<TH1F>("fraction_s_B", "subjet p_T", 120, 0, 600);
  fraction_c_B = book<TH1F>("fraction_c_B", "subjet p_T", 120, 0, 600);
  fraction_b_B = book<TH1F>("fraction_b_B", "subjet p_T", 120, 0, 600);
  fraction_g_B = book<TH1F>("fraction_g_B", "subjet p_T", 120, 0, 600);

  fraction_ud_C = book<TH1F>("fraction_ud_C", "subjet p_T", 120, 0, 600);
  fraction_s_C = book<TH1F>("fraction_s_C", "subjet p_T", 120, 0, 600);
  fraction_c_C = book<TH1F>("fraction_c_C", "subjet p_T", 120, 0, 600);
  fraction_b_C = book<TH1F>("fraction_b_C", "subjet p_T", 120, 0, 600);
  fraction_g_C = book<TH1F>("fraction_g_C", "subjet p_T", 120, 0, 600);

  N_matches = book<TH1F>("N_matches", "number of particles matched to subjet", 21, -0.5, 20.5);

  // handle for ttbargen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  // handle for jets
  h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");
  h_hadjets_gen=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");

}



void GenHists_particles::fill(const Event & event){

  const auto & ttbargen = event.get(h_ttbargen);
  // get weight
  double weight = event.weight;

  // get jets
  std::vector<TopJet> hadjets = event.get(h_hadjets);
  std::vector<Jet> subjets = hadjets[0].subjets();
  std::vector<GenTopJet> hadjets_gen = event.get(h_hadjets_gen);

  // cout tops
  int n_top = 0;
  std::vector<GenParticle>* genparts = event.genparticles;
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle p = genparts->at(i);
    if(abs(p.pdgId()) == 6){
      ++n_top;
    }
  }

  number_top->Fill(n_top, weight);

  // get particles from ttbargen class
  GenParticle tophad, toplep, bot, q1, q2, bot_lep, lep1, lep2, lepton, neutrino;
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
    return;
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

  hadtop_pt->Fill(tophad.pt(), weight);
  leptop_pt->Fill(toplep.pt(), weight);
  pthadtop_ptleptop->Fill(tophad.pt(), toplep.pt(), weight);
  hadtop_mass->Fill(tophad.v4().M(), weight);
  leptop_mass->Fill(toplep.v4().M(), weight);
  lepton_pt->Fill(lepton.pt(), weight);
  hadtop_phi->Fill(tophad.phi(), weight);
  leptop_phi->Fill(toplep.phi(), weight);
  deltaR_hadtop_b->Fill(deltaR(tophad, bot), weight);
  deltaR_leptop_b->Fill(deltaR(toplep, bot_lep), weight);
  deltaR_lep_b->Fill(deltaR(lepton, bot_lep), weight);
  deltaR_lep_neu->Fill(deltaR(lepton, neutrino), weight);
  deltaR_hadtop_leptop->Fill(deltaR(tophad, toplep), weight);
  deltaR_bot_jet1->Fill(deltaR(bot, hadjets.at(0)), weight);

  if(hadjets.size()>0){
    deltaR_hadtop_jet1->Fill(deltaR(tophad, hadjets.at(0)), weight);
    deltaPT_hadtop_jet1->Fill(hadjets.at(0).pt() - tophad.pt(), weight);
  }
  if(hadjets_gen.size()>0){
    deltaR_hadtop_genjet1->Fill(deltaR(tophad, hadjets_gen.at(0)), weight);
    deltaPT_hadtop_genjet1->Fill(hadjets_gen.at(0).pt() - tophad.pt(), weight);
  }

  // check if the b-quark ends up in the pT-leading subjet
  double dR_min = 1000;
  double dR_temp;
  int b_index = -1;
  // sort jets by pT and get index of b-jet
  sort_by_pt<Jet> (subjets);
  for(unsigned int i=0; i<subjets.size(); i++){
    dR_temp = deltaR(bot, subjets[i]);
    if(dR_temp < dR_min){
      dR_min = dR_temp;
      if(dR_min < 0.4) b_index = i; // b is only matched if minimum dR is smaller 0.4
    }
  }
  MinDeltaR_bot_subjet->Fill(dR_min, weight);
  bjet_index->Fill(b_index, weight);
  if(b_index == -1) bjet_index_nomatch->Fill(0.0, weight);
  ////

  // this is for measuring the flavor fraction
  // match genparts to reco subjets and fill subjet pT if matched
  // A: just count jet to flavor if dR(parton, jet) < 0.4
  for(unsigned int j=0; j<subjets.size(); j++){
    Jet subjet = subjets[j];
    double pt = subjet.pt();
    int nMatches = 0;
    for (unsigned int i=0; i < genparts->size(); ++i){
      GenParticle p = genparts->at(i);
      if(i==0 || i==1) continue;        // skip initial particles
      if(abs(p.pdgId()) == 6) continue; // skip top quarks
      if(deltaR(subjet, p) < 0.4){
        if(abs(p.pdgId()) == 1 || abs(p.pdgId()) == 2){
          fraction_ud_A->Fill(pt, weight);
          nMatches ++;
        }
        else if(abs(p.pdgId()) == 3){
          fraction_s_A->Fill(pt, weight);
          nMatches ++;
        }
        else if(abs(p.pdgId()) == 4){
          fraction_c_A->Fill(pt, weight);
          nMatches ++;
        }
        else if(abs(p.pdgId()) == 5){
          fraction_b_A->Fill(pt, weight);
          nMatches ++;
        }
        else if(abs(p.pdgId()) == 9 || abs(p.pdgId()) == 21){
          fraction_g_A->Fill(pt, weight);
          nMatches ++;
        }
      }
    }
    N_matches->Fill(nMatches, weight);
    if(nMatches == 0){
      fraction_nomatch_A->Fill(pt, weight);
    }
  }
  ////

  // B: just count flavor with highest pT
  // C: zähle parton nur mit pT Anteil an allen gematchten partonen
  for(unsigned int j=0; j<subjets.size(); j++){
    Jet subjet = subjets[j];
    double pt = subjet.pt();
    double E_ud = 0;
    double E_s = 0;
    double E_c = 0;
    double E_b = 0;
    double E_g = 0;
    double E_ges = 0;
    for (unsigned int i=0; i < genparts->size(); ++i){
      GenParticle p = genparts->at(i);
      if(i==0 || i==1) continue;        // skip initial particles
      if(abs(p.pdgId()) == 6) continue; // skip top quarks
      if(deltaR(subjet, p) < 0.4){
        if(abs(p.pdgId()) == 1 || abs(p.pdgId()) == 2){
          E_ud += p.energy();
          E_ges += p.energy();
        }
        else if(abs(p.pdgId()) == 3){
          E_s += p.energy();
          E_ges += p.energy();
        }
        else if(abs(p.pdgId()) == 4){
          E_c += p.energy();
          E_ges += p.energy();
        }
        else if(abs(p.pdgId()) == 5){
          E_b += p.energy();
          E_ges += p.energy();
        }
        else if(abs(p.pdgId()) == 9 || abs(p.pdgId()) == 21){
          E_g += p.energy();
          E_ges += p.energy();
        }
      }
    }
    if(E_ud > E_s  && E_ud > E_c  && E_ud > E_b  && E_ud > E_g) fraction_ud_B->Fill(pt, weight);
    if(E_s  > E_ud && E_s  > E_c  && E_s  > E_b  && E_s  > E_g) fraction_s_B->Fill(pt, weight);
    if(E_c  > E_ud && E_c  > E_s  && E_c  > E_b  && E_c  > E_g) fraction_c_B->Fill(pt, weight);
    if(E_b  > E_ud && E_b  > E_s  && E_b  > E_c  && E_b  > E_g) fraction_b_B->Fill(pt, weight);
    if(E_g  > E_ud && E_g  > E_s  && E_g  > E_c  && E_g  > E_b) fraction_g_B->Fill(pt, weight);

    if(E_ges != 0){
      fraction_ud_C->Fill(pt, weight*E_ud/E_ges);
      fraction_s_C->Fill(pt, weight*E_s/E_ges);
      fraction_c_C->Fill(pt, weight*E_c/E_ges);
      fraction_b_C->Fill(pt, weight*E_b/E_ges);
      fraction_g_C->Fill(pt, weight*E_g/E_ges);
    }
  }
  ////

  // plot dR of quarks to subjets vs subjet pT to get flavor fraktions depending on pT
  // first get the gluon from the radiation
  // This is always the entry with index 4
  GenParticle gluon = genparts->at(4);

  for(unsigned int i=0; i<subjets.size(); i++){
    double pt = subjets[i].pt();
    double dR_b = deltaR(subjets[i], bot);
    dR_pt_b_subjet->Fill(dR_b, pt, weight);

    double dR_ud1=100, dR_ud2=100, dR_ud=100;
    if(abs(q1.pdgId()) == 1 || abs(q1.pdgId()) == 2 ){
      dR_ud1 = deltaR(subjets[i], q1);
      dR_pt_ud_subjet->Fill(dR_ud1, pt, weight);
    }
    if(abs(q2.pdgId()) == 1 || abs(q2.pdgId()) == 2 ){
      dR_ud2 = deltaR(subjets[i], q2);
      dR_pt_ud_subjet->Fill(dR_ud2, pt, weight);
    }
    double dR_g = deltaR(subjets[i], gluon);
    dR_pt_g_subjet->Fill(dR_g, pt, weight);


    // get closest particle to subjet
    if(dR_ud1 < dR_ud2) dR_ud = dR_ud1;
    else                dR_ud = dR_ud2;
    int index_closest_particle = 0; // 1:gluon, 2:light, 3:b
    if(dR_g < dR_ud && dR_g < dR_b && dR_g < 0.4) index_closest_particle = 1;
    if(dR_ud < dR_g && dR_ud < dR_b && dR_ud < 0.4) index_closest_particle = 2;
    if(dR_b < dR_g && dR_b < dR_ud && dR_b < 0.4) index_closest_particle = 3;
    if(pt >   0 && pt <  50) closest_particle_0->Fill(index_closest_particle, weight);
    if(pt >  50 && pt < 100) closest_particle_50->Fill(index_closest_particle, weight);
    if(pt > 100 && pt < 200) closest_particle_100->Fill(index_closest_particle, weight);
    if(pt > 200 && pt < 300) closest_particle_200->Fill(index_closest_particle, weight);
    if(pt > 300 && pt < 400) closest_particle_300->Fill(index_closest_particle, weight);
    if(pt > 400 && pt < 500) closest_particle_400->Fill(index_closest_particle, weight);
    if(pt > 500 && pt < 600) closest_particle_500->Fill(index_closest_particle, weight);

  }

  // is b quark inside fatjet?
  bool b_in_fat = false;
  if(deltaR(bot, hadjets.at(0)) < 1.2) b_in_fat = true;
  if(b_in_fat) b_inside_fatjet->Fill(1.0, weight);
  else         b_inside_fatjet->Fill(0.0, weight);

  // is b quark inside fatjet but still not matched
  if(b_in_fat) InsideFat_MinDeltaR_bot_subjet->Fill(dR_min, weight);
  ////


  double had_phi, lep_phi, delta_phi;
  had_phi = tophad.phi() + M_PI;
  lep_phi = toplep.phi() + M_PI;
  delta_phi = fabs(had_phi - lep_phi);
  deltaPhi_hadtop_leptop->Fill(delta_phi, weight);

  return;
}
