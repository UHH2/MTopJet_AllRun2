#include "UHH2/MTopJet/include/GenSelections.h"
#include <iostream>


uhh2::TTbarSemilep::TTbarSemilep(uhh2::Context& ctx):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}

bool uhh2::TTbarSemilep::passes(const uhh2::Event& event){

    const auto & ttbargen = event.get(h_ttbargen);
    bool semilep = false;

    if(ttbargen.DecayChannel() == TTbarGen::e_muhad || ttbargen.DecayChannel() == TTbarGen::e_ehad){
      semilep = true;
    }
    return semilep;
}


uhh2::TopHadPT::TopHadPT(uhh2::Context& ctx, float ptmin):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")), ptmin_(ptmin) {}


bool uhh2::TopHadPT::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    bool pass_pt = true;
    GenParticle tophad = ttbargen.TopHad();
    float toppt = tophad.pt();
    if(toppt < ptmin_) pass_pt=false;
    return pass_pt;
}


uhh2::Matching::Matching(uhh2::Context& ctx, const std::string & name, float jetradius):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}


bool uhh2::Matching::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    std::vector<Jet> jets = event.get(h_jets);
    Jet jet1;
    if(jets.size()>0) jet1 = jets.at(0);
    bool matched = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, q1, q2; 
    if(jets.size() > 0){
      if(ttbargen.IsTopHadronicDecay()){
	bot = ttbargen.bTop();
	q1 = ttbargen.Wdecay1();
	q2 = ttbargen.Wdecay2();
      }
      else if(ttbargen.IsAntiTopHadronicDecay()){
	bot = ttbargen.bAntitop();
	q1 = ttbargen.WMinusdecay1();
	q2 = ttbargen.WMinusdecay2();
      }
    }
      //check if particles from hadronic top are clustered into jet
      if((deltaR(bot, jet1)<=jetradius_) && (deltaR(q1, jet1)<=jetradius_) && (deltaR(q2, jet1)<=jetradius_)) matched = true;
      return matched;
}

uhh2::Matching_HOTVR::Matching_HOTVR(uhh2::Context& ctx, const std::string & name, double rho):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  rho_(rho) {}


bool uhh2::Matching_HOTVR::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    std::vector<Jet> jets = event.get(h_jets);
    Jet jet1;
    TLorentzVector jet1_v4;
    if(jets.size()>0) jet1 = jets.at(0);
    bool matched = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, q1, q2; 
    if(jets.size() > 0){
      if(ttbargen.IsTopHadronicDecay()){
	bot = ttbargen.bTop();
	q1 = ttbargen.Wdecay1();
	q2 = ttbargen.Wdecay2();
      }
      else if(ttbargen.IsAntiTopHadronicDecay()){
	bot = ttbargen.bAntitop();
	q1 = ttbargen.WMinusdecay1();
	q2 = ttbargen.WMinusdecay2();
      }
    }

    jet1_v4.SetPxPyPzE(jet1.v4().Px(), jet1.v4().Py(), jet1.v4().Pz(), jet1.v4().E());
    double pt = jet1_v4.Pt();
    double Reff = rho_/pt;
    if(Reff > 1.5) Reff = 1.5;
    if(Reff < 0.1) Reff = 0.1;

    //check if particles from hadronic top are clustered into jet
    if((deltaR(bot, jet1)<=Reff) && (deltaR(q1, jet1)<=Reff) && (deltaR(q2, jet1)<=Reff)) matched = true;
    return matched;
}

uhh2::Matching_XCone::Matching_XCone(uhh2::Context& ctx, const std::string & name):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}


bool uhh2::Matching_XCone::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    std::vector<Jet> jets = event.get(h_jets);
    Jet jet1, jet2, jet3, jet4;
    if(jets.size()>0) jet1 = jets.at(0);
    if(jets.size()>1) jet2 = jets.at(1);
    if(jets.size()>2) jet3 = jets.at(2);
    if(jets.size()>3) jet4 = jets.at(3);
    bool matched = false;
    bool matched_q1 = false;
    bool matched_q2 = false;
    bool matched_bot = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, q1, q2, lep1, lep2, lepton; 
    if(jets.size() > 0){
      if(ttbargen.IsTopHadronicDecay()){
	bot = ttbargen.bTop();
	q1 = ttbargen.Wdecay1();
	q2 = ttbargen.Wdecay2();
	lep1 = ttbargen.WMinusdecay1();
	lep2 = ttbargen.WMinusdecay2();
      }
      else if(ttbargen.IsAntiTopHadronicDecay()){
	bot = ttbargen.bAntitop();
	q1 = ttbargen.WMinusdecay1();
	q2 = ttbargen.WMinusdecay2();
	lep1 = ttbargen.Wdecay1();
	lep2 = ttbargen.Wdecay2();
      }
      if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
	lepton = lep1;
      }
      else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
	lepton = lep2;
      }
    }

    // claculate distance to Lepton
    double dR1, dR2, dR3, dR4;
    dR1 = deltaR(jet1, lepton);
    dR2 = deltaR(jet2, lepton);
    dR3 = deltaR(jet3, lepton);
    dR4 = deltaR(jet4, lepton);
    if(dR1 < dR2 && dR1 < dR3 && dR1 < dR4){
      if(deltaR(q1, jet2) < 0.4 || deltaR(q1, jet3) < 0.4 || deltaR(q1, jet4) < 0.4) matched_q1 = true;
      if(deltaR(q2, jet2) < 0.4 || deltaR(q2, jet3) < 0.4 || deltaR(q2, jet4) < 0.4) matched_q2 = true;
      if(deltaR(bot, jet2) < 0.4 || deltaR(bot, jet3) < 0.4 || deltaR(bot, jet4) < 0.4) matched_bot = true;
    }
    if(dR2 < dR1 && dR2 < dR3 && dR2 < dR4){
      if(deltaR(q1, jet1) < 0.4 || deltaR(q1, jet3) < 0.4 || deltaR(q1, jet4) < 0.4) matched_q1 = true;
      if(deltaR(q2, jet1) < 0.4 || deltaR(q2, jet3) < 0.4 || deltaR(q2, jet4) < 0.4) matched_q2 = true;
      if(deltaR(bot, jet1) < 0.4 || deltaR(bot, jet3) < 0.4 || deltaR(bot, jet4) < 0.4) matched_bot = true;
    }    
    if(dR3 < dR1 && dR3 < dR2 && dR3 < dR4){
      if(deltaR(q1, jet1) < 0.4 || deltaR(q1, jet2) < 0.4 || deltaR(q1, jet4) < 0.4) matched_q1 = true;
      if(deltaR(q2, jet1) < 0.4 || deltaR(q2, jet2) < 0.4 || deltaR(q2, jet4) < 0.4) matched_q2 = true;
      if(deltaR(bot, jet1) < 0.4 || deltaR(bot, jet2) < 0.4 || deltaR(bot, jet4) < 0.4) matched_bot = true;
    }
    if(dR4 < dR1 && dR4 < dR2 && dR4 < dR3){
      if(deltaR(q1, jet1) < 0.4 || deltaR(q1, jet2) < 0.4 || deltaR(q1, jet3) < 0.4) matched_q1 = true;
      if(deltaR(q2, jet1) < 0.4 || deltaR(q2, jet2) < 0.4 || deltaR(q2, jet3) < 0.4) matched_q2 = true;
      if(deltaR(bot, jet1) < 0.4 || deltaR(bot, jet2) < 0.4 || deltaR(bot, jet3) < 0.4) matched_bot = true;
    }

    if(matched_q1 && matched_q2 && matched_bot) matched = true;
    return matched;
}

uhh2::Matching_XCone23::Matching_XCone23(uhh2::Context& ctx, const std::string & name):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}


bool uhh2::Matching_XCone23::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    std::vector<Jet> jets = event.get(h_jets);
    Jet jet1, jet2, jet3;
    if(jets.size()>0) jet1 = jets.at(0);
    if(jets.size()>1) jet2 = jets.at(1);
    if(jets.size()>2) jet3 = jets.at(2);
    bool matched = false;
    bool matched_q1 = false;
    bool matched_q2 = false;
    bool matched_bot = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, q1, q2, lep1, lep2, lepton; 
    if(ttbargen.IsTopHadronicDecay()){
      bot = ttbargen.bTop();
      q1 = ttbargen.Wdecay1();
      q2 = ttbargen.Wdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      bot = ttbargen.bAntitop();
      q1 = ttbargen.WMinusdecay1();
      q2 = ttbargen.WMinusdecay2();
    }
    

    if(deltaR(q1, jet1) < 0.4 || deltaR(q1, jet2) < 0.4 || deltaR(q1, jet3) < 0.4) matched_q1 = true;
    if(deltaR(q2, jet1) < 0.4 || deltaR(q2, jet2) < 0.4 || deltaR(q2, jet3) < 0.4) matched_q2 = true;
    if(deltaR(bot, jet1) < 0.4 || deltaR(bot, jet2) < 0.4 || deltaR(bot, jet3) < 0.4) matched_bot = true;
    
    if(matched_q1 && matched_q2 && matched_bot) matched = true;
    return matched;
}

uhh2::Matching_XCone_botlep_lep::Matching_XCone_botlep_lep(uhh2::Context& ctx, const std::string & name):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}


bool uhh2::Matching_XCone_botlep_lep::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    bool botlep_lep = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, bot_lep, lep1, lep2, lepton, q1, q2; 
    if(ttbargen.IsTopHadronicDecay()){
      bot = ttbargen.bTop();
      bot_lep = ttbargen.bAntitop();
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
      q1 = ttbargen.Wdecay1();
      q2 = ttbargen.Wdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      bot = ttbargen.bAntitop();
      bot_lep = ttbargen.bTop();
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
      q1 = ttbargen.WMinusdecay1();
      q2 = ttbargen.WMinusdecay2();
    }
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    

    // claculate distance to Lepton
    double dR_lep_bot, dR_lep_botlep, dR_lep_q1, dR_lep_q2;
    dR_lep_bot = deltaR(bot, lepton);
    dR_lep_botlep = deltaR(bot_lep, lepton);
    dR_lep_q1 = deltaR(q1, lepton);
    dR_lep_q2 = deltaR(q2, lepton);
    if(dR_lep_botlep < dR_lep_bot && dR_lep_botlep < dR_lep_q1 && dR_lep_botlep < dR_lep_q1) botlep_lep = true;
    if(dR_lep_bot < dR_lep_botlep && dR_lep_bot < dR_lep_q1 && dR_lep_bot < dR_lep_q1) botlep_lep = false;
    if(dR_lep_q1 < dR_lep_botlep && dR_lep_q1 < dR_lep_bot && dR_lep_q1 < dR_lep_q2) botlep_lep = false;
    if(dR_lep_q2 < dR_lep_botlep && dR_lep_q2 < dR_lep_bot && dR_lep_q2 < dR_lep_q1) botlep_lep = false;

    return botlep_lep;
}

uhh2::Matching_top::Matching_top(uhh2::Context& ctx, float jetradius):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}


bool uhh2::Matching_top::passes(const uhh2::Event& event){
    const auto & ttbargen = event.get(h_ttbargen);
    std::vector<GenTopJet> jets = event.get(h_gentopjets);
    GenTopJet jet1;
    if(jets.size()>0) jet1 = jets.at(0);
    bool matched = false;
    // get stable particles from ttbar decay and sort them into leptonic and hadronic
    GenParticle bot, q1, q2; 
    if(jets.size() > 0){
      if(ttbargen.IsTopHadronicDecay()){
	bot = ttbargen.bTop();
	q1 = ttbargen.Wdecay1();
	q2 = ttbargen.Wdecay2();
      }
      else if(ttbargen.IsAntiTopHadronicDecay()){
	bot = ttbargen.bAntitop();
	q1 = ttbargen.WMinusdecay1();
	q2 = ttbargen.WMinusdecay2();
      }
    }
      //check if particles from hadronic top are clustered into jet
      if((deltaR(bot, jet1)<=jetradius_) && (deltaR(q1, jet1)<=jetradius_) && (deltaR(q2, jet1)<=jetradius_)) matched = true;
      return matched;
}

uhh2::LeadingJetPT::LeadingJetPT(uhh2::Context& ctx, const std::string & name, float ptcut):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  ptcut_(ptcut) {}

bool uhh2::LeadingJetPT::passes(const uhh2::Event& event){
  bool pass_jetpt = false;
  std::vector<Jet> jets = event.get(h_jets);
  Jet jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float pt = jet1.pt();
    if(pt > ptcut_) pass_jetpt = true;
  }
  return pass_jetpt;
}

uhh2::LeadingJetPT_gen::LeadingJetPT_gen(uhh2::Context& ctx, const std::string & name, float ptcut):
  h_jets(ctx.get_handle<std::vector<Particle>>(name)),
  ptcut_(ptcut) {}

bool uhh2::LeadingJetPT_gen::passes(const uhh2::Event& event){
  bool pass_jetpt = false;
  std::vector<Particle> jets = event.get(h_jets);
  Particle jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float pt = jet1.pt();
    if(pt > ptcut_) pass_jetpt = true;
  }
  return pass_jetpt;
}

uhh2::LeadingTopJetPT::LeadingTopJetPT(uhh2::Context& ctx, float ptcut):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  ptcut_(ptcut) {}

bool uhh2::LeadingTopJetPT::passes(const uhh2::Event& event){
  bool pass_jetpt = false;
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  GenTopJet jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float pt = jet1.pt();
    if(pt > ptcut_) pass_jetpt = true;
  }
  return pass_jetpt;
}

uhh2::MassCut::MassCut(uhh2::Context& ctx, const std::string & name):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}

bool uhh2::MassCut::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  const auto & ttbargen = event.get(h_ttbargen);
  Jet jet1, jet2;
  GenParticle lepton, lep1, lep2;
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  bool pass_masscut = false;
  if(jets.size()>1){
    jet1 = jets.at(0);
    jet1 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(jets.at(1).pt(),jets.at(1).eta(),jets.at(1).phi(),jets.at(1).energy()); //v4 of first jet
    lepton1_v4.SetPtEtaPhiE(lepton.pt(),lepton.eta(),lepton.phi(),lepton.energy()); //v4 of lepton
    jet2_lep_v4 = jet2_v4 + lepton1_v4; // v4 of (lepton+2nd jet)

    if(jet1_v4.M() > jet2_lep_v4.M()) pass_masscut = true;
  }
  return pass_masscut;
}

uhh2::MassCut_gen::MassCut_gen(uhh2::Context& ctx, const std::string & hadname, const std::string & lepname):
  h_hadjets(ctx.get_handle<std::vector<Particle>>(hadname)),
  h_lepjets(ctx.get_handle<std::vector<Particle>>(lepname)){}

bool uhh2::MassCut_gen::passes(const uhh2::Event& event){
  std::vector<Particle> hadjets = event.get(h_hadjets);
  std::vector<Particle> lepjets = event.get(h_lepjets);
  Particle jet1, jet2;
  TLorentzVector jet1_v4, jet2_v4;
  bool pass_masscut = false;
  if(hadjets.size()>0 && lepjets.size()>0){
    jet1_v4.SetPtEtaPhiE(hadjets.at(0).pt(),hadjets.at(0).eta(),hadjets.at(0).phi(),hadjets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(lepjets.at(0).pt(),lepjets.at(0).eta(),lepjets.at(0).phi(),lepjets.at(0).energy()); //v4 of first jet
    if(jet1_v4.M() > jet2_v4.M()) pass_masscut = true;
  }
  return pass_masscut;
}


uhh2::MassCut_top::MassCut_top(uhh2::Context& ctx):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}

bool uhh2::MassCut_top::passes(const uhh2::Event& event){
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  const auto & ttbargen = event.get(h_ttbargen);
  GenTopJet jet1, jet2;
  GenParticle lepton, lep1, lep2;
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  bool pass_masscut = false;
  if(jets.size()>1){
    jet1 = jets.at(0);
    jet1 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(jets.at(1).pt(),jets.at(1).eta(),jets.at(1).phi(),jets.at(1).energy()); //v4 of first jet
    lepton1_v4.SetPtEtaPhiE(lepton.pt(),lepton.eta(),lepton.phi(),lepton.energy()); //v4 of lepton
    jet2_lep_v4 = jet2_v4 + lepton1_v4; // v4 of (lepton+2nd jet)

    if(jet1_v4.M() > jet2_lep_v4.M()) pass_masscut = true;
  }
  return pass_masscut;
}

uhh2::DeltaPhiCut::DeltaPhiCut(uhh2::Context& ctx, const std::string & name, float jetradius):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}

bool uhh2::DeltaPhiCut::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  const auto & ttbargen = event.get(h_ttbargen);
  Jet jet2;
  GenParticle lepton, lep1, lep2;
  bool pass_deltaR = false;
  if(jets.size()>1){
    jet2 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    
    if(abs(lepton.phi() - jet2.phi()) < jetradius_) pass_deltaR = true;
  }
  return pass_deltaR;
}

uhh2::DeltaRCut::DeltaRCut(uhh2::Context& ctx, const std::string & name, float jetradius):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}

bool uhh2::DeltaRCut::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  const auto & ttbargen = event.get(h_ttbargen);
  Jet jet2;
  GenParticle lepton, lep1, lep2;
  bool pass_deltaR = false;
  if(jets.size()>1){
    jet2 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    
    if(deltaR(lepton, jet2) < jetradius_) pass_deltaR = true;
  }
  return pass_deltaR;
}

uhh2::DeltaRCut_HOTVR::DeltaRCut_HOTVR(uhh2::Context& ctx, const std::string & name, double rho):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  rho_(rho) {}

bool uhh2::DeltaRCut_HOTVR::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  const auto & ttbargen = event.get(h_ttbargen);
  Jet jet2;
  TLorentzVector jet2_v4;
  GenParticle lepton, lep1, lep2;
  bool pass_deltaR = false;
  if(jets.size()>1){
    jet2 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }

    jet2_v4.SetPxPyPzE(jet2.v4().Px(), jet2.v4().Py(), jet2.v4().Pz(), jet2.v4().E());
    double pt = jet2_v4.Pt();
    double Reff = rho_/pt;
    if(Reff > 1.5) Reff = 1.5;
    if(Reff < 0.1) Reff = 0.1;

    if(deltaR(lepton, jet2) < Reff) pass_deltaR = true;
  }
  return pass_deltaR;
}

uhh2::DeltaRCut_top::DeltaRCut_top(uhh2::Context& ctx, float jetradius):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}

bool uhh2::DeltaRCut_top::passes(const uhh2::Event& event){
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  const auto & ttbargen = event.get(h_ttbargen);
  GenTopJet jet2;
  GenParticle lepton, lep1, lep2;
  bool pass_deltaR = false;
  if(jets.size()>1){
    jet2 = jets.at(1);
    if(ttbargen.IsTopHadronicDecay()){
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
    
    if(deltaR(lepton, jet2) < jetradius_) pass_deltaR = true;
  }
  return pass_deltaR;
}

uhh2::NGenJets::NGenJets(uhh2::Context& ctx, const std::string & name, float min_pt, float min, float max):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  min_pt_(min_pt),
  min_(min),
  max_(max) {}

bool uhh2::NGenJets::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  bool pass_NGen = false;
  int counter = 0;
  for(unsigned int i = 0; i<jets.size(); ++i){
    if(jets.at(i).pt() >= min_pt_){
      ++counter;
    }
  }
  if(counter >= min_ && counter <= max_) pass_NGen = true;
  return pass_NGen;
}

uhh2::NGenTopJets::NGenTopJets(uhh2::Context& ctx, float min_pt, float min, float max):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  min_pt_(min_pt),
  min_(min),
  max_(max) {}

bool uhh2::NGenTopJets::passes(const uhh2::Event& event){
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  bool pass_NGen = false;
  int counter = 0;
  for(unsigned int i = 0; i<jets.size(); ++i){
    if(jets.at(i).pt() >= min_pt_){
      ++counter;
    }
  }
  if(counter >= min_ && counter <= max_) pass_NGen = true;
  return pass_NGen;
}

GenJetLeptonCleaner::GenJetLeptonCleaner(uhh2::Context& ctx, const std::string & name, float jetradius):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}

bool GenJetLeptonCleaner::process(uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  const auto & ttbargen = event.get(h_ttbargen);
  std::vector<Jet> cleaned_jets;
  Jet jet;
  GenParticle lepton, lep1, lep2;

  // check if lepton is in jet1 or jet2
  if(ttbargen.IsTopHadronicDecay()){
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  //check which lep is neutrino and which is elec/muon
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
  }


  // perform cleaning
  TLorentzVector jet_v4, lepton_v4, jetlep_v4;
  for(unsigned int i = 0; i < jets.size(); ++i){
    jet = jets.at(i);

    if(deltaR(lepton, jet) < jetradius_){

      // std::cout<<"pt before: "<< jet.pt()<<std::endl;

     jet_v4.SetPxPyPzE(jet.v4().Px(), jet.v4().Py(), jet.v4().Pz(), jet.v4().E());
     lepton_v4.SetPxPyPzE(lepton.v4().Px(), lepton.v4().Py(), lepton.v4().Pz(), lepton.v4().E()); 
     jetlep_v4 = jet_v4 - lepton_v4;

     jet.set_pt(jetlep_v4.Pt());
     jet.set_eta(jetlep_v4.Eta());
     jet.set_phi(jetlep_v4.Phi());
     jet.set_energy(jetlep_v4.E());

      // std::cout<<"pt after: "<< jet.pt()<<std::endl;

    }
  cleaned_jets.push_back(jet);

  }
  sort_by_pt<Jet>(cleaned_jets); // Sort Jets by pT
  event.set(h_jets, cleaned_jets);

  return true;
}

GenTopJetLeptonCleaner::GenTopJetLeptonCleaner(uhh2::Context& ctx, float jetradius):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  jetradius_(jetradius) {}

bool GenTopJetLeptonCleaner::process(uhh2::Event& event){
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  const auto & ttbargen = event.get(h_ttbargen);
  std::vector<GenTopJet> cleaned_jets;
  GenTopJet jet;
  GenParticle lepton, lep1, lep2;

  // check if lepton is in jet1 or jet2
  if(ttbargen.IsTopHadronicDecay()){
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  //check which lep is neutrino and which is elec/muon
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
  }


  // perform cleaning
  TLorentzVector jet_v4, lepton_v4, jetlep_v4;
  for(unsigned int i = 0; i < jets.size(); ++i){
    jet = jets.at(i);

    if(deltaR(lepton, jet) < jetradius_){

      // std::cout<<"pt before: "<< jet.pt()<<std::endl;

     jet_v4.SetPxPyPzE(jet.v4().Px(), jet.v4().Py(), jet.v4().Pz(), jet.v4().E());
     lepton_v4.SetPxPyPzE(lepton.v4().Px(), lepton.v4().Py(), lepton.v4().Pz(), lepton.v4().E()); 
     jetlep_v4 = jet_v4 - lepton_v4;

     jet.set_pt(jetlep_v4.Pt());
     jet.set_eta(jetlep_v4.Eta());
     jet.set_phi(jetlep_v4.Phi());
     jet.set_energy(jetlep_v4.E());

      // std::cout<<"pt after: "<< jet.pt()<<std::endl;

    }
  cleaned_jets.push_back(jet);

  }
  sort_by_pt<GenTopJet>(cleaned_jets); // Sort Jets by pT
  event.set(h_gentopjets, cleaned_jets);

  return true;
}


uhh2::MassCutGen1::MassCutGen1(uhh2::Context& ctx, const std::string & name, float M_min, float M_max):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  M_min_(M_min),
  M_max_(M_max) {}

bool MassCutGen1::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  TLorentzVector jet1_v4;
  jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
  return (M_min_ < jet1_v4.M() && jet1_v4.M() < M_max_);
}


uhh2::MassCutGen1_top::MassCutGen1_top(uhh2::Context& ctx, float M_min, float M_max):
  h_gentopjets(ctx.get_handle<std::vector<GenTopJet>>("gentopjets")),
  M_min_(M_min),
  M_max_(M_max) {}

bool MassCutGen1_top::passes(const uhh2::Event& event){
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  TLorentzVector jet1_v4;
  jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
  return (M_min_ < jet1_v4.M() && jet1_v4.M() < M_max_);
}
