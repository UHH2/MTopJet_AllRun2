#include "UHH2/MTopJet/include/MTopJetGenSelections.h"
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


uhh2::MassCutGen1::MassCutGen1(uhh2::Context& ctx, const std::string & name, float M_min, float M_max):
  h_jets(ctx.get_handle<std::vector<Jet>>(name)),
  M_min_(M_min),
  M_max_(M_max) {}

bool MassCutGen1::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  TLorentzVector jet1_v4;
  jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy());
  return (M_min_ < jet1_v4.M() && jet1_v4.M() < M_max_);
}
