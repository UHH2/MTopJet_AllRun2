#include <UHH2/MTopJet/include/RecoSelections_topjet.h>


uhh2::LeadingRecoJetPT_topjet::LeadingRecoJetPT_topjet(uhh2::Context& ctx, float ptcut):
  h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets")),
  ptcut_(ptcut) {}

bool uhh2::LeadingRecoJetPT_topjet::passes(const uhh2::Event& event){
  bool pass_jetpt = false;
  std::vector<TopJet> jets = event.get(h_topjets);
  Jet jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float pt = jet1.pt();
    if(pt > ptcut_) pass_jetpt = true;
  }
  return pass_jetpt;
}
////////////////////////////////////////////////////////

uhh2::MassCutReco_topjet::MassCutReco_topjet(uhh2::Context& ctx):
  h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets")) {}

bool uhh2::MassCutReco_topjet::passes(const uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_topjets);
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }
  else return false;

  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  bool pass_masscut = false;

  if(jets.size()>1){
    jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(jets.at(1).pt(),jets.at(1).eta(),jets.at(1).phi(),jets.at(1).energy()); //v4 of first jet
    lepton1_v4.SetPtEtaPhiE(lepton.pt(),lepton.eta(),lepton.phi(),lepton.energy()); //v4 of lepton
    jet2_lep_v4 = jet2_v4 + lepton1_v4; // v4 of (lepton+2nd jet)

    if(jet1_v4.M() > jet2_lep_v4.M()) pass_masscut = true;
  }
  return pass_masscut;
}
////////////////////////////////////////////////////////

uhh2::DeltaRCutReco_topjet::DeltaRCutReco_topjet(uhh2::Context& ctx, float jetradius):
  h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets")),
  jetradius_(jetradius) {}

bool uhh2::DeltaRCutReco_topjet::passes(const uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_topjets);
  Jet jet2;
  bool pass_deltaR = false;
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }
  else return false;

  if(jets.size()>1){
    jet2 = jets.at(1);
    if(deltaR(lepton, jet2) < jetradius_) pass_deltaR = true;
  }
  return pass_deltaR;
}
////////////////////////////////////////////////////////

uhh2::NRecoJets_topjet::NRecoJets_topjet(uhh2::Context& ctx, float min_pt, float min, float max):
  h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets")),
  min_pt_(min_pt),
  min_(min),
  max_(max) {}

bool uhh2::NRecoJets_topjet::passes(const uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_topjets);
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

////////////////////////////////////////////////////////

RecoJetLeptonCleaner_topjet::RecoJetLeptonCleaner_topjet(uhh2::Context& ctx):
  h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets"))
   {}

bool RecoJetLeptonCleaner_topjet::process(uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_topjets);
  std::vector<TopJet> cleaned_jets;
  TopJet jet;

  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }


  // perform cleaning
  double jetradius_ = 0.8;
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
  sort_by_pt<TopJet>(cleaned_jets); // Sort Jets by pT
  event.set(h_topjets, cleaned_jets);

  return true;
}

////////////////////////////////////////////////////////
uhh2::deltaRCut_topjet::deltaRCut_topjet(float max_deltaR):
  max_deltaR_(max_deltaR) {}

bool uhh2::deltaRCut_topjet::passes(const uhh2::Event& event){

  assert(event.muons && event.electrons && event.topjets);
  if(event.topjets->size()<2) return false;
  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- deltaRCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  float dR=100;
  if(event.muons->size()!=0){
    dR = deltaR(event.muons->at(0), event.topjets->at(1));
  }
  else dR = deltaR(event.electrons->at(0), event.topjets->at(1));

  return (dR < max_deltaR_);
}

////////////////////////////////////////////////////////

uhh2::Jet2Cut_topjet::Jet2Cut_topjet(float min_pt):
  min_pt_(min_pt) {}

bool uhh2::Jet2Cut_topjet::passes(const uhh2::Event& event){
    const Particle* jet2 = &event.topjets->at(1);
    float jet_pt = jet2->pt();

    return jet_pt > min_pt_;

}

////////////////////////////////////////////////////////
