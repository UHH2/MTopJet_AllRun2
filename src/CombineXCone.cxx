#include "UHH2/MTopJet/include/CombineXCone.h"

// Get Primary Lepton
Particle CombineXCone::GetLepton(uhh2::Event & event){

  Particle lepton;
  bool muonislep = false;
  bool elecislep = false;

  if(event.muons->size() > 0 && event.electrons->size() > 0){
    if(event.muons->at(0).pt() > event.electrons->at(0).pt()) muonislep = true;
    else elecislep = true;
  }
  else if(event.muons->size() > 0) muonislep = true;
  else if(event.electrons->size() > 0) elecislep = true;

  if(muonislep) lepton = event.muons->at(0);
  if(elecislep) lepton = event.electrons->at(0);

  return lepton;
}

// Combine Subjets to jet (do not use subjet if pt < ptmin)
Jet CombineXCone::AddSubjets(vector<Jet> subjets, double ptmin){
  double px=0, py=0, pz=0, E=0;
  TLorentzVector jet_v4;
  Jet jet;
  for(unsigned int i=0; i < subjets.size(); ++i){
    if(subjets.at(i).pt() < ptmin) continue;
    px += subjets.at(i).v4().Px();
    py += subjets.at(i).v4().Py();
    pz += subjets.at(i).v4().Pz();
    E += subjets.at(i).v4().E();
  }
  jet_v4.SetPxPyPzE(px, py, pz, E);
  jet.set_pt(jet_v4.Pt());
  jet.set_eta(jet_v4.Eta());
  jet.set_phi(jet_v4.Phi());
  jet.set_energy(jet_v4.E());

  return jet;
}


// Get final Jets from 3+3 Method on Reco level
CombineXCone33::CombineXCone33(uhh2::Context & ctx):
  h_xcone33hadjets(ctx.declare_event_output<std::vector<Jet>>("XCone33_had_Combined")),
  h_xcone33lepjets(ctx.declare_event_output<std::vector<Jet>>("XCone33_lep_Combined")),
  h_fatjets(ctx.get_handle<std::vector<TopJet>>("XConeTopJets")) {}

bool CombineXCone33::process(uhh2::Event & event){
  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<TopJet> jets = event.get(h_fatjets);
  CombineXCone combine;
  Particle lepton = combine.GetLepton(event);

  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet (deltaR) -----------------------------------------
  //---------------------------------------------------------------------------------------
  TopJet fathadjet, fatlepjet;
  float dR1 = deltaR(lepton, jets.at(0));
  float dR2 = deltaR(lepton, jets.at(1));
  if(dR1 < dR2){
    fatlepjet = jets.at(0);
    fathadjet = jets.at(1);
  }
  else{
    fatlepjet = jets.at(1);
    fathadjet = jets.at(0);
  }

  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<Jet> subjets_lep = fatlepjet.subjets();
  std::vector<Jet> subjets_had = fathadjet.subjets();
  Jet lepjet = combine.AddSubjets(subjets_lep, 30);
  Jet hadjet = combine.AddSubjets(subjets_had, 30);
  vector<Jet> hadjets;
  vector<Jet> lepjets;
  hadjets.push_back(hadjet);
  lepjets.push_back(lepjet);

  //---------------------------------------------------------------------------------------
  //--------------------------------- Write Jets ------------------------------------------
  //---------------------------------------------------------------------------------------
  event.set(h_xcone33hadjets, hadjets);
  event.set(h_xcone33lepjets, lepjets);

  return true;
}
