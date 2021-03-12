#include <UHH2/MTopJet/include/CombineXCone.h>

// TODO: remove ptcut from combine

bool CombineXCone::FindLepton(uhh2::Event & event){
  bool lepton_in_event = false;
  if(event.muons->size() > 0 || event.electrons->size() > 0) lepton_in_event = true;
  return lepton_in_event;
}

bool CombineXCone::FindLepton_gen(uhh2::Event & event){
  bool lepton_in_event = false;
  int lep_counter = 0;
  std::vector<GenParticle>* genparts = event.genparticles;
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle p = genparts->at(i);
    if(abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13 ) ++lep_counter;
  }

  if(lep_counter > 0) lepton_in_event = true;
  return lepton_in_event;
}


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

GenParticle CombineXCone::GetLepton_gen(uhh2::Event & event){

  GenParticle lepton;
  std::vector<GenParticle>* genparts = event.genparticles;

  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle p = genparts->at(i);
    if(abs(p.pdgId()) == 13) lepton = p;
    else if(abs(p.pdgId()) == 11) lepton = p;
  }
  return lepton;
}

// Combine Subjets to jet (do not use subjet if pt < ptmin)
TopJet CombineXCone::CreateTopJetFromSubjets(vector<Jet> subjets, double ptmin, double etamax){
  double px=0, py=0, pz=0, E=0;
  TLorentzVector jet_v4;
  TopJet jet;
  // first store subjets
  jet.set_subjets(subjets);
  // only combine if all subjets have pt > ptmin
  bool good_jet = true;
  for(unsigned int i=0; i < subjets.size(); ++i){
    if(subjets.at(i).pt() < ptmin) good_jet = false;
    if(subjets.at(i).eta() > etamax) good_jet = false;
    if(subjets.at(i).eta() < -etamax) good_jet = false;
  }
  if(good_jet && subjets.size() >= 3){
    for(unsigned int i=0; i < subjets.size(); ++i){
      if(subjets.at(i).pt() < ptmin) continue;
      px += subjets.at(i).v4().Px();
      py += subjets.at(i).v4().Py();
      pz += subjets.at(i).v4().Pz();
      E += subjets.at(i).v4().E();
    }
  }
  jet_v4.SetPxPyPzE(px, py, pz, E);
  jet.set_pt(jet_v4.Pt());
  jet.set_eta(jet_v4.Eta());
  jet.set_phi(jet_v4.Phi());
  jet.set_energy(jet_v4.E());

  return jet;
}

GenTopJet CombineXCone::CreateTopJetFromSubjets_gen(vector<GenJet> subjets, double ptmin, double etamax){
  double px=0, py=0, pz=0, E=0;
  TLorentzVector jet_v4;
  GenTopJet jet;
  // first store subjets
  for(auto sj: subjets) jet.add_subjet(sj);
  // only combine if all subjets have pt > ptmin
  bool good_jet = true;
  for(unsigned int i=0; i < subjets.size(); ++i){
    if(subjets.at(i).pt() < ptmin) good_jet = false;
    if(subjets.at(i).eta() > etamax) good_jet = false;
    if(subjets.at(i).eta() < -etamax) good_jet = false;
  }
  if(good_jet){
    for(unsigned int i=0; i < subjets.size(); ++i){
      px += subjets.at(i).v4().Px();
      py += subjets.at(i).v4().Py();
      pz += subjets.at(i).v4().Pz();
      E += subjets.at(i).v4().E();
    }
  }
  jet_v4.SetPxPyPzE(px, py, pz, E);
  jet.set_pt(jet_v4.Pt());
  jet.set_eta(jet_v4.Eta());
  jet.set_phi(jet_v4.Phi());
  jet.set_energy(jet_v4.E());

  return jet;
}

TopJet CombineXCone::FindHadjet(uhh2::Event &event, const vector<TopJet>& fatjets){
  CombineXCone* combine = new CombineXCone();
  bool lepton_in_event = combine->FindLepton(event);
  TopJet fathadjet, fatlepjet;
  if(lepton_in_event){
    Particle lepton = combine->GetLepton(event);
    float dR1 = deltaR(lepton, fatjets.at(0));
    float dR2 = deltaR(lepton, fatjets.at(1));
    if(dR1 < dR2){
      fatlepjet = fatjets.at(0);
      fathadjet = fatjets.at(1);
    }
    else{
      fatlepjet = fatjets.at(1);
      fathadjet = fatjets.at(0);
    }
  }
  else{
    fatlepjet = fatjets.at(1);
    fathadjet = fatjets.at(0);
  }

  return fathadjet;
}

/*
██████  ███████  ██████
██   ██ ██      ██
██████  █████   ██
██   ██ ██      ██
██   ██ ███████  ██████
*/


// Get final Jets from 3+3 Method on Reco level
CombineXCone33::CombineXCone33(uhh2::Context & ctx, const std::string & name_had, const std::string & name_lep, const std::string & name_fat):
h_xcone33hadjets(ctx.declare_event_output<std::vector<TopJet>>(name_had)),
h_xcone33lepjets(ctx.declare_event_output<std::vector<TopJet>>(name_lep)),
h_fatjets(ctx.get_handle<std::vector<TopJet>>(name_fat)){}

bool CombineXCone33::process(uhh2::Event & event){
  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<TopJet> jets = event.get(h_fatjets);
  if(jets.size() < 2) return false;
  CombineXCone* combine = new CombineXCone();
  bool lepton_in_event = combine->FindLepton(event);

  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet (deltaR) -----------------------------------------
  //---------------------------------------------------------------------------------------
  TopJet fathadjet, fatlepjet;
  if(lepton_in_event){
    Particle lepton = combine->GetLepton(event);
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
  }
  else{
    fatlepjet = jets.at(1);
    fathadjet = jets.at(0);
  }


  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  TopJet lepjet, hadjet;
  std::vector<Jet> subjets_lep = fatlepjet.subjets();
  std::vector<Jet> subjets_had = fathadjet.subjets();
  lepjet = combine->CreateTopJetFromSubjets(subjets_lep, 0, 2.5);
  hadjet = combine->CreateTopJetFromSubjets(subjets_had, 0, 2.5);

  vector<TopJet> hadjets;
  vector<TopJet> lepjets;
  hadjets.push_back(hadjet);
  lepjets.push_back(lepjet);
  //---------------------------------------------------------------------------------------
  //--------------------------------- Write Jets ------------------------------------------
  //---------------------------------------------------------------------------------------
  event.set(h_xcone33hadjets, hadjets);
  event.set(h_xcone33lepjets, lepjets);

  return true;
}

//---------------------------------------------------------------------------------------
//--------------------------------- W Subjets -------------------------------------------
//---------------------------------------------------------------------------------------

vector<int> CombineXCone33::GetWSubjetsIndices(uhh2::Event & event){
  bool debug = false;
  // Select Subjets for Wjet ---------------------------------------------------------------------------------------
  // Select only ak4 jets with are close to the hadronic xcone jet
  std::vector<TopJet> fatjets = event.get(h_fatjets);

  //------------- define had and lep jet (deltaR) -----------------------------------------
  CombineXCone* combine = new CombineXCone();
  bool lepton_in_event = combine->FindLepton(event);
  TopJet fathadjet, fatlepjet;
  if(lepton_in_event){
    Particle lepton = combine->GetLepton(event);
    float dR1 = deltaR(lepton, fatjets.at(0));
    float dR2 = deltaR(lepton, fatjets.at(1));
    if(dR1 < dR2){
      fatlepjet = fatjets.at(0);
      fathadjet = fatjets.at(1);
    }
    else{
      fatlepjet = fatjets.at(1);
      fathadjet = fatjets.at(0);
    }
  }
  else{
    fatlepjet = fatjets.at(1);
    fathadjet = fatjets.at(0);
  }

  //------------- Collect AK4 jets --------------------------------------------------------
  std::vector<Jet> ak4_matched_xcone_had;
  for(const auto& j : *event.jets){
    if(deltaR(j, fathadjet)<1.2) ak4_matched_xcone_had.push_back(j);
  }

  if(debug) std::cout << " - search for AK4 jet with higest btag" << '\n';
  std::vector<double> btag_v; // Using a vector to identifie events with ambigous btag results (two high btags i.e. through a c-quark)
  double btag_high = 0;
  Jet ak4_highest_btag;

  if(ak4_matched_xcone_had.size()==0) return {0,1};
  else{
    for(unsigned int i=0; i<ak4_matched_xcone_had.size(); i++){
      btag_v.push_back(ak4_matched_xcone_had[i].btag_DeepJet());
      if(i==0){
        btag_high = btag_v[i];
        ak4_highest_btag = ak4_matched_xcone_had[i];
      }
      if(i>0 && btag_high<btag_v[i]){
        btag_high = btag_v[i];
        ak4_highest_btag = ak4_matched_xcone_had[i];
      }
    }
  }

  if(debug) std::cout << " - Identify Wjet: Distance AK4 to subjets" << '\n';
  std::vector<Jet> xcone_had_subjets = fathadjet.subjets();
  // Matching AK4 jet with highest btag to the three XCone subjets
  double dR_subjet1_ak4 = deltaR(xcone_had_subjets[0], ak4_highest_btag);
  double dR_subjet2_ak4 = deltaR(xcone_had_subjets[1], ak4_highest_btag);
  double dR_subjet3_ak4 = deltaR(xcone_had_subjets[2], ak4_highest_btag);

  if(debug) std::cout << " - Identify Wjet: Compare Distance" << '\n';
  int index_Wjet1 = -1; int index_Wjet2 = -1;
  if(dR_subjet1_ak4<dR_subjet2_ak4 && dR_subjet1_ak4<dR_subjet3_ak4){ // bjet is 0
    index_Wjet1 = 1;
    index_Wjet2 = 2;
  }
  else if(dR_subjet2_ak4<dR_subjet1_ak4 && dR_subjet2_ak4<dR_subjet3_ak4){ // bjet is 1
    index_Wjet1 = 0;
    index_Wjet2 = 2;
  }
  else{ // bjet is 2
    index_Wjet1 = 0;
    index_Wjet2 = 1;
  }

  vector<int> WSubjetsIndex = {index_Wjet1, index_Wjet2};
  return WSubjetsIndex;
}

TopJet CombineXCone33::GetHadronicWJet(uhh2::Event & event, const vector<int>& Indices){
  bool debug = false;
  // Select Subjets for Wjet ---------------------------------------------------------------------------------------
  // Select only ak4 jets with are close to the hadronic xcone jet
  std::vector<TopJet> fatjets = event.get(h_fatjets);

  //------------- define had and lep jet (deltaR) -----------------------------------------
  if(debug) cout << "GetHadronicWJet: FindHadjet" << endl;
  CombineXCone* combine = new CombineXCone();
  TopJet fathadjet = combine->FindHadjet(event, fatjets);

  std::vector<Jet> xcone_had_subjets = fathadjet.subjets();
  Jet WSubjet1 = xcone_had_subjets[Indices[0]];
  Jet WSubjet2 = xcone_had_subjets[Indices[1]];
  vector<Jet> WSubjets = {WSubjet1, WSubjet2};

  TLorentzVector jet_v4;
  TopJet WJet;

  // first store subjets
  if(debug) cout << "GetHadronicWJet: Creat WJet from Subjets" << endl;
  WJet.set_subjets(WSubjets);

  double px=0, py=0, pz=0, E=0;
  for(unsigned int i=0; i < WSubjets.size(); ++i){
    px += WSubjets.at(i).v4().Px();
    py += WSubjets.at(i).v4().Py();
    pz += WSubjets.at(i).v4().Pz();
    E += WSubjets.at(i).v4().E();
  }
  jet_v4.SetPxPyPzE(px, py, pz, E);
  WJet.set_pt(jet_v4.Pt());
  WJet.set_eta(jet_v4.Eta());
  WJet.set_phi(jet_v4.Phi());
  WJet.set_energy(jet_v4.E());

  return WJet;
}

/*
.██████  ███████ ███    ██
██       ██      ████   ██
██   ███ █████   ██ ██  ██
██    ██ ██      ██  ██ ██
.██████  ███████ ██   ████
*/


// Get final Jets from 3+3 Method on Gen level
CombineXCone33_gen::CombineXCone33_gen(uhh2::Context & ctx, bool isTTbar):
h_GENxcone33hadjets(ctx.declare_event_output<std::vector<GenTopJet>>("GEN_XCone33_had_Combined")),
h_GENxcone33lepjets(ctx.declare_event_output<std::vector<GenTopJet>>("GEN_XCone33_lep_Combined")),
h_GENfatjets(ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets")),
h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
isTTbar_ (isTTbar) {}

bool CombineXCone33_gen::process(uhh2::Event & event){
  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  GenParticle lepton;
  std::vector<GenTopJet> jets = event.get(h_GENfatjets);
  if(jets.size() < 2) return false;
  CombineXCone* combine = new CombineXCone();
  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet (deltaR) -----------------------------------------
  //---------------------------------------------------------------------------------------
  GenTopJet fathadjet, fatlepjet;
  if(isTTbar_) {
    TTbarGen ttbargen = event.get(h_ttbargen);
    if(ttbargen.IsSemiLeptonicDecay()){
      lepton = ttbargen.ChargedLepton();
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
    }
    else{
      fatlepjet = jets.at(1);
      fathadjet = jets.at(0);
    }
  }
  else{
    fatlepjet = jets.at(1);
    fathadjet = jets.at(0);
  }
  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<GenJet> subjets_lep = fatlepjet.subjets();
  std::vector<GenJet> subjets_had = fathadjet.subjets();
  GenTopJet lepjet = combine->CreateTopJetFromSubjets_gen(subjets_lep, 0, 2.5);
  GenTopJet hadjet = combine->CreateTopJetFromSubjets_gen(subjets_had, 0, 2.5);
  vector<GenTopJet> hadjets;
  vector<GenTopJet> lepjets;
  hadjets.push_back(hadjet);
  lepjets.push_back(lepjet);

  //---------------------------------------------------------------------------------------
  //--------------------------------- Write Jets ------------------------------------------
  //---------------------------------------------------------------------------------------
  event.set(h_GENxcone33hadjets, hadjets);
  event.set(h_GENxcone33lepjets, lepjets);

  // delete combine;
  return true;
}

/*
.█████  ██      ██          ██   ██  █████  ██████
██   ██ ██      ██          ██   ██ ██   ██ ██   ██
███████ ██      ██          ███████ ███████ ██   ██
██   ██ ██      ██          ██   ██ ██   ██ ██   ██
██   ██ ███████ ███████     ██   ██ ██   ██ ██████
*/


// Get final Jets from 3+3 Method on Reco level for ALL HAD Selection
CombineXConeAllHad::CombineXConeAllHad(uhh2::Context & ctx, const std::string & name_had, const std::string & name_lep, const std::string & name_fat):
h_xcone33hadjets(ctx.declare_event_output<std::vector<TopJet>>(name_had)),
h_xcone33lepjets(ctx.declare_event_output<std::vector<TopJet>>(name_lep)),
h_fatjets(ctx.get_handle<std::vector<TopJet>>(name_fat)){}

bool CombineXConeAllHad::process(uhh2::Event & event){
  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<TopJet> jets = event.get(h_fatjets);
  if(jets.size() < 2) return false;
  CombineXCone* combine = new CombineXCone();

  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet (deltaR) -----------------------------------------
  //---------------------------------------------------------------------------------------
  TopJet fathadjet = jets.at(0);
  TopJet fatlepjet = jets.at(1);
  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  TopJet lepjet, hadjet;
  std::vector<Jet> subjets_lep = fatlepjet.subjets();
  std::vector<Jet> subjets_had = fathadjet.subjets();
  lepjet = combine->CreateTopJetFromSubjets(subjets_lep, 0, 2.5);
  hadjet = combine->CreateTopJetFromSubjets(subjets_had, 0, 2.5);

  vector<TopJet> hadjets;
  vector<TopJet> lepjets;
  hadjets.push_back(hadjet);
  lepjets.push_back(lepjet);
  //---------------------------------------------------------------------------------------
  //--------------------------------- Write Jets ------------------------------------------
  //---------------------------------------------------------------------------------------
  event.set(h_xcone33hadjets, hadjets);
  event.set(h_xcone33lepjets, lepjets);

  return true;
}

// Get final Jets from 3+3 Method on Gen level
CombineXConeAllHad_gen::CombineXConeAllHad_gen(uhh2::Context & ctx):
h_GENxcone33hadjets(ctx.declare_event_output<std::vector<GenTopJet>>("GEN_XCone33_had_Combined")),
h_GENxcone33lepjets(ctx.declare_event_output<std::vector<GenTopJet>>("GEN_XCone33_lep_Combined")),
h_GENfatjets(ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets")) {}

bool CombineXConeAllHad_gen::process(uhh2::Event & event){
  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<GenTopJet> jets = event.get(h_GENfatjets);
  if(jets.size() < 2) return false;
  CombineXCone* combine = new CombineXCone();

  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet  -------------------------------------------------
  //---------------------------------------------------------------------------------------
  GenTopJet fathadjet, fatlepjet;
  fatlepjet = jets.at(1);
  fathadjet = jets.at(0);

  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<GenJet> subjets_lep = fatlepjet.subjets();
  std::vector<GenJet> subjets_had = fathadjet.subjets();
  GenTopJet lepjet = combine->CreateTopJetFromSubjets_gen(subjets_lep, 0, 2.5);
  GenTopJet hadjet = combine->CreateTopJetFromSubjets_gen(subjets_had, 0, 2.5);
  vector<GenTopJet> hadjets;
  vector<GenTopJet> lepjets;
  hadjets.push_back(hadjet);
  lepjets.push_back(lepjet);

  //---------------------------------------------------------------------------------------
  //--------------------------------- Write Jets ------------------------------------------
  //---------------------------------------------------------------------------------------
  event.set(h_GENxcone33hadjets, hadjets);
  event.set(h_GENxcone33lepjets, lepjets);

  // delete combine;
  return true;
}



CopyJets::CopyJets(uhh2::Context & ctx, const std::string & old_name, const std::string & new_name):
h_new(ctx.declare_event_output<std::vector<TopJet>>(new_name)),
h_old(ctx.get_handle<std::vector<TopJet>>(old_name)) {}


bool CopyJets::process(uhh2::Event & event){
  std::vector<TopJet> jets = event.get(h_old);
  event.set(h_new, jets);

  return true;
}
