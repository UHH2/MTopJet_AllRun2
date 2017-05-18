#include "UHH2/MTopJet/include/JetCorrections_xcone.h"

// JEC_factor_raw has to be set to a non 0 value for JEC
std::vector<TopJet> set_JEC_factor(std::vector<TopJet> jets){
  Jet new_subjet;
  vector<Jet> new_subjets;
  TopJet new_fatjet;
  vector<TopJet> new_fatjets;
  for(unsigned int i=0; i<jets.size(); i++){
    new_fatjet = jets.at(i);
    new_subjets.clear();
    for(unsigned int j=0; j<new_fatjet.subjets().size(); j++){
      new_subjet = jets.at(i).subjets().at(j);
      new_subjet.set_JEC_factor_raw(1.);
      new_subjets.push_back(new_subjet);
    }
    new_fatjet.set_subjets(new_subjets);
    new_fatjets.push_back(new_fatjet);
  }
  return new_fatjets;
}

JetCorrections_xcone::JetCorrections_xcone(){}

void JetCorrections_xcone::init(uhh2::Context & ctx, const std::string& jet_collection){
  h_topjets = ctx.get_handle<std::vector<TopJet>>("XConeTopJets");
  isMC = (ctx.get("dataset_type") == "MC");
  jet_corrector_MC.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, jet_collection));
  jet_corrector_BCD.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA, jet_collection));
  jet_corrector_EFearly.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA, jet_collection));
  jet_corrector_FlateG.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA, jet_collection));
  jet_corrector_H.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA, jet_collection));
}

bool JetCorrections_xcone::process(uhh2::Event & event){
  std::vector<TopJet> jets = event.get(h_topjets);
  event.set(h_topjets, set_JEC_factor(jets)); // first set JEC_factor_raw to a non-0 value
  if(isMC)jet_corrector_MC->process(event);
  else{
    if(event.run <= runnr_BCD)         jet_corrector_BCD->process(event);
    else if(event.run < runnr_EFearly) jet_corrector_EFearly->process(event);
    else if(event.run <= runnr_FlateG) jet_corrector_FlateG->process(event);
    else if(event.run > runnr_FlateG)  jet_corrector_H->process(event);
    else throw runtime_error("Jet Correction: run number not covered by if-statements in process-routine.");
  }
}

// JER SMEARING for XCone Jets

void JER_Smearer::apply_smearing(std::vector<Jet>& rec_jets, const std::vector<Particle>& gen_jets, LorentzVector& met){

  for(unsigned int i=0; i<rec_jets.size(); ++i){

    auto& jet = rec_jets.at(i);

    // find next genjet:
    auto closest_genjet = closestParticle(jet, gen_jets);
    // ignore unmatched jets (=no genjets at all or large DeltaR), or jets with very low genjet pt:
    if(closest_genjet == nullptr || uhh2::deltaR(*closest_genjet, jet) > 0.3) continue;
    const float genpt = closest_genjet->pt();
    if(genpt < 15.) continue;

    LorentzVector jet_v4 = jet.v4();

    const float recopt = jet_v4.pt();
    const float abseta = fabs(jet_v4.eta());

    int ieta(-1);
    for(unsigned int idx=0; idx<JER_SFs_.size(); ++idx){

      const float min_eta = idx ? JER_SFs_.at(idx-1).at(0) : 0.;
      const float max_eta =       JER_SFs_.at(idx)  .at(0);

      if(min_eta <= abseta && abseta < max_eta){ ieta = idx; break;}
    }
    if(ieta < 0) {
      std::cout << "WARNING: JetResolutionSmearer: index for JER-smearing SF not found for jet with |eta| = " << abseta << std::endl;
      std::cout << "         no JER smearing is applied." << std::endl;
      continue;
    }

    float c;
    if     (direction ==  0) c = JER_SFs_.at(ieta).at(1);
    else if(direction ==  1) c = JER_SFs_.at(ieta).at(2);
    else if(direction == -1) c = JER_SFs_.at(ieta).at(3);
    else throw std::runtime_error("JERSmearer::process -- invalid value for JER 'direction' (must be 0, +1 or -1): "+std::to_string(direction));

    const float new_pt = std::max(0.0f, genpt + c * (recopt - genpt));
    jet_v4 *= new_pt / recopt;

    float factor_raw = jet.JEC_factor_raw();
    factor_raw *= recopt/new_pt;

    jet.set_JEC_factor_raw(factor_raw);
    jet.set_v4(jet_v4);
  }

  return;
}

JER_Smearer::JER_Smearer(uhh2::Context& ctx, const std::string& recjet_label, const std::string& genjet_label, const bool allow_met_smearing, const JERSmearing::SFtype1& JER_sf){

  if(ctx.get("meta_jer_applied__"+recjet_label, "") != "true") ctx.set_metadata("jer_applied__"+recjet_label, "true");
  else throw std::runtime_error("JERSmearer::JERSmearer -- JER smearing already applied to this RECO-jets collection: "+recjet_label);

  const std::string& dir = ctx.get("jersmear_direction", "nominal");
  if     (dir == "nominal") direction =  0;
  else if(dir == "up")      direction =  1;
  else if(dir == "down")    direction = -1;
  else throw std::runtime_error("JERSmearer::JERSmearer -- invalid value jersmear_direction='"+dir+"' (valid: 'nominal', 'up', 'down')");

  h_rectopjets_ = ctx.get_handle<std::vector<TopJet> >   (recjet_label);
  h_gentopjets_ = ctx.get_handle<std::vector<GenTopJet> >(genjet_label);

  JER_SFs_ = JER_sf;
}

bool JER_Smearer::process(uhh2::Event& evt){

  if(evt.isRealData) return true;

  LorentzVector met;
  if(evt.met) met = evt.met->v4();
  const std::vector<Particle>*  gen_jets(0);
  std::vector<GenTopJet> gen_topjets = evt.get(h_gentopjets_);
  std::vector<Jet> rec_jets(0);

  for(unsigned int i=0; i<gen_topjets.size(); i++){
    gen_jets = &(evt.get(h_gentopjets_).at(i).subjets());
    rec_jets = evt.get(h_rectopjets_).at(i).subjets();
    apply_smearing(rec_jets   , *gen_jets   , met);
  }

  return true;
}
