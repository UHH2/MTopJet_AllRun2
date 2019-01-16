#include "UHH2/MTopJet/include/NonClosureUncertainty.h"

double NonClosureUncertainty::get_factor(double pt){
  // factor is zero for pt > 400 anyways
  if(pt > 400) pt = 400;
  if(pt < 0)   pt = 0;
  // the entry in the file gives the value f = (rec - gen) / gen
  // in order to bring f to a perfect 0, one has to multiply rec with 1 / (f+1)
  // this is exactly the factor that is applyed for the variation.
  double factor = 1;
  double entry = NonClosureFunction->Eval(pt);
  factor = 1.0 / (entry+1);
  return factor;
}


void NonClosureUncertainty::get_function(){
  TString filename = "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/NonClosure/FlavorJEC_nonClosure.root";
  TFile *file = new TFile(filename);
  TString histname = "nonclosure";
  NonClosureFunction = (TGraph*)file->Get(histname);
}

Particle NonClosureUncertainty::GetLepton(uhh2::Event & event){

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


NonClosureUncertainty::NonClosureUncertainty(uhh2::Context & ctx):
h_oldjets(ctx.get_handle<std::vector<TopJet>>("xconeCHS")){
  get_function(); // read file here, not in every event
}

bool NonClosureUncertainty::process(uhh2::Event & event){
  std::vector<TopJet> oldjets = event.get(h_oldjets);
  std::vector<TopJet> newjets;
  for(unsigned int i=0; i < oldjets.size(); i++) newjets.push_back(oldjets.at(i));

  // get had jet, to only correct had subjets (in l+jets case)
  int had_nr = 0;
  int lep_nr = 1;
  Particle lepton = GetLepton(event);
  if(deltaR(lepton, oldjets.at(0)) < deltaR(lepton, oldjets.at(1))){
    had_nr = 1;
    lep_nr = 0;
  }

  // leave lep subjets unchanged
  newjets.at(lep_nr).set_subjets(oldjets.at(lep_nr).subjets());

  // now correct had subjets
  std::vector<Jet> oldsubjets = oldjets.at(had_nr).subjets();
  std::vector<Jet> newsubjets;
  LorentzVector newjet_v4;
  Jet newjet;
  float factor;
  for(unsigned int i=0; i < oldsubjets.size(); i++){
    factor = get_factor(oldsubjets.at(i).pt());
    newjet_v4 = oldsubjets.at(i).v4() * factor;
    newjet = oldsubjets.at(i); // this acutally is not required
    newjet.set_v4(newjet_v4);  // Because of this, all JEC factors are set to the defaul value of 1
    newsubjets.push_back(newjet);
  }
  newjets.at(had_nr).set_subjets(newsubjets);

  event.set(h_oldjets, newjets);

  return true;
}
