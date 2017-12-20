#include "UHH2/MTopJet/include/CorrectionFactor.h"

double CorrectionFactor::get_factor(double pt, double eta){

  // get eta bin
  int etabin = 0;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(unsigned int i = 0; i < N_eta_bins; i++){
    if(eta > eta_bins[i]) etabin = i;
  }

  // get factor from function
  if(pt > 400) pt = 400;
  double factor = par[etabin][0] + par[etabin][1] * pt + par[etabin][2] * pt * pt;

  // cout << "Eta Bin     = " << etabin << endl;
  // cout << "Parameter 0 = " << par[etabin][0] << endl;
  // cout << "Parameter 1 = " << par[etabin][1] << endl;
  // cout << "Parameter 2 = " << par[etabin][2] << endl;
  // cout << "Parameter 3 = " << par[etabin][3] << endl;
  // cout << "Factor      = " << factor << endl;

  return factor;
}

void CorrectionFactor::get_values(){
  // get factors in pt and eta bins from file
  std::ifstream file("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/CorrectionFile/CorrectionFactors_new.txt");
  int eta;
  double p0, p1, p2, p3;
  while(file >> eta >> p0 >> p1 >> p2){
    cout<< "eta = " << eta << "  p0 = " << p0 << "  p1 = " << p1 << "  p2 = " << p2 <<endl; 
    par[eta][0] = p0;
    par[eta][1] = p1;
    par[eta][2] = p2;
  }
  file.close();

  return;
}

Particle GetLepton(uhh2::Event & event){

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


CorrectionFactor::CorrectionFactor(uhh2::Context & ctx, const std::string & name):
  h_oldjets(ctx.get_handle<std::vector<TopJet>>("XConeTopJets")),
  h_newjets(ctx.declare_event_output<std::vector<TopJet>>(name))
{
  get_values(); // read file here, not in every event
}

bool CorrectionFactor::process(uhh2::Event & event){
  std::vector<TopJet> oldjets = event.get(h_oldjets);
  std::vector<TopJet> newjets;
  for(unsigned int i=0; i < oldjets.size(); i++) newjets.push_back(oldjets.at(i));

  // get had jet, to only correct had subjets
  Particle lepton = GetLepton(event);
  int had_nr = 0;
  int lep_nr = 1;
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
    factor = get_factor(oldsubjets.at(i).pt(), oldsubjets.at(i).eta());
    newjet_v4 = oldsubjets.at(i).v4() * factor;
    newjet.set_v4(newjet_v4);
    newsubjets.push_back(newjet);
    // cout << "factor: " << factor << endl;
  }

  newjets.at(had_nr).set_subjets(newsubjets);
  event.set(h_newjets, newjets);

  return true;
}
