#include "UHH2/MTopJet/include/CorrectionFactor.h"

double CorrectionFactor::get_factor(double pt, double eta){

  // get eta bin
  int etabin = 0;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    if(eta > eta_bins[i]) etabin = i;
  }

  // get factor from function
  if(pt > 400) pt = 400;
  if(pt < 0)   pt = 0;


  double factor;
  if(CorUp || CorDown) factor = CorrectionGraphs[etabin]->Eval(pt);
  else factor = CorrectionFunctions[etabin]->Eval(pt);

  // cout << "Eta Bin = " << etabin << endl;
  // cout << "PT      = " << pt << endl;
  // cout << "Factor  = " << factor << endl;

  return factor;
}

void CorrectionFactor::get_function(){

  TString dir = "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  // filename = dir + "Correction_noCHS.root";
  filename = dir + "Correction.root";
  TFile *file = new TFile(filename);

  TString histname;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    stringstream BinString;
    BinString << i;
    histname = BinString.str();
    if(CorUp){
      histname += "/Up";
      CorrectionGraphs.push_back((TGraph*)file->Get(histname));
    }
    else if(CorDown){
      histname += "/Down";
      CorrectionGraphs.push_back((TGraph*)file->Get(histname));
    }
    else{
      histname += "/Central";
      CorrectionFunctions.push_back((TF1*)file->Get(histname));
    }
  }
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


CorrectionFactor::CorrectionFactor(uhh2::Context & ctx, const std::string & name, std::string corvar):
  h_oldjets(ctx.get_handle<std::vector<TopJet>>("XConeTopJets")),
  h_newjets(ctx.declare_event_output<std::vector<TopJet>>(name))
{
  CorUp = false;
  CorDown=false;

  if(corvar == "up")   CorUp = true;
  else if(corvar == "down") CorDown = true;

  get_function(); // read file here, not in every event
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
