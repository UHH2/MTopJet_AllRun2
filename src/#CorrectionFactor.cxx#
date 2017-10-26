#include "UHH2/MTopJet/include/CorrectionFactor.h"

// double CorrectionFactor::get_factor_table(double pt, double eta){

//   // get pt bin
//   int ptbin = 0;
//   for(unsigned int i = 0; i < 6; i++){
//     if(pt > pt_bins[i]) ptbin = i;
//   }

//   // get eta bin
//   int etabin = 0;
//   for(unsigned int i = 0; i < 6; i++){
//     if(eta > eta_bins[i]) etabin = i;
//   }

//   // get factor from array
//   double factor = arr[ptbin][etabin];

//   return factor;
// }

double CorrectionFactor::get_factor(double pt, double eta){

  // get eta bin
  int etabin = 0;
  for(unsigned int i = 0; i < 6; i++){
    if(eta > eta_bins[i]) etabin = i;
  }

  // get factor from function
  if(pt > 425) pt = 425;
  double factor = par[etabin][0] + par[etabin][1] * pt + par[etabin][2] * pt * pt;
  return factor;
}

// void CorrectionFactor::get_values_table(){

//   // get factors in pt and eta bins from file
//   std::ifstream file("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/CorrectionFile/CorrectionFactors_table.txt");
//   int p, e;
//   double f;
//   while (file >> p >> e >> f) {
//     arr[p][e] = f;
//   }
//   file.close();

//   return;
// }

void CorrectionFactor::get_values(){
  // get factors in pt and eta bins from file
  // std::ifstream file("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/CorrectionFile/CorrectionFactors.txt");
  // int eta, p0, p1, p2, p3;
  // while (file >> eta >> p0 >> p1 >> p2) {
  //   cout<< "eta = " << eta << "  p0 = " << p0 << "  p1 = " << p1 << "  p2 = " << p2 << endl; 
  //   par[eta][0] = p0;
  //   par[eta][1] = p1;
  //   par[eta][2] = p2;
  // }
  // file.close();

  par[0][0] = 0.9004245128;
  par[0][1] = 0.0003551003;
  par[0][2] = -0.0000004440;

  par[1][0] = 0.9143051034;
  par[1][1] = 0.0002610321;
  par[1][2] = -0.0000003132;

  par[2][0] = 0.9112724866;
  par[2][1] = 0.0003547596;
  par[2][2] = -0.0000005067;

  par[3][0] = 0.9342017103;
  par[3][1] = 0.0001811427;
  par[3][2] = -0.0000001937;

  par[4][0] = 0.9357155465;
  par[4][1] = 0.0001761840;
  par[4][2] = -0.0000001986;

  par[5][0] = 0.9346114769;
  par[5][1] = 0.0001799118;
  par[5][2] = -0.0000002205;

  par[6][0] = 0.9384138337;
  par[6][1] = 0.0001473672;
  par[6][2] = -0.0000001366;

  par[7][0] = 0.9325012830;
  par[7][1] = 0.0002224927;
  par[7][2] = -0.0000002989;

  par[8][0] = 0.9333083751;
  par[8][1] = 0.0002215052;
  par[8][2] = -0.0000003100;

  par[9][0] = 0.9258555506;
  par[9][1] = 0.0002659339;
  par[9][2] = -0.0000003721;

  par[10][0] = 0.9212095065;
  par[10][1] = 0.0002461865;
  par[10][2] = -0.0000003176;

  par[11][0] = 0.9021580514;
  par[11][1] = 0.0003397708;
  par[11][2] = -0.0000004066;

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
  bool use_table = false;
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
