#include "UHH2/MTopJet/include/JetCorrections_xcone.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"

#include <string>

using namespace std;
using namespace uhh2;


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

void JetCorrections_xcone::init(uhh2::Context & ctx, const std::string& jet_collection_rec){
  h_topjets = ctx.get_handle<std::vector<TopJet>>(jet_collection_rec);

  // setup JEC for XCone
  isMC = (ctx.get("dataset_type") == "MC");
  jet_corrector_MC.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, jet_collection_rec));

  // jet_corrector_MC_b.reset(new GenericSubJetCorrector_flavor(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorB, jet_collection_rec, "b"));
  // jet_corrector_MC_ud.reset(new GenericSubJetCorrector_flavor(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorUD, jet_collection_rec, "ud"));

  jet_corrector_BCD.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA, jet_collection_rec));
  jet_corrector_EFearly.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA, jet_collection_rec));
  jet_corrector_FlateG.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA, jet_collection_rec));
  jet_corrector_H.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA, jet_collection_rec));
}

bool JetCorrections_xcone::process(uhh2::Event & event){
  std::vector<TopJet> jets = event.get(h_topjets);
  // event.set(h_topjets, set_JEC_factor(jets)); // first set JEC_factor_raw to a non-0 value
  if(isMC){
    //jet_corrector_MC_b->process(event);
    //jet_corrector_MC_ud->process(event);
    jet_corrector_MC->process(event);
  }
  else{
    if(event.run <= runnr_BCD)         jet_corrector_BCD->process(event);
    else if(event.run < runnr_EFearly) jet_corrector_EFearly->process(event);
    else if(event.run <= runnr_FlateG) jet_corrector_FlateG->process(event);
    else if(event.run > runnr_FlateG)  jet_corrector_H->process(event);
    else throw runtime_error("Jet Correction: run number not covered by if-statements in process-routine.");
  }

  return true;
}


JER_Smearer_xcone::JER_Smearer_xcone(){}

void JER_Smearer_xcone::init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen){
  JER_Smearer.reset(new GenericJetResolutionSmearer(ctx, jet_collection_rec, jet_collection_gen, false,  JERSmearing::SF_13TeV_2016_07Aug2017_v1, "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"));
}

bool JER_Smearer_xcone::process(uhh2::Event & event){
  JER_Smearer->process(event);
  return true;
}
