#include "UHH2/MTopJet/include/JetCorrections_xcone.h"

//update SF from Summer16 07Aug measurment
extern const JERSmearing::SFtype1 SF_13TeV_updated = {
  // 0 = upper jet-eta limit
  // 1 = JER SF
  // 2 = JER SF + 1sigma
  // 3 = JER SF - 1sigma
  {{0.5, 1.1595, 1.095, 1.224}},
  {{0.8, 1.1948, 1.1296, 1.26}},
  {{1.1, 1.1464, 1.0832, 1.2096}},
  {{1.3, 1.1609, 1.0584, 1.2634}},
  {{1.7, 1.1278, 1.0292, 1.2264}},
  {{1.9, 1.1000, 0.9921, 1.2079}},
  {{2.1, 1.1426, 1.0212, 1.264}},
  {{2.3, 1.1512, 1.0372, 1.2652}},
  {{2.5, 1.2963, 1.0592, 1.5334}},
  {{2.8, 1.3418, 1.1327, 1.5509}},
  {{3.0, 1.7788, 1.578, 1.9796}},
  {{3.2, 1.1869, 1.0626, 1.3112}},
  {{5.0, 1.1922, 1.0434, 1.341}},
};

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

void JetCorrections_xcone::init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen){
  h_topjets = ctx.get_handle<std::vector<TopJet>>(jet_collection_rec);
  h_GENtopjets = ctx.get_handle<std::vector<GenTopJet>>(jet_collection_gen);

  // setup JEC for XCone
  isMC = (ctx.get("dataset_type") == "MC");
  jet_corrector_MC.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PF_MC, jet_collection_rec));
  jet_corrector_BCD.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PF_DATA, jet_collection_rec));
  jet_corrector_EFearly.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PF_DATA, jet_collection_rec));
  jet_corrector_FlateG.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PF_DATA, jet_collection_rec));
  jet_corrector_H.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PF_DATA, jet_collection_rec));
  //   jet_corrector_MC.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, jet_collection));
  //   jet_corrector_BCD.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA, jet_collection));
  //   jet_corrector_EFearly.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA, jet_collection));
  //   jet_corrector_FlateG.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA, jet_collection));
  //   jet_corrector_H.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA, jet_collection));

  // setup JER for XCone
  TString ResolutionFileName = "Summer16_25nsV1_MC_SF_AK4PF.txt";
  JER_Smearer.reset(new GenericJetResolutionSmearer(ctx, jet_collection_rec, jet_collection_gen, true, SF_13TeV_updated, ResolutionFileName));
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
  JER_Smearer->process(event);
  return true;
}
