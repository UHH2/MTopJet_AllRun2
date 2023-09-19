#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/core/include/Utils.h>

#include <UHH2/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <UHH2/JetMETObjects/interface/JetCorrectorParameters.h>

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

  string jec_tag_2016 = "Summer16_07Aug2017";
  string jec_ver_2016 = "11"; // 11

  string jec_tag_2017 = "Fall17_17Nov2017";
  string jec_ver_2017 = "32";

  string jec_tag_2018 = "Autumn18";
  string jec_ver_2018 = "19";

  string jec_jet_coll = "AK4PFchs";

  if (isMC)
  {
    jet_corrector_MC.reset(new YearSwitcher(ctx));
    jet_corrector_MC->setup2016v3(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2017v2(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2018(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll),jet_collection_rec));
  }

  else
  {
    jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      jec_switcher_16->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr),jet_collection_rec));
    }

    jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    for (const auto & runItr : runPeriods2017) {
      jec_switcher_17->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr),jet_collection_rec));
    }

    jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    for (const auto & runItr : runPeriods2018) {
      jec_switcher_18->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr),jet_collection_rec));
    }

    jet_corrector_data.reset(new YearSwitcher(ctx));
    jet_corrector_data->setup2016v3(jec_switcher_16);
    jet_corrector_data->setup2017v2(jec_switcher_17);
    jet_corrector_data->setup2018(jec_switcher_18);
  }

  // The new SetUp chooses the right Modul by itself. The old one has to be done by telling the program which Module to choose.

}

bool JetCorrections_xcone::process(uhh2::Event & event){
  std::vector<TopJet> jets = event.get(h_topjets);
  event.set(h_topjets, set_JEC_factor(jets)); // first set JEC_factor_raw to a non-0 value
  if(isMC) jet_corrector_MC->process(event);
  else jet_corrector_data->process(event);

  return true;
}


JER_Smearer_xcone::JER_Smearer_xcone(){}

void JER_Smearer_xcone::init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen, const std::string& fat_sub){
  h_rectopjets_ = ctx.get_handle<std::vector<TopJet> >   (jet_collection_rec);
  h_gentopjets_ = ctx.get_handle<std::vector<GenTopJet> >(jet_collection_gen);
  if(fat_sub == "sub") use_subjets = true;
  else                 use_subjets = false;


  JERSmearing::SFtype1 JER_sf;

  std::string jer_tag = "";
  std::string jetCollection = "AK4PFchs";
  const Year & year = extract_year(ctx);

  if (year == Year::is2016v3) {
    jer_tag = "Summer16_25nsV1";
  } else if (year == Year::is2017v2) {
    jer_tag = "Fall17_V3";
  } else if (year == Year::is2018) {
    jer_tag = "Autumn18_V7";
  } else {
    throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
  }

  jer_factors_wo = {};
  jer_factors_with = {};

  JER_Smearer.reset(new GenericJetResolutionSmearer(ctx,  "xconeCHS", jet_collection_gen, JERFiles::JERPathStringMC(jer_tag,jetCollection,"SF"), JERFiles::JERPathStringMC(jer_tag,jetCollection,"PtResolution")));

}

std::vector<double> JER_Smearer_xcone::JER_factors(int index){
  if(index!=0 && index!=1){
    throw runtime_error("Wrong index while obtain JER factors; " + to_string(index));
  }
  vector<double> f_jer = {};
  for(unsigned int i=0; i<jer_factors_wo[index].size(); i++){
    f_jer.push_back(jer_factors_wo[index][i]/jer_factors_with[index][i]);
  }
  return f_jer;
}


bool JER_Smearer_xcone::process(uhh2::Event & event){
  if(use_subjets){
    if(event.isRealData) return true;
    LorentzVector met;
    if(event.met) met = event.met->v4();
    std::vector<TopJet>*    rec_topjets;
    std::vector<GenTopJet>* gen_topjets;
    rec_topjets = &event.get(h_rectopjets_);
    gen_topjets = &event.get(h_gentopjets_);

    // match gentopjet to rectopjet
    double dR0 = deltaR(rec_topjets->at(0), gen_topjets->at(0));
    double dR1 = deltaR(rec_topjets->at(0), gen_topjets->at(1));
    unsigned int i0, i1;
    if(dR0 < dR1){
      i0 = 0;
      i1 = 1;
    }
    else {
      i0 = 1;
      i1 = 0;
    }

    // first smear subjets of topjet with index 0
    std::vector<Jet> rec_subjets0;
    for(unsigned int j=0; j<rec_topjets->at(0).subjets().size(); j++){
      rec_subjets0.push_back((rec_topjets->at(0).subjets().at(j)));
    }
    std::vector<GenJet> gen_subjets0;
    for(unsigned int j=0; j<gen_topjets->at(i0).subjets().size(); j++){
      gen_subjets0.push_back((gen_topjets->at(i0).subjets().at(j)) );
    }
    JER_Smearer->apply_JER_smearing(rec_subjets0, gen_subjets0, 0.4, event.rho);
    rec_topjets->at(0).set_subjets(rec_subjets0);

    // then smear subjets of topjet with index 1
    std::vector<Jet> rec_subjets1;
    for(unsigned int j=0; j<rec_topjets->at(1).subjets().size(); j++){
      rec_subjets1.push_back((rec_topjets->at(1).subjets().at(j)));
    }
    std::vector<GenJet> gen_subjets1;
    for(unsigned int j=0; j<gen_topjets->at(i1).subjets().size(); j++){
      gen_subjets1.push_back((gen_topjets->at(i1).subjets().at(j)) );
    }
    JER_Smearer->apply_JER_smearing(rec_subjets1, gen_subjets1, 0.4, event.rho);
    rec_topjets->at(1).set_subjets(rec_subjets1);
  }
  else

  JER_Smearer->process(event); // Zu dem else ???????????????????????????????????????????

  return true;
}
