#include "UHH2/MTopJet/include/GenericSubJetCorrector_flavor.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"

#include <string>

using namespace std;
using namespace uhh2;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//            replace L2Relative with flavor
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorUD = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorUD = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorUD = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorUD = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorUD = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_ud_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorB = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorB = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorB = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorB = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorB = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_b_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V10_Flavor_Pythia8_MC/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

void correct_jet_flavor(FactorizedJetCorrector & corrector, Jet & jet, const Event & event, JetCorrectionUncertainty* jec_unc, int jec_unc_direction){
  auto factor_raw = jet.JEC_factor_raw();
  corrector.setJetPt(jet.pt() * factor_raw);
  corrector.setJetEta(jet.eta());
  corrector.setJetE(jet.energy() * factor_raw);
  corrector.setJetA(jet.jetArea());
  corrector.setRho(event.rho);
  auto correctionfactors = corrector.getSubCorrections();
  auto correctionfactor_L1  = correctionfactors.front();
  auto correctionfactor = correctionfactors.back();

  LorentzVector jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);

  if(jec_unc_direction!=0){
    if (jec_unc==NULL){
      std::cerr << "JEC variation should be applied, but JEC uncertainty object is NULL! Abort." << std::endl;
      exit(EXIT_FAILURE);
    }
    // ignore jets with very low pt or high eta, avoiding a crash from the JESUncertainty tool
    double pt = jet_v4_corrected.Pt();
    double eta = jet_v4_corrected.Eta();
    if (!(pt<5. || fabs(eta)>5.)) {

      jec_unc->setJetEta(eta);
      jec_unc->setJetPt(pt);

      double unc = 0.;
      if (jec_unc_direction == 1){
        unc = jec_unc->getUncertainty(1);
        correctionfactor *= (1 + fabs(unc));
      } else if (jec_unc_direction == -1){
        unc = jec_unc->getUncertainty(-1);
        correctionfactor *= (1 - fabs(unc));
      }
      jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);
    }
  }


  jet.set_v4(jet_v4_corrected);
  jet.set_JEC_factor_raw(1. / correctionfactor);
  jet.set_JEC_L1factor_raw(correctionfactor_L1);

}


namespace {

  // to share some code between JetCorrector and JetLeptonCleaner, provide some methods
  // dealing with jet energy corrections here:
  std::unique_ptr<FactorizedJetCorrector> build_corrector_flavor(const std::vector<std::string> & filenames){
    std::vector<JetCorrectorParameters> pars;
    for(const auto & filename : filenames){
      pars.emplace_back(locate_file(filename));
    }
    return uhh2::make_unique<FactorizedJetCorrector>(pars);
  }

  JetCorrectionUncertainty* corrector_uncertainty_flavor(uhh2::Context & ctx, const std::vector<std::string> & filenames, int &direction){

    auto dir = ctx.get("jecsmear_direction", "nominal");
    if(dir == "up"){
      direction = 1;
    }
    else if(dir == "down"){
      direction = -1;
    }
    else if(dir != "nominal"){
      // direction = 0 is default
      throw runtime_error("JetCorrector: invalid value jecsmear_direction='" + dir + "' (valid: 'nominal', 'up', 'down')");
    }

    //initialize JetCorrectionUncertainty if shift direction is not "nominal", else return NULL pointer
    if(direction!=0){
      //take name from the L1FastJet correction (0th element of filenames) and replace "L1FastJet" by "Uncertainty" to get the proper name of the uncertainty file
      TString unc_file = locate_file(filenames[0]);
      if (unc_file.Contains("L1FastJet")) {
        unc_file.ReplaceAll("L1FastJet","Uncertainty");
      }
      else if (unc_file.Contains("L2Relative")) {
        unc_file.ReplaceAll("L2Relative","Uncertainty");
      }
      else {
        throw runtime_error("WARNING No JEC Uncertainty File found!");
      }
      JetCorrectionUncertainty* jec_uncertainty = new JetCorrectionUncertainty(unc_file.Data());
      return jec_uncertainty;
    }
    return NULL;

  }

}


GenericSubJetCorrector_flavor::GenericSubJetCorrector_flavor(uhh2::Context & ctx, const std::vector<std::string> & filenames, const std::string & collectionname, std::string flavor_){
  corrector = build_corrector_flavor(filenames);
  direction = 0;
  flavor = flavor_;
  jec_uncertainty = corrector_uncertainty_flavor(ctx, filenames, direction) ;
  h_jets = ctx.get_handle<std::vector<TopJet> >(collectionname);
}

bool GenericSubJetCorrector_flavor::process(uhh2::Event & event){

  const auto topjets = &event.get(h_jets);
  assert(topjets);
  for(auto & topjet : *topjets){
    auto subjets = topjet.subjets();
    int index_bjet = 0;
    double ptmax = 0;
    for(unsigned int i=0; i<subjets.size(); i++){
      if(subjets[i].pt() > ptmax) ptmax = subjets[i].pt();
      index_bjet = i;
    }
    for (unsigned int i=0; i<subjets.size(); i++) {
      if(flavor == "ud"){
        if(i != index_bjet) correct_jet_flavor(*corrector, subjets[i], event, jec_uncertainty, direction);
      }
      else if(flavor == "b"){
        if(i == index_bjet) correct_jet_flavor(*corrector, subjets[i], event, jec_uncertainty, direction);
      }
      else{
        correct_jet_flavor(*corrector, subjets[i], event, jec_uncertainty, direction);
      }
    }
    topjet.set_subjets(move(subjets));
  }
  return true;
}

// note: implement here because only here (and not in the header file), the destructor of FactorizedJetCorrector is known
GenericSubJetCorrector_flavor::~GenericSubJetCorrector_flavor(){}
