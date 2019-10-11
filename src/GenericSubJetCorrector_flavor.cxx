#include <UHH2/MTopJet/include/GenericSubJetCorrector_flavor.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/core/include/Utils.h>

#include <UHH2/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <UHH2/JetMETObjects/interface/JetCorrectorParameters.h>

#include <string>

using namespace std;
using namespace uhh2;

// TODO: replace data JER files with v15

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//            replace L2Relative with flavor
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/*
██████   █████  ████████  █████      ██    ██ ██████
██   ██ ██   ██    ██    ██   ██     ██    ██ ██   ██
██   ██ ███████    ██    ███████     ██    ██ ██   ██
██   ██ ██   ██    ██    ██   ██     ██    ██ ██   ██
██████  ██   ██    ██    ██   ██      ██████  ██████
*/

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

/*
██████   █████  ████████  █████      ███████
██   ██ ██   ██    ██    ██   ██     ██
██   ██ ███████    ██    ███████     ███████
██   ██ ██   ██    ██    ██   ██          ██
██████  ██   ██    ██    ██   ██     ███████
*/
const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorS = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorS = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorS = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorS = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};

/*
██████   █████  ████████  █████       ██████
██   ██ ██   ██    ██    ██   ██     ██
██   ██ ███████    ██    ███████     ██
██   ██ ██   ██    ██    ██   ██     ██
██████  ██   ██    ██    ██   ██      ██████
*/
const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorC = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorC = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorC = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorC = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};

/*
██████   █████  ████████  █████      ██████
██   ██ ██   ██    ██    ██   ██     ██   ██
██   ██ ███████    ██    ███████     ██████
██   ██ ██   ██    ██    ██   ██     ██   ██
██████  ██   ██    ██    ██   ██     ██████
*/

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


/*
██████   █████  ████████  █████       ██████
██   ██ ██   ██    ██    ██   ██     ██
██   ██ ███████    ██    ███████     ██   ███
██   ██ ██   ██    ██    ██   ██     ██    ██
██████  ██   ██    ██    ██   ██      ██████
*/

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorG = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorG = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorG = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorG = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};


/*
███    ███  ██████     ██    ██ ██████
████  ████ ██          ██    ██ ██   ██
██ ████ ██ ██          ██    ██ ██   ██
██  ██  ██ ██          ██    ██ ██   ██
██      ██  ██████      ██████  ██████
*/



const std::vector<std::string> JERFiles::Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorUD = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_ud_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

/*
███    ███  ██████     ███████
████  ████ ██          ██
██ ████ ██ ██          ███████
██  ██  ██ ██               ██
██      ██  ██████     ███████
*/

const std::vector<std::string> JERFiles::Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorS = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_s_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

/*
███    ███  ██████      ██████
████  ████ ██          ██
██ ████ ██ ██          ██
██  ██  ██ ██          ██
██      ██  ██████      ██████
*/

const std::vector<std::string> JERFiles::Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorC = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_c_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

/*
███    ███  ██████     ██████
████  ████ ██          ██   ██
██ ████ ██ ██          ██████
██  ██  ██ ██          ██   ██
██      ██  ██████     ██████
*/

const std::vector<std::string> JERFiles::Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorB = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_b_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
};

/*
███    ███  ██████      ██████
████  ████ ██          ██
██ ████ ██ ██          ██   ███
██  ██  ██ ██          ██    ██
██      ██  ██████      ██████
*/

const std::vector<std::string> JERFiles::Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorG = {
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_g_L5Flavor_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_07Aug2017_V15_Flavor_Pythia8_MC/Summer16_07Aug2017_V15_Flavor_Pythia8_MC_L3Absolute_AK4PFchs.txt",
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

/*
███████ ██       █████  ██    ██  ██████  ██████       ██████  ██████  ██████  ██████  ███████  ██████ ████████  ██████  ██████
██      ██      ██   ██ ██    ██ ██    ██ ██   ██     ██      ██    ██ ██   ██ ██   ██ ██      ██         ██    ██    ██ ██   ██
█████   ██      ███████ ██    ██ ██    ██ ██████      ██      ██    ██ ██████  ██████  █████   ██         ██    ██    ██ ██████
██      ██      ██   ██  ██  ██  ██    ██ ██   ██     ██      ██    ██ ██   ██ ██   ██ ██      ██         ██    ██    ██ ██   ██
██      ███████ ██   ██   ████    ██████  ██   ██      ██████  ██████  ██   ██ ██   ██ ███████  ██████    ██     ██████  ██   ██
*/




GenericSubJetCorrector_flavor::GenericSubJetCorrector_flavor(uhh2::Context & ctx, const std::vector<std::vector<std::string>>&  filenames, const std::string & collectionname){
  for(unsigned int i=0; i<filenames.size(); i++){
    correctors.push_back(build_corrector_flavor(filenames[i]));
    direction = 0;
    jec_uncertainty.push_back(corrector_uncertainty_flavor(ctx, filenames[i], direction));
  }
  h_jets = ctx.get_handle<std::vector<TopJet> >(collectionname);
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  ReadFlavorFractionsFile();
}

bool GenericSubJetCorrector_flavor::process(uhh2::Event & event){

  const auto topjets = &event.get(h_jets);
  assert(topjets);
  for(auto & topjet : *topjets){
    auto subjets = topjet.subjets();

    for (unsigned int i=0; i<subjets.size(); i++) {
      // first do a copy of every subjet for every flavor corrector
      // also do one copy that becomes the final Jet in the end
      vector<Jet> raw_jets;
      for(unsigned int j=0; j<correctors.size(); j++){
        raw_jets.push_back(subjets[i]);
      }
      // also store raw pt to get correct flavor fraction
      // set a maximum value for pT here because flavor fraction lacks of statistics at very high pT
      double pt_raw = subjets[i].pt();
      double pt_max = 500.;
      if(pt_raw > pt_max) pt_raw = pt_max;
      // now put raw jet in every corrector and extract L2L3 factor
      // also store L1_factors (all of them should be the same)
      vector<double> L1_factors;
      vector<double> L2L3_factors;
      for(unsigned int j=0; j<correctors.size(); j++){
        correct_jet_flavor(*correctors[j], raw_jets[j], event, jec_uncertainty[j], direction);
        // attention: factor_raw brings corrected jet back to raw
        // but L1factor_raw brings raw back to L1
        // since I only need to extract L2L3, it has to be like this:
        double L2L3_f = 1./(raw_jets[j].JEC_L1factor_raw() * raw_jets[j].JEC_factor_raw() );
        L2L3_factors.push_back(L2L3_f);
        L1_factors.push_back(raw_jets[j].JEC_L1factor_raw());
      }
      // now calculate final L2L3 factor from flavor_fractions
      double final_L2L3 = 0;
      for(unsigned int k=0; k<L2L3_factors.size(); k++){
        double flavor_fraction = FlavorFractions[k]->Eval(pt_raw);
        final_L2L3 += L2L3_factors[k] * flavor_fraction;
      }
      // multiply with L1 for total correction factor
      double final_L1 = L1_factors[0];
      double final_factor = final_L1 * final_L2L3;
      // apply total correction to jet, also set correct factor_raw
      LorentzVector jet_v4_corrected = subjets[i].v4() * final_factor;
      subjets[i].set_v4(jet_v4_corrected);
      subjets[i].set_JEC_factor_raw(1. / final_L2L3);
      subjets[i].set_JEC_L1factor_raw(final_L1);

    }
    topjet.set_subjets(move(subjets));
  }
  return true;
}

// note: implement here because only here (and not in the header file), the destructor of FactorizedJetCorrector is known
GenericSubJetCorrector_flavor::~GenericSubJetCorrector_flavor(){}

void GenericSubJetCorrector_flavor::ReadFlavorFractionsFile(){
  TString dir = "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/FlavorFractions/";
  TString filename;
  filename = dir + "FlavorFractions.root";
  TFile *file = new TFile(filename);
  FlavorFractions.push_back((TGraph*)file->Get("ud_fraction"));
  FlavorFractions.push_back((TGraph*)file->Get("s_fraction"));
  FlavorFractions.push_back((TGraph*)file->Get("c_fraction"));
  FlavorFractions.push_back((TGraph*)file->Get("b_fraction"));
  FlavorFractions.push_back((TGraph*)file->Get("g_fraction"));
}

//
// // function to find b-quark in jet
// int GenericSubJetCorrector_flavor::IndexBJet(const uhh2::Event & event, TopJet jet){
//   int index = -1;
//   vector<Jet> subjets = jet.subjets();
//   const auto & ttbargen = event.get(h_ttbargen);
//   GenParticle b1 = ttbargen.bTop();
//   GenParticle b2 = ttbargen.bAntitop();
//   GenParticle bot;
//   // first match b to TopJet
//   if(deltaR(b1, jet) < deltaR(b2, jet)) bot = b1;
//   else bot = b2;
//
//   // now find subjet closest to bottom quark
//   double dR = 100;
//   for(unsigned int i=0; i<subjets.size(); i++){
//     if(deltaR(bot, subjets[i]) < dR){
//       dR = deltaR(bot, subjets[i]);
//       if(dR < 0.4) index = i;
//     }
//   }
//   return index;
// }
//
// // function to find b-quark in jet
// int GenericSubJetCorrector_flavor::IndexBJet_match_btag(const uhh2::Event & event, TopJet jet){
//   int index = -1;
//   vector<Jet> subjets = jet.subjets();
//   vector<Jet> btags;
//   // first, create a vector of btag that are inside fatjet
//   for(const auto& ak4 : *event.jets){
//     if(CSVBTag(CSVBTag::WP_TIGHT)(ak4, event)){
//       if( deltaR(ak4, jet) < 1.2){
//         btags.push_back(ak4);
//       }
//     }
//   }
//   if(btags.size() == 0 ) return -1;
//   // if there are still multiple jets, take leading in pT
//   sort_by_pt<Jet>(btags);
//   // now find subjet closest to btag ak4
//   double dR_min = 100;
//   for(unsigned int i=0; i<subjets.size(); i++){
//     if(deltaR(btags[0], subjets[i]) < dR_min){
//       dR_min = deltaR(btags[0], subjets[i]);
//       if(dR_min < 0.4) index = i;
//     }
//   }
//
//   return index;
// }
