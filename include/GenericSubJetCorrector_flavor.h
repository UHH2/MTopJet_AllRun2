#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TRandom.h"
#include "TFormula.h"
#include "TFile.h"
#include "TF1.h"

#include <iostream>
#include <fstream>
#include <string>

class FactorizedJetCorrector;

namespace JERFiles {
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK4PFchs_MC_flavorB;
}

void correct_jet_flavor(FactorizedJetCorrector & corrector, Jet & jet, const uhh2::Event & event, JetCorrectionUncertainty* jec_unc = NULL, int jec_unc_direction=0);

class GenericSubJetCorrector_flavor: public uhh2::AnalysisModule {
public:
  explicit GenericSubJetCorrector_flavor(uhh2::Context & ctx, const std::vector<std::string> & filenames, const std::string & collectionname, std::string flavor_);
  virtual bool process(uhh2::Event & event) override;
  virtual ~GenericSubJetCorrector_flavor();

private:
  std::unique_ptr<FactorizedJetCorrector> corrector;
  uhh2::Event::Handle<std::vector<TopJet> > h_jets;
  JetCorrectionUncertainty* jec_uncertainty;
  int direction = 0; // -1 = down, +1 = up, 0 = nominal
  std::string flavor;
};
