#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <UHH2/common/include/JetIds.h>
#include "TRandom.h"
#include "TFormula.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>

class FactorizedJetCorrector;

namespace JERFiles {
  // DATA
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorUD;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA_flavorB;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA_flavorB;
  // MC
  extern const std::vector<std::string> Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorUD;
  extern const std::vector<std::string> Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorC;
  extern const std::vector<std::string> Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorS;
  extern const std::vector<std::string> Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorB;
  extern const std::vector<std::string> Summer16_07Aug2017_V15_L123_AK4PFchs_MC_flavorG;
}

void correct_jet_flavor(FactorizedJetCorrector & corrector, Jet & jet, const uhh2::Event & event, JetCorrectionUncertainty* jec_unc = NULL, int jec_unc_direction=0);

class GenericSubJetCorrector_flavor: public uhh2::AnalysisModule {
public:
  explicit GenericSubJetCorrector_flavor(uhh2::Context & ctx, const std::vector<std::vector<std::string>> & filenames, const std::string & collectionname);
  virtual bool process(uhh2::Event & event) override;
  virtual ~GenericSubJetCorrector_flavor();

private:
  void ReadFlavorFractionsFile();
  int IndexBJet(const uhh2::Event & event, TopJet jet);
  int IndexBJet_match_btag(const uhh2::Event & event, TopJet jet);
  std::vector<std::unique_ptr<FactorizedJetCorrector>> correctors;
  std::vector<JetCorrectionUncertainty*> jec_uncertainty;
  uhh2::Event::Handle<std::vector<TopJet> > h_jets;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  int direction = 0; // -1 = down, +1 = up, 0 = nominal
  std::vector<TGraph*> FlavorFractions;

};
