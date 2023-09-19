#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TRandom.h"
#include "TFormula.h"
#include "TFile.h"
#include "TF1.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/MTopJet/include/GenericSubJetCorrector_flavor.h"

#include <iostream>
#include <fstream>

std::vector<TopJet> set_JEC_factor(std::vector<TopJet> jets);

class JetCorrections_xcone : public uhh2::AnalysisModule {

public:
  JetCorrections_xcone();
  void init(uhh2::Context & ctx, const std::string& jet_collection_rec);
  virtual bool process(uhh2::Event & event) override;

private:
  std::unique_ptr<YearSwitcher> jet_corrector_MC, jet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18;
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;
  bool isMC;
  uhh2::Event::Handle<std::vector<TopJet>>h_topjets;

};


class JER_Smearer_xcone : public uhh2::AnalysisModule {

public:
  JER_Smearer_xcone();
  void init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen, const std::string& fat_sub);
  std::vector<double> JER_factors(int index);
  virtual bool process(uhh2::Event & event) override;


private:
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;
  bool use_subjets;
  uhh2::Event::Handle<std::vector<TopJet>>    h_rectopjets_;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets_;
  std::vector<std::vector<double>> jer_factors_wo, jer_factors_with;
};
