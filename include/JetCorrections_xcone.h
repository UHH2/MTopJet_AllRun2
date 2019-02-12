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
  std::unique_ptr<GenericSubJetCorrector_flavor> jet_corrector_MC_flavor;
  std::unique_ptr<GenericSubJetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;
  bool isMC;
  uhh2::Event::Handle<std::vector<TopJet>>h_topjets;

};


class JER_Smearer_xcone : public uhh2::AnalysisModule {

public:
  JER_Smearer_xcone();
  void init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen, const std::string& fat_sub);
  virtual bool process(uhh2::Event & event) override;

private:
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;
  bool use_subjets;
  uhh2::Event::Handle<std::vector<TopJet>>    h_rectopjets_;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets_;
};
