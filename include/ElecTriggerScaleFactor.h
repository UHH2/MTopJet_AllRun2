#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/core/include/LorentzVector.h>

#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>

#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/utils.h>
#include <UHH2/MTopJet/include/GenSelections.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>



namespace uhh2 {
  class ElectronTriggerSF: public uhh2::AnalysisModule{

  public:
    explicit ElectronTriggerSF(uhh2::Context & ctx);
    virtual bool process(uhh2::Event & event) override;

  private:
    double sf;
    TString sysdirection;

  };

}
