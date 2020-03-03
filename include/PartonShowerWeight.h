#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/JetCorrections.h"

class PartonShowerWeight : public uhh2::AnalysisModule {
 public:

  explicit PartonShowerWeight(uhh2::Context & ctx, std::string sys);
  virtual bool process(uhh2::Event & event) override;

 private:
  int weightIndex;

};
