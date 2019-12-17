#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/BTagCalibrationStandalone.h"

#include <math.h>
#include <vector>
#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

class ScaleFactorHists: public uhh2::Hists {
public:
  // use the same constructor arguments as Hists for forwarding:
  ScaleFactorHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

protected:

  TH1F *out_of_bounds;

  BTag btag_;
  std::unique_ptr<BTagCalibrationReader> calib_up_;
  std::unique_ptr<BTagCalibrationReader> calib_;
  std::unique_ptr<BTagCalibrationReader> calib_down_;
}
