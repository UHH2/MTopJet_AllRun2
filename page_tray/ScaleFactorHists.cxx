#include <UHH2/MTopJet/include/ScaleFactorHists.h>

using namespace uhh2;

ScaleFactorHists::ScaleFactorHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  out_of_bounds = book<TH1F>("Not in bounds", "true", 1, -2, 2);

  btag_(BTag(tagger, wp));
  BTagEntry::OperatingPoint op = (BTagEntry::OperatingPoint) btag_.GetWorkingPoint();
  calib_up_.reset(new BTagCalibrationReader(op, "up"));
  calib_.reset(new BTagCalibrationReader(op, "central"));
  calib_down_.reset(new BTagCalibrationReader(op, "down"));
}

void ScaleFactorHists::fill(const Event & event){



  auto btagentry_flav = flav == 5 ? BTagEntry::FLAV_B : (
                            flav == 4 ? BTagEntry::FLAV_C :
                                BTagEntry::FLAV_UDSG);

  auto sf_bounds = calib_->min_max_pt(btagentry_flav, abs_eta);

  float pt_for_eval = pt;
  bool is_out_of_bounds = false;
  if (pt < sf_bounds.first) {
    pt_for_eval = sf_bounds.first + 1e-5;
    is_out_of_bounds = true;
  } else if (pt > sf_bounds.second) {
    pt_for_eval = sf_bounds.second - 0.1;
    is_out_of_bounds = true;
  }

  out_of_bounds->Fill(is_out_of_bounds, weight);

  return;
}
