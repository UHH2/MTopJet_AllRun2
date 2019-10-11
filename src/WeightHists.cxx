#include <UHH2/MTopJet/include/WeightHists.h>

using namespace uhh2;

WeightHists::WeightHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  h_weight_w = book<TH1F>("weight_w", "weight (filled with weights)", 200, -1, 1);
  h_weight_1 = book<TH1F>("weight_1", "weight (filled without weights)", 200, -1, 1);
}



void WeightHists::fill(const Event & event){
  double weight = event.weight;
  h_weight_w->Fill(weight, weight);
  h_weight_1->Fill(weight, 1);
}
