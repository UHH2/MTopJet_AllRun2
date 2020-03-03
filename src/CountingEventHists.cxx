#include <UHH2/MTopJet/include/CountingEventHists.h>


CountingEventHists::CountingEventHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  Events = book<TH1F>("Events", "events", 1, 0, 1);
}

void CountingEventHists::fill(const Event & event){
  double weight = event.weight;
  Events->Fill(1, weight);
}
