#include <UHH2/MTopJet/include/ControlHists.h>

using namespace uhh2;

// ###########################################################################################################################
// ###########################################################################################################################
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


// ###########################################################################################################################
// ###########################################################################################################################
CountingEventHists::CountingEventHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  events_w = book<TH1F>("Events_weight", "events with weight", 2, 0, 2);
  events_nw = book<TH1F>("Events_noweight", "events with no weight", 2, 0, 2);
}

void CountingEventHists::fill(const Event & event){
  double weight = event.weight;
  events_w->Fill(1, weight);
  events_nw->Fill(1, 1);
}


// ###########################################################################################################################
// ###########################################################################################################################
WeightRangeHists::WeightRangeHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  h_weight_big_w = book<TH1F>("weight_big_w", "weight (filled with weights)", 1000, 0, 10000);
  h_weight_big_1 = book<TH1F>("weight_big_1", "weight (filled without weights)", 1000, 0, 10000);

  h_weight_medium_w = book<TH1F>("weight__medium_w", "weight (filled with weights)", 1000, 0, 1000);
  h_weight_medium_1 = book<TH1F>("weight__medium_1", "weight (filled without weights)", 1000, 0, 1000);

  h_weight_small_w = book<TH1F>("weight__small_w", "weight (filled with weights)", 100, 0, 10);
  h_weight_small_1 = book<TH1F>("weight__small_1", "weight (filled without weights)", 100, 0, 10);

  h_weight_negativ_big_w = book<TH1F>("weight_negativ_big_w", "weight (filled with weights)", 1000, -10000, 0);
  h_weight_negativ_big_1 = book<TH1F>("weight_negativ_big_1", "weight (filled without weights)", 1000, -10000, 0);

  h_weight_negativ_medium_w = book<TH1F>("weight_negativ_medium_w", "weight (filled with weights)", 1000, -1000, 0);
  h_weight_negativ_medium_1 = book<TH1F>("weight_negativ_medium_1", "weight (filled without weights)", 1000, -1000, 0);

  h_weight_negativ_small_w = book<TH1F>("weight_negativ_small_w", "weight (filled with weights)", 100, -10, 0);
  h_weight_negativ_small_1 = book<TH1F>("weight_negativ_small_1", "weight (filled without weights)", 100, -10, 0);

  events = book<TH1F>("Events", "events", 1, 0, 1);
  events_weight = book<TH1F>("Events", "events with weights", 1, 0, 1);
}

void WeightRangeHists::fill(const Event & event){
  double weight = event.weight;
  h_weight_big_w->Fill(weight, weight);
  h_weight_big_1->Fill(weight, 1);

  h_weight_medium_w->Fill(weight, weight);
  h_weight_medium_1->Fill(weight, 1);

  h_weight_small_w->Fill(weight, weight);
  h_weight_small_1->Fill(weight, 1);

  h_weight_negativ_big_w->Fill(weight, weight);
  h_weight_negativ_big_1->Fill(weight, 1);

  h_weight_negativ_medium_w->Fill(weight, weight);
  h_weight_negativ_medium_1->Fill(weight, 1);

  h_weight_negativ_small_w->Fill(weight, weight);
  h_weight_negativ_small_1->Fill(weight, 1);

  events->Fill(1, weight);
  events_weight->Fill(1, 1);
}

// ###########################################################################################################################
// ###########################################################################################################################
GenWeightRangeHists::GenWeightRangeHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  h_weight_big_w = book<TH1F>("gen_weight_big_w", "gen weight (filled with weights)", 1000, 0, 10000);
  h_weight_big_1 = book<TH1F>("gen_weight_big_1", "gen weight (filled without weights)", 1000, 0, 10000);

  h_weight_medium_w = book<TH1F>("gen_weight_medium_w", "gen weight (filled with weights)", 1000, 0, 1000);
  h_weight_medium_1 = book<TH1F>("gen_weight_medium_1", "gen weight (filled without weights)", 1000, 0, 1000);

  h_weight_small_w = book<TH1F>("gen_weight_small_w", "gen weight (filled with weights)", 100, 0, 10);
  h_weight_small_1 = book<TH1F>("gen_weight_small_1", "gen weight (filled without weights)", 100, 0, 10);

  h_weight_negativ_big_w = book<TH1F>("gen_weight_negativ_big_w", "gen weight (filled with weights)", 1000, -10000, 0);
  h_weight_negativ_big_1 = book<TH1F>("gen_weight_negativ_big_1", "gen weight (filled without weights)", 1000, -10000, 0);

  h_weight_negativ_medium_w = book<TH1F>("gen_weight_negativ_medium_w", "gen weight (filled with weights)", 1000, -1000, 0);
  h_weight_negativ_medium_1 = book<TH1F>("gen_weight_negativ_medium_1", "gen weight (filled without weights)", 1000, -1000, 0);

  h_weight_negativ_small_w = book<TH1F>("gen_weight_negativ_small_w", "gen weight (filled with weights)", 100, -10, 0);
  h_weight_negativ_small_1 = book<TH1F>("gen_weight_negativ_small_1", "gen weight (filled without weights)", 100, -10, 0);

  events = book<TH1F>("events", "events", 1, 0, 1);
  events_weight = book<TH1F>("events_weights", "events with gen weights", 1, 0, 1);

}

void GenWeightRangeHists::fill(const Event & event){
  double gen_weight = event.weight;
  h_weight_big_w->Fill(gen_weight, gen_weight);
  h_weight_big_1->Fill(gen_weight, 1);

  h_weight_medium_w->Fill(gen_weight, gen_weight);
  h_weight_medium_1->Fill(gen_weight, 1);

  h_weight_small_w->Fill(gen_weight, gen_weight);
  h_weight_small_1->Fill(gen_weight, 1);

  h_weight_negativ_big_w->Fill(gen_weight, gen_weight);
  h_weight_negativ_big_1->Fill(gen_weight, 1);

  h_weight_negativ_medium_w->Fill(gen_weight, gen_weight);
  h_weight_negativ_medium_1->Fill(gen_weight, 1);

  h_weight_negativ_small_w->Fill(gen_weight, gen_weight);
  h_weight_negativ_small_1->Fill(gen_weight, 1);

  events->Fill(1, gen_weight);
  events_weight->Fill(1, 1);
}

// ###########################################################################################################################
// ###########################################################################################################################
JetMassScaleHists::JetMassScaleHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  h_mass = book<TH1F>("hadjet_jms_mass", "m_{jet} [GeV]", 50, 0, 500);

  events = book<TH1F>("events", "events", 1, 0, 1);
  events_weight = book<TH1F>("events_weights", "events with gen weights", 1, 0, 1);
}

void JetMassScaleHists::fill(const Event & event){// dummy
  double weight = event.weight; // avoid warning
  weight +=0;                   // avoid warning
}

void JetMassScaleHists::fill_mass(const Event & event, const double mass){
  double weight = event.weight;
  h_mass->Fill(mass, weight);
  events->Fill(1, weight);
  events_weight->Fill(1, 1);
}
