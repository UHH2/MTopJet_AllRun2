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
  events_w = book<TH1F>("n_events", "events with weight", 1, 0, 1);
  events_nw = book<TH1F>("n_events_noweight", "events with no weight", 1, 0, 1);
}

void CountingEventHists::fill(const Event & event){
  double weight = event.weight;
  events_w->Fill(0.5, weight);
  events_nw->Fill(0.5, 1);
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
  h_wmass = book<TH1F>("wjet_jms_mass", "m_{W} [GeV]", 50, 0, 300);

  events = book<TH1F>("events", "events", 1, 0, 1);
  events_weight = book<TH1F>("events_weights", "events with gen weights", 1, 0, 1);
}

void JetMassScaleHists::fill(const Event & event){// dummy
  double weight = event.weight; // avoid warning
  weight +=0;                   // avoid warning
}

void JetMassScaleHists::fill_mass(const Event & event, const double mass, const double wmass){
  double weight = event.weight;
  h_mass->Fill(mass, weight);
  h_wmass->Fill(wmass, weight);
  events->Fill(0.5, weight);
  events_weight->Fill(0.5, 1);
}

// ###########################################################################################################################
// ###########################################################################################################################
MissingPtHist::MissingPtHist(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  MPT = book<TH1F>("MPT", "p_{T}^{miss}", 200,0,1000);
}

void MissingPtHist::fill(const Event & event){
  MPT->Fill(event.met->pt(), event.weight);
}

// ###########################################################################################################################
// ###########################################################################################################################
PositionBTagHists::PositionBTagHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  Position1stBtag = book<TH1F>("Position1stBtag", "Located", 3,0,3);
  Position2ndBtag = book<TH1F>("Position2ndBtag", "Located", 3,0,3);
  PositionBoth = book<TH1F>("PositionBoth", "Located", 9,-1,8);

  h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined");
  h_lepjets=ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined");
  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
}

void PositionBTagHists::fill(const Event & event){
  std::vector<TopJet> hadjets = event.get(h_hadjets);
  std::vector<TopJet> lepjets = event.get(h_lepjets);
  std::vector<TopJet> fatjets = event.get(h_fatjets);

  // get had jet from fat jets for softdrop mass
  int index_had = 0;
  int index_lep = 1;
  if(deltaR(lepjets.at(0), fatjets.at(0)) < deltaR(hadjets.at(0), fatjets.at(0))){
    index_had = 1;
    index_lep = 0;
  }

  std::vector<double> btags; // Using a vector to identifie events with ambigous btag results (two high btags i.e. through a c-quark)
  for(const auto& j : *event.jets){
    double btag = j.btag_DeepJet();
    btags.push_back(btag);
  }
  std::vector<double> sorted = btags; // Using a vector to identifie events with ambigous btag results (two high btags i.e. through a c-quark)
  sort(sorted.begin(), sorted.end(), greater<double>());

  int index_first = -1;
  int index_second = -1;

  auto it = find(btags.begin(), btags.end(), sorted[0]);
  if (it != btags.end()) index_first = it - btags.begin();
  else throw runtime_error("<E> Highest btag not in vector - check");

  auto ak4 = *event.jets;
  bool b1_in_had = deltaR(ak4[index_first], fatjets[index_had])<1.2;
  bool b1_in_lep = deltaR(ak4[index_first], fatjets[index_lep])<1.2;

  int match1 = b1_in_had?0:b1_in_lep?1:2;
  Position1stBtag->Fill(match1, event.weight);

  if(btags.size()>1){
    auto it = find(btags.begin(), btags.end(), sorted[1]);
    if (it != btags.end()) index_second = it - btags.begin();
    else throw runtime_error("<E> 2nd highest btag not in vector - check");

    bool b2_in_had = deltaR(ak4[index_second], fatjets[index_had])<1.2;
    bool b2_in_lep = deltaR(ak4[index_second], fatjets[index_lep])<1.2;

    bool noHad   =  !b1_in_had && !b2_in_had;
    bool noLep   =  !b1_in_lep && !b2_in_lep;
    bool bothHad =   b1_in_had &&  b2_in_had;
    bool bothLep =   b1_in_lep &&  b2_in_lep;
    bool oneHad  = (!b1_in_had &&  b2_in_had) || ( b1_in_had && !b2_in_had);
    bool oneLep  = (!b1_in_lep &&  b2_in_lep) || ( b1_in_lep && !b2_in_lep);


    int match=-1;
    if( noHad && noLep )        match=0;
    else if( bothHad )          match=1;
    else if( bothLep )          match=2;
    else if( oneHad && oneLep)  match=3;
    else if( oneHad )           match=4;
    else if( oneLep )           match=5;
    else if( oneHad && !oneLep) match=6;
    else if(!oneHad && oneLep)  match=7;

    int match2 = b2_in_had?0:b2_in_lep?1:2;
    Position2ndBtag->Fill(match2, event.weight);
    PositionBoth->Fill(match, event.weight);
  }

}
