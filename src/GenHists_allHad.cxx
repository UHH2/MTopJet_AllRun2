#include "UHH2/MTopJet/include/GenHists_allHad.h"


GenHists_allHad::GenHists_allHad(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // book all histograms here
  number_top = book<TH1F>("number_top", "number of tops", 10, 0, 10);
  top_pt = book<TH1F>("top_pt", "p_{T}", 100, 0, 1000);
  top_Dphi = book<TH1F>("top_Dphi", "#{Delta} #{Phi} (top, antitop)", 30, 0, 4);

  index_top_antitop = book<TH1F>("index_top_antitop", "top and antitop are matched to different jets", 2, -0.5, 1.5);

  bjet_index = book<TH1F>("bjet_index", "jet index including b-quark", 3, -0.5, 2.5);
  bjet_index_nomatch_sub = book<TH1F>("bjet_index_nomatch_sub", "no match found inside fatjet", 1, -0.5, 0.5);
  bjet_index_nomatch_fat = book<TH1F>("bjet_index_nomatch_fat", "b-quark not inside fatjet", 1, -0.5, 0.5);

  b_in_leading_sub = book<TH1F>("b_in_leading_sub", "b-quark inside leading subjet", 2, -0.5, 1.5);

  deltaR_b_bjet = book<TH1F>("deltaR_b_bjet", "#{Delta}R(b quark, b jet)", 20, 0.0, 0.5);

  // handle for ttbargen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  // handle for jets
  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
  h_fatjets_gen=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

}



void GenHists_allHad::fill(const Event & event){

  const auto & ttbargen = event.get(h_ttbargen);
  // get weight
  double weight = event.weight;

  // get jets
  std::vector<TopJet> fatjets = event.get(h_fatjets);
  std::vector<GenTopJet> fatjets_gen = event.get(h_fatjets_gen);

  // count tops
  int n_top = 0;
  std::vector<GenParticle>* genparts = event.genparticles;
  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle p = genparts->at(i);
    if(abs(p.pdgId()) == 6){
      ++n_top;
    }
  }

  number_top->Fill(n_top, weight);

  // get particles from ttbargen class
  GenParticle top, antitop, bot, antibot;
  if(ttbargen.IsTopHadronicDecay() && ttbargen.IsAntiTopHadronicDecay()){
    top = ttbargen.Top();
    antitop = ttbargen.Antitop();
    bot = ttbargen.bTop();
    antibot = ttbargen.bAntitop();
  }
  else return;


  double deltaPhi = fabs(top.phi() - antitop.phi());
  if(deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;
  top_Dphi->Fill(deltaPhi, weight);

  // match tops with jets
  if(fatjets.size() != 2) return;
  int i_top, i_antitop;
  if(deltaR(top, fatjets[0]) < deltaR(top, fatjets[1])){
    i_top = 0;
  }
  else{
    i_top = 1;
  }
  if(deltaR(antitop, fatjets[0]) < deltaR(antitop, fatjets[1])){
    i_antitop = 0;
  }
  else{
    i_antitop = 1;
  }

  bool good_jet_top = (fatjets[i_top].pt() > 400);
  bool good_jet_antitop = (fatjets[i_antitop].pt() > 400);

  if(good_jet_top && good_jet_antitop){
    if(i_top == i_antitop) index_top_antitop->Fill(0.0,weight);
    else                   index_top_antitop->Fill(1.0,weight);
  }

  if(good_jet_top)     top_pt->Fill(top.pt(), weight);
  if(good_jet_antitop) top_pt->Fill(antitop.pt(), weight);

  // check if the b-quark ends up in the pT-leading subjet
  // sort jets by pT and get index of b-jet
  double dR_bot = 1000;
  double dR_antibot = 1000;
  double dR_temp = 1000;
  int bot_index = -1;
  int antibot_index = -1;
  std::vector<Jet> subjets_top = fatjets[i_top].subjets();
  std::vector<Jet> subjets_antitop = fatjets[i_antitop].subjets();
  sort_by_pt<Jet> (subjets_top);
  sort_by_pt<Jet> (subjets_antitop);


  if(good_jet_top){
    if( deltaR(bot, fatjets[i_top]) < 1.2 ){
      for(unsigned int i=0; i<subjets_top.size(); i++){
        dR_temp = deltaR(bot, subjets_top[i]);
        if(dR_temp < dR_bot && dR_temp < 0.4){
          dR_bot = dR_temp;
          bot_index = i;
        }
      }
      if(bot_index != -1){
        bjet_index->Fill(bot_index, weight);
        deltaR_b_bjet->Fill(deltaR(bot, subjets_top[bot_index]), weight);
      }
      if(bot_index == -1) bjet_index_nomatch_sub->Fill(0.0, weight);
    }
    else bjet_index_nomatch_fat->Fill(0.0, weight);
  }

  if(good_jet_antitop){
    if( deltaR(antibot, fatjets[i_antitop]) < 1.2 ){
      for(unsigned int i=0; i<subjets_antitop.size(); i++){
        dR_temp = deltaR(antibot, subjets_antitop[i]);
        if(dR_temp < dR_antibot && dR_temp < 0.4){
          dR_antibot = dR_temp;
          antibot_index = i;
        }
      }
      if(antibot_index != -1){
        bjet_index->Fill(antibot_index, weight);
        deltaR_b_bjet->Fill(deltaR(antibot,subjets_antitop[antibot_index]), weight);
      }
      if(antibot_index == -1) bjet_index_nomatch_sub->Fill(0.0, weight);
    }
    else bjet_index_nomatch_fat->Fill(0.0, weight);
  }
  ////

  for(auto fatjet: fatjets){
    if(fatjet.pt() < 400) continue;
    bool binjet;
    std::vector<Jet> subjets = fatjet.subjets();
    int i_ptmax = -1;
    double ptmax = 0;
    for(unsigned int i=0; i<subjets.size(); i++){
      if(subjets[i].pt() > ptmax){
        ptmax = subjets[i].pt();
        i_ptmax = i;
      }
    }
    if(i_ptmax == -1) return;
    if(deltaR(bot, subjets[i_ptmax])<0.4 || deltaR(antibot, subjets[i_ptmax])<0.4) binjet = true;
    else binjet = false;

    if(binjet) b_in_leading_sub->Fill(1.0,weight);
    else       b_in_leading_sub->Fill(0.0,weight);
  }



  return;
}
