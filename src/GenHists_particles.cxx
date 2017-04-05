#include "UHH2/MTopJet/include/GenHists_particles.h"


GenHists_particles::GenHists_particles(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  number_top = book<TH1F>("number_top", "number of tops", 10, 0, 10);
  hadtop_pt = book<TH1F>("hadtop_pt", "p_{T}", 100, 0, 1000);
  leptop_pt = book<TH1F>("leptop_pt", "p_{T}", 100, 0, 1000);
  lepton_pt = book<TH1F>("lepton_pt", "p_{T}", 100, 0, 1000);
  deltaR_hadtop_b = book<TH1F>("deltaR_hadtop_b", "#Delta R(had. top, b)", 30, 0, 6);
  deltaR_leptop_b = book<TH1F>("deltaR_leptop_b", "#Delta R(lep. top, b)", 30, 0, 6);
  deltaR_lep_b = book<TH1F>("deltaR_lep_b", "#Delta R(lep, b)", 30, 0, 6);
  deltaR_lep_neu = book<TH1F>("deltaR_lep_neu", "#Delta R(lepton, neutrino)", 30, 0, 6);
  deltaR_hadtop_leptop = book<TH1F>("deltaR_hadtop_leptop", "#Delta R(had. top, lep. top)", 30, 0, 6);
  deltaPhi_hadtop_leptop = book<TH1F>("deltaPhi_hadtop_leptop", "#Delta #Phi(had. top, lep. top)", 30, 0, 6);


  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
}



void GenHists_particles::fill(const Event & event){

  const auto & ttbargen = event.get(h_ttbargen);
  // get weight
  double weight = event.weight;

  // cout tops
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
  GenParticle tophad, toplep, bot, q1, q2, bot_lep, lep1, lep2, lepton, neutrino;
  if(ttbargen.IsTopHadronicDecay()){
    tophad = ttbargen.Top();
    toplep = ttbargen.Antitop();
    bot = ttbargen.bTop();
    q1 = ttbargen.Wdecay1();
    q2 = ttbargen.Wdecay2();
    bot_lep = ttbargen.bAntitop();
    lep1 = ttbargen.WMinusdecay1();
    lep2 = ttbargen.WMinusdecay2();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    tophad = ttbargen.Antitop();
    toplep = ttbargen.Top();
    bot = ttbargen.bAntitop();
    q1 = ttbargen.WMinusdecay1();
    q2 = ttbargen.WMinusdecay2();
    bot_lep = ttbargen.bTop();
    lep1 = ttbargen.Wdecay1();
    lep2 = ttbargen.Wdecay2();
  }
  else if(!(ttbargen.IsTopHadronicDecay() || ttbargen.IsAntiTopHadronicDecay())){
    return;
  }

  //check which lep is neutrino and which is elec/muon
  if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
    lepton = lep1;
    neutrino = lep2;
  }
  else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
    lepton = lep2;
    neutrino = lep1;
  }

  hadtop_pt->Fill(tophad.pt(), weight);
  leptop_pt->Fill(toplep.pt(), weight);
  lepton_pt->Fill(lepton.pt(), weight);
  deltaR_hadtop_b->Fill(deltaR(tophad, bot), weight);
  deltaR_leptop_b->Fill(deltaR(toplep, bot_lep), weight);
  deltaR_lep_b->Fill(deltaR(lepton, bot_lep), weight);
  deltaR_lep_neu->Fill(deltaR(lepton, neutrino), weight);
  deltaR_hadtop_leptop->Fill(deltaR(tophad, toplep), weight);
  deltaPhi_hadtop_leptop->Fill(abs(tophad.phi() - toplep.phi()), weight);

  return;
}


