#include "UHH2/MTopJet/include/GenHists_particles.h"


GenHists_particles::GenHists_particles(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  W_pt = book<TH1F>("W pt", "p_{T}", 100, 0, 1000);
  top_pt = book<TH1F>("top pt", "p_{T}", 100, 0, 1000);
  elec_pt = book<TH1F>("elec pt", "p_{T}", 100, 0, 1000);
  muon_pt = book<TH1F>("muon pt", "p_{T}", 100, 0, 1000);




}



void GenHists_particles::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- get W  -------------------------------------------
  //---------------------------------------------------------------------------------------
  GenParticle W, top;
  std::vector<GenParticle>* genparts = event.genparticles;

  for (unsigned int i=0; i<(genparts->size()); ++i){
    GenParticle p = genparts->at(i);
    if(abs(p.pdgId()) == 24){
      W_pt->Fill(p.pt(), weight);
    }
    if(abs(p.pdgId()) == 6){
      top_pt->Fill(p.pt(), weight);
    }
    if(abs(p.pdgId()) == 11){
      elec_pt->Fill(p.pt(), weight);
    }
    if(abs(p.pdgId()) == 13){
      muon_pt->Fill(p.pt(), weight);
    }
  }

}


