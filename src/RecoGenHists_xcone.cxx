#include "UHH2/MTopJet/include/RecoGenHists_xcone.h"


RecoGenHists_xcone::RecoGenHists_xcone(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  MassReso = book<TH1F>("MassResolution", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  PtReso = book<TH1F>("PtResolution", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  DeltaR_Rec1_Gen1 = book<TH1F>("DeltaR_Rec1-Gen1", "#Delta R (rec1, gen1) ", 40, 0, 4);
  DeltaR_Rec2_Gen2 = book<TH1F>("DeltaR_Rec2-Gen2", "#Delta R (rec2, gen2) ", 40, 0, 4);
  pt_rec1_gen1_beforeMatching = book<TH1F>("pt_rec1_gen1_beforeMatching", "p^{rec}_{T, jet1} - p^{gen}_{T, jet1}", 30, -300, 300);
  pt_rec2_gen2_beforeMatching = book<TH1F>("pt_rec2_gen2_beforeMatching", "p^{rec}_{T, jet2} - p^{gen}_{T, jet2}", 30, -300, 300);
  pt_rec1_gen1_afterMatching = book<TH1F>("pt_rec1_gen1_afterMatching", "p^{rec}_{T, jet1} - p^{gen}_{T, jet1}", 30, -300, 300);
  pt_rec2_gen2_afterMatching = book<TH1F>("pt_rec2_gen2_afterMatching", "p^{rec}_{T, jet2} - p^{gen}_{T, jet2}", 30, -300, 300);
  GenJetMass = book<TH1F>("GenJetMass", "M^{gen}_{jet1}", 50, 0, 500);

  RecGenMass = book<TH2F>("RecGenMass", "x=M_Gen y=M_Rec", 50, 0, 500., 50, 0, 500.);
  RecGenPT = book<TH2F>("RecGenPT", "x=PT_Gen y=PT_Rec", 50, 0, 2000., 50, 0, 2000.);

  // handle for clustered jets
  h_recjets_had=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
  h_recjets_lep=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined");
  h_genjets_had=ctx.get_handle<std::vector<Particle>>("GEN_XCone33_had_Combined");
  h_genjets_lep=ctx.get_handle<std::vector<Particle>>("GEN_XCone33_lep_Combined");

}



void RecoGenHists_xcone::fill(const Event & event){


  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
 // define all objects needed
  std::vector<Jet> rec_hadjets = event.get(h_recjets_had);
  std::vector<Jet> rec_lepjets = event.get(h_recjets_lep);
  std::vector<Particle> gen_hadjets = event.get(h_genjets_had);
  std::vector<Particle> gen_lepjets = event.get(h_genjets_lep);


  TLorentzVector rec_jet1_v4, rec_jet2_v4, gen_jet1_v4, gen_jet2_v4;
  Jet rec_jet1, rec_jet2;
  Particle gen_jet1, gen_jet2;
  if(rec_hadjets.size()>0){
    rec_jet1 = rec_hadjets.at(0);
    rec_jet1_v4.SetPxPyPzE(rec_jet1.v4().Px(), rec_jet1.v4().Py(), rec_jet1.v4().Pz(), rec_jet1.v4().E());
  }
  if(gen_hadjets.size()>0){
    gen_jet1 = gen_hadjets.at(0);
    gen_jet1_v4.SetPxPyPzE(gen_jet1.v4().Px(), gen_jet1.v4().Py(), gen_jet1.v4().Pz(), gen_jet1.v4().E());
  }

  if(rec_lepjets.size()>0){
    rec_jet2 = rec_lepjets.at(0);
    rec_jet2_v4.SetPxPyPzE(rec_jet2.v4().Px(), rec_jet2.v4().Py(), rec_jet2.v4().Pz(), rec_jet2.v4().E()); 
  }
  if(gen_lepjets.size()>0){
    gen_jet2 = gen_lepjets.at(0);
    gen_jet2_v4.SetPxPyPzE(gen_jet2.v4().Px(), gen_jet2.v4().Py(), gen_jet2.v4().Pz(), gen_jet2.v4().E()); 
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------





  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;

  if(gen_hadjets.size() > 0 && gen_lepjets.size() > 0 && rec_hadjets.size() > 0&& rec_lepjets.size() > 0 ){
    DeltaR_Rec1_Gen1->Fill(deltaR(rec_jet1, gen_jet1), weight);
    DeltaR_Rec2_Gen2->Fill(deltaR(rec_jet2, gen_jet2), weight);
    pt_rec1_gen1_beforeMatching->Fill( rec_jet1_v4.Pt() - gen_jet1_v4.Pt(), weight);
    pt_rec2_gen2_beforeMatching->Fill( rec_jet2_v4.Pt() - gen_jet2_v4.Pt(), weight);
    GenJetMass->Fill(gen_jet1_v4.M(), weight);
    MassReso->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    PtReso->Fill( (rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
    pt_rec1_gen1_afterMatching->Fill( rec_jet1_v4.Pt() - gen_jet1_v4.Pt(), weight);
    pt_rec2_gen2_afterMatching->Fill( rec_jet2_v4.Pt() - gen_jet2_v4.Pt(), weight);

    RecGenMass->Fill(gen_jet1_v4.M(), rec_jet1_v4.M(), weight);
    RecGenPT->Fill(gen_jet1_v4.Pt(), rec_jet1_v4.Pt(), weight);
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 

}


