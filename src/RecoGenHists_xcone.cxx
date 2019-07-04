#include "UHH2/MTopJet/include/RecoGenHists_xcone.h"


RecoGenHists_xcone::RecoGenHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type): Hists(ctx, dirname){
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

  PtReso_mass1 = book<TH1F>("PtResolution_mass1", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_mass2 = book<TH1F>("PtResolution_mass2", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_mass3 = book<TH1F>("PtResolution_mass3", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_mass4 = book<TH1F>("PtResolution_mass4", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_mass5 = book<TH1F>("PtResolution_mass5", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_mass6 = book<TH1F>("PtResolution_mass6", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  MassReso_mass1 = book<TH1F>("MassResolution_mass1", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_mass2 = book<TH1F>("MassResolution_mass2", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_mass3 = book<TH1F>("MassResolution_mass3", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_mass4 = book<TH1F>("MassResolution_mass4", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_mass5 = book<TH1F>("MassResolution_mass5", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_mass6 = book<TH1F>("MassResolution_mass6", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);

  PtReso_pt1 = book<TH1F>("PtResolution_pt1", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt2 = book<TH1F>("PtResolution_pt2", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt3 = book<TH1F>("PtResolution_pt3", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt4 = book<TH1F>("PtResolution_pt4", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt5 = book<TH1F>("PtResolution_pt5", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt6 = book<TH1F>("PtResolution_pt6", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt1_rec = book<TH1F>("PtResolution_pt1_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt2_rec = book<TH1F>("PtResolution_pt2_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt3_rec = book<TH1F>("PtResolution_pt3_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt4_rec = book<TH1F>("PtResolution_pt4_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt5_rec = book<TH1F>("PtResolution_pt5_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);
  PtReso_pt6_rec = book<TH1F>("PtResolution_pt6_rec", "(p^{rec}_{T, jet1} - p^{gen}_{T, jet1}) / p^{gen}_{T, jet1}) ", 90, -1.5, 1.5);

  MassReso_pt1 = book<TH1F>("MassResolution_pt1", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt2 = book<TH1F>("MassResolution_pt2", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt3 = book<TH1F>("MassResolution_pt3", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt4 = book<TH1F>("MassResolution_pt4", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt5 = book<TH1F>("MassResolution_pt5", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt6 = book<TH1F>("MassResolution_pt6", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 180, -1.5, 1.5);
  MassReso_pt1_rec = book<TH1F>("MassResolution_pt1_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_pt2_rec = book<TH1F>("MassResolution_pt2_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_pt3_rec = book<TH1F>("MassResolution_pt3_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_pt4_rec = book<TH1F>("MassResolution_pt4_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_pt5_rec = book<TH1F>("MassResolution_pt5_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);
  MassReso_pt6_rec = book<TH1F>("MassResolution_pt6_rec", "(M^{rec}_{jet1} - M^{gen}_{jet1}) / M^{gen}_{jet1}) ", 90, -1.5, 1.5);

  // handle for clustered jets
  if(type == "jec"){
    h_recjets_had=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined");
    h_recjets_lep=ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined");
  }
  else if(type == "raw"){
    h_recjets_had=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");
    h_recjets_lep=ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined_noJEC");
  }
  else if(type == "cor"){
    h_recjets_had=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
    h_recjets_lep=ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined_Corrected");
  }

  h_genjets_had=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");
  h_genjets_lep=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_lep_Combined");

}



void RecoGenHists_xcone::fill(const Event & event){


  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
 // define all objects needed
  std::vector<TopJet> rec_hadjets = event.get(h_recjets_had);
  std::vector<TopJet> rec_lepjets = event.get(h_recjets_lep);
  std::vector<GenTopJet> gen_hadjets = event.get(h_genjets_had);
  std::vector<GenTopJet> gen_lepjets = event.get(h_genjets_lep);


  TLorentzVector rec_jet1_v4, rec_jet2_v4, gen_jet1_v4, gen_jet2_v4;
  TopJet rec_jet1, rec_jet2;
  GenTopJet gen_jet1, gen_jet2;
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

    // resolution binned in Mass
    if(gen_jet1_v4.M() < 150) {
      PtReso_mass1->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass1->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.M() >= 150 && gen_jet1_v4.M() < 175){
      PtReso_mass2->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass2->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.M() >= 175 && gen_jet1_v4.M() < 200){
      PtReso_mass3->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass3->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.M() >= 200 && gen_jet1_v4.M() < 250){
      PtReso_mass4->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass4->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.M() >= 250 && gen_jet1_v4.M() < 300){
      PtReso_mass5->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass5->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.M() >= 300 && gen_jet1_v4.M() < 500){
      PtReso_mass6->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_mass6->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    // resolution binned in PT gen
    if(gen_jet1_v4.Pt() < 450) {
      PtReso_pt1->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt1->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.Pt() >= 450 && gen_jet1_v4.Pt() < 500){
      PtReso_pt2->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt2->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.Pt() >= 500 && gen_jet1_v4.Pt() < 600){
      PtReso_pt3->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt3->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.Pt() >= 600 && gen_jet1_v4.Pt() < 700){
      PtReso_pt4->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt4->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.Pt() >= 700 && gen_jet1_v4.Pt() < 800){
      PtReso_pt5->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt5->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(gen_jet1_v4.Pt() >= 800){
      PtReso_pt6->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt6->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    // resolution binned in PT rec
    if(rec_jet1_v4.Pt() < 450) {
      PtReso_pt1_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt1_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(rec_jet1_v4.Pt() >= 450 && rec_jet1_v4.Pt() < 500){
      PtReso_pt2_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt2_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(rec_jet1_v4.Pt() >= 500 && rec_jet1_v4.Pt() < 600){
      PtReso_pt3_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt3_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(rec_jet1_v4.Pt() >= 600 && rec_jet1_v4.Pt() < 700){
      PtReso_pt4_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt4_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(rec_jet1_v4.Pt() >= 700 && rec_jet1_v4.Pt() < 800){
      PtReso_pt5_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt5_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }
    if(rec_jet1_v4.Pt() >= 800){
      PtReso_pt6_rec->Fill((rec_jet1_v4.Pt() - gen_jet1_v4.Pt())/gen_jet1_v4.Pt(), weight );
      MassReso_pt6_rec->Fill( (rec_jet1_v4.M() - gen_jet1_v4.M())/gen_jet1_v4.M(), weight );
    }

  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


}
