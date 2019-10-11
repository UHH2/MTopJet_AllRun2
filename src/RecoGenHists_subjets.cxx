#include <UHH2/MTopJet/include/RecoGenHists_subjets.h>


RecoGenHists_subjets::RecoGenHists_subjets(uhh2::Context & ctx, const std::string & dirname,  const std::string & type): Hists(ctx, dirname){
  // book all histograms here
  MassReso = book<TH1F>("MassResolution", "(M^{rec}_{jet} - M^{gen}_{jet}) / M^{gen}_{jet}) ", 90, -1.5, 1.5);
  PtReso = book<TH1F>("PtResolution", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  area_all = book<TH1F>("area_all", "area all jets", 50, 0., 1.0);

  PtReso_1 = book<TH1F>("PtResolution_1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_2 = book<TH1F>("PtResolution_2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_3 = book<TH1F>("PtResolution_3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_4 = book<TH1F>("PtResolution_4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_5 = book<TH1F>("PtResolution_5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_6 = book<TH1F>("PtResolution_6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_7 = book<TH1F>("PtResolution_7", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_8 = book<TH1F>("PtResolution_8", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_9 = book<TH1F>("PtResolution_9", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_10 = book<TH1F>("PtResolution_10", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);

  PtReso_rec1 = book<TH1F>("PtResolution_rec1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec2 = book<TH1F>("PtResolution_rec2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec3 = book<TH1F>("PtResolution_rec3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec4 = book<TH1F>("PtResolution_rec4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec5 = book<TH1F>("PtResolution_rec5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec6 = book<TH1F>("PtResolution_rec6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec7 = book<TH1F>("PtResolution_rec7", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec8 = book<TH1F>("PtResolution_rec8", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec9 = book<TH1F>("PtResolution_rec9", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_rec10 = book<TH1F>("PtResolution_rec10", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);

  PtRec_1 = book<TH1F>("PtRec_1", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_2 = book<TH1F>("PtRec_2", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_3 = book<TH1F>("PtRec_3", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_4 = book<TH1F>("PtRec_4", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_5 = book<TH1F>("PtRec_5", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_6 = book<TH1F>("PtRec_6", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_7 = book<TH1F>("PtRec_7", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_8 = book<TH1F>("PtRec_8", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_9 = book<TH1F>("PtRec_9", "p^{rec}_{T, jet}", 200, 0, 1000);
  PtRec_10 = book<TH1F>("PtRec_10", "p^{rec}_{T, jet}", 200, 0, 1000);

  MatchedGen_1 = book<TH1F>("MatchedGen_1", "count", 1, 0, 2);
  MatchedGen_2 = book<TH1F>("MatchedGen_2", "count", 1, 0, 2);
  MatchedGen_3 = book<TH1F>("MatchedGen_3", "count", 1, 0, 2);
  MatchedGen_4 = book<TH1F>("MatchedGen_4", "count", 1, 0, 2);
  MatchedGen_5 = book<TH1F>("MatchedGen_5", "count", 1, 0, 2);
  MatchedGen_6 = book<TH1F>("MatchedGen_6", "count", 1, 0, 2);
  MatchedGen_7 = book<TH1F>("MatchedGen_7", "count", 1, 0, 2);
  MatchedGen_8 = book<TH1F>("MatchedGen_8", "count", 1, 0, 2);
  MatchedGen_9 = book<TH1F>("MatchedGen_9", "count", 1, 0, 2);
  MatchedGen_10 = book<TH1F>("MatchedGen_10", "count", 1, 0, 2);

  MatchedRec_1 = book<TH1F>("MatchedRec_1", "count", 1, 0, 2);
  MatchedRec_2 = book<TH1F>("MatchedRec_2", "count", 1, 0, 2);
  MatchedRec_3 = book<TH1F>("MatchedRec_3", "count", 1, 0, 2);
  MatchedRec_4 = book<TH1F>("MatchedRec_4", "count", 1, 0, 2);
  MatchedRec_5 = book<TH1F>("MatchedRec_5", "count", 1, 0, 2);
  MatchedRec_6 = book<TH1F>("MatchedRec_6", "count", 1, 0, 2);
  MatchedRec_7 = book<TH1F>("MatchedRec_7", "count", 1, 0, 2);
  MatchedRec_8 = book<TH1F>("MatchedRec_8", "count", 1, 0, 2);
  MatchedRec_9 = book<TH1F>("MatchedRec_9", "count", 1, 0, 2);
  MatchedRec_10 = book<TH1F>("MatchedRec_10", "count", 1, 0, 2);

  TotalGen_1 = book<TH1F>("TotalGen_1", "count", 1, 0, 2);
  TotalGen_2 = book<TH1F>("TotalGen_2", "count", 1, 0, 2);
  TotalGen_3 = book<TH1F>("TotalGen_3", "count", 1, 0, 2);
  TotalGen_4 = book<TH1F>("TotalGen_4", "count", 1, 0, 2);
  TotalGen_5 = book<TH1F>("TotalGen_5", "count", 1, 0, 2);
  TotalGen_6 = book<TH1F>("TotalGen_6", "count", 1, 0, 2);
  TotalGen_7 = book<TH1F>("TotalGen_7", "count", 1, 0, 2);
  TotalGen_8 = book<TH1F>("TotalGen_8", "count", 1, 0, 2);
  TotalGen_9 = book<TH1F>("TotalGen_9", "count", 1, 0, 2);
  TotalGen_10 = book<TH1F>("TotalGen_10", "count", 1, 0, 2);

  TotalRec_1 = book<TH1F>("TotalRec_1", "count", 1, 0, 2);
  TotalRec_2 = book<TH1F>("TotalRec_2", "count", 1, 0, 2);
  TotalRec_3 = book<TH1F>("TotalRec_3", "count", 1, 0, 2);
  TotalRec_4 = book<TH1F>("TotalRec_4", "count", 1, 0, 2);
  TotalRec_5 = book<TH1F>("TotalRec_5", "count", 1, 0, 2);
  TotalRec_6 = book<TH1F>("TotalRec_6", "count", 1, 0, 2);
  TotalRec_7 = book<TH1F>("TotalRec_7", "count", 1, 0, 2);
  TotalRec_8 = book<TH1F>("TotalRec_8", "count", 1, 0, 2);
  TotalRec_9 = book<TH1F>("TotalRec_9", "count", 1, 0, 2);
  TotalRec_10 = book<TH1F>("TotalRec_10", "count", 1, 0, 2);

  min_mass_Wjet_rec = book<TH1F>("min_mass_Wjet_rec", "min M_{ij}", 60, 60, 120);
  min_mass_Wjet_gen = book<TH1F>("min_mass_Wjet_gen", "min M_{ij}", 60, 60, 120);
  WMassReso = book<TH1F>("WMassResolution", "(M^{rec}_{W} - M^{gen}_{W}) / M^{gen}_{W}) ", 90, -1.5, 1.5);

  // handle for clustered jets
  if(type == "jec"){
    h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined");
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  }
  else if(type == "raw"){
    h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_noJEC");
  }
  else if(type == "cor"){
    h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
  }
  h_genjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

}



void RecoGenHists_subjets::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> rec = event.get(h_recjets);
  std::vector<TopJet> had = event.get(h_hadjets);
  std::vector<GenTopJet> gen = event.get(h_genjets);
  std::vector<Jet> rec_sub;
  std::vector<Particle> gen_sub;

  // only look in had subjets, to not correct for lepton
  int had_nr = 0;
  if(deltaR(had.at(0), rec.at(0)) > deltaR(had.at(0), rec.at(1))) had_nr = 1;

  for(unsigned int j=0; j<rec.at(had_nr).subjets().size(); j++){
    rec_sub.push_back(rec.at(had_nr).subjets().at(j));
  }

  // select which gen fatjet is closer to had jet and get those gen subjets
  int gen_index = 0;
  if(deltaR(had[0], gen[0]) > deltaR(had[0], gen[1]) ){
    gen_index = 1;
  }

  // get all subjets on gen level
  for(unsigned int j=0; j<gen.at(gen_index).subjets().size(); j++){
    gen_sub.push_back(gen.at(gen_index).subjets().at(j));
  }


  if(rec_sub.size() < 1) return;
  if(gen_sub.size() < 1) return;

  /* ******************************************************************************
  matching between gen and reco jets:
  - a rec jet is called isolated if the next jet is not within 2R. An isolated jet should be spherical and more simelar to an ak4 jet
  - for each reco jet, calc distance to all other gen jets
  - gen jet with lowest distance is a match if distance is < 0.2
  - then calculate resolution with reco jet and matched gen jet
  - to do: account for double counting
  ********************************************************************************* */

  double dR;
  double dR_temp;
  int nearest_j;
  double gen_pt, rec_pt, R;

  // do matching from rec to gen
  for(unsigned int i=0; i<rec_sub.size(); i++){
    dR = 1000;
    nearest_j = 100;
    area_all->Fill(rec_sub.at(i).jetArea(), weight);
    for(unsigned int j=0; j<gen_sub.size(); j++){
      dR_temp = uhh2::deltaR(rec_sub.at(i), gen_sub.at(j));
      if(dR_temp < dR){
        dR = dR_temp;
        nearest_j = j;
      }
    }
    if(nearest_j == 100) return;
    rec_pt = rec_sub.at(i).v4().Pt();
    if(nearest_j != 100 && dR <= 0.2){
      // count all matched subjets
      if(rec_pt <= 50) MatchedRec_1->Fill(1, weight );
      if(rec_pt > 50 && rec_pt <= 80) MatchedRec_2->Fill(1, weight );
      if(rec_pt > 80 && rec_pt <= 120) MatchedRec_3->Fill(1, weight );
      if(rec_pt > 120 && rec_pt <= 170) MatchedRec_4->Fill(1, weight );
      if(rec_pt > 170 && rec_pt <= 220) MatchedRec_5->Fill(1, weight );
      if(rec_pt > 220 && rec_pt <= 270) MatchedRec_6->Fill(1, weight );
      if(rec_pt > 270 && rec_pt <= 320) MatchedRec_7->Fill(1, weight );
      if(rec_pt > 320 && rec_pt <= 370) MatchedRec_8->Fill(1, weight );
      if(rec_pt > 370 && rec_pt <= 420) MatchedRec_9->Fill(1, weight );
      if(rec_pt > 420) MatchedRec_10->Fill(1, weight );
      ////
    }

  }

  // do matching from gen to rec
  for(unsigned int i=0; i<gen_sub.size(); i++){
    dR = 1000;
    nearest_j = 100;
    for(unsigned int j=0; j<rec_sub.size(); j++){
      dR_temp = uhh2::deltaR(rec_sub.at(j), gen_sub.at(i));
      if(dR_temp < dR){
        dR = dR_temp;
        nearest_j = j;
      }
    }

    if(nearest_j == 100) return;
    gen_pt = gen_sub.at(i).v4().Pt();
    rec_pt = rec_sub.at(nearest_j).v4().Pt();
    R = (rec_pt - gen_pt) / gen_pt;

    if(nearest_j != 100 && dR <= 0.2){
      // if( !RecoGenHists_subjets::findMatch(event, rec_sub.at(nearest_j)) ) continue;
      PtReso->Fill( R, weight );
      MassReso->Fill( (rec_sub.at(nearest_j).v4().M() - gen_sub.at(i).v4().M())/gen_sub.at(i).v4().M() , weight );
      if(gen_pt <= 50) PtReso_1->Fill( R, weight );
      if(gen_pt > 50 && gen_pt <= 80) PtReso_2->Fill( R, weight );
      if(gen_pt > 80 && gen_pt <= 120) PtReso_3->Fill( R, weight );
      if(gen_pt > 120 && gen_pt <= 170) PtReso_4->Fill( R, weight );
      if(gen_pt > 170 && gen_pt <= 220) PtReso_5->Fill( R, weight );
      if(gen_pt > 220 && gen_pt <= 270) PtReso_6->Fill( R, weight );
      if(gen_pt > 270 && gen_pt <= 320) PtReso_7->Fill( R, weight );
      if(gen_pt > 320 && gen_pt <= 370) PtReso_8->Fill( R, weight );
      if(gen_pt > 370 && gen_pt <= 420) PtReso_9->Fill( R, weight );
      if(gen_pt > 420) PtReso_10->Fill( R, weight );

      if(gen_pt <= 50) PtRec_1->Fill( rec_pt, weight );
      if(gen_pt > 50 && gen_pt <= 80) PtRec_2->Fill( rec_pt, weight );
      if(gen_pt > 80 && gen_pt <= 120) PtRec_3->Fill( rec_pt, weight );
      if(gen_pt > 120 && gen_pt <= 170) PtRec_4->Fill( rec_pt, weight );
      if(gen_pt > 170 && gen_pt <= 220) PtRec_5->Fill( rec_pt, weight );
      if(gen_pt > 220 && gen_pt <= 270) PtRec_6->Fill( rec_pt, weight );
      if(gen_pt > 270 && gen_pt <= 320) PtRec_7->Fill( rec_pt, weight );
      if(gen_pt > 320 && gen_pt <= 370) PtRec_8->Fill( rec_pt, weight );
      if(gen_pt > 370 && gen_pt <= 420) PtRec_9->Fill( rec_pt, weight );
      if(gen_pt > 420) PtRec_10->Fill( rec_pt, weight );

      if(rec_pt <= 50) PtReso_rec1->Fill( R, weight );
      if(rec_pt > 50 && rec_pt <= 80) PtReso_rec2->Fill( R, weight );
      if(rec_pt > 80 && rec_pt <= 120) PtReso_rec3->Fill( R, weight );
      if(rec_pt > 120 && rec_pt <= 170) PtReso_rec4->Fill( R, weight );
      if(rec_pt > 170 && rec_pt <= 220) PtReso_rec5->Fill( R, weight );
      if(rec_pt > 220 && rec_pt <= 270) PtReso_rec6->Fill( R, weight );
      if(rec_pt > 270 && rec_pt <= 320) PtReso_rec7->Fill( R, weight );
      if(rec_pt > 320 && rec_pt <= 370) PtReso_rec8->Fill( R, weight );
      if(rec_pt > 370 && rec_pt <= 420) PtReso_rec9->Fill( R, weight );
      if(rec_pt > 420) PtReso_rec10->Fill( R, weight );

      // count all matched subjets
      if(gen_pt <= 50) MatchedGen_1->Fill(1, weight );
      if(gen_pt > 50 && gen_pt <= 80) MatchedGen_2->Fill(1, weight );
      if(gen_pt > 80 && gen_pt <= 120) MatchedGen_3->Fill(1, weight );
      if(gen_pt > 120 && gen_pt <= 170) MatchedGen_4->Fill(1, weight );
      if(gen_pt > 170 && gen_pt <= 220) MatchedGen_5->Fill(1, weight );
      if(gen_pt > 220 && gen_pt <= 270) MatchedGen_6->Fill(1, weight );
      if(gen_pt > 270 && gen_pt <= 320) MatchedGen_7->Fill(1, weight );
      if(gen_pt > 320 && gen_pt <= 370) MatchedGen_8->Fill(1, weight );
      if(gen_pt > 370 && gen_pt <= 420) MatchedGen_9->Fill(1, weight );
      if(gen_pt > 420) MatchedGen_10->Fill(1, weight );
      ////
    }

  }

  //// Here: Count total subjets in bins
  for(unsigned int j=0; j<gen_sub.size(); j++){
    double gen_pt = gen_sub[j].pt();
    if(gen_pt <= 50) TotalGen_1->Fill(1, weight );
    if(gen_pt > 50 && gen_pt <= 80) TotalGen_2->Fill(1, weight );
    if(gen_pt > 80 && gen_pt <= 120) TotalGen_3->Fill(1, weight );
    if(gen_pt > 120 && gen_pt <= 170) TotalGen_4->Fill(1, weight );
    if(gen_pt > 170 && gen_pt <= 220) TotalGen_5->Fill(1, weight );
    if(gen_pt > 220 && gen_pt <= 270) TotalGen_6->Fill(1, weight );
    if(gen_pt > 270 && gen_pt <= 320) TotalGen_7->Fill(1, weight );
    if(gen_pt > 320 && gen_pt <= 370) TotalGen_8->Fill(1, weight );
    if(gen_pt > 370 && gen_pt <= 420) TotalGen_9->Fill(1, weight );
    if(gen_pt > 420) TotalGen_10->Fill(1, weight );
  }

  for(unsigned int i=0; i<rec_sub.size(); i++){
    double rec_pt = rec_sub[i].pt();
    if(rec_pt <= 50) TotalRec_1->Fill(1, weight );
    if(rec_pt > 50 && rec_pt <= 80) TotalRec_2->Fill(1, weight );
    if(rec_pt > 80 && rec_pt <= 120) TotalRec_3->Fill(1, weight );
    if(rec_pt > 120 && rec_pt <= 170) TotalRec_4->Fill(1, weight );
    if(rec_pt > 170 && rec_pt <= 220) TotalRec_5->Fill(1, weight );
    if(rec_pt > 220 && rec_pt <= 270) TotalRec_6->Fill(1, weight );
    if(rec_pt > 270 && rec_pt <= 320) TotalRec_7->Fill(1, weight );
    if(rec_pt > 320 && rec_pt <= 370) TotalRec_8->Fill(1, weight );
    if(rec_pt > 370 && rec_pt <= 420) TotalRec_9->Fill(1, weight );
    if(rec_pt > 420) TotalRec_10->Fill(1, weight );
  }
  //---------------------------------------------------------------------------------------
  //--------------------------------- add subjets to reconstruct W ------------------------
  //---------------------------------------------------------------------------------------
  TLorentzVector Wjet_rec;
  double M_min_rec = 1000;
  double M_temp = 0;
  double px, py, pz, E;
  for(unsigned int i=0; i<rec.at(0).subjets().size(); i++){
    for(unsigned int j=0; j<rec.at(0).subjets().size(); j++){
      if(i != j){
        px = rec.at(0).subjets().at(i).v4().Px() + rec.at(0).subjets().at(j).v4().Px();
        py = rec.at(0).subjets().at(i).v4().Py() + rec.at(0).subjets().at(j).v4().Py();
        pz = rec.at(0).subjets().at(i).v4().Pz() + rec.at(0).subjets().at(j).v4().Pz();
        E = rec.at(0).subjets().at(i).v4().E() + rec.at(0).subjets().at(j).v4().E();
        Wjet_rec.SetPxPyPzE(px, py, pz, E);
        M_temp = Wjet_rec.M();
        if(M_temp < M_min_rec) M_min_rec = M_temp;
      }
    }
  }
  TLorentzVector Wjet_gen;
  double M_min_gen = 1000;
  for(unsigned int i=0; i<gen.at(0).subjets().size(); i++){
    for(unsigned int j=0; j<gen.at(0).subjets().size(); j++){
      if(i != j){
        px = gen.at(0).subjets().at(i).v4().Px() + gen.at(0).subjets().at(j).v4().Px();
        py = gen.at(0).subjets().at(i).v4().Py() + gen.at(0).subjets().at(j).v4().Py();
        pz = gen.at(0).subjets().at(i).v4().Pz() + gen.at(0).subjets().at(j).v4().Pz();
        E = gen.at(0).subjets().at(i).v4().E() + gen.at(0).subjets().at(j).v4().E();
        Wjet_gen.SetPxPyPzE(px, py, pz, E);
        M_temp = Wjet_gen.M();
        if(M_temp < M_min_gen) M_min_gen = M_temp;
      }
    }
  }
  min_mass_Wjet_rec->Fill(M_min_rec, weight);
  min_mass_Wjet_gen->Fill(M_min_gen, weight);
  WMassReso->Fill( (M_min_rec - M_min_gen) / M_min_gen, weight);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


}

// only returns true if there is exactly one gen particle that is matched to subjet
bool RecoGenHists_subjets::findMatch(const Event & event, Jet jet){
  bool match = false;
  int nMatches = 0;
  std::vector<GenParticle>* genparts = event.genparticles;
  for (unsigned int i=0; i < genparts->size(); ++i){
    GenParticle p = genparts->at(i);
    if(i==0 || i==1) continue;        // skip initial particles
    int id = abs(p.pdgId());
    if(id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 9){
      double dR = deltaR(jet, p);
      if(dR < 0.2){
        nMatches++;
      }
    }
  }
  if(nMatches == 1) match = true;
  return match;
}
