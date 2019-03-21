#include "UHH2/MTopJet/include/RecoGenHists_allHad_flavor.h"


RecoGenHists_allHad_flavor::RecoGenHists_allHad_flavor(uhh2::Context & ctx, const std::string & dirname,  const std::string & type, double ptreccut_, double ptgencut_, TString flavor_): Hists(ctx, dirname){
  // book all histograms here
  pt_rec = book<TH1F>("pt_rec", "p_{T}", 100, 0.0, 500);
  bIndex = book<TH1F>("bIndex", "subjet index matched to b", 4, -1.5, 2.5);
  nsub_before_erase = book<TH1F>("nsub_before_erase", "number of subjets before erase", 5, -0.5, 4.5);
  nsub_after_erase = book<TH1F>("nsub_after_erase", "number of subjets after erase", 5, -0.5, 4.5);

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

  PtReso_area1 = book<TH1F>("PtResolution_area1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area2 = book<TH1F>("PtResolution_area2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area3 = book<TH1F>("PtResolution_area3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area4 = book<TH1F>("PtResolution_area4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area5 = book<TH1F>("PtResolution_area5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area6 = book<TH1F>("PtResolution_area6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area7 = book<TH1F>("PtResolution_area7", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_area8 = book<TH1F>("PtResolution_area8", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);


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

  min_mass_Wjet_rec = book<TH1F>("min_mass_Wjet_rec", "min M_{ij}", 60, 60, 120);
  min_mass_Wjet_gen = book<TH1F>("min_mass_Wjet_gen", "min M_{ij}", 60, 60, 120);
  WMassReso = book<TH1F>("WMassResolution", "(M^{rec}_{W} - M^{gen}_{W}) / M^{gen}_{W}) ", 90, -1.5, 1.5);

  // handle for clustered jets
  if(type == "jec"){
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  }
  else if(type == "raw"){
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_noJEC");
  }
  else if(type == "cor"){
    h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
  }
  h_genjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  // handle for ttbargen
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  // set min ptgen and min ptrec
  ptreccut = ptreccut_;
  ptgencut = ptgencut_;

  if(flavor_ == "b"){
    fill_b = true;
    fill_light = false;
  }
  else if(flavor_ == "light"){
    fill_b = false;
    fill_light = true;
  }
  else{
    fill_b = true;
    fill_light = true;
  }

}



void RecoGenHists_allHad_flavor::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> recjets = event.get(h_recjets);
  std::vector<GenTopJet> genjets = event.get(h_genjets);

  // do not continus if a jet collection is empty
  if(recjets.size() < 2) return;
  for(unsigned int i=0; i<recjets.size(); i++){
    if(recjets.at(i).subjets().size() < 3) return;
  }
  if(genjets.size() < 2) return;
  for(unsigned int i=0; i<genjets.size(); i++){
    if(genjets.at(i).subjets().size() < 3) return;
  }


  // first match fatjets in gen and reco (find index of genjet to reco jet at i=0)
  double dRfat = 100;
  double dRfat_temp = 100;
  int i1 = 0;
  int i2 = 1;
  for(unsigned int i=0; i<genjets.size(); i++){
    dRfat_temp = deltaR(recjets.at(0), genjets.at(i));
    if(dRfat_temp < dRfat){
      dRfat = dRfat_temp;
      i1 = i;
    }
  }
  if(i1 == 0) i2=1;
  else        i2=0;

  bool use_jet1 = false;
  bool use_jet2 = false;
  if( deltaR(recjets.at(0), genjets.at(i1)) < 1.2 ) use_jet1 = true;
  if( deltaR(recjets.at(1), genjets.at(i2)) < 1.2 ) use_jet2 = true;

  std::vector<Jet> rec_sub1 = recjets.at(0).subjets();
  std::vector<Jet> rec_sub2 = recjets.at(1).subjets();
  std::vector<Particle> gen_sub1 = genjets.at(i1).subjets();
  std::vector<Particle> gen_sub2 = genjets.at(i2).subjets();

  sort_by_pt<Jet> (rec_sub1);
  sort_by_pt<Jet> (rec_sub2);

  double ptsubmax = 0;
  for(auto subjet: rec_sub1){
    if(subjet.v4().Pt() < ptsubmax) use_jet1 = false;
  }
  for(auto subjet: rec_sub2){
    if(subjet.v4().Pt() < ptsubmax) use_jet2 = false;
  }
  if(recjets.at(0).v4().Pt() < ptreccut || genjets.at(i1).v4().Pt() < ptgencut) use_jet1 = false;
  if(recjets.at(1).v4().Pt() < ptreccut || genjets.at(i2).v4().Pt() < ptgencut) use_jet2 = false;

  int b_index1 = IndexBJet(event, recjets.at(0));
  int b_index2 = IndexBJet(event, recjets.at(1));

  if(use_jet1) bIndex->Fill(b_index1, weight);
  if(use_jet2) bIndex->Fill(b_index2, weight);

  if(b_index1 == -1) use_jet1 = false;
  if(b_index2 == -1) use_jet2 = false;

  if(!use_jet1 && !use_jet2) return;

  if(use_jet1) nsub_before_erase->Fill(rec_sub1.size(), weight);
  if(use_jet2) nsub_before_erase->Fill(rec_sub2.size(), weight);
  // now delete subjet(s) to only maintain 1 b-jet or 2 light jets
  vector<Jet> sub_new1;
  for(unsigned int i=0; i<rec_sub1.size(); i++){
    if(fill_b     && (int)i == b_index1) sub_new1.push_back(rec_sub1[i]);
    if(fill_light && (int)i != b_index1) sub_new1.push_back(rec_sub1[i]);
  }
  rec_sub1 = sub_new1;

  vector<Jet> sub_new2;
  for(unsigned int i=0; i<rec_sub2.size(); i++){
    if(fill_b     && (int)i == b_index2) sub_new2.push_back(rec_sub2[i]);
    if(fill_light && (int)i != b_index2) sub_new2.push_back(rec_sub2[i]);
  }
  rec_sub2 = sub_new2;
  ////

  if(use_jet1) nsub_after_erase->Fill(rec_sub1.size(), weight);
  if(use_jet2) nsub_after_erase->Fill(rec_sub2.size(), weight);
  /* ******************************************************************************
  matching between gen and reco jets:
  - for each reco jet, calc distance to all other gen jets
  - gen jet with lowest distance is a match if distance is < 0.2
  - then calculate resolution with reco jet and matched gen jet
  - to do: account for double counting
  ********************************************************************************* */

  double dR;
  double dR_temp;
  int nearest_j;
  double gen_pt;
  double rec_pt;
  double R;
  double area;
  // do matching
  if(use_jet1){
    for(unsigned int i=0; i<rec_sub1.size(); i++){
      dR = 1000;
      nearest_j = 100;
      for(unsigned int j=0; j<gen_sub1.size(); j++){
        dR_temp = uhh2::deltaR(rec_sub1.at(i), gen_sub1.at(j));
        if(dR_temp < dR){
          dR = dR_temp;
          nearest_j = j;
        }
      }
      gen_pt=gen_sub1.at(nearest_j).v4().Pt();
      rec_pt=rec_sub1.at(i).v4().Pt();
      area = rec_sub1.at(i).jetArea();
      R = (rec_pt - gen_pt) / gen_pt;
      pt_rec->Fill(rec_pt, weight);
      if(nearest_j != 100 && dR <= 0.2){
        PtReso->Fill( R, weight );
        MassReso->Fill( (rec_sub1.at(i).v4().M() - gen_sub1.at(nearest_j).v4().M())/gen_sub1.at(nearest_j).v4().M() , weight );
        // resolution vs genpt
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
        // resolution vs area
        if(area <= 0.30) PtReso_area1->Fill( R, weight );
        if(area > 0.30 && area <= 0.35) PtReso_area2->Fill( R, weight );
        if(area > 0.35 && area <= 0.38) PtReso_area3->Fill( R, weight );
        if(area > 0.38 && area <= 0.41) PtReso_area4->Fill( R, weight );
        if(area > 0.41 && area <= 0.44) PtReso_area5->Fill( R, weight );
        if(area > 0.44 && area <= 0.47) PtReso_area6->Fill( R, weight );
        if(area > 0.47 && area <= 0.5) PtReso_area7->Fill( R, weight );
        if(area > 0.5) PtReso_area8->Fill( R, weight );

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
      }
    }
  }
  if(use_jet2){
    for(unsigned int i=0; i<rec_sub2.size(); i++){
      dR = 1000;
      nearest_j = 100;
      for(unsigned int j=0; j<gen_sub2.size(); j++){
        dR_temp = uhh2::deltaR(rec_sub2.at(i), gen_sub2.at(j));
        if(dR_temp < dR){
          dR = dR_temp;
          nearest_j = j;
        }
      }
      gen_pt=gen_sub2.at(nearest_j).v4().Pt();
      rec_pt=rec_sub2.at(i).v4().Pt();
      area = rec_sub2.at(i).jetArea();
      R = (rec_pt - gen_pt) / gen_pt;
      pt_rec->Fill(rec_pt, weight);
      if(nearest_j != 100 && dR <= 0.2){
        PtReso->Fill( R, weight );
        MassReso->Fill( (rec_sub2.at(i).v4().M() - gen_sub2.at(nearest_j).v4().M())/gen_sub2.at(nearest_j).v4().M() , weight );
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
        // resolution vs area
        if(area <= 0.30) PtReso_area1->Fill( R, weight );
        if(area > 0.30 && area <= 0.35) PtReso_area2->Fill( R, weight );
        if(area > 0.35 && area <= 0.38) PtReso_area3->Fill( R, weight );
        if(area > 0.38 && area <= 0.41) PtReso_area4->Fill( R, weight );
        if(area > 0.41 && area <= 0.44) PtReso_area5->Fill( R, weight );
        if(area > 0.44 && area <= 0.47) PtReso_area6->Fill( R, weight );
        if(area > 0.47 && area <= 0.5) PtReso_area7->Fill( R, weight );
        if(area > 0.5) PtReso_area8->Fill( R, weight );

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
      }
    }
  }
}

int RecoGenHists_allHad_flavor::IndexBJet(const uhh2::Event & event, TopJet jet){
  int index = -1;
  vector<Jet> subjets = jet.subjets();
  const auto & ttbargen = event.get(h_ttbargen);
  GenParticle b1 = ttbargen.bTop();
  GenParticle b2 = ttbargen.bAntitop();
  GenParticle bot;
  // first match b to TopJet
  if(deltaR(b1, jet) < deltaR(b2, jet)) bot = b1;
  else bot = b2;

  // now find subjet closest to bottom quark
  double dR = 100;
  for(unsigned int i=0; i<subjets.size(); i++){
    if(deltaR(bot, subjets[i]) < dR){
      dR = deltaR(bot, subjets[i]);
      if(dR < 0.4) index = i;
    }
  }
  return index;
}
