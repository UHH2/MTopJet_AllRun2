#include "UHH2/MTopJet/include/RecoGenHists_subjets.h"


RecoGenHists_subjets::RecoGenHists_subjets(uhh2::Context & ctx, const std::string & dirname, bool use_JEC): Hists(ctx, dirname){
  // book all histograms here
  MassReso = book<TH1F>("MassResolution", "(M^{rec}_{jet} - M^{gen}_{jet}) / M^{gen}_{jet}) ", 90, -1.5, 1.5);
  PtReso = book<TH1F>("PtResolution", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  MassReso_iso = book<TH1F>("MassResolution_iso", "(M^{rec}_{jet} - M^{gen}_{jet}) / M^{gen}_{jet}) ", 90, -1.5, 1.5);
  PtReso_iso = book<TH1F>("PtResolution_iso", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  area_iso = book<TH1F>("area_iso", "area isolated jets", 50, 0., 1.0);
  area_all = book<TH1F>("area_all", "area all jets", 50, 0., 1.0);

  PtReso_1 = book<TH1F>("PtResolution_1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_2 = book<TH1F>("PtResolution_2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_3 = book<TH1F>("PtResolution_3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_4 = book<TH1F>("PtResolution_4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_5 = book<TH1F>("PtResolution_5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_6 = book<TH1F>("PtResolution_6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_1 = book<TH1F>("PtResolution_iso_1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_2 = book<TH1F>("PtResolution_iso_2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_3 = book<TH1F>("PtResolution_iso_3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_4 = book<TH1F>("PtResolution_iso_4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_5 = book<TH1F>("PtResolution_iso_5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_iso_6 = book<TH1F>("PtResolution_iso_6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);

  min_mass_Wjet_rec = book<TH1F>("min_mass_Wjet_rec", "min M_{ij}", 60, 60, 120);
  min_mass_Wjet_gen = book<TH1F>("min_mass_Wjet_gen", "min M_{ij}", 60, 60, 120);
  WMassReso = book<TH1F>("WMassResolution", "(M^{rec}_{W} - M^{gen}_{W}) / M^{gen}_{W}) ", 90, -1.5, 1.5);

  // handle for clustered jets
  if(use_JEC) h_recjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets");
  else h_recjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets_noJEC");

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
  std::vector<GenTopJet> gen = event.get(h_genjets);
  std::vector<Jet> rec_sub;
  std::vector<Particle> gen_sub;

  for(unsigned int i=0; i<rec.size(); i++){
    for(unsigned int j=0; j<rec.at(i).subjets().size(); j++){
      rec_sub.push_back(rec.at(i).subjets().at(j));
    }
  }
  for(unsigned int i=0; i<gen.size(); i++){
    for(unsigned int j=0; j<gen.at(i).subjets().size(); j++){
      gen_sub.push_back(gen.at(i).subjets().at(j));
    }
  }

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
  bool isolated = false;
  double gen_pt;

  // do matching
  for(unsigned int i=0; i<rec_sub.size(); i++){
    dR = 1000;
    nearest_j = 100;
    for(unsigned int i2=0; i2<rec_sub.size(); i2++){
      if(i != i2) isolated = (uhh2::deltaR(rec_sub.at(i), rec_sub.at(i2)) > 0.8);
    }
    area_all->Fill(rec_sub.at(i).jetArea(), weight);
    if(isolated) area_iso->Fill(rec_sub.at(i).jetArea(), weight);
    for(unsigned int j=0; j<gen_sub.size(); j++){
      dR_temp = uhh2::deltaR(rec_sub.at(i), gen_sub.at(j));
      if(dR_temp < dR){
	dR = dR_temp;
	nearest_j = j;
      }
    }
    gen_pt=gen_sub.at(nearest_j).v4().Pt();
    if(nearest_j != 100 && dR <= 0.2){
      PtReso->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      MassReso->Fill( (rec_sub.at(i).v4().M() - gen_sub.at(nearest_j).v4().M())/gen_sub.at(nearest_j).v4().M() , weight );
      if(gen_pt <= 50) PtReso_1->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(gen_pt > 50 && gen_pt <= 100) PtReso_2->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(gen_pt > 100 && gen_pt <= 200) PtReso_3->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(gen_pt > 200 && gen_pt <= 300) PtReso_4->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(gen_pt > 300 && gen_pt <= 400) PtReso_5->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(gen_pt > 400 && gen_pt <= 500) PtReso_6->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      if(isolated){
	PtReso_iso->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	MassReso_iso->Fill( (rec_sub.at(i).v4().M() - gen_sub.at(nearest_j).v4().M())/gen_sub.at(nearest_j).v4().M() , weight );
	if(gen_pt <= 50) PtReso_iso_1->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	if(gen_pt > 50 && gen_pt <= 100) PtReso_iso_2->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	if(gen_pt > 100 && gen_pt <= 200) PtReso_iso_3->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	if(gen_pt > 200 && gen_pt <= 300) PtReso_iso_4->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	if(gen_pt > 300 && gen_pt <= 400) PtReso_iso_5->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
	if(gen_pt > 400 && gen_pt <= 500) PtReso_iso_6->Fill( (rec_sub.at(i).v4().Pt() - gen_sub.at(nearest_j).v4().Pt())/gen_sub.at(nearest_j).v4().Pt(), weight );
      }

    }

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


