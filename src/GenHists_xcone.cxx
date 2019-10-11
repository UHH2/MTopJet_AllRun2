#include <UHH2/MTopJet/include/GenHists_xcone.h>


GenHists_xcone::GenHists_xcone(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Mass_HadJet33 = book<TH1F>("Mass_HadJet33", "m_{jet}", 50, 0, 500);
  Mass_HadJet33_B = book<TH1F>("Mass_HadJet33_B", "m_{jet}", 100, 0, 500);
  Mass_HadJet33_C = book<TH1F>("Mass_HadJet33_C", "m_{jet}", 500, 0, 500);
  Mass_HadJet33_rebin = book<TH1F>("Mass_HadJet33_rebin", "m_{jet}", 25, 0, 500);

  Float_t xbins[] = {0, 127, 142, 157, 172, 187, 203, 218, 240, 270, 300, 500};
  Mass_HadJet33_unfold = book<TH1F>("Mass_HadJet33_unfold", "M_{jet}", 11, xbins);

  Mass_LepJet33 = book<TH1F>("Mass_LepJet33", "m_{jet+lepton}", 50, 0, 500);
  MassDiff33 = book<TH1F>("MassDiff33", "m_{jet1} - m_{jet2 + lepton}", 40, -200, 200);

  PT_HadJet33 = book<TH1F>("PT_HadJet33", "p_{T}", 50, 0, 1000);
  PT_LepJet33 = book<TH1F>("PT_LepJet33", "p_{T}", 50, 0, 1000);

  number_HadJet33 = book<TH1F>("number_HadJet33", "number", 10, 0, 10);
  number_LepJet33 = book<TH1F>("number_LepJet33", "number", 10, 0, 10);

  SoftDropMass = book<TH1F>("SoftDropMass", "M_{jet1}^{SD}", 50, 0, 500);

  min_mass_Wjet = book<TH1F>("min_mass_Wjet", "min M_{ij}", 50, 50, 150);
  min_mass_Wjet_ptlow = book<TH1F>("min_mass_Wjet_ptlow", "min M_{ij}", 50, 50, 150);
  min_mass_Wjet_pthigh = book<TH1F>("min_mass_Wjet_pthigh", "min M_{ij}", 50, 50, 150);
  all_mass_Wjet = book<TH1F>("all_mass_Wjet", "all M_{ij}", 50, 50, 150);

  min_mass_Wjet_zoom = book<TH1F>("min_mass_Wjet_zoom", "min M_{ij}", 60, 60, 120);
  pt_Wjet = book<TH1F>("pt_Wjet", "p_{T,ij}", 100, 0, 500);
  pt_Wjet_i = book<TH1F>("pt_Wjet_i", "p_{T,i}", 100, 0, 500);
  pt_Wjet_j = book<TH1F>("pt_Wjet_j", "p_{T,j}", 100, 0, 500);
  pt_Wjet_diff = book<TH1F>("pt_Wjet_diff", "p_{T,i} - p_{T,j}", 200, -400, 400);
  pt_Wjets = book<TH2F>("pt_Wjets", "x=p_{T,i} y=p_{T,j}", 100, 0, 500, 100, 0, 500);
  Wjet_combination = book<TH1F>("Wjet_combination", "combination", 3, 0, 3);

  pt_subjets_had = book<TH1F>("pt_subjets_had", "had. subjets p_{T}", 100, 0, 500);
  pt_subjets_had1 = book<TH1F>("pt_subjets_had1", "1st had. subjet p_{T}", 100, 0, 500);
  pt_subjets_had2 = book<TH1F>("pt_subjets_had2", "2nd had. subjet p_{T}", 100, 0, 500);
  pt_subjets_had3 = book<TH1F>("pt_subjets_had3", "3rd had. subjet p_{T}", 100, 0, 500);

  // handle for clustered jets (2+3 und 3+3)
  h_fatjets33=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  h_hadjets33=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");
  h_lepjets33=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_lep_Combined");

  //h_softdrop=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets_softdrop");

}



void GenHists_xcone::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets  -------------------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<GenTopJet> hadjets33 = event.get(h_hadjets33);
  std::vector<GenTopJet> lepjets33 = event.get(h_lepjets33);
  std::vector<GenTopJet> fatjets33 = event.get(h_fatjets33);
  //std::vector<GenTopJet> softdrop = event.get(h_softdrop);

  if(hadjets33.size() == 0 || lepjets33.size() == 0) return;

  GenTopJet had33 = hadjets33.at(0);
  GenTopJet lep33 = lepjets33.at(0);

  // get had subjets
  int index_had = -1;
  if(deltaR(fatjets33[0], had33) < deltaR(fatjets33[0], lep33)) index_had = 0;
  else index_had = 1;
  std::vector<GenJet> had_subjets = fatjets33[index_had].subjets();

  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;

  LorentzVector lep33_v4;

  std::vector<GenParticle> leptons;
  bool is_semilep = false;
  for (unsigned int i=0; i<event.genparticles->size(); i++){
    if(abs(event.genparticles->at(i).pdgId()) == 11 || abs(event.genparticles->at(i).pdgId()) == 13 ){
      leptons.push_back(event.genparticles->at(i));
    }
  }
  if(leptons.size() == 1) is_semilep = true;
  if(is_semilep){
    lep33_v4 = lep33.v4() + leptons[0].v4();
  }
  else{
    lep33_v4 = lep33.v4();
  }

  Mass_HadJet33->Fill(had33.v4().M(), weight);
  Mass_HadJet33_rebin->Fill(had33.v4().M(), weight);
  Mass_HadJet33_unfold->Fill(had33.v4().M(), weight);
  Mass_HadJet33_B->Fill(had33.v4().M(), weight);
  Mass_HadJet33_C->Fill(had33.v4().M(), weight);


  Mass_LepJet33->Fill(lep33_v4.M(), weight);
  MassDiff33->Fill(had33.v4().M() - lep33_v4.M(), weight);

  PT_HadJet33->Fill(had33.pt(), weight);
  PT_LepJet33->Fill(lep33.pt(), weight);

  number_HadJet33->Fill(hadjets33.size(), weight);
  number_LepJet33->Fill(lepjets33.size(), weight);


  //SoftDropMass->Fill(softdrop.at(had_nr).v4().M(), weight);

  for(unsigned int i=0; i<had_subjets.size(); i++){
    double pt = had_subjets.at(i).pt();
    pt_subjets_had->Fill(pt, weight);
    if(i==0)      pt_subjets_had1->Fill(pt, weight);
    else if(i==1) pt_subjets_had2->Fill(pt, weight);
    else if(i==2) pt_subjets_had3->Fill(pt, weight);
  }
  //---------------------------------------------------------------------------------------
  //--------------------------------- add subjets to reconstruct W ------------------------
  //---------------------------------------------------------------------------------------
  if(had_subjets.size() == 3){
    double px, py, pz, E;
    TLorentzVector Wjet12, Wjet13, Wjet23;
    px = had_subjets.at(0).v4().Px() + had_subjets.at(1).v4().Px();
    py = had_subjets.at(0).v4().Py() + had_subjets.at(1).v4().Py();
    pz = had_subjets.at(0).v4().Pz() + had_subjets.at(1).v4().Pz();
    E = had_subjets.at(0).v4().E() + had_subjets.at(1).v4().E();
    Wjet12.SetPxPyPzE(px, py, pz, E);
    px = had_subjets.at(0).v4().Px() + had_subjets.at(2).v4().Px();
    py = had_subjets.at(0).v4().Py() + had_subjets.at(2).v4().Py();
    pz = had_subjets.at(0).v4().Pz() + had_subjets.at(2).v4().Pz();
    E = had_subjets.at(0).v4().E() + had_subjets.at(2).v4().E();
    Wjet13.SetPxPyPzE(px, py, pz, E);
    px = had_subjets.at(1).v4().Px() + had_subjets.at(2).v4().Px();
    py = had_subjets.at(1).v4().Py() + had_subjets.at(2).v4().Py();
    pz = had_subjets.at(1).v4().Pz() + had_subjets.at(2).v4().Pz();
    E = had_subjets.at(1).v4().E() + had_subjets.at(2).v4().E();
    Wjet23.SetPxPyPzE(px, py, pz, E);

    double M12, M13, M23, M_min = 1000;
    double pt1, pt2, pt3;
    pt1 = had_subjets.at(0).pt();
    pt2 = had_subjets.at(1).pt();
    pt3 = had_subjets.at(2).pt();
    M12 = Wjet12.M();
    M13 = Wjet13.M();
    M23 = Wjet23.M();
    if(M12 < M13 && M12 < M23){
      M_min = M12;
      if(pt1 <= 120 && pt2 <= 120) min_mass_Wjet_ptlow->Fill(M_min, weight);
      if(pt1 > 120 && pt2 > 120) min_mass_Wjet_pthigh->Fill(M_min, weight);
      pt_Wjet->Fill(Wjet12.Pt(), weight);
      pt_Wjet_i->Fill(pt1, weight);
      pt_Wjet_j->Fill(pt2, weight);
      pt_Wjet_diff->Fill(pt1 - pt2, weight);
      pt_Wjets->Fill(pt1, pt2, weight);
      Wjet_combination->Fill(0.5, weight);
    }
    if(M13 < M12 && M13 < M23){
      M_min = M13;
      if(pt1 <= 120 && pt3 <= 120) min_mass_Wjet_ptlow->Fill(M_min, weight);
      if(pt1 > 120 && pt3 > 120) min_mass_Wjet_pthigh->Fill(M_min, weight);
      pt_Wjet->Fill(Wjet13.Pt(), weight);
      pt_Wjet_i->Fill(pt1, weight);
      pt_Wjet_j->Fill(pt3, weight);
      pt_Wjet_diff->Fill(pt1 - pt3, weight);
      pt_Wjets->Fill(pt1, pt3, weight);
      Wjet_combination->Fill(1.5, weight);
    }
    if(M23 < M12 && M23 < M13){
      M_min = M23;
      if(pt2 <= 120 && pt3 <= 120) min_mass_Wjet_ptlow->Fill(M_min, weight);
      if(pt2 > 120 && pt3 > 120) min_mass_Wjet_pthigh->Fill(M_min, weight);
      pt_Wjet->Fill(Wjet23.Pt(), weight);
      pt_Wjet_i->Fill(pt2, weight);
      pt_Wjet_j->Fill(pt3, weight);
      pt_Wjet_diff->Fill(pt2 - pt3, weight);
      pt_Wjets->Fill(pt2, pt3, weight);
      Wjet_combination->Fill(2.5, weight);
    }

    if(M_min != 1000){
      min_mass_Wjet->Fill(M_min, weight);
      min_mass_Wjet_zoom->Fill(M_min, weight);
    }
  }
}
