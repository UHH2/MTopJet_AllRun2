#include <UHH2/MTopJet/include/RecoGenHists_xcone_topjet.h>


RecoGenHists_xcone_topjet::RecoGenHists_xcone_topjet(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // Reco Histograms -----------------------------------------------------------
  HadJetPt = book<TH1F>("pt_hadjet", "p_{T}", 50, 0, 1000);
  HadJetMass = book<TH1F>("M_hadjet", "m_{jet} [GeV]", 50, 0, 500);
  HadJetEta = book<TH1F>("eta_jet1", "#eta", 50, -5, 5);
  HadJetPhi = book<TH1F>("phi_jet1", "#phi", 50, -2*M_PI, 2*M_PI);

  distance_hadx_hadak = book<TH1F>("deltaR_hadx_hadak", "dR", 50, 0, 5);

  XCone_HadJetEta = book<TH1F>("eta_jet1", "#eta", 50, -5, 5);
  XCone_HadJetPhi = book<TH1F>("phi_jet1", "#phi", 50, -2*M_PI, 2*M_PI);

  HadJetTau1 = book<TH1F>("tau1_hadjet", "#tau_{1}", 50, 0, 1);
  HadJetTau2 = book<TH1F>("tau2_hadjet", "#tau_{2}", 50, 0, 1);
  HadJetTau3 = book<TH1F>("tau3_hadjet", "#tau_{3}", 50, 0, 1);

  HadJetTau32 = book<TH1F>("tau32_hadjet", "#tau_{3}/#tau_{2}", 50, 0, 1);
  HadJetTau23 = book<TH1F>("tau23_hadjet", "#tau_{2}/#tau_{3}", 50, 0, 1);

  number_ak8jets = book<TH1F>("number_ak8jets", "number AK8 jets", 10, 0, 10);
  number_xconejets = book<TH1F>("number_xconejets", "number XCone Jets", 10, 0, 10);

  // Gen Histograms ------------------------------------------------------------

  Number_matched_all = book<TH1F>("events_matched_all", "matched with all particles", 2, 0, 2);
  Number_matched_top = book<TH1F>("events_matched_top", "matched with top", 2, 0, 2);
  Number_matched_q1 = book<TH1F>("events_matched_q1", "matched with q1", 2, 0, 2);
  Number_matched_q2 = book<TH1F>("events_matched_q2", "matched with q2", 2, 0, 2);
  Number_matched_bottom = book<TH1F>("events_matched_bottom", "matched with bottom", 2, 0, 2);

  HadJetMass_fully_merged = book<TH1F>("M_hadjet_fully_merged", "m_{jet} [GeV]", 50, 0, 500);

  // handle for clustered jets -------------------------------------------------
  h_recjets=ctx.get_handle<std::vector<TopJet>>("jetsAk8CHSSubstructure_SoftDropCHS");
  h_xcone=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

}



void RecoGenHists_xcone_topjet::fill(const Event & event){


  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> rec_jets = event.get(h_recjets);
  std::vector<TopJet> xcone_jets = event.get(h_xcone);


  // Select AK8 Jet in range of XCone size.
  TLorentzVector rec_jet, xcone_jet, had_jet;
  // TLorentzVector rec_jet1_v4, rec_jet2_v4, gen_jet1_v4, gen_jet2_v4;
  std::vector<TLorentzVector> had_jets;
  std::vector<TopJet> had_jets_ak;
  TopJet had_jet_ak;

  xcone_jet.SetPxPyPzE(xcone_jets.at(0).v4().Px(), xcone_jets.at(0).v4().Py(), xcone_jets.at(0).v4().Pz(), xcone_jets.at(0).v4().E());

  for(unsigned int i=0; i<rec_jets.size(); i++){
    rec_jet.SetPxPyPzE(rec_jets.at(i).v4().Px(), rec_jets.at(i).v4().Py(), rec_jets.at(i).v4().Pz(), rec_jets.at(i).v4().E());
    if(deltaR(xcone_jets.at(0), rec_jets.at(i)) < 0.8){
      had_jets.push_back(rec_jet);
      had_jets_ak.push_back(rec_jets.at(i));
    }
  }

  double jet1_pt, jet2_pt;
  if(had_jets.size()==1){
    had_jet = had_jets.at(0);
    had_jet_ak = had_jets_ak.at(0);
  }
  else{
    jet1_pt = had_jets.at(0).Pt();
    jet2_pt = had_jets.at(1).Pt();
    if(jet1_pt>jet2_pt){
      had_jet=had_jets.at(0);
      had_jet_ak = had_jets_ak.at(0);
    }
    else{
      had_jet=had_jets.at(1);
      had_jet_ak = had_jets_ak.at(1);
    }
  }

  //---------------------------------------------------------------------------------------
  //---------------------------- Fill Reco Hists here -------------------------------------
  //---------------------------------------------------------------------------------------
  // get weight
  double weight = event.weight;

  HadJetPt->Fill(had_jet.Pt(), weight);
  HadJetMass->Fill(had_jet.M(), weight);
  HadJetEta->Fill(had_jet.Eta(), weight);
  HadJetPhi->Fill(had_jet.Phi(), weight);

  XCone_HadJetEta->Fill(xcone_jets.at(0).eta(), weight);
  XCone_HadJetPhi->Fill(xcone_jets.at(0).phi(), weight);

  HadJetTau1->Fill(had_jet_ak.tau1(), weight);
  HadJetTau2->Fill(had_jet_ak.tau2(), weight);
  HadJetTau3->Fill(had_jet_ak.tau3(), weight);
  HadJetTau32->Fill(had_jet_ak.tau3()/had_jet_ak.tau2(), weight);
  HadJetTau23->Fill(had_jet_ak.tau3()/had_jet_ak.tau2(), weight);

  distance_hadx_hadak->Fill(deltaR(xcone_jets.at(0), had_jet_ak), weight);

  number_ak8jets->Fill(rec_jets.size(), weight);
  number_xconejets->Fill(xcone_jets.size(), weight);
  //---------------------------------------------------------------------------------------
  //--------------------- Matching GenParticles -------------------------------------------
  //---------------------------------------------------------------------------------------

  TTbarGen ttbargen = event.get(h_ttbargen);
  GenParticle top, bottom, q1, q2;

  bool dilep = false;
  if(ttbargen.IsTopHadronicDecay()){
    bottom = ttbargen.bTop();
    q1 = ttbargen.Wdecay1();
    q2 = ttbargen.Wdecay2();
    top = ttbargen.Top();
  }
  else if(ttbargen.IsAntiTopHadronicDecay()){
    bottom = ttbargen.bAntitop();
    q1 = ttbargen.WMinusdecay1();
    q2 = ttbargen.WMinusdecay2();
    top = ttbargen.Antitop();
  }

  bool matched = false;
  bool matched_top = false;
  bool matched_q1 = false;
  bool matched_q2 = false;
  bool matched_bottom = false;
  if(deltaR(q1, had_jet_ak) < 0.8) matched_q1 = true;
  if(deltaR(q2, had_jet_ak) < 0.8) matched_q2 = true;
  if(deltaR(bottom, had_jet_ak) < 0.8) matched_bottom = true;

  if(matched_q1 && matched_q2 && matched_bottom) matched = true;

  //---------------------------------------------------------------------------------------
  //---------------------------- Fill Gen Hists here --------------------------------------
  //---------------------------------------------------------------------------------------

  if(matched_top) Number_matched_top->Fill(1, weight);
  else Number_matched_top->Fill(0.1, weight);

  if(matched_q1) Number_matched_q1->Fill(1, weight);
  else Number_matched_q1->Fill(0.1, weight);

  if(matched_q2) Number_matched_q2->Fill(1, weight);
  else Number_matched_q2->Fill(0.1, weight);

  if(matched_bottom) Number_matched_bottom->Fill(1, weight);
  else Number_matched_bottom->Fill(0.1, weight);

  if(matched){
    Number_matched_all->Fill(1, weight);
    HadJetMass_fully_merged->Fill(had_jet.M(), weight);
  }
  else Number_matched_all->Fill(0.1, weight);



  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


}
