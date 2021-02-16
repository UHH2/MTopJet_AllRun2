#include <UHH2/MTopJet/include/RecoGenHists_xcone_topjet.h>


RecoGenHists_xcone_topjet::RecoGenHists_xcone_topjet(uhh2::Context & ctx, const std::string & dirname, bool & isTTbar, const double & masscut): Hists(ctx, dirname), isTTbar_(isTTbar), masscut_(masscut){

  h_number_events = book<TH1F>("number_events_total", "events", 1, 0, 1);

  // -----------------------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------------------
  // Tau Histograms --------------------------------------------------------------------------------------------

  // Reco Histograms -------------------------------------------------------------------------------------------
  h_HadJetPt = book<TH1F>("ak8_hadjet_pt", "p_{T} [GeV]", 50, 0, 1000);
  h_HadJetMass = book<TH1F>("ak8_hadjet_mass", "m_{jet} [GeV]", 50, 0, 500);
  h_HadJetEta = book<TH1F>("ak8_hadjet_eta", "#eta", 50, -5, 5);
  h_HadJetPhi = book<TH1F>("ak8_hadjet_phi", "#phi", 50, -2*M_PI, 2*M_PI);

  h_distance_hadx_hadak = book<TH1F>("deltaR_hadjet_ak8_xcone", "dR", 60, 0, 3);

  h_XCone_HadJetEta = book<TH1F>("xcone_hadjet_eta", "#eta", 50, -5, 5);
  h_XCone_HadJetPhi = book<TH1F>("xcone_hadjet_phi", "#phi", 50, -2*M_PI, 2*M_PI);

  h_HadJetTau1 = book<TH1F>("ak8_hadjet_tau1", "#tau_{1}", 50, 0, 1);
  h_HadJetTau2 = book<TH1F>("ak8_hadjet_tau2", "#tau_{2}", 50, 0, 1);
  h_HadJetTau3 = book<TH1F>("ak8_hadjet_tau3", "#tau_{3}", 50, 0, 1);

  h_HadJetTau32 = book<TH1F>("ak8_hadjet_tau32", "#tau_{3}/#tau_{2}", 50, 0, 1);
  h_HadJetTau23 = book<TH1F>("ak8_hadjet_tau23", "#tau_{2}/#tau_{3}", 100, 0, 2);

  h_number_ak8_jets        = book<TH1F>("ak8_number_matched_hadjets", "number matched AK8 jets", 10, 0, 10);
  h_number_ak8_not_matched = book<TH1F>("ak8_number_unmatched_hadjets", "number unmatched AK8 jets", 1, 0, 1);
  h_number_ak8_had_jets    = book<TH1F>("ak8_number_hadjets", "number AK8 hadjets", 10, 0, 10);
  h_number_xcone_had_jet   = book<TH1F>("xcone_number_hadjets", "number XCone hadjets", 10, 0, 10);

  // Gen Histograms --------------------------------------------------------------------------------------------
  // if(isTTbar){
  //   h_Number_matched_all    = book<TH1F>("gen_events_matched_all", "matched with all particles", 2, 0, 2);
  //   h_Number_matched_top    = book<TH1F>("gen_events_matched_top", "matched with top", 2, 0, 2);
  //   h_Number_matched_q1     = book<TH1F>("gen_events_matched_q1", "matched with q1", 2, 0, 2);
  //   h_Number_matched_q2     = book<TH1F>("gen_events_matched_q2", "matched with q2", 2, 0, 2);
  //   h_Number_matched_bottom = book<TH1F>("gen_events_matched_bottom", "matched with bottom", 2, 0, 2);
  //
  //   h_HadJetMass_fullymerged = book<TH1F>("gen_mass_hadjet_fullymerged", "m_{jet} [GeV]", 50, 0, 500);
  //   h_HadJetMass_semimerged  = book<TH1F>("gen_mass_hadjet_semimerged", "m_{jet} [GeV]", 50, 0, 500);
  //   h_HadJetMass_notmerged   = book<TH1F>("gen_mass_hadjet_notmerged", "m_{jet} [GeV]", 50, 0, 500);
  //
  //   h_HadJetTau32_fullymerged = book<TH1F>("gen_tau32_hadjet_fullymerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
  //   h_HadJetTau32_semimerged  = book<TH1F>("gen_tau32_hadjet_semimerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
  //   h_HadJetTau32_notmerged   = book<TH1F>("gen_tau32_hadjet_notmerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
  // }

  // -----------------------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------------------
  // WJet Histograms -------------------------------------------------------------------------------------------

  // dummy hists
  h_number_ak4_jets = book<TH1F>("number_ak4_jets", "number AK4 jets", 20, 0, 20);
  h_number_ak4_matched_xcone_had = book<TH1F>("number_ak4_jets_matched_xcone_fathadjet", "number matched AK4 jets", 10, 0, 10);
  h_number_equal_distance = book<TH1F>("number_events_distance_ak4_xcone_equal", "events with equal distance", 1, 0, 1);

  h_ak4_btag = book<TH1F>("ak4_all_btag", "btag", 50, 0, 1);

  // distance ak4 subjets
  h_distance_ak4_xcone_subjet_btag = book<TH1F>("ak4_distance_xcone_btag", "dR_{small}(ak4, xcone)", 50, 0, 1);
  h_difference_distance_small_mid  = book<TH1F>("ak4_distance_subjets_small_mid", "dR_{small}-dR_{mid}", 40, 0, 4);
  h_difference_distance_small_high = book<TH1F>("ak4_distance_subjets_small_high", "dR_{small}-dR_{high}", 40, 0, 4);
  h_difference_distance_mid_high   = book<TH1F>("ak4_distance_subjets_mid_high", "dR_{mid}-dR_{high}", 40, 0, 4);

  // Btag
  h_btag_difference_high_mid              = book<TH1F>("btag_difference_high_mid", "btag_{high}-btag_{mid}", 50, 0, 1);
  h_btag_percentage_high_mid              = book<TH1F>("btag_percentage_high_mid", "btag_{mid}/btag_{high}", 50, 0, 1);
  h_btag_difference_high_mid_sel_btag_min = book<TH1F>("btag_difference_high_mid_sel_btag_min", "btag_{high}-btag_{mid}", 50, 0, 1);
  h_btag_percentage_high_mid_sel_btag_min = book<TH1F>("btag_percentage_high_mid_sel_btag_min", "btag_{mid}/btag_{high}", 50, 0, 1);

  // Method: Match highest btag --------------------------------------------------------------------------------

  //Wjet
  h_wjet_bigger_pt_subjet  = book<TH1F>("wjet_bigger_pt_subjet", "", 2, 0, 2);
  h_wjet_pt_match          = book<TH1F>("wjet_pt_match", "p_{T} [GeV]", 50, 0, 1000);
  h_wjet_pt_match_subjet1  = book<TH1F>("wjet_pt_match_subjet1", "p_{T} [GeV]", 100, 0, 1000);
  h_wjet_pt_match_subjet2  = book<TH1F>("wjet_pt_match_subjet2", "p_{T} [GeV]", 100, 0, 1000);
  h_wjet_pt_match_S1divS2  = book<TH1F>("wjet_pt_match_S1divS2", "p_{T} [GeV]", 100, 0, 10);
  h_wjet_pt_match_S1divW   = book<TH1F>("wjet_pt_match_S1divW", "p_{T} [GeV]", 20, 0, 1);
  h_wjet_pt_match_S1mS2    = book<TH1F>("wjet_pt_match_S1mS2", "p_{T} [GeV]", 51, -10, 500);
  h_wjet_pt_match_S1_S2    = book<TH2F>("wjet_pt_match_S1_S2", "W Subjet pt [GeV]", 100, 0, 1000, 100, 0, 1000);
  h_wjet_pt_match_S1divW_W = book<TH2F>("wjet_pt_match_S1divW_W", "W Subjet pt [GeV]", 100, 0, 1000, 100, 0, 1);
  h_wmass_match            = book<TH1F>("wmass_match", "m_{Wjet} [GeV]", 180, 0, 180);

  // pt bins
  h_wmass_match_ptbin_low    = book<TH1F>("wmass_match_ptbin_low", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wmass_match_ptbin_high   = book<TH1F>("wmass_match_ptbin_high", "m_{Wjet} [GeV]", 180, 0, 180);

  h_wmass_match_divbin_low   = book<TH1F>("wmass_match_divbin_low", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wmass_match_divbin_high  = book<TH1F>("wmass_match_divbin_high", "m_{Wjet} [GeV]", 180, 0, 180);

  h_wmass_match_ptdiv_hh     = book<TH1F>("wmass_match_ptdiv_hh", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wmass_match_ptdiv_hl     = book<TH1F>("wmass_match_ptdiv_hl", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wmass_match_ptdiv_lh     = book<TH1F>("wmass_match_ptdiv_lh", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wmass_match_ptdiv_ll     = book<TH1F>("wmass_match_ptdiv_ll", "m_{Wjet} [GeV]", 180, 0, 180);


  //Btag
  h_ak4_btag_high = book<TH1F>("ak4_highest_btag", "btag", 50, 0, 1);

  // Counting Events
  h_events_no_ak4 = book<TH1F>("number_events_no_ak4", "events with no matched ak4", 1, 0, 1);
  h_events_ak4    = book<TH1F>("number_events_with_ak4", "events with matched ak4", 1, 0, 1);

  // dummy
  h_number_no_close_subjets       = book<TH1F>("number_no_close_subjets", "no jet with dR<0.4", 1, 0, 1);
  h_number_one_close_subjets      = book<TH1F>("number_one_close_subjets", "one jet with dR<0.4", 1, 0, 1);
  h_number_multiple_close_subjets = book<TH1F>("number_multiple_close_subjets", "two jets with dR<0.4", 1, 0, 1);

  // BtagCut ----------------------------------------------------------------------------------------------------
  h_ak4_btag_high_btagcut = book<TH1F>("number_ak4_btag_high_btagcut", "matched AK4 jet highest btag", 50, 0, 1);

  // WJet
  h_wjet_pt_btagcut = book<TH1F>("wjet_pt_btagcut", "p_{T} [GeV]", 50, 0, 1000);
  h_wmass_btagcut   = book<TH1F>("wmass_btagcut", "m_{Wjet} [GeV]", 180, 0, 180);

  // Counting Events
  h_btagcut_missed = book<TH1F>("selection_miss_btagcut", "event missed btag selection", 1, 0, 1);
  h_btagcut_passed = book<TH1F>("selection_pass_btagcut", "event passed btag selection", 1, 0, 1);

  // Selection: One high btag ----------------------------------------------------------------------------------
  h_wmass_one_btag   = book<TH1F>("wmass_one_btag", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_one_btag = book<TH1F>("wjet_pt_one_btag", "p_{T} [GeV]", 50, 0, 1000);

  // with btagcut
  h_wmass_btagcut_one_btag   = book<TH1F>("wmass_btagcut_one_btag", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_btagcut_one_btag = book<TH1F>("wjet_pt_btagcut_one_btag", "p_{T} [GeV]", 50, 0, 1000);

  // Counting events
  h_pass_one_btag         = book<TH1F>("selection_pass_one_btag", "", 1, 0, 1);
  h_pass_btagcut_one_btag = book<TH1F>("selection_pass_btagcut_one_btag", "", 1, 0, 1);

  // Selection: One subjet -------------------------------------------------------------------------------------
  h_wmass_one_subjet   = book<TH1F>("wmass_one_subjet", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_one_subjet = book<TH1F>("wjet_pt_one_subjet", "p_{T} [GeV]", 50, 0, 1000);

  // with btagcut
  h_wmass_btagcut_one_subjet   = book<TH1F>("wmass_btagcut_one_subjet", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_btagcut_one_subjet = book<TH1F>("wjet_pt_btagcut_one_subjet", "p_{T} [GeV]", 50, 0, 1000);

  // Counting events
  h_pass_one_subjet         = book<TH1F>("selection_pass_one_subjet", "", 1, 0, 1);
  h_pass_btagcut_one_subjet = book<TH1F>("selection_pass_btagcut_one_subjet", "", 1, 0, 1);

  // Selection: One high btag and one subjet  ------------------------------------------------------------------
  h_wmass_one_btag_subjet   = book<TH1F>("wmass_one_btag_subjet", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_one_btag_subjet = book<TH1F>("wjet_pt_one_btag_subjet", "p_{T} [GeV]", 50, 0, 1000);

  // with btagcut
  h_wmass_btagcut_one_btag_subjet   = book<TH1F>("wmass_btagcut_one_btag_subjet", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_btagcut_one_btag_subjet = book<TH1F>("wjet_pt_btagcut_one_btag_subjet", "p_{T} [GeV]", 50, 0, 1000);

  // Counting events
  h_pass_one_btag_subjet         = book<TH1F>("selection_pass_one_btag_subjet", "", 1, 0, 1);
  h_pass_btagcut_one_btag_subjet = book<TH1F>("selection_pass_btagcut_one_btag_subjet", "", 1, 0, 1);

  // Selection btag>0.7 and minimal mass -----------------------------------------------------------------------
  h_wmass_combined_methods   = book<TH1F>("wmass_combined_methods", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_combined_methods = book<TH1F>("wjet_pt_combined", "p_{T} [GeV]", 50, 0, 1000);
  h_wjet_combination         = book<TH1F>("wjet_combination", "combination", 4, 0, 4);

  // Method - Min Mass -----------------------------------------------------------------------------------------
  h_wmass_min_mass   = book<TH1F>("wmass_min", "min m_{ij}", 180, 0, 180);
  h_wjet_pt_min_mass = book<TH1F>("wjet_pt_min", "p_{T} [GeV]", 50, 0, 1000);

  // Method - Compare ------------------------------------------------------------------------------------------
  h_wmass_compare   = book<TH1F>("wmass_compare", "m_{Wjet} [GeV]", 180, 0, 180);
  h_wjet_pt_compare = book<TH1F>("wjet_pt_compare", "p_{T} [GeV]", 50, 0, 1000);


  // -----------------------------------------------------------------------------------------------------------
  // handle for clustered jets ---------------------------------------------------------------------------------
  h_recjets  = ctx.get_handle<std::vector<TopJet>>("topjets"); // jetsAk8CHSSubstructure_SoftDropCHS
  h_hadjet   = ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  h_lepjet   = ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined_Corrected");
  h_fatjets  = ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

}

void RecoGenHists_xcone_topjet::fill(const Event & event){

  bool debug = false;
  double weight = event.weight;  // get weight
  h_number_events->Fill(0.5, weight);

  // ****************************************************************************************************************************************************
  // ****************************************************************************************************************************************************
  // ****************************************************************************************************************************************************
  /*
  .     ████████  █████  ██    ██
  .        ██    ██   ██ ██    ██
  .        ██    ███████ ██    ██
  .        ██    ██   ██ ██    ██
  .        ██    ██   ██  ██████
  */

  /*
  ██████  ███████  ██████
  ██   ██ ██      ██
  ██████  █████   ██
  ██   ██ ██      ██
  ██   ██ ███████  ██████
  */
  if(debug) std::cout << "Tau: Rec" << '\n';
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects--------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> ak8_jets      = event.get(h_recjets);
  std::vector<TopJet> xcone_had_jet = event.get(h_hadjet);
  std::vector<TopJet> xcone_lep_jet = event.get(h_lepjet);
  std::vector<TopJet> xcone_fatjets = event.get(h_fatjets);
  if(debug){
    std::cout << "Number HadJets: "<< xcone_had_jet.size() << '\n';
    std::cout << "Number LepJets: "<< xcone_lep_jet.size() << '\n';
    std::cout << "Number FatJets: "<< xcone_fatjets.size() << '\n';
  }
  // Select AK8 Jet in range of XCone size.
  // Looking at each jet seperatly by going through with a for-Loop
  TLorentzVector ak8_v4, xcone_had_v4;

  // TLorentzVector rec_jet1_v4, rec_jet2_v4, gen_jet1_v4, gen_jet2_v4;
  std::vector<TLorentzVector> ak8_had_v4_vector;
  std::vector<TopJet> ak8_had_vector;
  TopJet ak8_had;

  xcone_had_v4.SetPxPyPzE(xcone_had_jet.at(0).v4().Px(), xcone_had_jet.at(0).v4().Py(), xcone_had_jet.at(0).v4().Pz(), xcone_had_jet.at(0).v4().E());

  for(unsigned int i=0; i<ak8_jets.size(); i++){
    ak8_v4.SetPxPyPzE(ak8_jets.at(i).v4().Px(), ak8_jets.at(i).v4().Py(), ak8_jets.at(i).v4().Pz(), ak8_jets.at(i).v4().E());

    if(deltaR(xcone_had_jet.at(0), ak8_jets.at(i)) < 0.8){
      ak8_had_v4_vector.push_back(ak8_v4);
      ak8_had_vector.push_back(ak8_jets.at(i));
    }
  }

  // Selecting the hadronic AK8 jet
  TLorentzVector ak8_had_v4;
  double jet1_pt = 0;
  double jet2_pt = 0;
  if(ak8_had_v4_vector.size()>=1){
    if(ak8_had_v4_vector.size()==1){
      ak8_had_v4 = ak8_had_v4_vector.at(0);
      ak8_had    = ak8_had_vector.at(0);
    }
    else{
      jet1_pt = ak8_had_v4_vector.at(0).Pt();
      jet2_pt = ak8_had_v4_vector.at(1).Pt();
      if(jet1_pt>jet2_pt){
        ak8_had_v4 = ak8_had_v4_vector.at(0);
        ak8_had    = ak8_had_vector.at(0);
      }
      else{
        ak8_had_v4 = ak8_had_v4_vector.at(1);
        ak8_had    = ak8_had_vector.at(1);
      }
    }
  }
  else if(ak8_had_v4_vector.size()==0) h_number_ak8_not_matched->Fill(0.5, weight);

  // N-subjettiness -----------------------------------------------------------------------
  double tau_32_ak  = ak8_had.tau3()/ak8_had.tau2();
  mass = ak8_had_v4.M();
  tau32 = tau_32_ak;

  bool pass_masscut = false;
  if(ak8_had_v4.M()>masscut_) pass_masscut = true;

  /*
  ███████ ██ ██      ██          ██████  ███████  ██████
  ██      ██ ██      ██          ██   ██ ██      ██
  █████   ██ ██      ██          ██████  █████   ██
  ██      ██ ██      ██          ██   ██ ██      ██
  ██      ██ ███████ ███████     ██   ██ ███████  ██████
  */
  if(debug) std::cout << "Tau: Rec fill" << '\n';

  if(pass_masscut){

    h_HadJetPt->Fill(ak8_had_v4.Pt(), weight);
    h_HadJetMass->Fill(ak8_had_v4.M(), weight);
    h_HadJetEta->Fill(ak8_had_v4.Eta(), weight);
    h_HadJetPhi->Fill(ak8_had_v4.Phi(), weight);

    h_number_ak8_had_jets->Fill(ak8_had_v4_vector.size(), weight);

    h_XCone_HadJetEta->Fill(xcone_had_jet.at(0).eta(), weight);
    h_XCone_HadJetPhi->Fill(xcone_had_jet.at(0).phi(), weight);

    h_HadJetTau1->Fill(ak8_had.tau1(), weight);
    h_HadJetTau2->Fill(ak8_had.tau2(), weight);
    h_HadJetTau3->Fill(ak8_had.tau3(), weight);
    h_HadJetTau32->Fill(tau_32_ak, weight);
    h_HadJetTau23->Fill(1/tau_32_ak, weight);

    h_distance_hadx_hadak->Fill(deltaR(xcone_had_jet.at(0), ak8_had), weight);

    h_number_ak8_jets->Fill(ak8_jets.size(), weight);
    h_number_xcone_had_jet->Fill(xcone_had_jet.size(), weight);

    /*
    ███    ███  █████  ████████  ██████ ██   ██      ██████  ███████ ███    ██
    ████  ████ ██   ██    ██    ██      ██   ██     ██       ██      ████   ██
    ██ ████ ██ ███████    ██    ██      ███████     ██   ███ █████   ██ ██  ██
    ██  ██  ██ ██   ██    ██    ██      ██   ██     ██    ██ ██      ██  ██ ██
    ██      ██ ██   ██    ██     ██████ ██   ██      ██████  ███████ ██   ████
    */
    if(debug) std::cout << "Tau: Gen" << '\n';

    if(isTTbar_){
      TTbarGen ttbargen = event.get(h_ttbargen);
      GenParticle top, bottom, q1, q2;

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
      bool semimatched = false;
      bool notmatched = false;
      bool matched_top = false;
      bool matched_q1 = false;
      bool matched_q2 = false;
      bool matched_bottom = false;
      if(deltaR(q1, ak8_had) < 0.8) matched_q1 = true;
      if(deltaR(q2, ak8_had) < 0.8) matched_q2 = true;
      if(deltaR(bottom, ak8_had) < 0.8) matched_bottom = true;

      if(matched_q1 && matched_q2 && matched_bottom) matched = true;

      if((matched_q1 && !matched_q2 && matched_bottom) ||
      (!matched_q1 && matched_q2 && matched_bottom)) semimatched = true;


      if(!matched_q1 && !matched_q2 && !matched_bottom) notmatched = true;

      /*
      ███████ ██ ██      ██           ██████  ███████ ███    ██
      ██      ██ ██      ██          ██       ██      ████   ██
      █████   ██ ██      ██          ██   ███ █████   ██ ██  ██
      ██      ██ ██      ██          ██    ██ ██      ██  ██ ██
      ██      ██ ███████ ███████      ██████  ███████ ██   ████
      */
      if(debug) std::cout << "Tau: Gen fill" << '\n';


      // if(matched_top) h_Number_matched_top->Fill(1, weight);
      // else h_Number_matched_top->Fill(0.1, weight);

      // if(matched_q1) h_Number_matched_q1->Fill(1, weight);
      // else h_Number_matched_q1->Fill(0.1, weight);

      // if(matched_q2) h_Number_matched_q2->Fill(1, weight);
      // else h_Number_matched_q2->Fill(0.1, weight);

      // if(matched_bottom) h_Number_matched_bottom->Fill(1, weight);
      // else h_Number_matched_bottom->Fill(0.1, weight);

      if(semimatched){
        // h_HadJetTau32_semimerged->Fill(tau_32_ak, weight);
        // h_HadJetMass_semimerged->Fill(ak8_had_v4.M(), weight);
      }

      if(notmatched){
        // h_HadJetTau32_notmerged->Fill(tau_32_ak, weight);
        // h_HadJetMass_notmerged->Fill(ak8_had_v4.M(), weight);
      }

      if(matched){
        // h_Number_matched_all->Fill(1, weight);
        // h_HadJetMass_fullymerged->Fill(ak8_had_v4.M(), weight);
        // h_HadJetTau32_fullymerged->Fill(tau_32_ak, weight);
      }
      //else h_Number_matched_all->Fill(0.1, weight);
    }
  }

  // ****************************************************************************************************************************************************
  // ****************************************************************************************************************************************************
  // ****************************************************************************************************************************************************

  /*
  .  ██     ██      ██ ███████ ████████ ███████
  .  ██     ██      ██ ██         ██    ██
  .  ██  █  ██      ██ █████      ██    ███████
  .  ██ ███ ██ ██   ██ ██         ██         ██
  .   ███ ███   █████  ███████    ██    ███████
  */
  if(debug) std::cout << "Wjets" << '\n';

  /*
  At first all the different methods are implemented. Afterwards all the Selections are defined
  */

  int index_had = -1;
  if(deltaR(xcone_fatjets[0], xcone_had_jet[0]) < deltaR(xcone_fatjets[0], xcone_lep_jet[0])) index_had = 0;
  else index_had = 1;

  std::vector<Jet> xcone_had_subjets = xcone_fatjets[index_had].subjets();
  std::vector<Jet> ak4_jets = *event.jets;
  h_number_ak4_jets->Fill(ak4_jets.size(), weight);

  if(debug){
    std::cout << "Number HadSubjets: "<< xcone_had_jet[0].subjets().size() << '\n';
    std::cout << "Number LepSubjets: "<< xcone_lep_jet[0].subjets().size() << '\n';
    for(unsigned int jet=0;jet<xcone_fatjets.size(); jet ++){
      std::cout << "Number FatSubjets - FatJet"<<jet<<": "<< xcone_fatjets[jet].subjets().size() << '\n';
    }
  }

  /*
  ███    ███ ███████ ████████ ██   ██  ██████  ██████  ███████
  ████  ████ ██         ██    ██   ██ ██    ██ ██   ██ ██
  ██ ████ ██ █████      ██    ███████ ██    ██ ██   ██ ███████
  ██  ██  ██ ██         ██    ██   ██ ██    ██ ██   ██      ██
  ██      ██ ███████    ██    ██   ██  ██████  ██████  ███████
  */
  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Method: Btag ------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */
  if(debug) std::cout << "Wjets Method: Btag" << '\n';

  // Select only ak4 jets with are close to the hadronic xcone jet
  std::vector<Jet> ak4_matched_xcone_had;
  for(const auto& j : *event.jets){
    if(deltaR(j, xcone_fatjets[index_had])<1.2){
      ak4_matched_xcone_had.push_back(j);
      h_ak4_btag->Fill(j.btag_DeepJet(), weight); // FILL HIST
    }
  }
  h_number_ak4_matched_xcone_had->Fill(ak4_matched_xcone_had.size(), weight);

  /*
  Choosing the AK4 jet with the highest btag to later identify the one of
  the three XCone subjets as the b quark
  */
  if(debug) std::cout << " - search for AK4 jet with higest btag" << '\n';
  std::vector<double> btag_v; // Using a vector to identifie events with ambigous btag results (two high btags i.e. through a c-quark)
  double btag_high = 0;
  Jet ak4_highest_btag;
  if(ak4_matched_xcone_had.size()==0){
    h_events_no_ak4->Fill(0.5, weight);
    return;
  }
  else if(ak4_matched_xcone_had.size()>0){
    h_events_ak4->Fill(0.5, weight);
    for(unsigned int i=0; i<ak4_matched_xcone_had.size(); i++){
      btag_v.push_back(ak4_matched_xcone_had[i].btag_DeepJet());
      if(i==0){
        btag_high = btag_v[i];
        ak4_highest_btag = ak4_matched_xcone_had[i];
      }
      if(i>0 && btag_high<btag_v[i]){
        btag_high = btag_v[i];
        ak4_highest_btag = ak4_matched_xcone_had[i];
      }
    }
  }
  h_ak4_btag_high->Fill(ak4_highest_btag.btag_DeepJet(), weight);

  if(debug) std::cout << " - Identify Wjet: Distance AK4 to subjets" << '\n';
  // Matching AK4 jet with highest btag to the three XCone subjets
  if(xcone_had_subjets.size()<3) return;
  double dR_subjet1_ak4 = deltaR(xcone_had_subjets[0], ak4_highest_btag);
  double dR_subjet2_ak4 = deltaR(xcone_had_subjets[1], ak4_highest_btag);
  double dR_subjet3_ak4 = deltaR(xcone_had_subjets[2], ak4_highest_btag);
  std::vector<double> dR_subjets = {dR_subjet1_ak4, dR_subjet2_ak4, dR_subjet3_ak4};
  sort(dR_subjets.begin(), dR_subjets.end()); // vector orderd from small to big

  if(debug) std::cout << " - Identify Wjet: Compare Distance" << '\n';
  int index_bjet = -1; int index_Wjet1 = -1; int index_Wjet2 = -1;
  if(dR_subjet1_ak4<dR_subjet2_ak4 && dR_subjet1_ak4<dR_subjet3_ak4){
    index_bjet = 0;
    index_Wjet1 = 1;
    index_Wjet2 = 2;
  }
  else if(dR_subjet2_ak4<dR_subjet1_ak4 && dR_subjet2_ak4<dR_subjet3_ak4){
    index_bjet = 1;
    index_Wjet1 = 0;
    index_Wjet2 = 2;
  }
  else if(dR_subjet3_ak4<dR_subjet1_ak4 && dR_subjet3_ak4<dR_subjet2_ak4){
    index_bjet = 2;
    index_Wjet1 = 0;
    index_Wjet2 = 1;
  }
  else{
    h_number_equal_distance->Fill(0.5, weight);
    return;
  }

  // Selection for minimum btag ------------------------------------------------
  if(debug) std::cout << " - Btag Selection" << '\n';
  bool sel_btag_min = false;
  if(btag_high > 0.7) sel_btag_min = true;

  // Get WJet  -----------------------------------------------------------------
  if(debug) std::cout << " - build Wjet" << '\n';
  TLorentzVector Wjet_v4_match;
  double px_wjet = xcone_had_subjets[index_Wjet1].v4().Px() + xcone_had_subjets[index_Wjet2].v4().Px();
  double py_wjet = xcone_had_subjets[index_Wjet1].v4().Py() + xcone_had_subjets[index_Wjet2].v4().Py();
  double pz_wjet = xcone_had_subjets[index_Wjet1].v4().Pz() + xcone_had_subjets[index_Wjet2].v4().Pz();
  double E_wjet = xcone_had_subjets[index_Wjet1].v4().E() + xcone_had_subjets[index_Wjet2].v4().E();
  Wjet_v4_match.SetPxPyPzE(px_wjet, py_wjet, pz_wjet, E_wjet);

  TLorentzVector Wjet_v4_match_subjetA = lorentz_to_tlorentz(xcone_had_subjets[index_Wjet1].v4()); // in Vector_utils.h
  TLorentzVector Wjet_v4_match_subjetB = lorentz_to_tlorentz(xcone_had_subjets[index_Wjet2].v4());

  // Sort subjets --------------------------------------------------------------
  if(debug) std::cout << " - sort Wsubjets" << '\n';
  TLorentzVector Wjet_v4_match_subjet1;
  TLorentzVector Wjet_v4_match_subjet2;
  if(Wjet_v4_match_subjetA.Pt()<Wjet_v4_match_subjetB.Pt())
  {
    Wjet_v4_match_subjet1 = Wjet_v4_match_subjetB;
    Wjet_v4_match_subjet2 = Wjet_v4_match_subjetA;
  } else {
    Wjet_v4_match_subjet1 = Wjet_v4_match_subjetA;
    Wjet_v4_match_subjet2 = Wjet_v4_match_subjetB;
  }

  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Method: Min Mass --------------------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */
  if(debug) std::cout << "Wjets Method: Min Mass" << '\n';

  double M_min = 1000;
  double px=0, py=0, pz=0, E=0;
  TLorentzVector Wjet12, Wjet13, Wjet23;
  px = xcone_had_subjets.at(0).v4().Px() + xcone_had_subjets.at(1).v4().Px();
  py = xcone_had_subjets.at(0).v4().Py() + xcone_had_subjets.at(1).v4().Py();
  pz = xcone_had_subjets.at(0).v4().Pz() + xcone_had_subjets.at(1).v4().Pz();
  E = xcone_had_subjets.at(0).v4().E() + xcone_had_subjets.at(1).v4().E();
  Wjet12.SetPxPyPzE(px, py, pz, E);
  px = xcone_had_subjets.at(0).v4().Px() + xcone_had_subjets.at(2).v4().Px();
  py = xcone_had_subjets.at(0).v4().Py() + xcone_had_subjets.at(2).v4().Py();
  pz = xcone_had_subjets.at(0).v4().Pz() + xcone_had_subjets.at(2).v4().Pz();
  E = xcone_had_subjets.at(0).v4().E() + xcone_had_subjets.at(2).v4().E();
  Wjet13.SetPxPyPzE(px, py, pz, E);
  px = xcone_had_subjets.at(1).v4().Px() + xcone_had_subjets.at(2).v4().Px();
  py = xcone_had_subjets.at(1).v4().Py() + xcone_had_subjets.at(2).v4().Py();
  pz = xcone_had_subjets.at(1).v4().Pz() + xcone_had_subjets.at(2).v4().Pz();
  E = xcone_had_subjets.at(1).v4().E() + xcone_had_subjets.at(2).v4().E();
  Wjet23.SetPxPyPzE(px, py, pz, E);

  // Get WJet  -----------------------------------------------------------------
  TLorentzVector Wjet_v4_min;
  double M12, M13, M23;
  M12 = Wjet12.M();
  M13 = Wjet13.M();
  M23 = Wjet23.M();
  if(M12 < M13 && M12 < M23){
    Wjet_v4_min = Wjet12;
    if(!sel_btag_min) h_wjet_combination->Fill(0.5, weight);
  }
  if(M13 < M12 && M13 < M23){
    Wjet_v4_min = Wjet13;
    if(!sel_btag_min) h_wjet_combination->Fill(1.5, weight);
  }
  if(M23 < M12 && M23 < M13){
    Wjet_v4_min = Wjet23;
    if(!sel_btag_min) h_wjet_combination->Fill(2.5, weight);
  }

  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Method: Mass close to Wmass ---------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */
  if(debug) std::cout << "Wjets Method: Wmass" << '\n';

  double m_W = 80.4;
  double diff_W_M12 = abs(M12-m_W);
  double diff_W_M13 = abs(M13-m_W);
  double diff_W_M23 = abs(M23-m_W);

  // Get WJet  -----------------------------------------------------------------
  TLorentzVector Wjet_v4_compare;
  // double M_compare = 1000;
  if(diff_W_M12<diff_W_M13 && diff_W_M12<diff_W_M23) Wjet_v4_compare = Wjet12;
  if(diff_W_M13<diff_W_M23 && diff_W_M13<diff_W_M12) Wjet_v4_compare = Wjet13;
  if(diff_W_M23<diff_W_M13 && diff_W_M23<diff_W_M12) Wjet_v4_compare = Wjet23;

  /*
  ███████ ███████ ██      ███████  ██████ ████████ ██  ██████  ███    ██ ███████
  ██      ██      ██      ██      ██         ██    ██ ██    ██ ████   ██ ██
  ███████ █████   ██      █████   ██         ██    ██ ██    ██ ██ ██  ██ ███████
  .    ██ ██      ██      ██      ██         ██    ██ ██    ██ ██  ██ ██      ██
  ███████ ███████ ███████ ███████  ██████    ██    ██  ██████  ██   ████ ███████
  */
  if(debug) std::cout << "Wjets: Selection" << '\n';


  // Define Selections ---------------------------------------------------------
  bool sel_one_btag   = false; // Only events with a second btag smaller than 0.6
  bool sel_one_subjet = false;

  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Selection: Btag high ----------------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */

  /* Defined at the end of - Method: Btag - */

  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Selection: ambigous btag ------------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */
  if(debug) std::cout << "Wjets Selection: ambigious btag" << '\n';

  /*
  The Selection for a second high btag is only used, if the higest btag passes the
  minimum btag criteria. The vector btag_v can only have one btag entrie.
  */
  sort(btag_v.begin(), btag_v.end());
  int btag_v_size = btag_v.size();
  if(sel_btag_min){
    for(int btag=0; btag<btag_v_size; btag++){
      if(btag_v[btag]==btag_high) continue; // If only one entry, this condition is always true.
      if(btag_v[btag]<0.6) sel_one_btag = true;
    }
  }
  /*
  -----------------------------------------------------------------------------------------------------------
  --------------------------------  Selection: dR of subjets ------------------------------------------------
  -----------------------------------------------------------------------------------------------------------
  */
  if(debug) std::cout << "Wjets Selection: distance subjets" << '\n';

  double number_close_subjets=0;
  // Keep in mind that the vector dR_subjets is already sorted
  for(unsigned int subjet=0; subjet<dR_subjets.size(); subjet++){
    if(dR_subjets[subjet]<0.4) number_close_subjets++;
  }

  if(number_close_subjets==0) h_number_no_close_subjets->Fill(0.5, weight);
  if(number_close_subjets==1){
    h_number_one_close_subjets->Fill(0.5, weight);
    sel_one_subjet = true;
  }
  if(number_close_subjets>1) h_number_multiple_close_subjets->Fill(0.5, weight);


  /*
  ███████ ██ ██      ██
  ██      ██ ██      ██
  █████   ██ ██      ██
  ██      ██ ██      ██
  ██      ██ ███████ ███████
  */
  if(debug) std::cout << "Wjets fill" << '\n';

  /*
  . Method Combined:            New method - Mass from subjet closest to ak4 jet with btag > 0.7
  .                             Old Method - If highest btag of ak4 jet is smaller than 0.7,
  .                             the Wjet is reconstructed by the combination of subjets
  .                             Therefore, the Histogram is filled in sel_btag_min and !sel_btag_min.

  . Method Btag matched:        XCone subjet clostest to *ak4 jet with highest btag
  .                             (*all ak4 jets with dR(Xcone_fathadjet, ak4) < 1.2

  . Selection btag:             highest btag > 0.7

  . Selection discrepancy:      second highest btag is greater/smaller than btag_high*(0.9)

  . Selection many close jets:  many subjets with dR(subjet, ak4_higest_btag) < 0.4
  */

  // ---------------------------------------------------------------------------------------------------------
  // Counting Events -----------------------------------------------------------------------------------------
  if(sel_btag_min){
    h_btagcut_passed->Fill(0.5, weight);
    if(sel_one_btag) h_pass_btagcut_one_btag->Fill(0.5, weight);
    if(sel_one_subjet) h_pass_btagcut_one_subjet->Fill(0.5, weight);
    if(sel_one_btag && sel_one_subjet) h_pass_btagcut_one_btag_subjet->Fill(0.5, weight);
  }
  else h_btagcut_missed->Fill(0.5, weight);

  if(sel_one_btag) h_pass_one_btag->Fill(0.5, weight);
  if(sel_one_subjet) h_pass_one_subjet->Fill(0.5, weight);
  if(sel_one_btag && sel_one_subjet) h_pass_one_btag_subjet->Fill(0.5, weight);

  // ---------------------------------------------------------------------------------------------------------
  // Distance ak4 subjets ------------------------------------------------------------------------------------
  h_distance_ak4_xcone_subjet_btag->Fill(deltaR(xcone_had_subjets[index_bjet], ak4_highest_btag), weight);
  h_difference_distance_small_mid->Fill(abs(dR_subjets[0]-dR_subjets[1]), weight);
  h_difference_distance_small_high->Fill(abs(dR_subjets[0]-dR_subjets[2]), weight);
  h_difference_distance_mid_high->Fill(abs(dR_subjets[1]-dR_subjets[2]), weight);

  // ---------------------------------------------------------------------------------------------------------
  // Method - Min --------------------------------------------------------------------------------------------
  if(Wjet_v4_min.M() != 1000){ // dummy
    h_wmass_min_mass->Fill(Wjet_v4_min.M(), weight);
    h_wjet_pt_min_mass->Fill(Wjet_v4_min.Pt(), weight);
  }

  // ---------------------------------------------------------------------------------------------------------
  // Method - Compare ----------------------------------------------------------------------------------------
  if(Wjet_v4_compare.M() != 1000){ // dummy
    h_wmass_compare->Fill(Wjet_v4_compare.M(), weight);
    h_wjet_pt_compare->Fill(Wjet_v4_compare.Pt(), weight);
  }

  // ---------------------------------------------------------------------------------------------------------
  // Method - btag matched: ----------------------------------------------------------------------------------
  double wjet_matched_pt  = Wjet_v4_match.Pt();
  double wjet_match_s1_pt = Wjet_v4_match_subjet1.Pt();
  double wjet_match_s2_pt = Wjet_v4_match_subjet2.Pt();
  double S1divW           = wjet_match_s1_pt/wjet_matched_pt;
  h_wjet_pt_match->Fill(wjet_matched_pt, weight);
  h_wjet_pt_match_subjet1->Fill(wjet_match_s1_pt, weight);
  h_wjet_pt_match_subjet2->Fill(wjet_match_s2_pt, weight);
  h_wjet_pt_match_S1mS2->Fill(wjet_match_s1_pt-wjet_match_s2_pt, weight);
  h_wjet_pt_match_S1divW->Fill(S1divW, weight);

  // 2D ------------------------------------------------------------------------
  h_wjet_pt_match_S1_S2->Fill(wjet_match_s1_pt, wjet_match_s2_pt, weight);
  h_wjet_pt_match_S1divW_W->Fill( wjet_matched_pt, S1divW, weight);

  h_wmass_match->Fill(Wjet_v4_match.M(), weight);
  if(wjet_match_s1_pt<wjet_match_s2_pt) h_wjet_bigger_pt_subjet->Fill(0.5, weight);
  else                                  h_wjet_bigger_pt_subjet->Fill(1.5, weight);

  // Bins ----------------------------------------------------------------------
  bool pt_greater, div_greater;
  if(wjet_matched_pt<300) pt_greater = false;
  else                    pt_greater = true;
  if(S1divW<0.7)         div_greater = false;
  else                   div_greater = true;

  if(!pt_greater)   h_wmass_match_ptbin_low->Fill(Wjet_v4_match.M(), weight);
  if(pt_greater)    h_wmass_match_ptbin_high->Fill(Wjet_v4_match.M(), weight);

  if(!div_greater)  h_wmass_match_divbin_low->Fill(Wjet_v4_match.M(), weight);
  if(div_greater)   h_wmass_match_divbin_high->Fill(Wjet_v4_match.M(), weight);

  if(pt_greater && div_greater)   h_wmass_match_ptdiv_hh->Fill(Wjet_v4_match.M(), weight);
  if(pt_greater && !div_greater)  h_wmass_match_ptdiv_hl->Fill(Wjet_v4_match.M(), weight);
  if(!pt_greater && div_greater)  h_wmass_match_ptdiv_lh->Fill(Wjet_v4_match.M(), weight);
  if(!pt_greater && !div_greater) h_wmass_match_ptdiv_ll->Fill(Wjet_v4_match.M(), weight);

  // ---------------------------------------------------------------------------------------------------------
  // Method - match+cut --------------------------------------------------------------------------------------
  if(debug) std::cout << "Wjets fill: btagcut" << '\n';

  if(sel_btag_min){
    h_btagcut_passed->Fill(0.5, weight);
    h_ak4_btag_high_btagcut->Fill(ak4_highest_btag.btag_DeepJet(), weight);

    h_wjet_pt_btagcut->Fill(Wjet_v4_match.Pt(), weight);
    h_wmass_btagcut->Fill(Wjet_v4_match.M(), weight);

    // Method - Combined -------------------------------------------------------
    h_wjet_combination->Fill(3.5, weight);
    h_wmass_combined_methods->Fill(Wjet_v4_match.M(), weight);
    h_wjet_pt_combined_methods->Fill(Wjet_v4_match.Pt(), weight);

    // Additional Selections ---------------------------------------------------
    if(sel_one_btag) {
      h_wmass_btagcut_one_btag->Fill(Wjet_v4_match.M(), weight);
      h_wjet_pt_btagcut_one_btag->Fill(Wjet_v4_match.Pt(), weight);
    }
    if(sel_one_subjet) {
      h_wmass_btagcut_one_subjet->Fill(Wjet_v4_match.M(), weight);
      h_wjet_pt_btagcut_one_subjet->Fill(Wjet_v4_match.Pt(), weight);
    }
    if(sel_one_btag && sel_one_subjet){
      h_wmass_btagcut_one_btag_subjet->Fill(Wjet_v4_match.M(), weight); // btag discrepancy & many close jets
      h_wjet_pt_btagcut_one_btag_subjet->Fill(Wjet_v4_match.Pt(), weight); // btag discrepancy & many close jets
    }
  }

  if(!sel_btag_min && M_min != 1000){
    h_wmass_combined_methods->Fill(Wjet_v4_min.M(), weight);
    h_wjet_pt_combined_methods->Fill(Wjet_v4_min.Pt(), weight);
  }

  // ---------------------------------------------------------------------------------------------------------
  // Selections ----------------------------------------------------------------------------------------------
  if(debug) std::cout << "Wjets fill: additional selections" << '\n';

  if(sel_one_btag){ // one btag selection
    h_wmass_one_btag->Fill(Wjet_v4_match.M(), weight);
    h_wjet_pt_one_btag->Fill(Wjet_v4_match.Pt(), weight);
  }
  if(sel_one_subjet){ // one close subjet
    h_wmass_one_subjet->Fill(Wjet_v4_match.M(), weight);
    h_wjet_pt_one_subjet->Fill(Wjet_v4_match.Pt(), weight);
  }
  if(sel_one_btag && sel_one_subjet){
    h_wmass_one_btag_subjet->Fill(Wjet_v4_match.M(), weight); // one btag selection & one close jets
    h_wjet_pt_one_btag_subjet->Fill(Wjet_v4_match.Pt(), weight); // one btag selection & one close jets
  }

  // ---------------------------------------------------------------------------------------------------------
  // difference btag -----------------------------------------------------------------------------------------
  if(debug) std::cout << "Wjets fill: btag discrepancy" << '\n';

  if(btag_v_size>1){
    h_btag_difference_high_mid->Fill(btag_v[btag_v_size-1]-btag_v[btag_v_size-2], weight);
    h_btag_percentage_high_mid->Fill(btag_v[btag_v_size-2]/btag_v[btag_v_size-1], weight);
  }

  if(sel_btag_min) {
    if(btag_v.size()>1){
      h_btag_difference_high_mid_sel_btag_min->Fill(btag_v[btag_v_size-1]-btag_v[btag_v_size-2], weight);
      h_btag_percentage_high_mid_sel_btag_min->Fill(btag_v[btag_v_size-2]/btag_v[btag_v_size-1], weight);
    }
  }
}

double RecoGenHists_xcone_topjet::get_tau32(){
  return tau32;
}

double RecoGenHists_xcone_topjet::get_mass(){
  return mass;
}
