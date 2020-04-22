#include <UHH2/MTopJet/include/RecoGenHists_xcone_topjet.h>


RecoGenHists_xcone_topjet::RecoGenHists_xcone_topjet(uhh2::Context & ctx, const std::string & dirname, bool isTTbar, double masscut): Hists(ctx, dirname), isTTbar_(isTTbar), masscut_(masscut){
  // Reco Histograms -----------------------------------------------------------
  HadJetPt = book<TH1F>("pt_hadjet", "p_{T}", 50, 0, 1000);
  HadJetMass = book<TH1F>("M_hadjet", "m_{jet} [GeV]", 50, 0, 500);
  HadJetEta = book<TH1F>("eta_ak8_had", "#eta", 50, -5, 5);
  HadJetPhi = book<TH1F>("phi_ak8_had", "#phi", 50, -2*M_PI, 2*M_PI);

  distance_hadx_hadak = book<TH1F>("deltaR_hadx_hadak", "dR", 60, 0, 3);

  XCone_HadJetEta = book<TH1F>("eta_xcone_had", "#eta", 50, -5, 5);
  XCone_HadJetPhi = book<TH1F>("phi_xcone_had", "#phi", 50, -2*M_PI, 2*M_PI);

  HadJetTau1 = book<TH1F>("tau1_hadjet", "#tau_{1}", 50, 0, 1);
  HadJetTau2 = book<TH1F>("tau2_hadjet", "#tau_{2}", 50, 0, 1);
  HadJetTau3 = book<TH1F>("tau3_hadjet", "#tau_{3}", 50, 0, 1);

  HadJetTau32 = book<TH1F>("tau32_hadjet", "#tau_{3}/#tau_{2}", 50, 0, 1);
  HadJetTau23 = book<TH1F>("tau23_hadjet", "#tau_{2}/#tau_{3}", 100, 0, 2);

  number_ak8_jets = book<TH1F>("number_ak8_jets", "number AK8 jets", 10, 0, 10);
  number_xcone_jets = book<TH1F>("number_xconejets", "number XCone Jets", 10, 0, 10);
  number_ak8_had_jets = book<TH1F>("number_ak8_had_jets", "number_ak8_had_jets", 10, 0, 10);

  // Gen Histograms ------------------------------------------------------------
  if(isTTbar){
    Number_matched_all = book<TH1F>("events_matched_all", "matched with all particles", 2, 0, 2);
    Number_matched_top = book<TH1F>("events_matched_top", "matched with top", 2, 0, 2);
    Number_matched_q1 = book<TH1F>("events_matched_q1", "matched with q1", 2, 0, 2);
    Number_matched_q2 = book<TH1F>("events_matched_q2", "matched with q2", 2, 0, 2);
    Number_matched_bottom = book<TH1F>("events_matched_bottom", "matched with bottom", 2, 0, 2);

    HadJetMass_fullymerged = book<TH1F>("M_hadjet_fullymerged", "m_{jet} [GeV]", 50, 0, 500);
    HadJetMass_semimerged = book<TH1F>("M_hadjet_semimerged", "m_{jet} [GeV]", 50, 0, 500);
    HadJetMass_notmerged = book<TH1F>("M_hadjet_notmerged", "m_{jet} [GeV]", 50, 0, 500);

    HadJetTau32_fullymerged = book<TH1F>("tau32_hadjet_fullymerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
    HadJetTau32_semimerged = book<TH1F>("tau32_hadjet_semimerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
    HadJetTau32_notmerged = book<TH1F>("tau32_hadjet_notmerged", "#tau_{3}/#tau_{2}", 50, 0, 1);
  }

  // handle for clustered jets -------------------------------------------------
  h_recjets=ctx.get_handle<std::vector<TopJet>>("topjets"); // jetsAk8CHSSubstructure_SoftDropCHS
  h_xcone=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
}

void RecoGenHists_xcone_topjet::fill(const Event & event){
  /*
  ██████  ███████  ██████
  ██   ██ ██      ██
  ██████  █████   ██
  ██   ██ ██      ██
  ██   ██ ███████  ██████
  */

  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> rec_jets = event.get(h_recjets);
  std::vector<TopJet> xcone_jets = event.get(h_xcone);

  // Select AK8 Jet in range of XCone size.
  // Looking at each jet seperatly by going through with a for-Loop
  TLorentzVector rec_v4, xcone_v4;

  // TLorentzVector rec_jet1_v4, rec_jet2_v4, gen_jet1_v4, gen_jet2_v4;
  std::vector<TLorentzVector> had_v4_vector;
  std::vector<TopJet> had_ak8_vector;
  TopJet had_ak8;

  xcone_v4.SetPxPyPzE(xcone_jets.at(0).v4().Px(), xcone_jets.at(0).v4().Py(), xcone_jets.at(0).v4().Pz(), xcone_jets.at(0).v4().E());

  for(unsigned int i=0; i<rec_jets.size(); i++){
    rec_v4.SetPxPyPzE(rec_jets.at(i).v4().Px(), rec_jets.at(i).v4().Py(), rec_jets.at(i).v4().Pz(), rec_jets.at(i).v4().E());

    if(deltaR(xcone_jets.at(0), rec_jets.at(i)) < 0.8){
      had_v4_vector.push_back(rec_v4);
      had_ak8_vector.push_back(rec_jets.at(i));
    }
  }

  // Selecting the hadronic AK8 jet
  TLorentzVector had_v4;
  double jet1_pt = 0;
  double jet2_pt = 0;
  if(had_v4_vector.size()>=1){
    if(had_v4_vector.size()==1){
      had_v4 = had_v4_vector.at(0);
      had_ak8 = had_ak8_vector.at(0);
    }
    else{
      jet1_pt = had_v4_vector.at(0).Pt();
      jet2_pt = had_v4_vector.at(1).Pt();
      if(jet1_pt>jet2_pt){
        had_v4 = had_v4_vector.at(0);
        had_ak8 = had_ak8_vector.at(0);
      }
      else{
        had_v4 = had_v4_vector.at(1);
        had_ak8 = had_ak8_vector.at(1);
      }
    }

    double tau_32_ak = had_ak8.tau3()/had_ak8.tau2();
    bool pass_masscut = false;
    if(had_v4.M()>masscut_) pass_masscut = true;
    //---------------------------------------------------------------------------------------
    //---------------------------- Fill Reco Hists here -------------------------------------
    //---------------------------------------------------------------------------------------
    if(pass_masscut){
      // get weight
      double weight = event.weight;

      HadJetPt->Fill(had_v4.Pt(), weight);
      HadJetMass->Fill(had_v4.M(), weight);
      HadJetEta->Fill(had_v4.Eta(), weight);
      HadJetPhi->Fill(had_v4.Phi(), weight);

      number_ak8_had_jets->Fill(had_v4_vector.size(), 1);

      XCone_HadJetEta->Fill(xcone_jets.at(0).eta(), weight);
      XCone_HadJetPhi->Fill(xcone_jets.at(0).phi(), weight);

      HadJetTau1->Fill(had_ak8.tau1(), weight);
      HadJetTau2->Fill(had_ak8.tau2(), weight);
      HadJetTau3->Fill(had_ak8.tau3(), weight);
      HadJetTau32->Fill(tau_32_ak, weight);
      HadJetTau23->Fill(1/tau_32_ak, weight);

      distance_hadx_hadak->Fill(deltaR(xcone_jets.at(0), had_ak8), weight);

      number_ak8_jets->Fill(rec_jets.size(), weight);
      number_xcone_jets->Fill(xcone_jets.size(), weight);

      /*
      .██████  ███████ ███    ██
      ██       ██      ████   ██
      ██   ███ █████   ██ ██  ██
      ██    ██ ██      ██  ██ ██
      .██████  ███████ ██   ████
      */

      //---------------------------------------------------------------------------------------
      //--------------------- Matching GenParticles -------------------------------------------
      //---------------------------------------------------------------------------------------

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
        if(deltaR(q1, had_ak8) < 0.8) matched_q1 = true;
        if(deltaR(q2, had_ak8) < 0.8) matched_q2 = true;
        if(deltaR(bottom, had_ak8) < 0.8) matched_bottom = true;

        if(matched_q1 && matched_q2 && matched_bottom) matched = true;

        if((matched_q1 && !matched_q2 && matched_bottom) ||
        (!matched_q1 && matched_q2 && matched_bottom)) semimatched = true;


        if(!matched_q1 && !matched_q2 && !matched_bottom) notmatched = true;

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

        if(semimatched){
          HadJetTau32_semimerged->Fill(tau_32_ak, weight);
          HadJetMass_semimerged->Fill(had_v4.M(), weight);
        }

        if(notmatched){
          HadJetTau32_notmerged->Fill(tau_32_ak, weight);
          HadJetMass_notmerged->Fill(had_v4.M(), weight);
        }

        if(matched){
          Number_matched_all->Fill(1, weight);
          HadJetMass_fullymerged->Fill(had_v4.M(), weight);
          HadJetTau32_fullymerged->Fill(tau_32_ak, weight);
        }
        else Number_matched_all->Fill(0.1, weight);
      }
    }
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

}
