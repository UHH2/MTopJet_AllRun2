#include "UHH2/MTopJet/include/GenHists_topjet.h"


GenHists_topjet::GenHists_topjet(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here
  GenJetNumber = book<TH1F>("number_jets", "number", 21, 0, 20);

  GenJet1Mass = book<TH1F>("M_jet1", "M_{jet}", 50, 0, 500);
  GenJet2Mass = book<TH1F>("M_jet2", "M_{jet}", 50, 0, 500);
  GenJet2Mass_scale = book<TH1F>("M_jet2_scale", "M_{jet}", 150, -5, 10);
  GenJet2LepMass = book<TH1F>("M_jet2_lep", "M_{jet2 + lepton}", 50, 0, 500);
  GenJet2LepMass_scale = book<TH1F>("M_jet2_lep_scale", "M_{jet2 + lepton}", 150, -5, 10);
  GenJet2MassBool = book<TH1F>("Mass_jet2_0", "is M_{jet2} == 0", 2, 0, 1);
  GenJet2LepMassBool = book<TH1F>("Mass_jet2_lep_0", "is M_{jet2 + lep} == 0", 2, 0, 1);
  Mass1Mass2 = book<TH1F>("M_jet1-M_jet2+lep", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);

  GenJet1MassJet2LepMass = book<TH2F>("M_jet1_jet2", "x=M_{jet1} y=M_{jet2_lep}", 25, 0, 500, 25, 0, 500);
 
  GenJet1PT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  GenJet2PT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);
  GenJet1Jet2PT = book<TH1F>("pt_jet1-pt_jet2", "p_{T,jet1} - p_{T,jet2}", 80, -400, 400);
  GenJet3PT = book<TH1F>("pt_jet3", "p_{T}", 50, 0, 1000);
  GenJetPT = book<TH1F>("pt_all_jets", "p_{T}", 50, 0, 1000);
  LeptonPT = book<TH1F>("pt_lepton", "p_{T}", 50, 0, 1000);

  TopHadPT = book<TH1F>("pt_tophad", "p_{T}", 20, 0, 1000);
  TopLepPT = book<TH1F>("pt_toplep", "p_{T}", 20, 0, 1000);

  GenJet2Eta = book<TH1F>("eta_jet2", "#eta", 24, -3, 3);

  deltaR_lep1_jet1 = book<TH1F>("deltaR_lep1_jet1", "#Delta R(lep1,1st Jet)", 80, 0, 4.0);
  deltaR_lep2_jet1 = book<TH1F>("deltaR_lep2_jet1", "#Delta R(lep2,1st Jet)", 80, 0, 4.0);
  deltaR_botlep_jet1 = book<TH1F>("deltaR_botlep_jet1", "#Delta R(bot_lep,1st Jet)", 80, 0, 4.0);
  deltaR_bot_jet1 = book<TH1F>("deltaR_bot_jet1", "#Delta R(bot,1st Jet)", 80, 0, 4.0);
  deltaR_q1_jet1 = book<TH1F>("deltaR_q1_jet1", "#Delta R(q1,1st Jet)", 80, 0, 4.0);
  deltaR_q2_jet1 = book<TH1F>("deltaR_q2_jet1", "#Delta R(q2,1st Jet)", 80, 0, 4.0);
  deltaR_tophad_jet1 = book<TH1F>("deltaR_tophad_jet1", "#Delta R(had Top, 1st Jet)", 80, 0, 4.0);
  deltaR_toplep_jet1 = book<TH1F>("deltaR_toplep_jet1", "#Delta R(lep Top, 1st Jet)", 80, 0, 4.0);

  deltaR_lep1_jet2 = book<TH1F>("deltaR_lep1_jet2", "#Delta R(lep1,2nd Jet)", 80, 0, 4.0);
  deltaR_lep2_jet2 = book<TH1F>("deltaR_lep2_jet2", "#Delta R(lep2,2nd Jet)", 80, 0, 4.0);
  deltaR_botlep_jet2 = book<TH1F>("deltaR_botlep_jet2", "#Delta R(bot_lep,2nd Jet)", 80, 0, 4.0);
  deltaR_bot_jet2 = book<TH1F>("deltaR_bot_jet2", "#Delta R(bot,2nd Jet)", 80, 0, 4.0);
  deltaR_q1_jet2 = book<TH1F>("deltaR_q1_jet2", "#Delta R(q1,2nd Jet)", 80, 0, 4.0);
  deltaR_q2_jet2 = book<TH1F>("deltaR_q2_jet2", "#Delta R(q2,2nd Jet)", 80, 0, 4.0);
  deltaR_tophad_jet2 = book<TH1F>("deltaR_tophad_jet2", "#Delta R(had Top, 2nd Jet)", 80, 0, 4.0);
  deltaR_toplep_jet2 = book<TH1F>("deltaR_toplep_jet2", "#Delta R(lep Top, 2nd Jet)", 80, 0, 4.0);

  // deltaR_lep1_jet3 = book<TH1F>("deltaR_lep1_jet3", "#Delta R(lep1,3rd Jet)", 80, 0, 4.0);
  // deltaR_lep2_jet3 = book<TH1F>("deltaR_lep2_jet3", "#Delta R(lep2,3rd Jet)", 80, 0, 4.0);
  // deltaR_botlep_jet3 = book<TH1F>("deltaR_botlep_jet3", "#Delta R(bot_lep,3rd Jet)", 80, 0, 4.0);
  // deltaR_bot_jet3 = book<TH1F>("deltaR_bot_jet3", "#Delta R(bot,3rd Jet)", 80, 0, 4.0);
  // deltaR_q1_jet3 = book<TH1F>("deltaR_q1_jet3", "#Delta R(q1,3rd Jet)", 80, 0, 4.0);
  // deltaR_q2_jet3 = book<TH1F>("deltaR_q2_jet3", "#Delta R(q2,3rd Jet)", 80, 0, 4.0);
  // deltaR_tophad_jet3 = book<TH1F>("deltaR_tophad_jet3", "#Delta R(had Top, 3rd Jet)", 80, 0, 4.0);
  // deltaR_toplep_jet3 = book<TH1F>("deltaR_toplep_jet3", "#Delta R(lep Top, 3rd Jet)", 80, 0, 4.0);

  // deltaR_lep1_jet4 = book<TH1F>("deltaR_lep1_jet4", "#Delta R(lep1,4th Jet)", 80, 0, 4.0);
  // deltaR_lep2_jet4 = book<TH1F>("deltaR_lep2_jet4", "#Delta R(lep2,4th Jet)", 80, 0, 4.0);
  // deltaR_botlep_jet4 = book<TH1F>("deltaR_botlep_jet4", "#Delta R(bot_lep,4th Jet)", 80, 0, 4.0);
  // deltaR_bot_jet4 = book<TH1F>("deltaR_bot_jet4", "#Delta R(bot,4th Jet)", 80, 0, 4.0);
  // deltaR_q1_jet4 = book<TH1F>("deltaR_q1_jet4", "#Delta R(q1,4th Jet)", 80, 0, 4.0);
  // deltaR_q2_jet4 = book<TH1F>("deltaR_q2_jet4", "#Delta R(q2,4th Jet)", 80, 0, 4.0);
  // deltaR_tophad_jet4 = book<TH1F>("deltaR_tophad_jet4", "#Delta R(had Top, 4th Jet)", 80, 0, 4.0);
  // deltaR_toplep_jet4 = book<TH1F>("deltaR_toplep_jet4", "#Delta R(lep Top, 4th Jet)", 80, 0, 4.0);

  // deltaR_lep1_jet5 = book<TH1F>("deltaR_lep1_jet5", "#Delta R(lep1,5th Jet)", 80, 0, 4.0);
  // deltaR_lep2_jet5 = book<TH1F>("deltaR_lep2_jet5", "#Delta R(lep2,5th Jet)", 80, 0, 4.0);
  // deltaR_botlep_jet5 = book<TH1F>("deltaR_botlep_jet5", "#Delta R(bot_lep,5th Jet)", 80, 0, 4.0);
  // deltaR_bot_jet5 = book<TH1F>("deltaR_bot_jet5", "#Delta R(bot,5th Jet)", 80, 0, 4.0);
  // deltaR_q1_jet5 = book<TH1F>("deltaR_q1_jet5", "#Delta R(q1,5th Jet)", 80, 0, 4.0);
  // deltaR_q2_jet5 = book<TH1F>("deltaR_q2_jet5", "#Delta R(q2,5th Jet)", 80, 0, 4.0);
  // deltaR_tophad_jet5 = book<TH1F>("deltaR_tophad_jet5", "#Delta R(had Top, 5th Jet)", 80, 0, 4.0);
  // deltaR_toplep_jet5 = book<TH1F>("deltaR_toplep_jet5", "#Delta R(lep Top, 5th Jet)", 80, 0, 4.0);

  // deltaR_lep1_jet6 = book<TH1F>("deltaR_lep1_jet6", "#Delta R(lep1,6th Jet)", 80, 0, 4.0);
  // deltaR_lep2_jet6 = book<TH1F>("deltaR_lep2_jet6", "#Delta R(lep2,6th Jet)", 80, 0, 4.0);
  // deltaR_botlep_jet6 = book<TH1F>("deltaR_botlep_jet6", "#Delta R(bot_lep,6th Jet)", 80, 0, 4.0);
  // deltaR_bot_jet6 = book<TH1F>("deltaR_bot_jet6", "#Delta R(bot,6th Jet)", 80, 0, 4.0);
  // deltaR_q1_jet6 = book<TH1F>("deltaR_q1_jet6", "#Delta R(q1,6th Jet)", 80, 0, 4.0);
  // deltaR_q2_jet6 = book<TH1F>("deltaR_q2_jet6", "#Delta R(q2,6th Jet)", 80, 0, 4.0);
  // deltaR_tophad_jet6 = book<TH1F>("deltaR_tophad_jet6", "#Delta R(had Top, 6th Jet)", 80, 0, 4.0);
  // deltaR_toplep_jet6 = book<TH1F>("deltaR_toplep_jet6", "#Delta R(lep Top, 6th Jet)", 80, 0, 4.0);


  // handle for TTbarGen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
  // handle for clustered jets
  h_gentopjets=ctx.get_handle<std::vector<GenTopJet>>("gentopjets");
}



void GenHists_topjet::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get first 13 gen parts -------------------------------
  //---------------------------------------------------------------------------------------
  // GenParticle* genp;
  // int id;
  // for(unsigned int i=0; i<13; ++i){
  //    genp = &(genparts->at(i));
  //    id = genp->pdgId();
  //    cout<<"ID of genparticle nr. " << i <<" : "<< id << endl;
  // }
  //---------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------- 


  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  const auto & ttbargen = event.get(h_ttbargen);
  // define all objects needed
  std::vector<GenTopJet> jets = event.get(h_gentopjets);
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4, topjet1_v4, topjet2_v4;
  GenTopJet jet1,jet2,jet3,jet4,jet5,jet6;
  if(jets.size()>0) jet1 = jets.at(0);
  if(jets.size()>1) jet2 = jets.at(1);
  // if(jets.size()>2) jet3 = jets.at(2);
  // if(jets.size()>3) jet4 = jets.at(3);
  // if(jets.size()>4) jet5 = jets.at(4);
  // if(jets.size()>5) jet6 = jets.at(5);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


  //---------------------------------------------------------------------------------------
  //--------------------------------- Matching --------------------------------------------
  //---------------------------------------------------------------------------------------
  // bool matched = false;
  // get stable particles from ttbar decay and sort them into leptonic and hadronic
  GenParticle bot, q1, q2, bot_lep, lep1, lep2, lepton; //leptons already defined above
  if(jets.size() > 0){
    if(ttbargen.IsTopHadronicDecay()){
      bot = ttbargen.bTop();
      q1 = ttbargen.Wdecay1();
      q2 = ttbargen.Wdecay2();
      bot_lep = ttbargen.bAntitop();
      lep1 = ttbargen.WMinusdecay1();
      lep2 = ttbargen.WMinusdecay2();
    }
    else if(ttbargen.IsAntiTopHadronicDecay()){
      bot = ttbargen.bAntitop();
      q1 = ttbargen.WMinusdecay1();
      q2 = ttbargen.WMinusdecay2();
      bot_lep = ttbargen.bTop();
      lep1 = ttbargen.Wdecay1();
      lep2 = ttbargen.Wdecay2();
    }
    else if(!(ttbargen.IsTopHadronicDecay() || ttbargen.IsAntiTopHadronicDecay())){
	throw runtime_error("no hadronic Top found");
    }
   //check which lep is neutrino and which is elec/muon
    if(abs(lep1.pdgId()) == 11 || abs(lep1.pdgId()) == 13){
      lepton = lep1;
    }
    else if(abs(lep2.pdgId()) == 11 || abs(lep2.pdgId()) == 13){
      lepton = lep2;
    }
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
  


  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of 2 jets and lepton -------------------=-----------------
  //---------------------------------------------------------------------------------------
   if(jets.size() > 1){
     jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
     jet2_v4.SetPxPyPzE(jets.at(1).v4().Px(), jets.at(1).v4().Py(), jets.at(1).v4().Pz(), jets.at(1).v4().E()); //v4 of first jet
     lepton1_v4.SetPxPyPzE(lepton.v4().Px(), lepton.v4().Py(), lepton.v4().Pz(), lepton.v4().E()); //v4 of lepton
     jet2_lep_v4 = jet2_v4 + lepton1_v4;
   }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;
  ////

  GenJetNumber->Fill(jets.size(), weight);
  LeptonPT->Fill(lepton1_v4.Pt(),weight);

  float pt;
  for(unsigned int i = 0; i<jets.size(); ++i){
    pt = jets.at(i).pt();
    GenJetPT->Fill(pt, weight);
  }

  if((jets.size()) > 1){
    GenJet1Mass->Fill(jet1_v4.M(),weight);
    GenJet2Mass->Fill(jet2_v4.M(),weight);
    GenJet2Mass_scale->Fill(jet2_v4.M(),weight);
    GenJet2LepMass->Fill(jet2_lep_v4.M(),weight);
    GenJet2LepMass_scale->Fill(jet2_lep_v4.M(),weight);
    GenJet1PT->Fill(jet1_v4.Pt(),weight);
    GenJet2PT->Fill(jet2_v4.Pt(),weight);
    GenJet2Eta->Fill(jet2_v4.Eta(),weight);
    Mass1Mass2->Fill(jet1_v4.M() - jet2_lep_v4.M(), weight);
    GenJet1Jet2PT->Fill(jet1_v4.Pt() - jet2_v4.Pt(),weight);
    GenJet1MassJet2LepMass->Fill(jet1_v4.M(), jet2_lep_v4.M(), weight);
  }
  if((jets.size()) > 2){
    GenJet3PT->Fill(jets.at(2).pt(),weight);
  }
  // pT of had. top
  GenParticle tophad = ttbargen.TopHad();
  float tophadpt = tophad.pt();
  TopHadPT->Fill(tophadpt, weight);
  
  // pT of lep. top
  GenParticle toplep = ttbargen.TopLep();
  float topleppt = toplep.pt();
  TopLepPT->Fill(topleppt, weight);

  // check M(jet2) == 0
  if(jets.size() > 1){
    float mass0, mass0b;
    if(jet2_v4.M() == 0) mass0 = 0.9;
    else mass0 = 0.1;
    if(jet2_lep_v4.M() == 0) mass0b = 0.9;
    else mass0b = 0.1;
    GenJet2MassBool->Fill(mass0, weight);
    GenJet2LepMassBool->Fill(mass0b, weight);
  }

  

  // delta R Hists
  if(jets.size() > 0){
    deltaR_lep1_jet1->Fill(deltaR(jet1, lep1), weight);
    deltaR_lep2_jet1->Fill(deltaR(jet1, lep2), weight);
    deltaR_botlep_jet1->Fill(deltaR(jet1, bot_lep), weight);
    deltaR_q1_jet1->Fill(deltaR(jet1, q1), weight);
    deltaR_q2_jet1->Fill(deltaR(jet1, q2), weight);
    deltaR_bot_jet1->Fill(deltaR(jet1, bot), weight);
    deltaR_tophad_jet1->Fill(deltaR(jet1, tophad), weight);
    deltaR_toplep_jet1->Fill(deltaR(jet1, toplep), weight);
  }

  if(jets.size() > 1){
    deltaR_lep1_jet2->Fill(deltaR(jet2, lep1), weight);
    deltaR_lep2_jet2->Fill(deltaR(jet2, lep2), weight);
    deltaR_botlep_jet2->Fill(deltaR(jet2, bot_lep), weight);
    deltaR_q1_jet2->Fill(deltaR(jet2, q1), weight);
    deltaR_q2_jet2->Fill(deltaR(jet2, q2), weight);
    deltaR_bot_jet2->Fill(deltaR(jet2, bot), weight);
    deltaR_tophad_jet2->Fill(deltaR(jet2, tophad), weight);
    deltaR_toplep_jet2->Fill(deltaR(jet2, toplep), weight);
  }
  
  // if(jets.size() > 2){
  //   deltaR_lep1_jet3->Fill(deltaR(jet3, lep1), weight);
  //   deltaR_lep2_jet3->Fill(deltaR(jet3, lep2), weight);
  //   deltaR_botlep_jet3->Fill(deltaR(jet3, bot_lep), weight);
  //   deltaR_q1_jet3->Fill(deltaR(jet3, q1), weight);
  //   deltaR_q2_jet3->Fill(deltaR(jet3, q2), weight);
  //   deltaR_bot_jet3->Fill(deltaR(jet3, bot), weight);
  //   deltaR_tophad_jet3->Fill(deltaR(jet3, tophad), weight);
  //   deltaR_toplep_jet3->Fill(deltaR(jet3, toplep), weight);
  // }

  // if(jets.size() > 3){
  //   deltaR_lep1_jet4->Fill(deltaR(jet4, lep1), weight);
  //   deltaR_lep2_jet4->Fill(deltaR(jet4, lep2), weight);
  //   deltaR_botlep_jet4->Fill(deltaR(jet4, bot_lep), weight);
  //   deltaR_q1_jet4->Fill(deltaR(jet4, q1), weight);
  //   deltaR_q2_jet4->Fill(deltaR(jet4, q2), weight);
  //   deltaR_bot_jet4->Fill(deltaR(jet4, bot), weight);
  //   deltaR_tophad_jet4->Fill(deltaR(jet4, tophad), weight);
  //   deltaR_toplep_jet4->Fill(deltaR(jet4, toplep), weight);
  // }

  // if(jets.size() > 4){
  //   deltaR_lep1_jet5->Fill(deltaR(jet5, lep1), weight);
  //   deltaR_lep2_jet5->Fill(deltaR(jet5, lep2), weight);
  //   deltaR_botlep_jet5->Fill(deltaR(jet5, bot_lep), weight);
  //   deltaR_q1_jet5->Fill(deltaR(jet5, q1), weight);
  //   deltaR_q2_jet5->Fill(deltaR(jet5, q2), weight);
  //   deltaR_bot_jet5->Fill(deltaR(jet5, bot), weight);
  //   deltaR_tophad_jet5->Fill(deltaR(jet5, tophad), weight);
  //   deltaR_toplep_jet5->Fill(deltaR(jet5, toplep), weight);
  // }

  // if(jets.size() > 5){
  //   deltaR_lep1_jet6->Fill(deltaR(jet6, lep1), weight);
  //   deltaR_lep2_jet6->Fill(deltaR(jet6, lep2), weight);
  //   deltaR_botlep_jet6->Fill(deltaR(jet6, bot_lep), weight);
  //   deltaR_q1_jet6->Fill(deltaR(jet6, q1), weight);
  //   deltaR_q2_jet6->Fill(deltaR(jet6, q2), weight);
  //   deltaR_bot_jet6->Fill(deltaR(jet6, bot), weight);
  //   deltaR_tophad_jet6->Fill(deltaR(jet6, tophad), weight);
  //   deltaR_toplep_jet6->Fill(deltaR(jet6, toplep), weight);
  // }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------
  jet1_v4.Delete();
  jet2_v4.Delete();
  lepton1_v4.Delete();
  jet2_lep_v4.Delete();
  //---------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------- 

}


