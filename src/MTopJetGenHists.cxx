#include "UHH2/MTopJet/include/MTopJetGenHists.h"
#include "UHH2/MTopJet/include/JetCluster.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/PFParticle.h"

#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

MTopJetGenHists::MTopJetGenHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname): Hists(ctx, dirname){
  // book all histograms here

  GenJet1Mass = book<TH1F>("GenJet Mass", "M_{jet}", 50, 0, 500);
  GenJet1PT = book<TH1F>("p_{T} of 1st Jet", "p_{T}", 50, 0, 1000);
  GenJet2PT = book<TH1F>("p_{T} of 2nd Jet", "p_{T}", 50, 0, 1000);
  LeptonPT = book<TH1F>("p_{T} of Lepton", "p_{T}", 50, 0, 1000);

  TopHadPT = book<TH1F>("p_{T} of hadronic Top", "p_{T}", 20, 0, 1000);
  TopLepPT = book<TH1F>("p_{T} of leptonic Top", "p_{T}", 20, 0, 1000);

  deltaR_lep1_jet1 = book<TH1F>("deltaR_lep1_jet1", "#Delta R(lep1,1st Jet)", 80, 0, 4.0);
  deltaR_lep2_jet1 = book<TH1F>("deltaR_lep2_jet1", "#Delta R(lep2,1st Jet)", 80, 0, 4.0);
  deltaR_lep1_jet2 = book<TH1F>("deltaR_lep1_jet2", "#Delta R(lep1,2nd Jet)", 80, 0, 4.0);
  deltaR_lep2_jet2 = book<TH1F>("deltaR_lep2_jet2", "#Delta R(lep2,2nd Jet)", 80, 0, 4.0);
  deltaR_bot_lep_jet1 = book<TH1F>("deltaR_bot_lep_jet1", "#Delta R(bot_lep,1st Jet)", 80, 0, 4.0);
  deltaR_bot_jet1 = book<TH1F>("deltaR_bot_jet1", "#Delta R(bot,1st Jet)", 80, 0, 4.0);
  deltaR_q1_jet1 = book<TH1F>("deltaR_q1_jet1", "#Delta R(q1,1st Jet)", 80, 0, 4.0);
  deltaR_q2_jet1 = book<TH1F>("deltaR_q2_jet1", "#Delta R(q2,1st Jet)", 80, 0, 4.0);
  deltaR_tophad_jet1 = book<TH1F>("deltaR_tophad_jet1", "#Delta R(had Top, 1st Jet)", 80, 0, 4.0);
  deltaR_toplep_jet1 = book<TH1F>("deltaR_toplep_jet1", "#Delta R(lep Top, 1st Jet)", 80, 0, 4.0);

  // handle for PF particles
  // h_pfpart=ctx.get_handle<vector<PFParticle>>("PFParticles");

  // handle for TTbarGen class
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
  // handle for clustered jets
  h_jets=ctx.get_handle<std::vector<Jet>>(jetname);

}



void MTopJetGenHists::fill(const Event & event){

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
  //--------------------------------- Cluster Reco Jets -----------------------------------
  //---------------------------------------------------------------------------------------
  // std::vector<PFParticle> pfparts = event.get(h_pfpart);
  // JetCluster* jetc_reco=new JetCluster();
  // std::vector<fastjet::PseudoJet> reco_ak06; 
  // reco_ak06 = jetc_reco->get_recojets(&pfparts, JetCluster::e_akt, 0.8, 200);
  //---------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------- 
 

  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  const auto & ttbargen = event.get(h_ttbargen);
 // define all objects needed
  std::vector<Jet> jets = event.get(h_jets);
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  Jet jet1,jet2;
  if(jets.size()>0) jet1 = jets.at(0);
  if(jets.size()>1) jet2 = jets.at(1);
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
   if(jets.size() == 2){
    jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(jets.at(1).pt(),jets.at(1).eta(),jets.at(1).phi(),jets.at(1).energy()); //v4 of first jet
    lepton1_v4.SetPtEtaPhiE(lepton.pt(),lepton.eta(),lepton.phi(),lepton.energy()); //v4 of lepton
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;
  ////

  if((jets.size())==2){
    GenJet1Mass->Fill(jet1_v4.M(),weight);
    GenJet1PT->Fill(jet1_v4.Pt(),weight);
    GenJet2PT->Fill(jet2_v4.Pt(),weight);
  }
  LeptonPT->Fill(lepton1_v4.Pt(),weight);

  // pT of had. top
  GenParticle tophad = ttbargen.TopHad();
  float tophadpt = tophad.pt();
  TopHadPT->Fill(tophadpt, weight);

  // pT of lep. top
  GenParticle toplep = ttbargen.TopLep();
  float topleppt = toplep.pt();
  TopLepPT->Fill(topleppt, weight);

  // delta R Hists
  if(jets.size() == 2){
    deltaR_lep1_jet1->Fill(deltaR(jet1, lep1), weight);
    deltaR_lep2_jet1->Fill(deltaR(jet1, lep2), weight);
    deltaR_lep1_jet2->Fill(deltaR(jet2, lep1), weight);
    deltaR_lep2_jet2->Fill(deltaR(jet2, lep2), weight);
    deltaR_bot_lep_jet1->Fill(deltaR(jet1, bot_lep), weight);
    deltaR_q1_jet1->Fill(deltaR(jet1, q1), weight);
    deltaR_q2_jet1->Fill(deltaR(jet1, q2), weight);
    deltaR_bot_jet1->Fill(deltaR(jet1, bot), weight);
    deltaR_tophad_jet1->Fill(deltaR(jet1, tophad), weight);
    deltaR_toplep_jet1->Fill(deltaR(jet1, toplep), weight);
  }
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


