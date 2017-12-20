#include "UHH2/MTopJet/include/RecoHists_topjet.h"


RecoHists_topjet::RecoHists_topjet(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname): Hists(ctx, dirname){
  // book all histograms here
  RecoJetNumber = book<TH1F>("number_jets", "number", 21, 0, 20);

  RecoJet1Mass = book<TH1F>("M_jet1", "M_{jet}", 50, 0, 500);
  RecoJet1Mass_rebin = book<TH1F>("M_jet1_", "M_{jet}", 25, 0, 500);

  RecoJet2Mass = book<TH1F>("M_jet2", "M_{jet}", 50, 0, 500);
  Mass1Mass2 = book<TH1F>("M_jet1-M_jet2+lep", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);
 
  RecoJet1PT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  RecoJet2PT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);
  RecoJet1Jet2PT = book<TH1F>("pt_jet1-pt_jet2", "p_{T,jet1} - p_{T,jet2}", 80, -400, 400);
  RecoJet3PT = book<TH1F>("pt_jet3", "p_{T}", 50, 0, 1000);
  RecoJetPT = book<TH1F>("pt_all_jets", "p_{T}", 50, 0, 1000);
  LeptonPT = book<TH1F>("pt_lepton", "p_{T}", 50, 0, 1000);
  NLepton = book<TH1F>("number_lep", "number", 10, 0, 10); 
  RecoJet2Eta = book<TH1F>("eta_jet2", "#eta", 24, -3, 3);


  // handle for clustered jets
  h_jets=ctx.get_handle<std::vector<TopJet>>(jetname);

}



void RecoHists_topjet::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> jets = event.get(h_jets);
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  TopJet jet1,jet2;
  if(jets.size()>0) jet1 = jets.at(0);
  if(jets.size()>1) jet2 = jets.at(1);
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }
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

  RecoJetNumber->Fill(jets.size(), weight);
  LeptonPT->Fill(lepton1_v4.Pt(),weight);
  NLepton->Fill(event.muons->size()+event.electrons->size(), weight);
  float pt;
  for(unsigned int i = 0; i<jets.size(); ++i){
    pt = jets.at(i).pt();
    RecoJetPT->Fill(pt, weight);
  }

  if((jets.size()) > 1){
    RecoJet1Mass->Fill(jet1_v4.M(),weight);
    RecoJet1Mass_rebin->Fill(jet1_v4.M(),weight);
    RecoJet2Mass->Fill(jet2_v4.M(),weight);
    RecoJet1PT->Fill(jet1_v4.Pt(),weight);
    RecoJet2PT->Fill(jet2_v4.Pt(),weight);
    RecoJet2Eta->Fill(jet2_v4.Eta(),weight);
    Mass1Mass2->Fill(jet1_v4.M() - jet2_lep_v4.M(), weight);
    RecoJet1Jet2PT->Fill(jet1_v4.Pt() - jet2_v4.Pt(),weight);
  }
  if((jets.size()) > 2){
    RecoJet3PT->Fill(jets.at(2).pt(),weight);
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


