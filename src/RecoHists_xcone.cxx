#include "UHH2/MTopJet/include/RecoHists_xcone.h"


RecoHists_xcone::RecoHists_xcone(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  HadJetMass = book<TH1F>("M_jet1", "M_{jet}", 50, 0, 500);
  LepJetMass = book<TH1F>("M_jet2", "M_{jet}", 50, 0, 500);
  HadMassLepMass = book<TH1F>("M_jet1-M_jet2+lep", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);
 
  HadJetPT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  LepJetPT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);

  DeltaRDiff = book<TH1F>("dR1_dR2", "dR(lepton, hadjet) - dR(lepton, lepjet)", 60, -6, -6);

  // handle for clustered jets
  h_jets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets");

}



void RecoHists_xcone::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get subjets and lepton ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<TopJet> jets = event.get(h_jets);

  // identify lepton, there has to be ==1 muon OR ==1 elec 
  Particle lepton;
  if(event.muons->size() > 0 && event.electrons->size() > 0){
    return;
  }
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }

  //---------------------------------------------------------------------------------------
  //------------- define had and lep jet (deltaR) -----------------------------------------
  //---------------------------------------------------------------------------------------
  TopJet fathadjet, fatlepjet;
  float dR1 = deltaR(lepton, jets.at(0));
  float dR2 = deltaR(lepton, jets.at(1));
  float dR_lep, dR_had;
  if(dR1 < dR2){
    fatlepjet = jets.at(0);
    fathadjet = jets.at(1);
    dR_lep = dR1;
    dR_had = dR2;
  }
  else{
    fatlepjet = jets.at(1);
    fathadjet = jets.at(0);
    dR_lep = dR2;
    dR_had = dR1;
  }

  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of subjets and combine them ------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<Jet> subjets_lep = fatlepjet.subjets();
  std::vector<Jet> subjets_had = fathadjet.subjets();
  TLorentzVector hadjet_v4, lepjet_v4;
  double pxlep=0, pylep=0, pzlep=0, Elep=0;
  for(unsigned int i=0; i < subjets_lep.size(); ++i){
    pxlep += subjets_lep.at(i).v4().Px();
    pylep += subjets_lep.at(i).v4().Py();
    pzlep += subjets_lep.at(i).v4().Pz();
    Elep += subjets_lep.at(i).v4().E();
  }
  lepjet_v4.SetPxPyPzE(pxlep, pylep, pzlep, Elep);

  double pxhad=0, pyhad=0, pzhad=0, Ehad=0;
  for(unsigned int i=0; i < subjets_had.size(); ++i){
    pxhad += subjets_had.at(i).v4().Px();
    pyhad += subjets_had.at(i).v4().Py();
    pzhad += subjets_had.at(i).v4().Pz();
    Ehad += subjets_had.at(i).v4().E();
  }
  hadjet_v4.SetPxPyPzE(pxhad, pyhad, pzhad, Ehad);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;

  HadJetMass->Fill(hadjet_v4.M(), weight);
  LepJetMass->Fill(lepjet_v4.M(), weight);
  HadMassLepMass->Fill(hadjet_v4.M() - lepjet_v4.M(), weight);
 
  HadJetPT->Fill(hadjet_v4.Pt(), weight);
  LepJetPT->Fill(lepjet_v4.Pt(), weight);

  DeltaRDiff->Fill(dR_had - dR_lep, weight);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------
  hadjet_v4.Delete();
  lepjet_v4.Delete();
  //---------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------- 

}


