#include "UHH2/MTopJet/include/RecoHists_xcone.h"


RecoHists_xcone::RecoHists_xcone(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  HadJetMass = book<TH1F>("M_jet1", "M_{jet}", 50, 0, 500);
  HadJetMass_rebin = book<TH1F>("M_jet1_", "M_{jet}", 25, 0, 500);
  LepJetMass = book<TH1F>("M_jet2", "M_{jet}", 50, 0, 500);
  HadMassLepMass = book<TH1F>("M_jet1-M_jet2+lep", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);
 
  HadJetPT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  LepJetPT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);

  number_hadjet = book<TH1F>("number_hadjet", "number", 10, 0, 10);
  number_lepjet = book<TH1F>("number_lepjet", "number", 10, 0, 10);

  // DeltaRDiff = book<TH1F>("dR1_dR2", "dR(lepton, hadjet) - dR(lepton, lepjet)", 60, -6, -6);

  // handle for clustered jets
  h_hadjets=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
  h_lepjets=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined");

}



void RecoHists_xcone::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets and lepton ---------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<Jet> hadjets = event.get(h_hadjets);
  std::vector<Jet> lepjets = event.get(h_lepjets);

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
  if(hadjets.size() == 0 || lepjets.size() == 0) return;

  float dR_had = deltaR(lepton, hadjets.at(0));
  float dR_lep = deltaR(lepton, lepjets.at(0));
  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of Combined Jets ---------- ------------------------------
  //---------------------------------------------------------------------------------------
  TLorentzVector hadjet_v4, lepjet_v4;
  double pxlep, pylep, pzlep, Elep;
  pxlep = lepjets.at(0).v4().Px();
  pylep = lepjets.at(0).v4().Py();
  pzlep = lepjets.at(0).v4().Pz();
  Elep = lepjets.at(0).v4().E();
  lepjet_v4.SetPxPyPzE(pxlep, pylep, pzlep, Elep);

  double pxhad, pyhad, pzhad, Ehad;
  pxhad = hadjets.at(0).v4().Px();
  pyhad = hadjets.at(0).v4().Py();
  pzhad = hadjets.at(0).v4().Pz();
  Ehad = hadjets.at(0).v4().E();
  hadjet_v4.SetPxPyPzE(pxhad, pyhad, pzhad, Ehad);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;
  number_hadjet->Fill(hadjets.size(), weight); // just for checks
  number_lepjet->Fill(lepjets.size(), weight); // just for checks
  HadJetMass->Fill(hadjet_v4.M(), weight);
  HadJetMass_rebin->Fill(hadjet_v4.M(), weight);
  LepJetMass->Fill(lepjet_v4.M(), weight);
  HadMassLepMass->Fill(hadjet_v4.M() - lepjet_v4.M(), weight);
 
  HadJetPT->Fill(hadjet_v4.Pt(), weight);
  LepJetPT->Fill(lepjet_v4.Pt(), weight);

  // DeltaRDiff->Fill(dR_had - dR_lep, weight);
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


