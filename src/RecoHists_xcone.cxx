#include "UHH2/MTopJet/include/RecoHists_xcone.h"


RecoHists_xcone::RecoHists_xcone(uhh2::Context & ctx, const std::string & dirname, const std::string & type): Hists(ctx, dirname){
  // book all histograms here

  HadJetMass = book<TH1F>("M_jet1", "Leading Jet Mass [GeV]", 50, 0, 500);
  HadJetMass_rebin = book<TH1F>("M_jet1_", "Leading Jet Mass [GeV]", 25, 0, 500);
  LepJetMass = book<TH1F>("M_jet2", "m_{jet}", 50, 0, 500);
  HadMassLepMass = book<TH1F>("M_jet1-M_jet2+lep", "m_{jet1} - m_{jet2 + lepton}", 40, -200, 200);
 
  HadJetEta = book<TH1F>("eta_jet1", "#eta", 50, -5, 5);
  HadJetPhi = book<TH1F>("phi_jet1", "#phi", 50, -2*M_PI, 2*M_PI);
  LepJetEta = book<TH1F>("eta_jet2", "#eta", 50, -5, 5);
  LepJetPhi = book<TH1F>("phi_jet2", "#phi", 50, -2*M_PI, 2*M_PI);


  SoftdropMass_had = book<TH1F>("SoftdropMass_had", "Soft Drop Mass [GeV]", 25, 0, 500);
  SoftdropMass_Sel = book<TH1F>("SoftdropMass_Sel", "Soft Drop Mass [GeV]", 25, 0, 500);
  SoftdropMass_lep = book<TH1F>("SoftdropMass_lep", "m_{fat lep jet}", 25, 0, 500);

  HadJetPT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  LepJetPT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);

  RhoA_fat = book<TH1F>("RhoA_fat", "#rho A_{fat}", 50, 0, 100);
  RhoA = book<TH1F>("RhoA", "#rho A_{final}", 50, 0, 100);
  RhoA_diff = book<TH1F>("RhoA_diff", "#rho (A_{fat} - A_{final})", 50, -50, 50);
  E_diff = book<TH1F>("E_diff", "E_{fat} - E_{final}", 50, -50, 50);


  FatJetPT_had = book<TH1F>("FatJetPT_had", "p_{T, fat had jet}", 50, 0, 1000);
  FatJetPT_lep = book<TH1F>("FatJetPT_lep", "p_{T, fat lep jet}", 50, 0, 1000);

  FatJetPTDiff_had = book<TH1F>("FatJetPTDiff_had", "p_{T, fat} - p_{T, jet1}", 50, -50, 200);
  FatJetMassDiff_had = book<TH1F>("FatJetMassDiff_had", "M_{fat} - m_{jet1}", 50, -200, 200);
  FatJetPTDiff_lep = book<TH1F>("FatJetPTDiff_lep", "p_{T, fat} - p_{T, jet1}", 50, -50, 200);
  FatJetMassDiff_lep = book<TH1F>("FatJetMassDiff_lep", "M_{fat} - m_{jet1}", 50, -200, 200);

  number_hadjet = book<TH1F>("number_hadjet", "number", 10, 0, 10);
  number_lepjet = book<TH1F>("number_lepjet", "number", 10, 0, 10);

  Mass_Vertices = book<TH2F>("Mass_Vertices", "x=Pile-up y=m_{jet}", 50, 0, 50, 50, 0, 500);

  // DeltaRDiff = book<TH1F>("dR1_dR2", "dR(lepton, hadjet) - dR(lepton, lepjet)", 60, -6, -6);

  // handle for clustered jets
  if(type == "jec"){
    h_hadjets=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
    h_lepjets=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined");
    h_fatjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets");
  }
  else if(type == "raw"){
    h_hadjets=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined_noJEC");
    h_lepjets=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined_noJEC");
    h_fatjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets_noJEC");
  }
  else if(type == "cor"){
    h_hadjets=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined_Corrected");
    h_lepjets=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined_Corrected");
    h_fatjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets_Corrected");
  }
}



void RecoHists_xcone::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets and lepton ---------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<Jet> hadjets = event.get(h_hadjets);
  std::vector<Jet> lepjets = event.get(h_lepjets);
  std::vector<TopJet> fatjets = event.get(h_fatjets);
  double rho = event.rho;
  // Particle lepton;
  // if(event.muons->size() > 0 && event.electrons->size() > 0){
  //   return;
  // }
  // if(event.muons->size() > 0){
  //   lepton = event.muons->at(0);
  // }
  // else if(event.electrons->size() > 0){
  //   lepton = event.electrons->at(0);
  // }
  // if(hadjets.size() == 0 || lepjets.size() == 0) return;

  // float dR_had = deltaR(lepton, hadjets.at(0));
  // float dR_lep = deltaR(lepton, lepjets.at(0));

  // get had jet from fat jets for softdrop mass
  int nr_hadjet = 0;
  int nr_lepjet = 1;
  if(deltaR(lepjets.at(0), fatjets.at(0)) < deltaR(hadjets.at(0), fatjets.at(0))){
    nr_hadjet = 1;
    nr_lepjet = 0;
  }
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

  HadJetEta->Fill(hadjet_v4.Eta(), weight);
  HadJetPhi->Fill(hadjet_v4.Phi(), weight);
  LepJetEta->Fill(lepjet_v4.Eta(), weight);
  LepJetEta->Fill(lepjet_v4.Phi(), weight);

  Mass_Vertices->Fill(event.pvs->size(), hadjet_v4.M(), weight);

  SoftdropMass_had->Fill(fatjets.at(nr_hadjet).softdropmass(), weight);
  SoftdropMass_lep->Fill(fatjets.at(nr_lepjet).softdropmass(), weight);

  HadJetPT->Fill(hadjet_v4.Pt(), weight);
  LepJetPT->Fill(lepjet_v4.Pt(), weight);

  FatJetPT_had->Fill(fatjets.at(nr_hadjet).pt(), weight);
  FatJetPTDiff_had->Fill((fatjets.at(nr_hadjet).pt() - hadjet_v4.Pt()), weight);
  FatJetMassDiff_had->Fill((fatjets.at(nr_hadjet).softdropmass() - hadjet_v4.M()), weight);

  FatJetPT_lep->Fill(fatjets.at(nr_lepjet).pt(), weight);
  FatJetPTDiff_lep->Fill((fatjets.at(nr_lepjet).pt() - lepjet_v4.Pt()), weight);
  FatJetMassDiff_lep->Fill((fatjets.at(nr_lepjet).softdropmass() - lepjet_v4.M()), weight);
  if(fatjets.at(nr_hadjet).pt() > 400 && fatjets.at(nr_lepjet).softdropmass() < fatjets.at(nr_hadjet).softdropmass())SoftdropMass_Sel->Fill(fatjets.at(nr_hadjet).softdropmass(), weight);

  double rhoa_fat = rho * M_PI * 1.2 * 1.2;
  double rhoa, area = 0;
  for(unsigned int i = 0; i < fatjets.at(nr_hadjet).subjets().size(); i++){
    area += fatjets.at(nr_hadjet).subjets().at(i).jetArea();
  }
  rhoa = rho * area;

  RhoA_fat->Fill(rhoa_fat, weight);
  RhoA->Fill(rhoa, weight);
  RhoA_diff->Fill(rhoa_fat - rhoa, weight);
  E_diff->Fill((fatjets.at(nr_hadjet).v4().E() - hadjet_v4.E()), weight);

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


