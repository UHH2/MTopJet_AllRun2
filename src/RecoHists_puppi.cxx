#include <UHH2/MTopJet/include/RecoHists_puppi.h>


RecoHists_puppi::RecoHists_puppi(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  HadJetMass = book<TH1F>("M_jet1", "Leading Jet Mass [GeV]", 50, 0, 500);
  HadJetMass_B = book<TH1F>("M_jet1_B", "Leading Jet Mass [GeV]", 100, 0, 500);
  HadJetMass_rebin = book<TH1F>("M_jet1_", "Leading Jet Mass [GeV]", 25, 0, 500);

  LepJetMass = book<TH1F>("M_jet2", "m_{jet}", 50, 0, 500);
  HadMassLepMass = book<TH1F>("M_jet1-M_jet2+lep", "m_{jet1} - m_{jet2 + lepton}", 40, -200, 200);

  HadJetEta = book<TH1F>("eta_jet1", "#eta", 50, -5, 5);
  HadJetPhi = book<TH1F>("phi_jet1", "#phi", 50, -2*M_PI, 2*M_PI);
  LepJetEta = book<TH1F>("eta_jet2", "#eta", 50, -5, 5);
  LepJetPhi = book<TH1F>("phi_jet2", "#phi", 50, -2*M_PI, 2*M_PI);

  HadJetPT = book<TH1F>("pt_jet1", "p_{T}", 50, 0, 1000);
  LepJetPT = book<TH1F>("pt_jet2", "p_{T}", 50, 0, 1000);

  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconePuppi");

}



void RecoHists_puppi::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets and lepton ---------------------------------
  //---------------------------------------------------------------------------------------

  std::vector<TopJet> fatjets = event.get(h_fatjets);
  TLorentzVector hadjet_v4, lepjet_v4;
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
  else return;

  if(fatjets.size() < 2) return;
  if(fatjets.at(0).subjets().size() < 3) return;
  if(fatjets.at(1).subjets().size() < 3) return;

  float dR0 = deltaR(lepton, fatjets.at(0));
  float dR1 = deltaR(lepton, fatjets.at(1));


  // build final jets here
  if(dR0 > dR1){
    double px=0, py=0, pz=0, E=0;
    for(unsigned int i=0; i<fatjets.at(0).subjets().size(); i++){
      if(fatjets.at(0).subjets().at(i).pt() < 30) return;
      px +=  fatjets.at(0).subjets().at(i).v4().Px();
      py +=  fatjets.at(0).subjets().at(i).v4().Py();
      pz +=  fatjets.at(0).subjets().at(i).v4().Pz();
      E  +=  fatjets.at(0).subjets().at(i).v4().E();
    }
    hadjet_v4.SetPxPyPzE(px, py, pz, E);
    px=0;
    py=0;
    pz=0;
    E=0;
    for(unsigned int i=0; i<fatjets.at(1).subjets().size(); i++){
      px +=  fatjets.at(1).subjets().at(i).v4().Px();
      py +=  fatjets.at(1).subjets().at(i).v4().Py();
      pz +=  fatjets.at(1).subjets().at(i).v4().Pz();
      E  +=  fatjets.at(1).subjets().at(i).v4().E();
    }
    lepjet_v4.SetPxPyPzE(px, py, pz, E);
  }
  else{
    double px=0, py=0, pz=0, E=0;
    for(unsigned int i=0; i<fatjets.at(1).subjets().size(); i++){
      if(fatjets.at(1).subjets().at(i).pt() < 30) return;
      px +=  fatjets.at(1).subjets().at(i).v4().Px();
      py +=  fatjets.at(1).subjets().at(i).v4().Py();
      pz +=  fatjets.at(1).subjets().at(i).v4().Pz();
      E  +=  fatjets.at(1).subjets().at(i).v4().E();
    }
    hadjet_v4.SetPxPyPzE(px, py, pz, E);
    px=0;
    py=0;
    pz=0;
    E=0;
    for(unsigned int i=0; i<fatjets.at(0).subjets().size(); i++){
      px +=  fatjets.at(0).subjets().at(i).v4().Px();
      py +=  fatjets.at(0).subjets().at(i).v4().Py();
      pz +=  fatjets.at(0).subjets().at(i).v4().Pz();
      E  +=  fatjets.at(0).subjets().at(i).v4().E();
    }
    lepjet_v4.SetPxPyPzE(px, py, pz, E);
  }

  // do some selection to be comparable to chose
  if(hadjet_v4.Pt() < 400) return;
  if(hadjet_v4.M() < lepjet_v4.M()) return;
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;

  HadJetMass->Fill(hadjet_v4.M(), weight);
  HadJetMass_B->Fill(hadjet_v4.M(), weight);
  HadJetMass_rebin->Fill(hadjet_v4.M(), weight);
  LepJetMass->Fill(lepjet_v4.M(), weight);
  HadMassLepMass->Fill(hadjet_v4.M() - lepjet_v4.M(), weight);

  HadJetEta->Fill(hadjet_v4.Eta(), weight);
  HadJetPhi->Fill(hadjet_v4.Phi(), weight);
  LepJetEta->Fill(lepjet_v4.Eta(), weight);
  LepJetEta->Fill(lepjet_v4.Phi(), weight);


  HadJetPT->Fill(hadjet_v4.Pt(), weight);
  LepJetPT->Fill(lepjet_v4.Pt(), weight);


  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------
  hadjet_v4.Delete();
  lepjet_v4.Delete();
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

}
