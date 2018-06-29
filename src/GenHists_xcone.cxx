#include "UHH2/MTopJet/include/GenHists_xcone.h"


GenHists_xcone::GenHists_xcone(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // book all histograms here

  Mass_HadJet23 = book<TH1F>("Mass_HadJet23", "M_{jet}", 50, 0, 500);
  Mass_HadJet23_rebin = book<TH1F>("Mass_HadJet23_rebin", "M_{jet}", 25, 0, 500);
  Mass_LepJet23 = book<TH1F>("Mass_LepJet23", "M_{jet}", 50, 0, 500);
  MassDiff23 = book<TH1F>("MassDiff23", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);

  PT_HadJet23 = book<TH1F>("PT_HadJet23", "p_{T}", 50, 0, 1000);
  PT_LepJet23 = book<TH1F>("PT_LepJet23", "p_{T}", 50, 0, 1000);

  number_HadJet23 = book<TH1F>("number_HadJet23", "number", 10, 0, 10);
  number_LepJet23 = book<TH1F>("number_LepJet23", "number", 10, 0, 10);

  Mass_HadJet33 = book<TH1F>("Mass_HadJet33", "M_{jet}", 50, 0, 500);
  Mass_HadJet33_B = book<TH1F>("Mass_HadJet33_B", "M_{jet}", 100, 0, 500);
  Mass_HadJet33_C = book<TH1F>("Mass_HadJet33_C", "M_{jet}", 500, 0, 500);
  Mass_HadJet33_rebin = book<TH1F>("Mass_HadJet33_rebin", "M_{jet}", 25, 0, 500);

  Float_t xbins[] = {0, 127, 142, 157, 172, 187, 203, 218, 240, 270, 300, 500};
  Mass_HadJet33_unfold = book<TH1F>("Mass_HadJet33_unfold", "M_{jet}", 11, xbins);

  Mass_LepJet33 = book<TH1F>("Mass_LepJet33", "M_{jet}", 50, 0, 500);
  MassDiff33 = book<TH1F>("MassDiff33", "M_{jet1} - M_{jet2 + lepton}", 40, -200, 200);

  PT_HadJet33 = book<TH1F>("PT_HadJet33", "p_{T}", 50, 0, 1000);
  PT_LepJet33 = book<TH1F>("PT_LepJet33", "p_{T}", 50, 0, 1000);

  number_HadJet33 = book<TH1F>("number_HadJet33", "number", 10, 0, 10);
  number_LepJet33 = book<TH1F>("number_LepJet33", "number", 10, 0, 10);

  DeltaR_23_33 = book<TH1F>("DeltaR_23_33", "#Delta R(jet2^{2+3}, jet2^{33})", 60, -3, 3);
  DeltaMass_23_33 = book<TH1F>("DeltaMass_23_33", "#Delta Mass(jet2^{2+3}, jet2^{33})", 100, -100, 100);
  DeltaPT_23_33 = book<TH1F>("DeltaPT_23_33", "#Delta p_{T}(jet2^{2+3}, jet2^{33})", 100, -100, 100);

  DeltaR_jet2_23_33 = book<TH1F>("DeltaR_jet2_23_33", "#Delta R(jet2^{2+3}, jet2^{33})", 60, -3, 3);
  DeltaMass_jet2_23_33 = book<TH1F>("DeltaMass_jet2_23_33", "#Delta Mass(jet2^{2+3}, jet2^{33})", 100, -100, 100);
  DeltaPT_jet2_23_33 = book<TH1F>("DeltaPT_jet2_23_33", "#Delta p_{T}(jet2^{2+3}, jet2^{33})", 100, -100, 100);

  SoftDropMass = book<TH1F>("SoftDropMass", "M_{jet1}^{SD}", 50, 0, 500);

  // handle for clustered jets (2+3 und 3+3)
  h_hadjets33=ctx.get_handle<std::vector<Particle>>("GEN_XCone33_had_Combined");
  h_lepjets33=ctx.get_handle<std::vector<Particle>>("GEN_XCone33_lep_Combined");

  h_hadjets23=ctx.get_handle<std::vector<Particle>>("GEN_XCone23_had_Combined");
  h_lepjets23=ctx.get_handle<std::vector<Particle>>("GEN_XCone23_lep_Combined");

  //h_softdrop=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets_softdrop");

}



void GenHists_xcone::fill(const Event & event){

  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets  -------------------------------------------
  //---------------------------------------------------------------------------------------
  std::vector<Particle> hadjets23 = event.get(h_hadjets23);
  std::vector<Particle> lepjets23 = event.get(h_lepjets23);
  std::vector<Particle> hadjets33 = event.get(h_hadjets33);
  std::vector<Particle> lepjets33 = event.get(h_lepjets33);
  //std::vector<GenTopJet> softdrop = event.get(h_softdrop);

  if(hadjets23.size() == 0 || lepjets23.size() == 0) return;
  if(hadjets33.size() == 0 || lepjets33.size() == 0) return;

  Particle had23 = hadjets23.at(0);
  Particle lep23 = lepjets23.at(0);
  Particle had33 = hadjets33.at(0);
  Particle lep33 = lepjets33.at(0);

  int had_nr;
  //if(uhh2::deltaR(had33, softdrop.at(0)) < uhh2::deltaR(had33, softdrop.at(1))) had_nr = 0;
  //else had_nr = 1;
  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;

  Mass_HadJet23->Fill(had23.v4().M(), weight);
  Mass_HadJet23_rebin->Fill(had23.v4().M(), weight);
  Mass_LepJet23->Fill(lep23.v4().M(), weight);
  MassDiff23->Fill(had23.v4().M() - lep23.v4().M(), weight);

  PT_HadJet23->Fill(had23.pt(), weight);
  PT_LepJet23->Fill(lep23.pt(), weight);

  number_HadJet23->Fill(hadjets23.size(), weight);
  number_LepJet23->Fill(lepjets23.size(), weight);

  Mass_HadJet33->Fill(had33.v4().M(), weight);
  Mass_HadJet33_rebin->Fill(had33.v4().M(), weight);
  Mass_HadJet33_unfold->Fill(had33.v4().M(), weight);
  Mass_HadJet33_B->Fill(had33.v4().M(), weight);
  Mass_HadJet33_C->Fill(had33.v4().M(), weight);


  Mass_LepJet33->Fill(lep33.v4().M(), weight);
  MassDiff33->Fill(had33.v4().M() - lep33.v4().M(), weight);

  PT_HadJet33->Fill(had33.pt(), weight);
  PT_LepJet33->Fill(lep33.pt(), weight);

  number_HadJet33->Fill(hadjets33.size(), weight);
  number_LepJet33->Fill(lepjets33.size(), weight);


  DeltaR_23_33->Fill(abs(deltaR(had23, had33)), weight);
  DeltaMass_23_33->Fill(had23.v4().M() - had33.v4().M() ,weight);
  DeltaPT_23_33->Fill(had23.pt() - had33.pt() ,weight);

  DeltaR_jet2_23_33->Fill(abs(deltaR(lep23, lep33)), weight);
  DeltaMass_jet2_23_33->Fill(lep23.v4().M() - lep33.v4().M() ,weight);
  DeltaPT_jet2_23_33->Fill(lep23.pt() - lep33.pt() ,weight);

  //SoftDropMass->Fill(softdrop.at(had_nr).v4().M(), weight);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

}
