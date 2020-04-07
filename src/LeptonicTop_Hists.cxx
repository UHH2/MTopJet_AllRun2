#include <UHH2/MTopJet/include/LeptonicTop_Hists.h>

LeptonicTop_Hists::LeptonicTop_Hists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  muon_pt = book<TH1F>("muon_pt", "p_{T}^{muon}", 50, 0, 1000);
  muon_perp_full_momentum = book<TH1F>("muon_momentum_perp_full", "#||{#vec{p}_{#perp}^{muon}}", 200, 0, 200);
  muon_perp_lead_momentum = book<TH1F>("muon_momentum_perp_lead", "#||{#vec{p}_{#perp}^{muon}}", 200, 0, 200);
  muon_para_full_momentum = book<TH1F>("muon_momentum_para_full", "#||{#vec{p}_{#parallel}^{muon}}", 100, 0, 1000);
  muon_para_lead_momentum = book<TH1F>("muon_momentum_para_lead", "#||{#vec{p}_{#parallel}^{muon}}", 100, 0, 1000);

  muon_boost_perp_full_momentum = book<TH1F>("boosted_muon_momentum_perp_full", "#||{#vec{p}_{#perp}^{muon*}}", 200, 0, 200);
  muon_boost_perp_lead_momentum = book<TH1F>("boosted_muon_momentum_perp_lead", "#||{#vec{p}_{#perp}^{muon*}}", 200, 0, 200);
  muon_boost_para_full_momentum = book<TH1F>("boosted_muon_momentum_para_full", "#||{#vec{p}_{#parallel}^{muon*}}", 100, 0, 1000);
  muon_boost_para_lead_momentum = book<TH1F>("boosted_muon_momentum_para_lead", "#||{#vec{p}_{#parallel}^{muon*}}", 100, 0, 1000);

  top_boost_full_pt  = book<TH1F>("boosted_top_pt_full", "p_{T}^{top}", 200, 0, 200); // Should be 0
  top_boost_lead_pt  = book<TH1F>("boosted_top_pt_lead", "p_{T}^{top}", 200, 0, 200); // Should be 0
  muon_boost_full_pt = book<TH1F>("boosted_muon_pt_full", "p_{T}^{muon*}", 100, 0, 1000);
  muon_boost_lead_pt = book<TH1F>("boosted_muon_pt_lead", "p_{T}^{muon*}", 100, 0, 1000);

  Full_LepJetMass = book<TH1F>("M_full_lepjet", "m_{jet}^{had} [GeV]", 50, 0, 500);
  Lead_LepJetMass = book<TH1F>("M_lead_lepjet", "m_{jet}^{had} [GeV]", 50, 0, 500);
  leading_subjet_position = book<TH1F>("Leading_Subjet_Position", "Position", 5, 0, 5); // Controll Plot

  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
}

void LeptonicTop_Hists::fill(const Event & event){

  double weight = event.weight; // Get weight

  //=== Get Lepton =============================================================
  Particle lepton;
  if(event.muons->size() > 0) lepton = event.muons->at(0);
  else throw runtime_error("In LeptonicTop_Hists: Event should have one muon.");
  //=== Get Fatjets ============================================================
  std::vector<TopJet> fatjets = event.get(h_fatjets);
  if(fatjets.size() != 2) throw runtime_error("In LeptonicTop_Hists: fatjet should contain two jets.");

  //=== Find Lepjet ============================================================
  TopJet lepfatjet = fatjets.at(1);
  if(deltaR(lepton, fatjets.at(0)) < deltaR(lepton, fatjets.at(1))) lepfatjet = fatjets.at(0);
  std::vector<Jet> lep_subjets = lepfatjet.subjets();

  //=== Get Lepjet =============================================================
  Jet subjet1 = lep_subjets.at(0);
  Jet subjet2 = lep_subjets.at(1);
  Jet subjet3 = lep_subjets.at(2);
  LorentzVector lepjet_3subjets = subjet1.v4() + subjet2.v4() + subjet3.v4();

  //=== Get leading Subjet =====================================================
  Jet lead_subjet;
  if(subjet1.pt() > subjet2.pt()){
    if(subjet1.pt() > subjet3.pt()) lead_subjet = subjet1;
    else lead_subjet = subjet3;
  }
  else{
    if(subjet2.pt() > subjet3.pt()) lead_subjet = subjet2;
    else lead_subjet = subjet3;
  }

  if(lead_subjet == subjet1) leading_subjet_position->Fill(1, weight);
  else if(lead_subjet == subjet2) leading_subjet_position->Fill(2, weight);
  else if(lead_subjet == subjet3) leading_subjet_position->Fill(3, weight);
  else leading_subjet_position->Fill(4, weight);

  //=== Get Neutrinos ==========================================================
  // full: top is reconstructed by all three subjets, the lepton and the neutrino
  // lead: top is reconstructed by the leading subjet, the lepton and the neutrino
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction(lepton.v4(), event.met->v4());
  LorentzVector nu_lep, top_full, top_lead;

  // Neutrino Reconstruction returns a LorentzVector with one or two entries
  if(neutrinos.size() == 1) nu_lep = neutrinos.at(0);
  else{
    if(deltaR(neutrinos.at(0), lepton) < deltaR(neutrinos.at(1), lepton)) nu_lep = neutrinos.at(0);
    else nu_lep = neutrinos.at(1);
  }
  top_full = lepjet_3subjets + lepton.v4() + nu_lep;
  top_lead = lead_subjet.v4() + lepton.v4() + nu_lep;

  // ===========================================================================
  // ========================= Vector ==========================================
  // ===========================================================================
  // --- use vector to compress code -------------------------------------------
  vector<LorentzVector> top;
  top.push_back(top_full); // i == 0
  top.push_back(top_lead); // i == 1

  TVector3 lepton_v3 = toVector(lepton.v4());  // Get 3-Vector vom Lorentz
  TLorentzVector lepton_v4_T = lorentz_to_tlorentz(lepton.v4());  // Convert Lorentz to TLorentz

  for(unsigned int i=0; i<top.size(); i++){

    //=== Get Projection of Lepton perpendicular to Top ========================
    TVector3 top_v3 = toVector(top[i]);
    TVector3 lep_par_v3 = get_parallel_component(lepton_v3, top_v3);
    TVector3 lep_perp_v3 = get_perpendicular_component(lepton_v3, lep_par_v3);

    //=== Get LorentzBoost of Lepton into Top-Refrence Frame ====================
    TLorentzVector top_T = lorentz_to_tlorentz(top[i]);  // Convert Lorentz to TLorentz
    TVector3 boost_top = top_T.BoostVector();  // Get boost

    // Copy of TLV from lepton, since it will be boosted everytime in the for loop.
    TLorentzVector lepton_v4_T_boost = TLorentzVector(lepton_v4_T);
    lepton_v4_T_boost.Boost(boost_top); // Use boost of top on lepton

    // The Idea of the next step is the following:
    // Boost of muon in top-frame. Looking at the transformed pt (from CMS-frame) perpendicular to the top direction in the LAB-frame.
    TVector3 transformed_lep_v3 = lepton_v4_T_boost.Vect(); // Get 3-Vector from TLorentz
    TVector3 transformed_lep_par_v3 = get_parallel_component(transformed_lep_v3, top_v3);
    TVector3 transformed_lep_perp_v3 = get_perpendicular_component(transformed_lep_v3, transformed_lep_par_v3);

    //---------------------------------------------------------------------------------------
    //--------------------------------- Fill Hists here -------------------------------------
    //---------------------------------------------------------------------------------------

    if(i==0){
      muon_perp_full_momentum->Fill(lep_perp_v3.Mag(), weight);
      muon_para_full_momentum->Fill(lep_par_v3.Mag(), weight);
      muon_boost_perp_full_momentum->Fill(transformed_lep_perp_v3.Mag(), weight);
      muon_boost_para_full_momentum->Fill(transformed_lep_par_v3.Mag(), weight);
      muon_boost_full_pt->Fill(lepton_v4_T_boost.Pt(), weight);
      top_boost_full_pt->Fill(boost_top.Pt(), weight);
      Full_LepJetMass->Fill(top[i].M(), weight);
    }
    else if(i==1){
      muon_perp_lead_momentum->Fill(lep_perp_v3.Mag(), weight);
      muon_para_lead_momentum->Fill(lep_par_v3.Mag(), weight);
      muon_boost_perp_lead_momentum->Fill(transformed_lep_perp_v3.Mag(), weight);
      muon_boost_para_lead_momentum->Fill(transformed_lep_par_v3.Mag(), weight);
      muon_boost_lead_pt->Fill(lepton_v4_T_boost.Pt(), weight);
      top_boost_lead_pt->Fill(boost_top.Pt(), weight);
      Lead_LepJetMass->Fill(top[i].M(), weight);
    }
  }

  // Additinal Hists, not affected by the top
  muon_pt->Fill(lepton.pt(), weight);

}
