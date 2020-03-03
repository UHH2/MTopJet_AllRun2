#include <UHH2/MTopJet/include/LeptonicTop_Hists.h>

LeptonicTop_Hists::LeptonicTop_Hists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  MuonPt = book<TH1F>("MuonPt", "p_{T}^{muon}", 50, 0, 1000);
  MuonPerpendicularFullPt = book<TH1F>("MuonPtPerpFull", "p_{T,#perp}^{muon}", 200, 0, 200);
  MuonPerpendicularLeadPt = book<TH1F>("MuonPtPerpLead", "p_{T,#perp}^{muon}", 200, 0, 200);
  MuonParallelFullPt = book<TH1F>("MuonPtParFull", "p_{T,#parallel}^{muon}", 100, 0, 1000);
  MuonParallelLeadPt = book<TH1F>("MuonPtParLead", "p_{T,#parallel}^{muon}", 100, 0, 1000);

  TopBoostFullPt  = book<TH1F>("TopPtBoostFull", "p_{T}^{top}", 200, 0, 200); // Should be 0
  TopBoostLeadPt  = book<TH1F>("TopPtBoostLead", "p_{T}^{top}", 200, 0, 200); // Should be 0
  MuonBoostFullPt = book<TH1F>("MuonPtBoostFull", "p_{T}^{muon}", 100, 0, 1000);
  MuonBoostLeadPt = book<TH1F>("MuonPtBoostLead", "p_{T}^{muon}", 100, 0, 1000);

  Full_LepJetMass = book<TH1F>("M_full_lepjet", "m_{jet} [GeV]", 50, 0, 500);
  Lead_LepJetMass = book<TH1F>("M_lead_lepjet", "m_{jet} [GeV]", 50, 0, 500);
  leading_subjet_position = book<TH1F>("Leading Subjet Position", "Position", 5, 0, 5); // Controll Plot

  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");


}

void LeptonicTop_Hists::fill(const Event & event){

  double weight = event.weight; // Get weight

  // As you might notice, there is a confusion with TLorentzVectors and LorentzVectors.
  // This Code can surely be written more structured.

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
  else leading_subjet_position->Fill(5, weight);

  //=== Get Neutrinos ==========================================================
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction(lepton.v4(), event.met->v4());
  LorentzVector nu_lep, top_full, top_lead; // lepjet_full includes all three subjets, the lepton and the neutrino

  // Neutrino Reconstruction returns a LorentzVector with one or two entries
  if(neutrinos.size() == 1) nu_lep = neutrinos.at(0);
  else{
    if(deltaR(neutrinos.at(0), lepton) < deltaR(neutrinos.at(1), lepton)) nu_lep = neutrinos.at(0);
    else nu_lep = neutrinos.at(1);
  }
  top_full = lepjet_3subjets + lepton.v4() + nu_lep;
  top_lead = lead_subjet.v4() + lepton.v4() + nu_lep;

  //=== Get Projection of Lepton perpendicular to Top ==========================
  TVector3 lepton_v3 = toVector(lepton.v4()); // Momentum of lepton from 4-Vector
  TVector3 top_full_v3 = toVector(top_full);
  TVector3 top_lead_v3 = toVector(top_lead);

  // Magnitued of top 3-vector
  double mag_top_full = top_full_v3.Mag();
  double mag_top_lead = top_lead_v3.Mag();

  // lepton projection to tops
  top_full_v3.SetMag(lepton_v3 * top_full_v3);
  top_lead_v3.SetMag(lepton_v3 * top_lead_v3);

  // This part could also be included into SetMag
  TVector3 lep_par_full_v3 = top_full_v3 * (1/mag_top_full);
  TVector3 lep_par_lead_v3 = top_lead_v3 * (1/mag_top_lead);

  TVector3 lep_perp_full_v3 = lepton_v3 - lep_par_full_v3;
  TVector3 lep_perp_lead_v3 = lepton_v3 - lep_par_lead_v3;

  //=== Get LorentzBoost of Lepton into Top-Refrence Frame =====================

  // Convert Lorentz to TLorentz
  TLorentzVector lepton_v4_T_full = lorentz_to_tlorentz(lepton.v4());
  TLorentzVector top_full_T = lorentz_to_tlorentz(top_full);
  TLorentzVector top_lead_T = lorentz_to_tlorentz(top_lead);

  // so bekommst du den boost
  TVector3 boost_top_full = top_full_T.BoostVector();
  TVector3 boost_top_lead = top_lead_T.BoostVector();

  // und anwenden auf das Lepton
  // Here I am not sure, if boosting one Vector twice will lead to the same result.
  TLorentzVector lepton_v4_T_lead = lepton_v4_T_full;
  lepton_v4_T_full.Boost(boost_top_full);
  lepton_v4_T_lead.Boost(boost_top_lead);

  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  MuonPt->Fill(lepton.pt(), weight);
  MuonPerpendicularFullPt->Fill(lep_perp_full_v3.Mag(), weight);
  MuonPerpendicularLeadPt->Fill(lep_perp_lead_v3.Mag(), weight);
  MuonParallelFullPt->Fill(lep_par_full_v3.Mag(), weight);
  MuonParallelLeadPt->Fill(lep_par_lead_v3.Mag(), weight);

  MuonBoostFullPt->Fill(lepton_v4_T_full.Pt(), weight);
  MuonBoostLeadPt->Fill(lepton_v4_T_lead.Pt(), weight);
  TopBoostFullPt->Fill(boost_top_full.Pt(), weight);
  TopBoostLeadPt->Fill(boost_top_lead.Pt(), weight);

  Full_LepJetMass->Fill(top_full.M(), weight);
  Lead_LepJetMass->Fill(top_lead.M(), weight);

}


/*
.██████  ██████  ███    ██ ██    ██ ███████ ██████  ████████     ██       ██████  ██████  ███████ ███    ██ ████████ ███████
██      ██    ██ ████   ██ ██    ██ ██      ██   ██    ██        ██      ██    ██ ██   ██ ██      ████   ██    ██       ███
██      ██    ██ ██ ██  ██ ██    ██ █████   ██████     ██        ██      ██    ██ ██████  █████   ██ ██  ██    ██      ███
██      ██    ██ ██  ██ ██  ██  ██  ██      ██   ██    ██        ██      ██    ██ ██   ██ ██      ██  ██ ██    ██     ███
.██████  ██████  ██   ████   ████   ███████ ██   ██    ██        ███████  ██████  ██   ██ ███████ ██   ████    ██    ███████
*/

TLorentzVector lorentz_to_tlorentz(const LorentzVector v4){
  TLorentzVector lorentz_v4;
  double px, py, pz, Energy;
  px = v4.Px();
  py = v4.Py();
  pz = v4.Pz();
  Energy = v4.E();
  lorentz_v4.SetPxPyPzE(px, py, pz, Energy);

  return lorentz_v4;
}
