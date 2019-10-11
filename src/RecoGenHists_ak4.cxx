#include <UHH2/MTopJet/include/RecoGenHists_ak4.h>


RecoGenHists_ak4::RecoGenHists_ak4(uhh2::Context & ctx, const std::string & dirname, bool use_JEC_): Hists(ctx, dirname){
  // book all histograms here
  PtReso = book<TH1F>("PtResolution", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  area_all = book<TH1F>("area_all", "area all jets", 50, 0., 1.0);

  PtReso_1 = book<TH1F>("PtResolution_1", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_2 = book<TH1F>("PtResolution_2", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_3 = book<TH1F>("PtResolution_3", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_4 = book<TH1F>("PtResolution_4", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_5 = book<TH1F>("PtResolution_5", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_6 = book<TH1F>("PtResolution_6", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_7 = book<TH1F>("PtResolution_7", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_8 = book<TH1F>("PtResolution_8", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_9 = book<TH1F>("PtResolution_9", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);
  PtReso_10 = book<TH1F>("PtResolution_10", "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet}) ", 90, -1.5, 1.5);

  Energy_Check = book<TH1F>("Energy_Check", "E_{L1 corr. subjet} - (E_{raw subjet} - #rho A_{raw subjet})", 100, -40, 40);

  use_JEC = use_JEC_;
}
void RecoGenHists_ak4::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed

  std::vector<Jet> rec_jet = *event.jets;
  std::vector<GenJet> gen_jet = *event.genjets;

  if(rec_jet.size() < 1) return;
  if(gen_jet.size() < 1) return;

  /* ******************************************************************************
  matching between gen and reco jets:
  - a rec jet is called isolated if the next jet is not within 2R. An isolated jet should be spherical and more simelar to an ak4 jet
  - for each reco jet, calc distance to all other gen jets
  - gen jet with lowest distance is a match if distance is < 0.2
  - then calculate resolution with reco jet and matched gen jet
  - to do: account for double counting
  ********************************************************************************* */
  // if use_JEC is false, get JEC factors and recalculate pt
  double dR;
  double dR_temp;
  int nearest_j;
  double gen_pt, rec_pt, R;
  LorentzVector jet_v4;
  for(unsigned int i=0; i<rec_jet.size(); i++){
    dR = 1000;
    nearest_j = 100;
    area_all->Fill(rec_jet.at(i).jetArea(), weight);
    if(!use_JEC) jet_v4 = rec_jet.at(i).v4() * rec_jet.at(i).JEC_factor_raw();
    else jet_v4 = rec_jet.at(i).v4();
    rec_pt = jet_v4.Pt();
    for(unsigned int j=0; j<gen_jet.size(); j++){
      dR_temp = uhh2::deltaR(rec_jet.at(i), gen_jet.at(j));
      if(dR_temp < dR){
        dR = dR_temp;
        nearest_j = j;
      }
    }
    gen_pt=gen_jet.at(nearest_j).v4().Pt();
    R = (rec_pt - gen_pt) / gen_pt;
    if(nearest_j != 100 && dR <= 0.2){
      PtReso->Fill( (rec_pt - gen_pt)/gen_pt, weight );
      if(gen_pt <= 50) PtReso_1->Fill( R, weight );
      if(gen_pt > 50 && gen_pt <= 80) PtReso_2->Fill( R, weight );
      if(gen_pt > 80 && gen_pt <= 120) PtReso_3->Fill( R, weight );
      if(gen_pt > 120 && gen_pt <= 170) PtReso_4->Fill( R, weight );
      if(gen_pt > 170 && gen_pt <= 220) PtReso_5->Fill( R, weight );
      if(gen_pt > 220 && gen_pt <= 270) PtReso_6->Fill( R, weight );
      if(gen_pt > 270 && gen_pt <= 320) PtReso_7->Fill( R, weight );
      if(gen_pt > 320 && gen_pt <= 370) PtReso_8->Fill( R, weight );
      if(gen_pt > 370 && gen_pt <= 420) PtReso_9->Fill( R, weight );
      if(gen_pt > 420) PtReso_10->Fill( R, weight );
    }
  }

  // check: is Energy with Jec == E without JEC + rho A?
  double rho = event.rho;
  double area, add_E;
  LorentzVector JEC_jet, raw_jet;
  for(unsigned int i=0; i<rec_jet.size(); i++){
    area = rec_jet.at(i).jetArea();
    add_E = area * rho;
    JEC_jet = rec_jet.at(i).v4() * rec_jet.at(i).JEC_factor_raw() * rec_jet.at(i).JEC_L1factor_raw();
    raw_jet = rec_jet.at(i).v4() * rec_jet.at(i).JEC_factor_raw();
    Energy_Check->Fill( (JEC_jet.E() - (raw_jet.E() - add_E)), weight);
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


}
