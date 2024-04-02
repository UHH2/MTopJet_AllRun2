#include <UHH2/MTopJet/include/StoreKinematics.h>

StoreKinematics::StoreKinematics(uhh2::Context & ctx, std::vector<double> jmss, bool isTTbar, bool isMC, bool debug):
  isTTbar_(isTTbar),
  isMC_(isMC),
  debug_(debug),
  jmss_(jmss)
  {
    h_hadjets = ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected"); 
    h_lepjets = ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined_Corrected");  
    h_hadjets_raw = ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");

    h_sub1_pt_rec = ctx.declare_event_output<double>("sub1_pt_rec");
    h_sub1_px_rec = ctx.declare_event_output<double>("sub1_px_rec");
    h_sub1_py_rec = ctx.declare_event_output<double>("sub1_py_rec");
    h_sub1_pz_rec = ctx.declare_event_output<double>("sub1_pz_rec");
    h_sub1_E_rec = ctx.declare_event_output<double>("sub1_E_rec");
    h_sub1_mass_rec = ctx.declare_event_output<double>("sub1_mass_rec");
    h_sub1_eta_rec = ctx.declare_event_output<double>("sub1_eta_rec");
    h_sub1_phi_rec = ctx.declare_event_output<double>("sub1_phi_rec");
    h_sub1_pt_gen = ctx.declare_event_output<double>("sub1_pt_gen");
    h_sub1_px_gen = ctx.declare_event_output<double>("sub1_px_gen");
    h_sub1_py_gen = ctx.declare_event_output<double>("sub1_py_gen");
    h_sub1_pz_gen = ctx.declare_event_output<double>("sub1_pz_gen");
    h_sub1_E_gen = ctx.declare_event_output<double>("sub1_E_gen");
    h_sub1_mass_gen = ctx.declare_event_output<double>("sub1_mass_gen");
    h_sub1_eta_gen = ctx.declare_event_output<double>("sub1_eta_gen");
    h_sub1_phi_gen = ctx.declare_event_output<double>("sub1_phi_gen");
    h_sub1_uncor_pt_rec = ctx.declare_event_output<double>("sub1_uncor_pt_rec");
    h_sub1_uncor_px_rec = ctx.declare_event_output<double>("sub1_uncor_px_rec");
    h_sub1_uncor_py_rec = ctx.declare_event_output<double>("sub1_uncor_py_rec");
    h_sub1_uncor_pz_rec = ctx.declare_event_output<double>("sub1_uncor_pz_rec");
    h_sub1_uncor_E_rec = ctx.declare_event_output<double>("sub1_uncor_E_rec");
    h_sub1_uncor_mass_rec = ctx.declare_event_output<double>("sub1_uncor_mass_rec");
    h_sub1_uncor_eta_rec = ctx.declare_event_output<double>("sub1_uncor_eta_rec");
    h_sub1_uncor_phi_rec = ctx.declare_event_output<double>("sub1_uncor_phi_rec");
    h_sub1_factor_jec = ctx.declare_event_output<double>("sub1_factor_jec");
    h_sub1_sigma_jec = ctx.declare_event_output<double>("sub1_sigma_jec");
    h_sub1_factor_cor = ctx.declare_event_output<double>("sub1_factor_cor");
    h_sub1_sigma_cor = ctx.declare_event_output<double>("sub1_sigma_cor");
    h_sub1_factor_jer = ctx.declare_event_output<double>("sub1_factor_jer");

    h_sub2_pt_rec = ctx.declare_event_output<double>("sub2_pt_rec");
    h_sub2_px_rec = ctx.declare_event_output<double>("sub2_px_rec");
    h_sub2_py_rec = ctx.declare_event_output<double>("sub2_py_rec");
    h_sub2_pz_rec = ctx.declare_event_output<double>("sub2_pz_rec");
    h_sub2_E_rec = ctx.declare_event_output<double>("sub2_E_rec");
    h_sub2_mass_rec = ctx.declare_event_output<double>("sub2_mass_rec");
    h_sub2_eta_rec = ctx.declare_event_output<double>("sub2_eta_rec");
    h_sub2_phi_rec = ctx.declare_event_output<double>("sub2_phi_rec");
    h_sub2_pt_gen = ctx.declare_event_output<double>("sub2_pt_gen");
    h_sub2_px_gen = ctx.declare_event_output<double>("sub2_px_gen");
    h_sub2_py_gen = ctx.declare_event_output<double>("sub2_py_gen");
    h_sub2_pz_gen = ctx.declare_event_output<double>("sub2_pz_gen");
    h_sub2_E_gen = ctx.declare_event_output<double>("sub2_E_gen");
    h_sub2_mass_gen = ctx.declare_event_output<double>("sub2_mass_gen");
    h_sub2_eta_gen = ctx.declare_event_output<double>("sub2_eta_gen");
    h_sub2_phi_gen = ctx.declare_event_output<double>("sub2_phi_gen");
    h_sub2_uncor_pt_rec = ctx.declare_event_output<double>("sub2_uncor_pt_rec");
    h_sub2_uncor_px_rec = ctx.declare_event_output<double>("sub2_uncor_px_rec");
    h_sub2_uncor_py_rec = ctx.declare_event_output<double>("sub2_uncor_py_rec");
    h_sub2_uncor_pz_rec = ctx.declare_event_output<double>("sub2_uncor_pz_rec");
    h_sub2_uncor_E_rec = ctx.declare_event_output<double>("sub2_uncor_E_rec");
    h_sub2_uncor_mass_rec = ctx.declare_event_output<double>("sub2_uncor_mass_rec");
    h_sub2_uncor_eta_rec = ctx.declare_event_output<double>("sub2_uncor_eta_rec");
    h_sub2_uncor_phi_rec = ctx.declare_event_output<double>("sub2_uncor_phi_rec");
    h_sub2_factor_jec = ctx.declare_event_output<double>("sub2_factor_jec");
    h_sub2_sigma_jec = ctx.declare_event_output<double>("sub2_sigma_jec");
    h_sub2_factor_cor = ctx.declare_event_output<double>("sub2_factor_cor");
    h_sub2_sigma_cor = ctx.declare_event_output<double>("sub2_sigma_cor");
    h_sub2_factor_jer = ctx.declare_event_output<double>("sub2_factor_jer");
    h_sub3_pt_rec = ctx.declare_event_output<double>("sub3_pt_rec");
    h_sub3_px_rec = ctx.declare_event_output<double>("sub3_px_rec");
    h_sub3_py_rec = ctx.declare_event_output<double>("sub3_py_rec");
    h_sub3_pz_rec = ctx.declare_event_output<double>("sub3_pz_rec");
    h_sub3_E_rec = ctx.declare_event_output<double>("sub3_E_rec");
    h_sub3_mass_rec = ctx.declare_event_output<double>("sub3_mass_rec");
    h_sub3_eta_rec = ctx.declare_event_output<double>("sub3_eta_rec");
    h_sub3_phi_rec = ctx.declare_event_output<double>("sub3_phi_rec");
    h_sub3_pt_gen = ctx.declare_event_output<double>("sub3_pt_gen");
    h_sub3_px_gen = ctx.declare_event_output<double>("sub3_px_gen");
    h_sub3_py_gen = ctx.declare_event_output<double>("sub3_py_gen");
    h_sub3_pz_gen = ctx.declare_event_output<double>("sub3_pz_gen");
    h_sub3_E_gen = ctx.declare_event_output<double>("sub3_E_gen");
    h_sub3_mass_gen = ctx.declare_event_output<double>("sub3_mass_gen");
    h_sub3_eta_gen = ctx.declare_event_output<double>("sub3_eta_gen");
    h_sub3_phi_gen = ctx.declare_event_output<double>("sub3_phi_gen");
    h_sub3_uncor_pt_rec = ctx.declare_event_output<double>("sub3_uncor_pt_rec");
    h_sub3_uncor_px_rec = ctx.declare_event_output<double>("sub3_uncor_px_rec");
    h_sub3_uncor_py_rec = ctx.declare_event_output<double>("sub3_uncor_py_rec");
    h_sub3_uncor_pz_rec = ctx.declare_event_output<double>("sub3_uncor_pz_rec");
    h_sub3_uncor_E_rec = ctx.declare_event_output<double>("sub3_uncor_E_rec");
    h_sub3_uncor_mass_rec = ctx.declare_event_output<double>("sub3_uncor_mass_rec");
    h_sub3_uncor_eta_rec = ctx.declare_event_output<double>("sub3_uncor_eta_rec");
    h_sub3_uncor_phi_rec = ctx.declare_event_output<double>("sub3_uncor_phi_rec");
    h_sub3_factor_jec = ctx.declare_event_output<double>("sub3_factor_jec");
    h_sub3_sigma_jec = ctx.declare_event_output<double>("sub3_sigma_jec");
    h_sub3_factor_cor = ctx.declare_event_output<double>("sub3_factor_cor");
    h_sub3_sigma_cor = ctx.declare_event_output<double>("sub3_sigma_cor");
    h_sub3_factor_jer = ctx.declare_event_output<double>("sub3_factor_jer");
    
    h_factor_jms_jec = ctx.declare_event_output<double>("factor_jms_jec");
    h_factor_jms_cor = ctx.declare_event_output<double>("factor_jms_cor");

    h_sub12_dR_rec = ctx.declare_event_output<double>("sub12_dR_rec");
    h_sub23_dR_rec = ctx.declare_event_output<double>("sub23_dR_rec");
    h_sub13_dR_rec = ctx.declare_event_output<double>("sub13_dR_rec");
    h_sub12_dR_gen = ctx.declare_event_output<double>("sub12_dR_gen");
    h_sub23_dR_gen = ctx.declare_event_output<double>("sub23_dR_gen");
    h_sub13_dR_gen = ctx.declare_event_output<double>("sub13_dR_gen");

    h_sub12_mass_rec = ctx.declare_event_output<double>("sub12_mass_rec");
    h_sub23_mass_rec = ctx.declare_event_output<double>("sub23_mass_rec");
    h_sub13_mass_rec = ctx.declare_event_output<double>("sub13_mass_rec");
    h_sub12_mass_gen = ctx.declare_event_output<double>("sub12_mass_gen");
    h_sub23_mass_gen = ctx.declare_event_output<double>("sub23_mass_gen");
    h_sub13_mass_gen = ctx.declare_event_output<double>("sub13_mass_gen");

    h_minimal_dR_rec = ctx.declare_event_output<double>("minimal_dR_rec");
    h_minimal_dR_gen = ctx.declare_event_output<double>("minimal_dR_gen");

    h_matched_index1 = ctx.declare_event_output<int>("matched_index1");
    h_matched_index2 = ctx.declare_event_output<int>("matched_index2");
    h_matched_index3 = ctx.declare_event_output<int>("matched_index3");

    h_index_w1 = ctx.declare_event_output<int>("index_w1");
    h_index_w2 = ctx.declare_event_output<int>("index_w2");
    h_index_b = ctx.declare_event_output<int>("index_b");

    h_index_parton_b_rec = ctx.declare_event_output<int>("index_parton_b_rec");
    h_index_parton_b_gen = ctx.declare_event_output<int>("index_parton_b_gen");

    h_same_leading_order = ctx.declare_event_output<bool>("same_leading_order");
    h_all_subjets_matched_gen = ctx.declare_event_output<bool>("all_subjets_matched_gen");
    h_all_subjets_matched_rec = ctx.declare_event_output<bool>("all_subjets_matched_rec");
    h_mass_hadtop = ctx.declare_event_output<double>("mass_had_top");

    if(isTTbar){
      h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
      h_hadjets_gen = ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");

    }
  }

bool StoreKinematics::process(uhh2::Event & event){
  double dummy = event.weight; // avoid warning
  dummy += 0;
  return true;
}

void StoreKinematics::store(
  uhh2::Event & event, vector<Jet>& subjets_jms,
  vector<int>& w_index, vector<vector<double>>& jecs, vector<vector<double>>& cors, vector<double>& jers
){

  std::vector<TopJet>hadjets_raw = event.get(h_hadjets_raw);
  std::vector<TopJet>hadjets = event.get(h_hadjets);
  std::vector<TopJet>lepjets = event.get(h_lepjets);

  vector<Jet> subjets = hadjets[0].subjets();
  vector<Jet> subjets_raw = hadjets_raw[0].subjets();
  
  if(debug_) cout << printf("... Check #jer %2lu\n",jers.size());
  event.set(h_sub1_factor_jer, jers[0]);
  event.set(h_sub2_factor_jer, jers[1]);
  event.set(h_sub3_factor_jer, jers[2]);

  if(debug_) cout << printf("... Check #jms %2lu\n",jmss_.size());
  event.set(h_factor_jms_jec, jmss_[0]);
  event.set(h_factor_jms_cor, jmss_[1]);
  
  // ==================
  // =   Store JMS    =
  // ==================

  if(debug_) cout << printf("... Check #jms %2lu\n",subjets_jms.size());
  
  LorentzVector lv_sub1, lv_sub2, lv_sub3;
  // sort subjets on reco level
  lv_sub1 = subjets_jms[0].v4();
  lv_sub2 = subjets_jms[1].v4();
  lv_sub3 = subjets_jms[2].v4();

  // Store als TLorentzVector to add vectors
  TLorentzVector tlv_sub1(lv_sub1.Px(), lv_sub1.Py(), lv_sub1.Pz(), lv_sub1.E());
  TLorentzVector tlv_sub2(lv_sub2.Px(), lv_sub2.Py(), lv_sub2.Pz(), lv_sub2.E());
  TLorentzVector tlv_sub3(lv_sub3.Px(), lv_sub3.Py(), lv_sub3.Pz(), lv_sub3.E());
  TLorentzVector tlv_sub12(tlv_sub1+tlv_sub2);  
  TLorentzVector tlv_sub13(tlv_sub1+tlv_sub3);
  TLorentzVector tlv_sub23(tlv_sub2+tlv_sub3);  
  event.set(h_sub1_mass_rec, tlv_sub1.M());
  event.set(h_sub2_mass_rec, tlv_sub2.M());
  event.set(h_sub3_mass_rec, tlv_sub3.M());
  event.set(h_sub12_mass_rec, tlv_sub12.M());
  event.set(h_sub13_mass_rec, tlv_sub13.M());
  event.set(h_sub23_mass_rec, tlv_sub23.M());

  // ==================
  // =   Store JEC    =
  // ==================

  double dR_12 = deltaR(subjets.at(0).v4(), subjets.at(1).v4());
  double dR_13 = deltaR(subjets.at(0).v4(), subjets.at(2).v4());
  double dR_23 = deltaR(subjets.at(1).v4(), subjets.at(2).v4());
  double dR_min_rec = -1;
  if     (dR_12<dR_13 && dR_12<dR_23) dR_min_rec = dR_12;
  else if(dR_13<dR_12 && dR_13<dR_23) dR_min_rec = dR_13;
  else if(dR_23<dR_12 && dR_23<dR_13) dR_min_rec = dR_23;
  event.set(h_minimal_dR_rec, dR_min_rec);

  int lead1_rec = 0;
  int lead2_rec = 1;
  int lead3_rec = 2;
  double pt1=subjets.at(0).v4().Pt();
  double pt2=subjets.at(1).v4().Pt();
  double pt3=subjets.at(2).v4().Pt();
  if     (pt3<pt2 && pt2<pt1){lead1_rec=0;lead2_rec=1;lead3_rec=2;}
  else if(pt2<pt3 && pt3<pt1){lead1_rec=0;lead2_rec=2;lead3_rec=1;}
  else if(pt1<pt3 && pt3<pt2){lead1_rec=1;lead2_rec=2;lead3_rec=0;}
  else if(pt3<pt1 && pt1<pt2){lead1_rec=1;lead2_rec=0;lead3_rec=2;}
  else if(pt2<pt1 && pt1<pt3){lead1_rec=2;lead2_rec=0;lead3_rec=1;}
  else if(pt1<pt2 && pt2<pt3){lead1_rec=2;lead2_rec=1;lead3_rec=0;}

  Jet jet_w1 = subjets[w_index[0]];
  Jet jet_w2 = subjets[w_index[1]];
  Jet jet_b = subjets[w_index[2]];
  int index_w1 = w_index[0];
  int index_w2 = w_index[1];
  int index_b = w_index[2];
  double d_w1_1 = deltaR( jet_w1, subjets.at(lead1_rec) );
  double d_w1_2 = deltaR( jet_w1, subjets.at(lead2_rec) );
  double d_w1_3 = deltaR( jet_w1, subjets.at(lead3_rec) );
  double d_w2_1 = deltaR( jet_w2, subjets.at(lead1_rec) );
  double d_w2_2 = deltaR( jet_w2, subjets.at(lead2_rec) );
  double d_w2_3 = deltaR( jet_w2, subjets.at(lead3_rec) );
  double d_b_1 = deltaR( jet_b, subjets.at(lead1_rec) );
  double d_b_2 = deltaR( jet_b, subjets.at(lead2_rec) );
  double d_b_3 = deltaR( jet_b, subjets.at(lead3_rec) );
  if     (d_w1_1<d_w1_2 && d_w1_1<d_w1_3) index_w1 = 0;
  else if(d_w1_2<d_w1_1 && d_w1_2<d_w1_3) index_w1 = 1;
  else if(d_w1_3<d_w1_1 && d_w1_3<d_w1_2) index_w1 = 2;
  if     (d_w2_1<d_w2_2 && d_w2_1<d_w2_3) index_w2 = 0;
  else if(d_w2_2<d_w2_1 && d_w2_2<d_w2_3) index_w2 = 1;
  else if(d_w2_3<d_w2_1 && d_w2_3<d_w2_2) index_w2 = 2;
  if     (d_b_1<d_b_2 && d_b_1<d_b_3) index_b = 0;
  else if(d_b_2<d_b_1 && d_b_2<d_b_3) index_b = 1;
  else if(d_b_3<d_b_1 && d_b_3<d_b_2) index_b = 2;
  if(index_w1 == index_w2 || index_w1 == index_b || index_w2 == index_b){
    printf("Indices:\n");
    printf("Before: w1 %2d w2 %2d b %2d\n",w_index[0],w_index[1],w_index[2]);
    printf("After: w1 %2d w2 %2d b %2d\n",index_w1,index_w2,index_b);
    printf("\nPt:\n");
    printf("Before: pt1 %5.2f pt2 %5.2f pt3 %5.2f\n",subjets.at(0).v4().Pt(),subjets.at(1).v4().Pt(),subjets.at(2).v4().Pt());
    printf("After: pt1 %5.2f pt2 %5.2f pt3 %5.2f\n",subjets.at(lead1_rec).v4().Pt(),subjets.at(lead2_rec).v4().Pt(),subjets.at(lead3_rec).v4().Pt());
    printf("\nDistance:\n");
    printf("w1: %5.2f %5.2f %5.2f\n",d_w1_1,d_w1_2,d_w1_3);
    printf("w1: %5.2f %5.2f %5.2f\n",d_w2_1,d_w2_2,d_w2_3);
    printf("b : %5.2f %5.2f %5.2f\n",d_b_1,d_b_2,d_b_3);
    throw runtime_error("[ERROR] Ordering the w indices lead to equal indices!");
  }

  if(debug_){
    printf("Check sorted jet correction factors:\n");
    printf(" --- JEC (unsorted / sorted) nominal:\n");
    printf("\t%6.2f %6.2f %6.2f\n",jecs[0][0],jecs[1][0],jecs[2][0]);
    printf("\t%6.2f %6.2f %6.2f\n",jecs[lead1_rec][0],jecs[lead2_rec][0],jecs[lead3_rec][0]);
    printf(" --- JEC (unsorted / sorted) sigma:\n");
    printf("\t%6.2f %6.2f %6.2f\n",jecs[0][1],jecs[1][1],jecs[2][1]);
    printf("\t%6.2f %6.2f %6.2f\n",jecs[lead1_rec][1],jecs[lead2_rec][1],jecs[lead3_rec][1]);
    printf("\n --- XCone (unsorted / sorted) nominal:\n");
    printf("\t%6.2f %6.2f %6.2f\n",cors[0][0],cors[1][0],cors[2][0]);
    printf("\t%6.2f %6.2f %6.2f\n",cors[lead1_rec][0],cors[lead2_rec][0],cors[lead3_rec][0]);
    printf(" --- XCone (unsorted / sorted) sigma:\n");
    printf("\t%6.2f %6.2f %6.2f\n",cors[0][1],cors[1][1],cors[2][1]);
    printf("\t%6.2f %6.2f %6.2f\n",cors[lead1_rec][1],cors[lead2_rec][1],cors[lead3_rec][1]);
  }

  if(debug_) printf("... Start with Subjet 1\n");
  event.set(h_sub1_pt_rec, subjets.at(lead1_rec).v4().Pt());
  event.set(h_sub1_eta_rec, subjets.at(lead1_rec).eta());
  event.set(h_sub1_phi_rec, subjets.at(lead1_rec).phi());
  event.set(h_sub1_px_rec, subjets.at(lead1_rec).v4().Px());
  event.set(h_sub1_py_rec, subjets.at(lead1_rec).v4().Py());
  event.set(h_sub1_pz_rec, subjets.at(lead1_rec).v4().Pz());
  event.set(h_sub1_E_rec,  subjets.at(lead1_rec).v4().E());
  event.set(h_sub1_uncor_pt_rec, subjets_raw.at(lead1_rec).v4().Pt());
  event.set(h_sub1_uncor_eta_rec, subjets_raw.at(lead1_rec).eta());
  event.set(h_sub1_uncor_phi_rec, subjets_raw.at(lead1_rec).phi());
  event.set(h_sub1_uncor_px_rec, subjets_raw.at(lead1_rec).v4().Px());
  event.set(h_sub1_uncor_py_rec, subjets_raw.at(lead1_rec).v4().Py());
  event.set(h_sub1_uncor_pz_rec, subjets_raw.at(lead1_rec).v4().Pz());
  event.set(h_sub1_uncor_E_rec,  subjets_raw.at(lead1_rec).v4().E());
  event.set(h_sub1_uncor_mass_rec,  subjets_raw.at(lead1_rec).v4().M());
  event.set(h_sub1_factor_jec, jecs[lead1_rec][0]);
  event.set(h_sub1_sigma_jec, jecs[lead1_rec][1]);
  event.set(h_sub1_factor_cor, cors[lead1_rec][0]);
  event.set(h_sub1_sigma_cor, cors[lead1_rec][1]);

  if(debug_) printf("... Start with Subjet 2\n");
  event.set(h_sub2_pt_rec, subjets.at(lead2_rec).v4().Pt());
  event.set(h_sub2_eta_rec, subjets.at(lead2_rec).eta());
  event.set(h_sub2_phi_rec, subjets.at(lead2_rec).phi());
  event.set(h_sub2_px_rec, subjets.at(lead2_rec).v4().Px());
  event.set(h_sub2_py_rec, subjets.at(lead2_rec).v4().Py());
  event.set(h_sub2_pz_rec, subjets.at(lead2_rec).v4().Pz());
  event.set(h_sub2_E_rec,  subjets.at(lead2_rec).v4().E());
  event.set(h_sub2_uncor_pt_rec, subjets_raw.at(lead2_rec).v4().Pt());
  event.set(h_sub2_uncor_eta_rec, subjets_raw.at(lead2_rec).eta());
  event.set(h_sub2_uncor_phi_rec, subjets_raw.at(lead2_rec).phi());
  event.set(h_sub2_uncor_px_rec, subjets_raw.at(lead2_rec).v4().Px());
  event.set(h_sub2_uncor_py_rec, subjets_raw.at(lead2_rec).v4().Py());
  event.set(h_sub2_uncor_pz_rec, subjets_raw.at(lead2_rec).v4().Pz());
  event.set(h_sub2_uncor_E_rec,  subjets_raw.at(lead2_rec).v4().E());
  event.set(h_sub2_uncor_mass_rec,  subjets_raw.at(lead2_rec).v4().M());
  event.set(h_sub2_factor_jec, jecs[lead2_rec][0]);
  event.set(h_sub2_sigma_jec, jecs[lead2_rec][1]);
  event.set(h_sub2_factor_cor, cors[lead2_rec][0]);
  event.set(h_sub2_sigma_cor, cors[lead2_rec][1]);

  if(debug_) printf("... Start with Subjet 3\n");
  event.set(h_sub3_pt_rec, subjets.at(lead3_rec).v4().Pt());
  event.set(h_sub3_eta_rec, subjets.at(lead3_rec).eta());
  event.set(h_sub3_phi_rec, subjets.at(lead3_rec).phi());
  event.set(h_sub3_px_rec, subjets.at(lead3_rec).v4().Px());
  event.set(h_sub3_py_rec, subjets.at(lead3_rec).v4().Py());
  event.set(h_sub3_pz_rec, subjets.at(lead3_rec).v4().Pz());
  event.set(h_sub3_E_rec,  subjets.at(lead3_rec).v4().E());
  event.set(h_sub3_uncor_pt_rec, subjets_raw.at(lead3_rec).v4().Pt());
  event.set(h_sub3_uncor_eta_rec, subjets_raw.at(lead3_rec).eta());
  event.set(h_sub3_uncor_phi_rec, subjets_raw.at(lead3_rec).phi());
  event.set(h_sub3_uncor_px_rec, subjets_raw.at(lead3_rec).v4().Px());
  event.set(h_sub3_uncor_py_rec, subjets_raw.at(lead3_rec).v4().Py());
  event.set(h_sub3_uncor_pz_rec, subjets_raw.at(lead3_rec).v4().Pz());
  event.set(h_sub3_uncor_E_rec,  subjets_raw.at(lead3_rec).v4().E());
  event.set(h_sub3_uncor_mass_rec,  subjets_raw.at(lead3_rec).v4().M());
  event.set(h_sub3_factor_jec, jecs[lead3_rec][0]);
  event.set(h_sub3_sigma_jec, jecs[lead3_rec][1]);
  event.set(h_sub3_factor_cor, cors[lead3_rec][0]);
  event.set(h_sub3_sigma_cor, cors[lead3_rec][1]);

  if(debug_) printf("... Start with indices\n");
  event.set(h_index_w1,  index_w1+1);
  event.set(h_index_w2,  index_w2+1);
  event.set(h_index_b,  index_b+1);
  vector<Jet> subjets_pt_sorted = {subjets[lead1_rec], subjets[lead2_rec], subjets[lead3_rec]};

  map<int, double> rec_sub_pt, rec_sub_dR;
  // select which gen fatjet is closer to had jet and get those gen subjets
  if(debug_) printf("... Start with distances\n");
  rec_sub_dR[12] = deltaR(subjets.at(lead1_rec), subjets.at(lead2_rec));
  rec_sub_dR[23] = deltaR(subjets.at(lead2_rec), subjets.at(lead3_rec));
  rec_sub_dR[13] = deltaR(subjets.at(lead1_rec), subjets.at(lead3_rec));
  if(debug_) printf("... Start with setting distances\n");
  event.set(h_sub12_dR_rec, rec_sub_dR[12]);
  event.set(h_sub23_dR_rec, rec_sub_dR[23]);
  event.set(h_sub13_dR_rec, rec_sub_dR[13]);

  // ==================
  // =   Store GEN    =
  // ==================

  std::vector<GenTopJet> hadjets_gen; // used later in the code as well
  vector<TLorentzVector> tlv_subjets_gen;
  vector<GenJet> subjets_gen;
  vector<GenJet> subjets_gen_pt_sorted;

  int lead1_gen=0;
  int lead2_gen=1;
  int lead3_gen=2;
  bool mismatch = false;
  bool fill_gen = false; // to fill non TTbar and mismatched at once

  if(isTTbar_){ // before is TTbar
    if(debug_) printf("... Start with gen jets\n");
    hadjets_gen = event.get(h_hadjets_gen);
    subjets_gen = hadjets_gen[0].subjets();

    if(debug_) printf("... Start with mismatch\n");
    if( deltaR(lepjets.at(0), hadjets_gen.at(0)) < deltaR(hadjets.at(0), hadjets_gen.at(0)) ) mismatch=true;

    // select which gen fatjet is closer to had jet and get those gen subjets
    if(debug_) cout << "Store pt and dR of subjets on gen level with #subjets=" << hadjets_gen.size() << endl;
    if( !(subjets_gen.size()<3 || mismatch) ) fill_gen = true;

    if(fill_gen){
      LorentzVector lv_sub1_gen = hadjets_gen[0].subjets().at(0).v4();
      LorentzVector lv_sub2_gen = hadjets_gen[0].subjets().at(1).v4();
      LorentzVector lv_sub3_gen = hadjets_gen[0].subjets().at(2).v4();
      TLorentzVector tlv_sub1_gen(lv_sub1_gen.Px(), lv_sub1_gen.Py(), lv_sub1_gen.Pz(), lv_sub1_gen.E());
      TLorentzVector tlv_sub2_gen(lv_sub2_gen.Px(), lv_sub2_gen.Py(), lv_sub2_gen.Pz(), lv_sub2_gen.E());
      TLorentzVector tlv_sub3_gen(lv_sub3_gen.Px(), lv_sub3_gen.Py(), lv_sub3_gen.Pz(), lv_sub3_gen.E());
      tlv_subjets_gen.push_back(tlv_sub1_gen);
      tlv_subjets_gen.push_back(tlv_sub2_gen);
      tlv_subjets_gen.push_back(tlv_sub3_gen);
      double pt1=tlv_subjets_gen[0].Pt();
      double pt2=tlv_subjets_gen[1].Pt();
      double pt3=tlv_subjets_gen[2].Pt();

      // ///////////////////////////////////////////////////
      // In one event two subjets have the same pt, hence <=
      //       74.30       74.30      338.70
      //        0           1           2
      // dR01   0.63 dR02   0.47 dR12   0.92
      if(pt3<=pt1 && pt2<=pt1){
        lead1_gen=0;
        if(pt3<=pt2) {lead2_gen=1; lead3_gen=2;}
        else         {lead2_gen=2; lead3_gen=1;}
      }
      else if(pt1<=pt2 && pt3<=pt2){
        lead1_gen=1;
        if(pt3<=pt1) {lead2_gen=0; lead3_gen=2;}
        else         {lead2_gen=2; lead3_gen=0;}
      }
      else if(pt1<=pt3 && pt2<=pt3){
        lead1_gen=2;
        if(pt2<=pt1) {lead2_gen=0; lead3_gen=1;}
        else         {lead2_gen=1; lead3_gen=0;}
      }
      else{
        printf("pt1   %6.2f pt2   %6.2f pt3   %6.2f\n",pt1,pt2,pt3);
        printf("lead1 %6d   lead2 %6d   lead2 %6d\n",lead1_gen,lead2_gen,lead3_gen);
        printf("dR01  %6.2f dR02  %6.2f dR12  %6.2f\n",
          deltaR(hadjets_gen[0].subjets().at(0), hadjets_gen[0].subjets().at(1)),
          deltaR(hadjets_gen[0].subjets().at(0), hadjets_gen[0].subjets().at(2)),
          deltaR(hadjets_gen[0].subjets().at(1), hadjets_gen[0].subjets().at(2))
        );
        throw runtime_error("[ERROR] in gen subjet ordering; Who is the largest now?");
      }

      if(lead1_gen==lead2_gen || lead1_gen==lead3_gen || lead2_gen==lead3_gen){
        printf("lead1 %3d lead2 %3d lead3 %3d\n",lead1_gen,lead2_gen,lead3_gen);
        throw runtime_error("[ERROR] leading indices are equal!");
      }

      subjets_gen_pt_sorted = {subjets_gen[lead1_gen], subjets_gen[lead2_gen], subjets_gen[lead3_gen]};
      TLorentzVector tlv_lead1_gen = tlv_subjets_gen[lead1_gen];
      TLorentzVector tlv_lead2_gen = tlv_subjets_gen[lead2_gen];
      TLorentzVector tlv_lead3_gen = tlv_subjets_gen[lead3_gen];
      map<int, double> gen_sub_dR;
      gen_sub_dR[12] = deltaR(hadjets_gen[0].subjets().at(lead1_gen), hadjets_gen[0].subjets().at(lead2_gen));
      gen_sub_dR[23] = deltaR(hadjets_gen[0].subjets().at(lead2_gen), hadjets_gen[0].subjets().at(lead3_gen));
      gen_sub_dR[13] = deltaR(hadjets_gen[0].subjets().at(lead1_gen), hadjets_gen[0].subjets().at(lead3_gen));
      double dR_min_gen=-1;
      if     (gen_sub_dR[12]<gen_sub_dR[13] && gen_sub_dR[12]<gen_sub_dR[23]) dR_min_gen = gen_sub_dR[12];
      else if(gen_sub_dR[13]<gen_sub_dR[12] && gen_sub_dR[13]<gen_sub_dR[23]) dR_min_gen = gen_sub_dR[13];
      else if(gen_sub_dR[23]<gen_sub_dR[12] && gen_sub_dR[23]<gen_sub_dR[13]) dR_min_gen = gen_sub_dR[23];
      event.set(h_minimal_dR_gen, dR_min_gen);
      event.set(h_sub12_dR_gen, gen_sub_dR[12]);
      event.set(h_sub23_dR_gen, gen_sub_dR[23]);
      event.set(h_sub13_dR_gen, gen_sub_dR[13]);

      event.set(h_sub1_pt_gen, tlv_lead1_gen.Pt());
      event.set(h_sub2_pt_gen, tlv_lead2_gen.Pt());
      event.set(h_sub3_pt_gen, tlv_lead3_gen.Pt());
      event.set(h_sub1_mass_gen, tlv_lead1_gen.M());
      event.set(h_sub2_mass_gen, tlv_lead2_gen.M());
      event.set(h_sub3_mass_gen, tlv_lead3_gen.M());
      event.set(h_sub12_mass_gen, tlv_sub12.M());
      event.set(h_sub13_mass_gen, tlv_sub13.M());
      event.set(h_sub23_mass_gen, tlv_sub23.M());
      event.set(h_sub1_px_gen, tlv_lead1_gen.Px());
      event.set(h_sub1_py_gen, tlv_lead1_gen.Py());
      event.set(h_sub1_pz_gen, tlv_lead1_gen.Pz());
      event.set(h_sub1_E_gen,  tlv_lead1_gen.E());
      event.set(h_sub1_eta_gen,  subjets_gen[lead1_gen].eta());
      event.set(h_sub1_phi_gen,  subjets_gen[lead1_gen].phi());
      event.set(h_sub2_px_gen, tlv_lead2_gen.Px());
      event.set(h_sub2_py_gen, tlv_lead2_gen.Py());
      event.set(h_sub2_pz_gen, tlv_lead2_gen.Pz());
      event.set(h_sub2_E_gen,  tlv_lead2_gen.E());
      event.set(h_sub2_eta_gen,  subjets_gen[lead2_gen].eta());
      event.set(h_sub2_phi_gen,  subjets_gen[lead2_gen].phi());
      event.set(h_sub3_px_gen, tlv_lead3_gen.Px());
      event.set(h_sub3_py_gen, tlv_lead3_gen.Py());
      event.set(h_sub3_pz_gen, tlv_lead3_gen.Pz());
      event.set(h_sub3_E_gen,  tlv_lead3_gen.E());
      event.set(h_sub3_eta_gen,  subjets_gen[lead3_gen].eta());
      event.set(h_sub3_phi_gen,  subjets_gen[lead3_gen].phi());
    }
  }
  if(!fill_gen){
    event.set(h_minimal_dR_gen, -1);
    event.set(h_sub1_pt_gen, -1);
    event.set(h_sub2_pt_gen, -1);
    event.set(h_sub3_pt_gen, -1);
    event.set(h_sub12_dR_gen, -1);
    event.set(h_sub23_dR_gen, -1);
    event.set(h_sub13_dR_gen, -1);
    event.set(h_sub1_mass_gen, -1);
    event.set(h_sub2_mass_gen, -1);
    event.set(h_sub3_mass_gen, -1);
    event.set(h_sub12_mass_gen, -1);
    event.set(h_sub23_mass_gen, -1);
    event.set(h_sub13_mass_gen, -1);
    event.set(h_sub1_px_gen, -1);
    event.set(h_sub1_py_gen, -1);
    event.set(h_sub1_pz_gen, -1);
    event.set(h_sub1_E_gen, -1);
    event.set(h_sub1_eta_gen, -1);
    event.set(h_sub1_phi_gen, -1);
    event.set(h_sub2_px_gen, -1);
    event.set(h_sub2_py_gen, -1);
    event.set(h_sub2_pz_gen, -1);
    event.set(h_sub2_E_gen, -1);
    event.set(h_sub2_eta_gen, -1);
    event.set(h_sub2_phi_gen, -1);
    event.set(h_sub3_px_gen, -1);
    event.set(h_sub3_py_gen, -1);
    event.set(h_sub3_pz_gen, -1);
    event.set(h_sub3_E_gen, -1);
    event.set(h_sub3_eta_gen, -1);
    event.set(h_sub3_phi_gen, -1);
  }

  // Match reco jets to gen jets to get the same ordering --------------------------------------------------------------
  // Index of matched reco jet to gen jet
  int matched_index1 = -1;
  int matched_index2 = -1;
  int matched_index3 = -1;
  if(isTTbar_){
    // dR_<REC CLUSTERED><GEN LEAD>
    if(hadjets_gen[0].subjets().size()>=3 && !mismatch){
      double dR_11 = deltaR(subjets_pt_sorted.at(0), subjets_gen_pt_sorted.at(0));
      double dR_12 = deltaR(subjets_pt_sorted.at(0), subjets_gen_pt_sorted.at(1));
      double dR_13 = deltaR(subjets_pt_sorted.at(0), subjets_gen_pt_sorted.at(2));
      double dR_21 = deltaR(subjets_pt_sorted.at(1), subjets_gen_pt_sorted.at(0));
      double dR_22 = deltaR(subjets_pt_sorted.at(1), subjets_gen_pt_sorted.at(1));
      double dR_23 = deltaR(subjets_pt_sorted.at(1), subjets_gen_pt_sorted.at(2));
      double dR_31 = deltaR(subjets_pt_sorted.at(2), subjets_gen_pt_sorted.at(0));
      double dR_32 = deltaR(subjets_pt_sorted.at(2), subjets_gen_pt_sorted.at(1));
      double dR_33 = deltaR(subjets_pt_sorted.at(2), subjets_gen_pt_sorted.at(2));
      if (dR_11<=dR_31 && dR_11<=dR_21){
        if(dR_11 < 0.2){ // if lowest distance >0.2, all distances are
          matched_index1=1;
          // remove first rec and gen jet from list
          if(dR_22<=dR_32) {
            if(dR_22<0.2) matched_index2=2;
            if(dR_33<0.2) matched_index3=3;
          }
          else{
            if(dR_32<0.2) matched_index2=3;
            if(dR_23<0.2) matched_index3=2;
          }
        }
      }
      else if(dR_21<=dR_31 && dR_21<=dR_11){
        if(dR_21 < 0.2){
          matched_index1=2;
          // remove first rec and gen jet from list
          if(dR_12<=dR_32) {
            if(dR_12<0.2) matched_index2=1;
            if(dR_33<0.2) matched_index3=3;
          }
          else{
            if(dR_32<0.2) matched_index2=3;
            if(dR_23<0.2) matched_index3=1;
          }
        }
      }
      else if(dR_31<=dR_21 && dR_31<=dR_11){
        if(dR_31 < 0.2){
          matched_index1=3;
          // remove first rec and gen jet from list
          if(dR_12<=dR_22) {
            if(dR_12<0.2) matched_index2=1;
            if(dR_23<0.2) matched_index3=2;
          }
          else{
            if(dR_22<0.2) matched_index2=2;
            if(dR_13<0.2) matched_index3=1;
          }
        }
      }
      else{
        printf("dR11 %3.2f dR12 %3.2f dR13 %3.2f\n",dR_11,dR_21,dR_31);
        throw runtime_error("No jet is the clostest to leading!");
      }
      // if((matched_index1==matched_index2 || matched_index1==matched_index3 || matched_index2==matched_index3)){
      //   printf("matched_index1 %3d matched_index2 %3d matched_index3 %3d\n",matched_index1,matched_index2,matched_index3);
      //   throw runtime_error("[ERROR] more then one leading reco jet!");
      // }
    }
  }
  if(debug_) printf("... Start with matched indices\n");
  event.set(h_matched_index1, matched_index1);
  event.set(h_matched_index2, matched_index2);
  event.set(h_matched_index3, matched_index3);

  if(debug_) printf("... Start with same leading order\n");
  bool same_leading_order = matched_index1==1 && matched_index2==2 && matched_index3==3;
  event.set(h_same_leading_order,same_leading_order);

  // ==================
  // = Store Matching =
  // ==================

  bool all_subjets_matched_rec = false;
  bool all_subjets_matched_gen = false;
  int index_parton_b_rec = -1;
  int index_parton_b_gen = -1;
  double top_mass = 172.5;
  if(isTTbar_){
    const auto & ttbargen = event.get(h_ttbargen);

    GenParticle top, bot, w1, w2;
    double dR_top_had = deltaR(ttbargen.Top(),hadjets.at(0));
    double dR_antitop_had = deltaR(ttbargen.Antitop(),hadjets.at(0));
    if(dR_top_had < dR_antitop_had){
      top = ttbargen.Top();
      bot = ttbargen.bTop();
      w1 = ttbargen.Wdecay1();
      w2 = ttbargen.Wdecay2();
    }
    else{
      top = ttbargen.Antitop();
      bot = ttbargen.bAntitop();
      w1 = ttbargen.WMinusdecay1();
      w2 = ttbargen.WMinusdecay2();
    }

    LorentzVector lv_top = top.v4();
    TLorentzVector tlv_top(lv_top.Px(), lv_top.Py(), lv_top.Pz(), lv_top.E());
    top_mass = tlv_top.M();

    bool b_matched_rec = false;
    bool w1_matched_rec = false;
    bool w2_matched_rec = false;
    bool b_matched_gen = false;
    bool w1_matched_gen = false;
    bool w2_matched_gen = false;
    for(unsigned int i=0; i<subjets.size(); i++){
      // for(auto& s: subjets){
      Jet s = subjets[i];
      double dR_b = deltaR(s, bot); 
      double dR_w1 = deltaR(s, w1); 
      double dR_w2 = deltaR(s, w2);
      if(dR_b<0.4 && dR_w1>0.4 && dR_w2>0.4){
        b_matched_rec = true;
        index_parton_b_rec = i+1;
      }
      else if(dR_b>0.4 && dR_w1<0.4 && dR_w2>0.4) w1_matched_rec = true;
      else if(dR_b>0.4 && dR_w1>0.4 && dR_w2<0.4) w2_matched_rec = true;
    }
    for(unsigned int i=0; i<hadjets_gen[0].subjets().size(); i++){
      // for(auto& s: hadjets_gen[0].subjets()){
      GenJet s = hadjets_gen[0].subjets()[i];
      double dR_b = deltaR(s, bot); 
      double dR_w1 = deltaR(s, w1); 
      double dR_w2 = deltaR(s, w2);
      if(dR_b<0.4 && dR_w1>0.4 && dR_w2>0.4){
        b_matched_gen = true;
        index_parton_b_gen = i+1;
      }
      else if(dR_b>0.4 && dR_w1<0.4 && dR_w2>0.4) w1_matched_gen = true;
      else if(dR_b>0.4 && dR_w1>0.4 && dR_w2<0.4) w2_matched_gen = true;
    }
    all_subjets_matched_rec = b_matched_rec && w1_matched_rec && w2_matched_rec;
    all_subjets_matched_gen = b_matched_gen && w1_matched_gen && w2_matched_gen;
  }
  if(debug_) printf("... set up matched handesl\n");
  event.set(h_all_subjets_matched_rec, all_subjets_matched_rec);
  event.set(h_all_subjets_matched_gen, all_subjets_matched_gen);
  event.set(h_index_parton_b_rec, index_parton_b_rec);
  event.set(h_index_parton_b_gen, index_parton_b_gen);
  event.set(h_mass_hadtop, top_mass);

}
