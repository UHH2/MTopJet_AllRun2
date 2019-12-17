#include <UHH2/MTopJet/include/GenHists_process.h>

// STIL INCLUDE NUMBER OF T W B AND NUMBER OF LEPTON

GenHists_process::GenHists_process(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  number_each_particle = book<TH1F>("number_each_particle", "number of each particle", 25, 0, 25);
  number_each_particle_cut = book<TH1F>("number_each_particle_cut", "number of each particle after cut", 25, 0, 25);

  mother_muon_1 = book<TH1F>("mother1_id_muon", "muon mother_1 Id", 25, 0, 25);
  mother_muon_1_cut = book<TH1F>("mother1_id_muon_cut", "muon mother Id after cut", 25, 0, 25);
  mother_muon_1_outside = book<TH1F>("mother1_id_muon_outside", "muon mother Id outside of cut", 25, 0, 25);
  mother_muon_2 = book<TH1F>("mother2_id_muon", "muon mother_2 Id", 25, 0, 25);

  daughter_muon = book<TH1F>("daughter_id_muon", "muon daughter Id", 25, 0, 25);
  daughter_muon_cut = book<TH1F>("daughter_id_muon_cut", "muon daughter Id after cut", 25, 0, 25);

  number_muon = book<TH1F>("number_muon", "number of muon", 10, 0, 10);
  number_muon_cut = book<TH1F>("number_muon_cut", "number of muon after cut", 10, 0, 10);

  number_top = book<TH1F>("number_top", "number of top", 10, 0, 10);
  number_top_cut = book<TH1F>("number_top_cut", "number of top after cut", 10, 0, 10);

  number_bottom = book<TH1F>("number_bottom", "number of bottom", 10, 0, 10);
  number_bottom_cut = book<TH1F>("number_bottom_cut", "number of bottom after cut", 10, 0, 10);

  number_W = book<TH1F>("number_W", "number of W", 10, 0, 10);
  number_W_cut = book<TH1F>("number_W_cut", "number of W after cut", 10, 0, 10);

  number_genparticles = book<TH1F>("number_genparticles", "number of genparticles", 300, 0, 300);

  only_one_lepton_cut = book<TH1F>("only_one_lepton_cut", "only one muon or elec after cut", 2, 0, 2);
  only_one_lepton = book<TH1F>("only_one_lepton", "only one muon or elec", 2, 0, 2);

  all_status  = book<TH1F>("all_status", "each status", 30, 0, 30);
  status_muon = book<TH1F>("muon_status", "muon status", 30, 0, 30);
  status_muon_cut = book<TH1F>("muon_status_cut", "muon status after cut", 30, 0, 30);

  index_status = book<TH2F>("index_status", "x - index and y - status", 30, 0, 30, 70, 0, 70);

  deltaR_muon1_muon2 = book<TH1F>("deltaR_muon1_muon2", "deltaR m1 m2", 30, 0, 6);
  deltaR_muon1_muon3 = book<TH1F>("deltaR_muon1_muon3", "deltaR m1 m3", 30, 0, 6);
  deltaR_muon1_muon4 = book<TH1F>("deltaR_muon1_muon4", "deltaR m1 m4", 30, 0, 6);
  deltaR_muon1_muon5 = book<TH1F>("deltaR_muon1_muon5", "deltaR m1 m5", 30, 0, 6);
  deltaR_muon1_muon6 = book<TH1F>("deltaR_muon1_muon6", "deltaR m1 m6", 30, 0, 6);
  deltaR_muon1_muon7 = book<TH1F>("deltaR_muon1_muon7", "deltaR m1 m7", 30, 0, 6);
}

//---------------------------------------------------------------------------------------
//--------------------------------- Fill Hists here -------------------------------------
//---------------------------------------------------------------------------------------

void GenHists_process::fill(const Event & event){

  double weight = event.weight;  // get weight

  int count_muon = 0;
  int count_elec = 0;
  int count_top = 0;
  int count_bottom = 0;
  int count_w = 0;
  int gen_status = 0;
  int gen_index = 0;
  int gen_id= 0;

  bool is_muon = false;
  bool one_lepton = false;
  bool muon_1_in = false;

  int m1_id, m2_id, d1_id;

  GenParticle muon_1;

  //std::vector<GenParticle>* genparticle_after_cut;
  std::vector<GenParticle>* genparts = event.genparticles;

  unsigned int number_gen = genparts->size();
  number_genparticles->Fill(number_gen, weight);

  for (unsigned int i=0; i<number_gen; ++i){
    GenParticle genp = genparts->at(i);

    // ================================ VARIABLES =============================================================
    gen_status = genp.status(); // Get status
    gen_index = genp.index(); // Get index
    gen_id = abs(genp.pdgId()); // Get absolute Id

    //Get mother particles
    // x.mother(y, z) gives the GenParticle and x.mother1/2 gives the position in GenParticle list of mother
    // Initial Particles AND Particles after Index 13 (All of them ???) have no mother
    auto m1 = genp.mother(genparts, 1);
    auto m2 = genp.mother(genparts, 2);
    auto d1 = genp.daughter(genparts, 1);
    // set id to 99999 to be out of range of Hists
    if(genp.mother1() == 65535) m1_id = 0;
    else m1_id = abs(m1->pdgId());
    if(genp.mother2() == 65535) m2_id = 0;
    else m2_id = abs(m2->pdgId());

    //daughter1
    if(genp.daughter1() == 65535) d1_id = 0;
    else d1_id = abs(d1->pdgId());

    // count muon
    if(gen_id == 13){
      is_muon = true;
      count_muon += 1;

      // ================================ deltaR ==============================================================

      if((count_muon == 1) && (i < 13)){
        muon_1 = genp;
        muon_1_in = true;
      }
      if(muon_1_in){
        if(count_muon == 2) deltaR_muon1_muon2->Fill(deltaR(muon_1, genp), weight);
        if(count_muon == 3) deltaR_muon1_muon3->Fill(deltaR(muon_1, genp), weight);
        if(count_muon == 4) deltaR_muon1_muon4->Fill(deltaR(muon_1, genp), weight);
        if(count_muon == 5) deltaR_muon1_muon5->Fill(deltaR(muon_1, genp), weight);
        if(count_muon == 6) deltaR_muon1_muon6->Fill(deltaR(muon_1, genp), weight);
        if(count_muon == 7) deltaR_muon1_muon7->Fill(deltaR(muon_1, genp), weight);
      }
    }
    else is_muon = false;

    if(gen_id == 11) count_elec += 1;     // count elec
    if(gen_id == 6) count_top += 1;       // count top
    if(gen_id == 5) count_bottom += 1;    // count bottom
    if(gen_id == 24) count_w += 1;        // count W

    // Only one muon or elec
    if((count_elec + count_muon) == 1) one_lepton = true;
    else one_lepton = false;

    // ================================ CUT ===================================================================

    if(i < 13){
      number_each_particle_cut->Fill(gen_id, weight);
      if(is_muon && (m1_id != 0)){
        mother_muon_1_cut->Fill(m1_id, weight);
        daughter_muon_cut->Fill(d1_id, weight);
      }

      if(i == 12){
        only_one_lepton_cut->Fill(one_lepton, weight);
        number_muon_cut->Fill(count_muon, weight);
        number_top_cut->Fill(count_top, weight);
        number_bottom_cut->Fill(count_bottom, weight);
        number_W_cut->Fill(count_w, weight);
      }
    }
    // ================================ NO CUT ================================================================

    number_each_particle->Fill(gen_id, weight);
    if(is_muon){
      mother_muon_1->Fill(m1_id, weight);
      mother_muon_2->Fill(m2_id, weight);
      daughter_muon->Fill(d1_id, weight);
    }
    if(i>12) mother_muon_1_outside->Fill(m1_id, weight);

    // ================================ STATUS ================================================================

    all_status->Fill(gen_status, weight);
    index_status->Fill(gen_index, gen_status, weight);

    if(gen_id == 13){
      if(i < 13) status_muon_cut->Fill(gen_status, weight);
      status_muon->Fill(gen_status, weight);
    }
  }

  // ================================ END OF LOOP ============================================================

  only_one_lepton->Fill(one_lepton, weight);
  number_muon->Fill(count_muon, weight);
  number_top->Fill(count_top, weight);
  number_bottom->Fill(count_bottom, weight);
  number_W->Fill(count_w, weight);

  // determine deltaR

}
