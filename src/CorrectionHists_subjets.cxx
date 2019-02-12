#include "UHH2/MTopJet/include/CorrectionHists_subjets.h"
#include <vector>

CorrectionHists_subjets::CorrectionHists_subjets(uhh2::Context & ctx, const std::string & dirname, const std::string & type): Hists(ctx, dirname){

  // setup pt and eta bins
  ptgen_binning = {0, 80, 130, 180, 250, 350, 10000};
  ptrec_binning = {0, 80, 130, 180, 250, 350, 10000};
  etarec_binning = {-10, -1.5, -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 10};

  // setup hists for every pt-eta bin
  unsigned int no_ptbins = 6; // rows
  unsigned int no_etabins = 12; // columns

  TH1F* initial_value;
  TH2F* initial_value_2d;
  pt_reso.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  pt_rec.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  pt_gen.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  pt_eta.resize(no_ptbins, std::vector<TH2F*>(no_etabins, initial_value_2d));
  event_count.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));

  for(unsigned int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(unsigned int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      std::string pt_string = std::to_string(pt_bin);
      std::string eta_string = std::to_string(eta_bin);
      pt_reso[pt_bin][eta_bin] = book<TH1F>("PtReso_"+pt_string +eta_string, "p^{rec}_{T, jet} / p^{gen}_{T, jet}", 90, -1, 3);
      pt_rec[pt_bin][eta_bin] = book<TH1F>("PtRec_"+pt_string +eta_string, "p^{rec}_{T, jet}", 100, 0, 1000);
      pt_gen[pt_bin][eta_bin] = book<TH1F>("PtGen_"+pt_string +eta_string, "p^{gen}_{T, jet}", 100, 0, 1000);
      pt_eta[pt_bin][eta_bin] = book<TH2F>("pt_eta_"+pt_string +eta_string, "x=pt y=eta",  50, 0, 500, 20, -4, 4);
      event_count[pt_bin][eta_bin] = book<TH1F>("Count_"+pt_string +eta_string, "a.u.", 1, 0.5, 1.5);
    }
  }

  book<TH2F>("TopJetMass1_TopJetMass2", "x=M_Top1 y=M_Top2", 30, 0, 300., 30, 0, 300.);


  // handle for clustered jets
  h_recjets_noJEC=ctx.get_handle<std::vector<TopJet>>("xconeCHS_noJEC");
  h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");
  h_genjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");
  h_hadgenjets=ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");

  if(type == "jec")     h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  else if(type == "raw")h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_noJEC");
  else if(type == "cor")h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");{

  }


}



void CorrectionHists_subjets::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> recjets = event.get(h_recjets);
  std::vector<TopJet> recjets_noJEC = event.get(h_recjets_noJEC);
  std::vector<TopJet> hadjets = event.get(h_hadjets);

  std::vector<GenTopJet> genjets = event.get(h_genjets);
  std::vector<GenTopJet> hadgenjets = event.get(h_hadgenjets);

  // do not continus if a jet collection is empty
  if(recjets.size() < 1) return;
  for(unsigned int i=0; i<recjets.size(); i++){
    if(recjets.at(i).subjets().size() < 1) return;
  }
  if(recjets_noJEC.size() < 1) return;
  for(unsigned int i=0; i<recjets_noJEC.size(); i++){
    if(recjets_noJEC.at(i).subjets().size() < 1) return;
  }
  if(hadjets.size() < 1) return;
  if(genjets.size() < 1) return;
  for(unsigned int i=0; i<genjets.size(); i++){
    if(genjets.at(i).subjets().size() < 1) return;
  }
  if(hadgenjets.size() < 1) return;




  // first find had. jet and only use those subjets (without using JEC)
  double dRrec = 100;
  double dRrec_temp = 100;
  int i_had = 0;
  for(unsigned int i=0; i<recjets_noJEC.size(); i++){
    dRrec_temp = deltaR(recjets_noJEC.at(i), hadjets.at(0));
    if(dRrec_temp < dRrec){
      dRrec = dRrec_temp;
      i_had = i;
    }
  }

  std::vector<Jet> rec_sub_noJEC = recjets_noJEC.at(i_had).subjets();
  std::vector<Jet> rec_sub = recjets.at(i_had).subjets();

  // also get had subjets on gen level
  double dRgen = 100;
  double dRgen_temp = 100;
  int i_genhad = 0;
  for(unsigned int i=0; i<genjets.size(); i++){
    dRgen_temp = deltaR(genjets.at(i), hadgenjets.at(0));
    if(dRgen_temp < dRgen){
      dRgen = dRgen_temp;
      i_genhad = i;
    }
  }

  std::vector<Particle> gen_sub = genjets.at(i_genhad).subjets();

  // do pt cut on subjets without JEC, only continue for pt > 30
  bool pt_valid = true;
  for(unsigned int i=0; i<rec_sub_noJEC.size(); i++){
    if(rec_sub_noJEC.at(i).v4().Pt() < 30) pt_valid = false;
  }
  if(!pt_valid) return;


  /* ******************************************************************************
  matching between gen and reco jets:
  - for each reco jet, calc distance to all other gen jets
  - gen jet with lowest distance is a match if distance is < 0.2
  - then calculate resolution with reco jet and matched gen jet
  - to do: account for double counting
  ********************************************************************************* */

  double dR;
  double dR_temp;
  int nearest_j;
  double gen_pt;
  //double gen_eta;
  double rec_pt;
  double rec_eta;
  double R;
  int ptrec_bin = 100;
  int ptgen_bin = 100;
  int etarec_bin = 100;
  // do matching
  for(unsigned int i=0; i<rec_sub.size(); i++){
    dR = 1000;
    nearest_j = 100;
    for(unsigned int j=0; j<gen_sub.size(); j++){
      dR_temp = uhh2::deltaR(rec_sub.at(i), gen_sub.at(j));
      if(dR_temp < dR){
        dR = dR_temp;
        nearest_j = j;
      }
    }
    gen_pt=gen_sub.at(nearest_j).v4().Pt();
    rec_pt=rec_sub.at(i).v4().Pt();
    rec_eta=rec_sub.at(i).v4().Eta();
    if(gen_pt == 0) continue;
    R = rec_pt/gen_pt;
    if(nearest_j != 100 && dR <= 0.2){
      // bins in pt_rec
      for(unsigned int i = 0; i < (ptrec_binning.size() - 1); i++){
        if(rec_pt > ptrec_binning[i] && rec_pt <= ptrec_binning[i+1]) ptrec_bin = i;
      }
      for(unsigned int i = 0; i < (ptgen_binning.size() - 1); i++){
        if(gen_pt > ptgen_binning[i] && gen_pt <= ptgen_binning[i+1]) ptgen_bin = i;
      }

      for(unsigned int i = 0; i < (etarec_binning.size() - 1); i++){
        if(rec_eta > etarec_binning[i] && rec_eta< etarec_binning[i+1]) etarec_bin = i;
      }

      if(ptrec_bin != 100 && etarec_bin != 100){
        pt_reso[ptgen_bin][etarec_bin]->Fill(R, weight);
        pt_rec[ptgen_bin][etarec_bin]->Fill(rec_pt, weight);
        pt_gen[ptgen_bin][etarec_bin]->Fill(gen_pt, weight);
        pt_eta[ptgen_bin][etarec_bin]->Fill(rec_pt, rec_eta, weight);
        event_count[ptgen_bin][etarec_bin]->Fill(1, weight);

      }
    }

  }

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

}
