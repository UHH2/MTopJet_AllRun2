#include "UHH2/MTopJet/include/CorrectionHists_allHad.h"
#include <vector>

CorrectionHists_allHad::CorrectionHists_allHad(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // setup pt and eta bins
  ptgen_binning = {0, 50, 70, 90, 120, 150, 180, 250, 350, 10000};
  ptrec_binning = {0, 50, 70, 90, 120, 150, 180, 250, 350, 10000};
  etarec_binning = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191};


  // setup hists for every pt-eta bin
  unsigned int no_ptbins = 9; // rows
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
      pt_rec[pt_bin][eta_bin] = book<TH1F>("PtRec_"+pt_string +eta_string, "p^{rec}_{T, jet}", 200, 0, 1000);
      pt_gen[pt_bin][eta_bin] = book<TH1F>("PtGen_"+pt_string +eta_string, "p^{gen}_{T, jet}", 200, 0, 1000);
      pt_eta[pt_bin][eta_bin] = book<TH2F>("pt_eta_"+pt_string +eta_string, "x=pt y=eta",  50, 0, 500, 20, -6, 6);
      event_count[pt_bin][eta_bin] = book<TH1F>("Count_"+pt_string +eta_string, "a.u.", 1, 0.5, 1.5);
    }
  }

  book<TH2F>("TopJetMass1_TopJetMass2", "x=M_Top1 y=M_Top2", 30, 0, 300., 30, 0, 300.);


  // handle for clustered jets
  h_genjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");
  h_recjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
}



void CorrectionHists_allHad::fill(const Event & event){

  // get weight
  double weight = event.weight;
  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> recjets = event.get(h_recjets);
  std::vector<GenTopJet> genjets = event.get(h_genjets);

  // do not continus if a jet collection is empty
  if(recjets.size() < 2) return;
  for(unsigned int i=0; i<recjets.size(); i++){
    if(recjets.at(i).subjets().size() < 3) return;
  }
  if(genjets.size() < 2) return;
  for(unsigned int i=0; i<genjets.size(); i++){
    if(genjets.at(i).subjets().size() < 3) return;
  }

  // first match fatjets in gen and reco (find index of genjet to reco jet at i=0)
  double dRfat = 100;
  double dRfat_temp = 100;
  int i1 = 0;
  int i2 = 1;
  for(unsigned int i=0; i<genjets.size(); i++){
    dRfat_temp = deltaR(recjets.at(0), genjets.at(i));
    if(dRfat_temp < dRfat){
      dRfat = dRfat_temp;
      i1 = i;
    }
  }
  if(i1 == 0) i2=1;
  else        i2=0;

  bool use_jet1 = false;
  bool use_jet2 = false;
  if( deltaR(recjets.at(0), genjets.at(i1)) < 1.2 ) use_jet1 = true;
  if( deltaR(recjets.at(1), genjets.at(i2)) < 1.2 ) use_jet2 = true;

  std::vector<Jet> rec_sub1 = recjets.at(0).subjets();
  std::vector<Jet> rec_sub2 = recjets.at(1).subjets();
  std::vector<Particle> gen_sub1 = genjets.at(i1).subjets();
  std::vector<Particle> gen_sub2 = genjets.at(i2).subjets();


  // do pt cut on subjets and fatjets
  // for(unsigned int i=0; i<rec_sub1.size(); i++){
  //   if(rec_sub1.at(i).v4().Pt() < 30) use_jet1 = false;
  // }
  // for(unsigned int i=0; i<rec_sub2.size(); i++){
  //   if(rec_sub2.at(i).v4().Pt() < 30) use_jet2 = false;
  // }
  if(recjets.at(0).v4().Pt() < 350) use_jet1 = false;
  if(recjets.at(1).v4().Pt() < 350) use_jet2 = false;


  if(!use_jet1 && !use_jet2) return;


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
  double rec_pt;
  double rec_eta;
  double R;
  int ptrec_bin = 100;
  int ptgen_bin = 100;
  int etarec_bin = 100;
  // do matching
  if(use_jet1){
    for(unsigned int i=0; i<rec_sub1.size(); i++){
      dR = 1000;
      nearest_j = 100;
      for(unsigned int j=0; j<gen_sub1.size(); j++){
        dR_temp = uhh2::deltaR(rec_sub1.at(i), gen_sub1.at(j));
        if(dR_temp < dR){
          dR = dR_temp;
          nearest_j = j;
        }
      }
      gen_pt=gen_sub1.at(nearest_j).v4().Pt();
      rec_pt=rec_sub1.at(i).v4().Pt();
      rec_eta=rec_sub1.at(i).v4().Eta();
      if(rec_eta < 0) rec_eta *= -1; // nur Betrag von eta nutzen
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
  }
  if(use_jet2){
    for(unsigned int i=0; i<rec_sub2.size(); i++){
      dR = 1000;
      nearest_j = 100;
      for(unsigned int j=0; j<gen_sub2.size(); j++){
        dR_temp = uhh2::deltaR(rec_sub2.at(i), gen_sub2.at(j));
        if(dR_temp < dR){
          dR = dR_temp;
          nearest_j = j;
        }
      }
      gen_pt=gen_sub2.at(nearest_j).v4().Pt();
      rec_pt=rec_sub2.at(i).v4().Pt();
      rec_eta=rec_sub2.at(i).v4().Eta();
      if(rec_eta < 0) rec_eta *= -1; // nur Betrag von eta nutzen
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
  }

}
