#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/TTbarGenHists.h>

#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/MTopJetSelections.h>
#include <UHH2/MTopJet/include/MTopJetGenSelections.h>
#include <UHH2/MTopJet/include/MTopJetGenHists.h>
#include <UHH2/MTopJet/include/MTopJetRecoHists.h>
#include <UHH2/MTopJet/include/MTopJetRecoGenHists.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/JetCluster.h>
class MTopJetGenPostSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetGenPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:

  // cleaners
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;
  // std::unique_ptr<JetLeptonCleaner_by_KEYmatching> topjetlepton_cleaner;
  // selections
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_gen;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_rec;

  std::unique_ptr<uhh2::AnalysisModule> jetprod; // to produce Jets
  std::unique_ptr<uhh2::Selection> matching;
  std::unique_ptr<uhh2::Selection> topjetpt;
  std::unique_ptr<uhh2::Selection> masscut;
  std::unique_ptr<uhh2::Selection> n_genjets;
  std::unique_ptr<uhh2::Selection> deltaR;
  std::unique_ptr<uhh2::Selection> n_recjets;
  std::unique_ptr<uhh2::Selection> topjetpt_rec;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> deltaR_rec;
  std::unique_ptr<uhh2::Selection> masscut_rec;
  std::unique_ptr<uhh2::Selection> massbin1;
  std::unique_ptr<uhh2::Selection> massbin2;
  std::unique_ptr<uhh2::Selection> massbin3;
  std::unique_ptr<uhh2::Selection> massbin4;
  std::unique_ptr<uhh2::Selection> massbin5;
  std::unique_ptr<uhh2::Selection> massbin6;
  std::unique_ptr<uhh2::Selection> massbin7;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_GenHists0a, h_GenHists0b, h_GenHists1, h_GenHists2, h_GenHists3, h_GenHists4, h_GenHists1_m, h_GenHists2_m, h_GenHists3_m, h_GenHists4_m, h_GenHists1_u, h_GenHists2_u, h_GenHists3_u, h_GenHists4_u;
  std::unique_ptr<Hists> h_RecHists0a, h_RecHists0b, h_RecHists0c, h_RecHists0d, h_RecHists1, h_RecHists2, h_RecHists3, h_RecHists4;
  std::unique_ptr<Hists> h_RecGenHists0, h_RecGenHists1, h_RecGenHists2, h_RecGenHists3, h_RecGenHists4, h_RecGenHists5, h_RecGenHists6, h_RecGenHists7;

};

MTopJetGenPostSelectionModule::MTopJetGenPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  // const bool isMC = (ctx.get("dataset_type") == "MC");
  
  // set up which jet to use
  const std::string jet_label_gen("ak08_gen");
  const std::string jet_label_rec("ak08_rec");
  float jet_radius = 0.8;

  ////

  //// COMMON MODULES

  const std::string ttbar_gen_label("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  // To produce Jets:
  jetprod.reset(new RecoJetProducer(ctx, jet_label_rec, 150, jet_radius)); // set to akt

  ////

  //// OBJ CLEANING
  const     MuonId muoSR(AndId<Muon>    (PtEtaCut  (45., 2.1), MuonIDMedium()));
  const ElectronId eleSR(AndId<Electron>(PtEtaSCCut(45., 2.5), ElectronID_Spring15_25ns_tight_noIso));
  muoSR_cleaner.reset(new     MuonCleaner(muoSR));
  eleSR_cleaner.reset(new ElectronCleaner(eleSR));
  cleaner_gen.reset(new GenJetLeptonCleaner(ctx, jet_label_gen, jet_radius));
  cleaner_rec.reset(new RecoJetLeptonCleaner(ctx, jet_label_rec, jet_radius));

  // set up Hists classes:
  // GEN
  h_GenHists0a.reset(new MTopJetGenHists(ctx, "00a_GenHists_before_Cleaner", jet_label_gen));
  h_GenHists0b.reset(new MTopJetGenHists(ctx, "00b_GenHists_after_Cleaner", jet_label_gen));

  h_GenHists1.reset(new MTopJetGenHists(ctx, "01_GenHists_JetN", jet_label_gen));
  h_GenHists1_m.reset(new MTopJetGenHists(ctx, "01_GenHists_JetN_matched", jet_label_gen));
  h_GenHists1_u.reset(new MTopJetGenHists(ctx, "01_GenHists_JetN_unmatched", jet_label_gen));
  h_GenHists2.reset(new MTopJetGenHists(ctx, "02_GenHists_pT_Cut", jet_label_gen));
  h_GenHists2_m.reset(new MTopJetGenHists(ctx, "02_GenHists_pT_Cut_matched", jet_label_gen));
  h_GenHists2_u.reset(new MTopJetGenHists(ctx, "02_GenHists_pT_Cut_unmatched", jet_label_gen));
  h_GenHists3.reset(new MTopJetGenHists(ctx, "03_GenHists_deltaR", jet_label_gen));
  h_GenHists3_m.reset(new MTopJetGenHists(ctx, "03_GenHists_deltaR_matched", jet_label_gen));
  h_GenHists3_u.reset(new MTopJetGenHists(ctx, "03_GenHists_deltaR_unmatched", jet_label_gen));
  h_GenHists4.reset(new MTopJetGenHists(ctx, "04_GenHists_masscut", jet_label_gen));
  h_GenHists4_m.reset(new MTopJetGenHists(ctx, "04_GenHists_masscut_matched", jet_label_gen));
  h_GenHists4_u.reset(new MTopJetGenHists(ctx, "04_GenHists_masscut_unmatched", jet_label_gen));

  //RECO
  h_RecHists0a.reset(new MTopJetRecoHists(ctx, "00a_RecHists_before_Cleaner", jet_label_rec));
  h_RecHists0b.reset(new MTopJetRecoHists(ctx, "00b_RecHists_after_Cleaner", jet_label_rec));
  h_RecHists0c.reset(new MTopJetRecoHists(ctx, "00c_RecHists_LepSel", jet_label_rec));
  h_RecHists0d.reset(new MTopJetRecoHists(ctx, "00d_RecHists_background_reduction", jet_label_rec));

  h_RecHists1.reset(new MTopJetRecoHists(ctx, "01_RecHists_JetN", jet_label_rec));
  h_RecHists2.reset(new MTopJetRecoHists(ctx, "02_RecHists_pT_Cut", jet_label_rec));
  h_RecHists3.reset(new MTopJetRecoHists(ctx, "03_RecHists_deltaR", jet_label_rec));
  h_RecHists4.reset(new MTopJetRecoHists(ctx, "04_RecHists_masscut", jet_label_rec));

  //GEN+RECO (resolution)
  h_RecGenHists0.reset(new MTopJetRecoGenHists(ctx, "00_RecGenHists_all", jet_label_rec, jet_label_gen));
  h_RecGenHists1.reset(new MTopJetRecoGenHists(ctx, "01_RecGenHists_000to050", jet_label_rec, jet_label_gen));
  h_RecGenHists2.reset(new MTopJetRecoGenHists(ctx, "02_RecGenHists_050to100", jet_label_rec, jet_label_gen));
  h_RecGenHists3.reset(new MTopJetRecoGenHists(ctx, "03_RecGenHists_100to150", jet_label_rec, jet_label_gen));
  h_RecGenHists4.reset(new MTopJetRecoGenHists(ctx, "04_RecGenHists_150to200", jet_label_rec, jet_label_gen));
  h_RecGenHists5.reset(new MTopJetRecoGenHists(ctx, "05_RecGenHists_200to250", jet_label_rec, jet_label_gen));
  h_RecGenHists6.reset(new MTopJetRecoGenHists(ctx, "06_RecGenHists_250to300", jet_label_rec, jet_label_gen));
  h_RecGenHists7.reset(new MTopJetRecoGenHists(ctx, "07_RecGenHists_300to500", jet_label_rec, jet_label_gen));

  // EVENT SELECTION
  // GEN
  n_genjets.reset(new NGenJets(ctx, jet_label_gen, 150, 2, 2));  // ==2 jets with pt > 150
  topjetpt.reset(new LeadingJetPT(ctx, jet_label_gen, 400)); // leading jet pt > 400
  deltaR.reset(new DeltaRCut(ctx, jet_label_gen, jet_radius)); 
  masscut.reset(new MassCut(ctx, jet_label_gen));
  matching.reset(new Matching(ctx, jet_label_gen, jet_radius));
  ////

  // RECO
  n_recjets.reset(new NRecoJets(ctx, jet_label_rec, 150, 2, 2));  // ==2 jets with pt > 150
  topjetpt_rec.reset(new LeadingRecoJetPT(ctx, jet_label_rec, 400)); // leading jet pt > 400
  met_sel.reset(new METCut  (20, uhh2::infinity));
  htlep_sel.reset(new HTlepCut(100, uhh2::infinity));
  twodcut_sel.reset(new TwoDCut1(.4, 40.));
  deltaR_rec.reset(new DeltaRCutReco(ctx, jet_label_rec, jet_radius)); 
  masscut_rec.reset(new MassCutReco(ctx, jet_label_rec));
  ////

  // Selection for Mass binning
  massbin1.reset(new MassCutGen1(ctx, jet_label_gen, 0, 50));
  massbin2.reset(new MassCutGen1(ctx, jet_label_gen, 50, 100));
  massbin3.reset(new MassCutGen1(ctx, jet_label_gen, 100, 150));
  massbin4.reset(new MassCutGen1(ctx, jet_label_gen, 150, 200));
  massbin5.reset(new MassCutGen1(ctx, jet_label_gen, 200, 250));
  massbin6.reset(new MassCutGen1(ctx, jet_label_gen, 250, 300));
  massbin7.reset(new MassCutGen1(ctx, jet_label_gen, 300, 500));
}

bool MTopJetGenPostSelectionModule::process(uhh2::Event& event){

  //  COMMON MODULES
  // set to true / false to run analysis
  bool produce_jet = false;
  bool do_gensel = true;
  bool do_recsel = false;

  ttgenprod->process(event);
 

  if(produce_jet){
    jetprod->process(event); // to produce jets
    return true;
  }

  h_GenHists0a->fill(event);
  // ---------------------------------------------------------------------------------------
  // ---------------------------- apply Gen Selection --------------------------------------
  // ---------------------------------------------------------------------------------------
  // Cleaner
  if(do_gensel){
    cleaner_gen->process(event);
    cleaner_rec->process(event);
  }
  h_GenHists0b->fill(event);

  // ==2 Jets with pT > 150
  if(do_gensel && !(n_genjets->passes(event))) return false;
  h_GenHists1->fill(event);

  if(do_gensel && matching->passes(event)){
    h_GenHists1_m->fill(event);
  }
  if(do_gensel && !(matching->passes(event))){
    h_GenHists1_u->fill(event);
  }
  // pT cut on leading Jet
  if(do_gensel && !(topjetpt->passes(event))) return false;
  h_GenHists2->fill(event);

  if(do_gensel && matching->passes(event)){
    h_GenHists2_m->fill(event);
  }
  if(do_gensel && !(matching->passes(event))){
    h_GenHists2_u->fill(event);
  }

  // deltaR(lepton, 2nd jet) < jet radius
  if(do_gensel && !(deltaR->passes(event))) return false;
  h_GenHists3->fill(event);

  if(do_gensel && matching->passes(event)){
    h_GenHists3_m->fill(event);
  }
  if(do_gensel && !(matching->passes(event))){
    h_GenHists3_u->fill(event);
  }

  // m(1st jet) > m(2nd jet + lepton)
  if(do_gensel && !(masscut->passes(event))) return false;
  h_GenHists4->fill(event);

  if(do_gensel && matching->passes(event)){
    h_GenHists4_m->fill(event);
  }
  if(do_gensel && !(matching->passes(event))){
    h_GenHists4_u->fill(event);
  }
  // ---------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------

  h_RecGenHists0->fill(event);
  if(massbin1->passes(event)) h_RecGenHists1->fill(event);
  if(massbin2->passes(event)) h_RecGenHists2->fill(event);
  if(massbin3->passes(event)) h_RecGenHists3->fill(event);
  if(massbin4->passes(event)) h_RecGenHists4->fill(event);
  if(massbin5->passes(event)) h_RecGenHists5->fill(event);
  if(massbin6->passes(event)) h_RecGenHists6->fill(event);
  if(massbin7->passes(event)) h_RecGenHists7->fill(event);

  // ---------------------------------------------------------------------------------------
  // ---------------------------- apply Reco Selection -------------------------------------
  // ---------------------------------------------------------------------------------------
  h_RecHists0a->fill(event);
  // Cleaner
  if(do_recsel){
    muoSR_cleaner->process(event);
    sort_by_pt<Muon>(*event.muons);

    eleSR_cleaner->process(event);
    sort_by_pt<Electron>(*event.electrons);

    cleaner_rec->process(event);
  }
  h_RecHists0b->fill(event);

  //select one elec/muon
  if(do_recsel){
    const bool pass_lep1 = (((event.muons->size() == 1)&&(event.electrons->size() == 0)) ||((event.muons->size() == 0)&&(event.electrons->size() == 1)));
    if(!pass_lep1) return false;
  }
  h_RecHists0c->fill(event);

  // reduction of background (MET, 2D, Ht_lep)
  if(do_recsel && !(met_sel->passes(event))) return false;
  if(do_recsel && !(twodcut_sel->passes(event))) return false;
  if(do_recsel && !(htlep_sel->passes(event))) return false;
  h_RecHists0d->fill(event);

  // ==2 Jets with pT > 150
  if(do_recsel && !(n_recjets->passes(event))) return false;
  h_RecHists1->fill(event);
  h_RecGenHists1->fill(event);

  // pT cut on leading Jet
  if(do_recsel && !(topjetpt_rec->passes(event))) return false;
  h_RecHists2->fill(event);

  // deltaR(lepton, 2nd jet) < jet radius
  if(do_recsel && !(deltaR_rec->passes(event))) return false;
  h_RecHists3->fill(event);

  // m(1st jet) > m(2nd jet + lepton)
  if(do_recsel && !(masscut_rec->passes(event))) return false;
  h_RecHists4->fill(event);
  // ---------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------


  return true; //false to delete all collections and keep hists
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetGenPostSelectionModule)
