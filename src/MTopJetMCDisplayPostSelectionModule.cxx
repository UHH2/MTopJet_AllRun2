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
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>

#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/RecoSelections_topjet.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/GenHists.h>
#include <UHH2/MTopJet/include/RecoHists.h>
#include <UHH2/MTopJet/include/RecoGenHists.h>
#include <UHH2/MTopJet/include/GenHists_topjet.h>
#include <UHH2/MTopJet/include/RecoHists_topjet.h>
#include <UHH2/MTopJet/include/RecoGenHists_topjet.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/ClusteringHists.h>
#include <UHH2/MTopJet/include/JetCluster.h>

class MTopJetMCDisplayPostSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetMCDisplayPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
   uhh2::Event::Handle<std::vector<Jet>>h_jets;


  // cleaners
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;
  std::unique_ptr<TopJetCorrector> topjet_corrector;
  std::unique_ptr<SubJetCorrector> topjet_subjet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer> topjetER_smearer;
  std::unique_ptr<TopJetCleaner> topjet_cleaner;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> topjetlepton_cleaner;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> jetcluster;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_gen;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_topgen;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_rec;
  std::unique_ptr<uhh2::AnalysisModule> cleaner_toprec;

  std::unique_ptr<uhh2::Selection> matching;
  std::unique_ptr<uhh2::Selection> matching_botlep_lep;
  std::unique_ptr<uhh2::Selection> matching_HOTVR;
  std::unique_ptr<uhh2::Selection> matching_XCone;
  // gen sel
  std::unique_ptr<uhh2::Selection> topjetpt;
  std::unique_ptr<uhh2::Selection> masscut;
  std::unique_ptr<uhh2::Selection> n_genjets;
  std::unique_ptr<uhh2::Selection> deltaR;
  std::unique_ptr<uhh2::Selection> deltaPhi;
  std::unique_ptr<uhh2::Selection> deltaR_HOTVR;
  std::unique_ptr<uhh2::Selection> matching_top;
  std::unique_ptr<uhh2::Selection> topjetpt_top;
  std::unique_ptr<uhh2::Selection> masscut_top;
  std::unique_ptr<uhh2::Selection> n_genjets_top;
  std::unique_ptr<uhh2::Selection> deltaR_top;
  // reco sel for XCone, HOTVR, ...
  std::unique_ptr<uhh2::Selection> n_recjets;
  std::unique_ptr<uhh2::Selection> topjetpt_rec;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> deltaR_rec;
  std::unique_ptr<uhh2::Selection> masscut_rec;
  // reco sel for topjets
  std::unique_ptr<uhh2::Selection> n_recjets_top;
  std::unique_ptr<uhh2::Selection> topjetpt_rec_top;
  std::unique_ptr<uhh2::Selection> deltaR_rec_top;
  std::unique_ptr<uhh2::Selection> masscut_rec_top;
  // seperate Mass bins for resolution studies
  std::unique_ptr<uhh2::Selection> massbin1;
  std::unique_ptr<uhh2::Selection> massbin2;
  std::unique_ptr<uhh2::Selection> massbin3;
  std::unique_ptr<uhh2::Selection> massbin4;
  std::unique_ptr<uhh2::Selection> massbin5;
  std::unique_ptr<uhh2::Selection> massbin6;
  std::unique_ptr<uhh2::Selection> massbin_HOTVR_low;
  std::unique_ptr<uhh2::Selection> massbin_HOTVR_peak;
  std::unique_ptr<uhh2::Selection> massbin_HOTVR_high;
  std::unique_ptr<uhh2::Selection> massbin1_top;
  std::unique_ptr<uhh2::Selection> massbin2_top;
  std::unique_ptr<uhh2::Selection> massbin3_top;
  std::unique_ptr<uhh2::Selection> massbin4_top;
  std::unique_ptr<uhh2::Selection> massbin5_top;
  std::unique_ptr<uhh2::Selection> massbin6_top;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_GenHists0a, h_GenHists0b, h_GenHists1, h_GenHists2, h_GenHists3, h_GenHists4, h_GenHists1_m, h_GenHists2_m, h_GenHists3_m, h_GenHists4_m,  h_GenHists4_m_XCone, h_GenHists1_u, h_GenHists2_u, h_GenHists3_u, h_GenHists4_u,  h_GenHists4_u_XCone, h_GenHists5_LowMass, h_GenHists5_PeakMass, h_GenHists5_HighMass;
  std::unique_ptr<Hists> h_GenHists0a_top, h_GenHists0b_top, h_GenHists1_top, h_GenHists2_top, h_GenHists3_top, h_GenHists4_top, h_GenHists1_m_top, h_GenHists2_m_top, h_GenHists3_m_top, h_GenHists4_m_top, h_GenHists1_u_top, h_GenHists2_u_top, h_GenHists3_u_top, h_GenHists4_u_top;
  std::unique_ptr<Hists> h_RecHists0a, h_RecHists0b, h_RecHists0c, h_RecHists0d, h_RecHists1, h_RecHists2, h_RecHists3, h_RecHists4;
  std::unique_ptr<Hists> h_RecHists0a_top, h_RecHists0b_top, h_RecHists0c_top, h_RecHists0d_top, h_RecHists1_top, h_RecHists2_top, h_RecHists3_top, h_RecHists4_top;
  std::unique_ptr<Hists> h_RecGenHists0, h_RecGenHists1, h_RecGenHists2, h_RecGenHists3, h_RecGenHists4, h_RecGenHists5, h_RecGenHists6;
  std::unique_ptr<Hists> h_RecGenHists0_top, h_RecGenHists1_top, h_RecGenHists2_top, h_RecGenHists3_top, h_RecGenHists4_top, h_RecGenHists5_top, h_RecGenHists6_top;
  std::unique_ptr<Hists> h_Elec, h_Muon, h_Elec2, h_Muon2;
  std::unique_ptr<Hists> Cluster01, Cluster02, Cluster03, Cluster04, Cluster05, Cluster06, Cluster07,  Cluster08, Cluster09, Cluster10, Cluster11, Cluster12, Cluster13, Cluster14, Cluster15, Cluster16, Cluster17,  Cluster18, Cluster19, Cluster20;
  int eventcounter = 0;
  int histcounter = 0;
  bool foundallhists = false;
};

MTopJetMCDisplayPostSelectionModule::MTopJetMCDisplayPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  const bool isMC = (ctx.get("dataset_type") == "MC");

  const std::string jet_label_gen("xcone23_gen_fatjets");
  const std::string jet_label_rec("reco_xcone23jets");
  float jet_radius = 0.4;
  float jet_radius_dR_Cut = 2.;
  ////

  //// COMMON MODULES

  const std::string ttbar_gen_label("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  h_jets=ctx.get_handle<std::vector<Jet>>("xcone23_gen_fatjets");

  jetcluster.reset(new GenXCONE23JetProducer(ctx, "xcone23_gen_fatjets", "xcone23_gen_subjets", 0.0, 0.0, 0.0));


  //// OBJ CLEANING
  // const     MuonId muoSR(AndId<Muon>    (PtEtaCut  (45., 2.1), MuonIDMedium()));
  // const ElectronId eleSR(AndId<Electron>(PtEtaSCCut(45., 2.5), ElectronID_Spring15_25ns_tight_noIso));
  // muoSR_cleaner.reset(new     MuonCleaner(muoSR));
  // eleSR_cleaner.reset(new ElectronCleaner(eleSR));
  // cleaner_gen.reset(new GenJetLeptonCleaner(ctx, jet_label_gen, 1));
  // cleaner_topgen.reset(new GenTopJetLeptonCleaner(ctx, 0.8));
  // cleaner_rec.reset(new RecoJetLeptonCleaner(ctx, jet_label_rec, jet_radius));
  // cleaner_toprec.reset(new RecoTopJetLeptonCleaner(ctx, 0.8));

  // JEC
  // std::vector<std::string> JEC_AK4, JEC_AK8;
  // if(isMC){
  //   JEC_AK4 = JERFiles::Fall15_25ns_L123_AK4PFchs_MC;
  //   JEC_AK8 = JERFiles::Fall15_25ns_L123_AK8PFchs_MC;
  // }
  // else {
  //   JEC_AK4 = JERFiles::Fall15_25ns_L123_AK4PFchs_DATA;
  //   JEC_AK8 = JERFiles::Fall15_25ns_L123_AK8PFchs_DATA;
  // }

  // topjet_corrector.reset(new TopJetCorrector(ctx, JEC_AK8));
  // topjet_subjet_corrector.reset(new SubJetCorrector(ctx, JEC_AK4));
  // if(isMC){
  //   ctx.declare_event_input<std::vector<Particle> >(ctx.get("GenTopJetCollection"), "gentopjets");
  //   topjetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "gentopjets", false));
  // }
  topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.4))));
  // topjetlepton_cleaner.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK8, "topjets"));

  //// set up Hists classes:
  // GEN
  // h_GenHists0a.reset(new GenHists(ctx, "00a_GenHists_before_Cleaner", jet_label_gen));
  // h_GenHists0b.reset(new GenHists(ctx, "00b_GenHists_after_Cleaner", jet_label_gen));
  // h_GenHists1.reset(new GenHists(ctx, "01_GenHists_JetN", jet_label_gen));
  // h_GenHists1_m.reset(new GenHists(ctx, "01_GenHists_JetN_matched", jet_label_gen));
  // h_GenHists1_u.reset(new GenHists(ctx, "01_GenHists_JetN_unmatched", jet_label_gen));
  // h_GenHists2.reset(new GenHists(ctx, "02_GenHists_pT_Cut", jet_label_gen));
  // h_GenHists2_m.reset(new GenHists(ctx, "02_GenHists_pT_Cut_matched", jet_label_gen));
  // h_GenHists2_u.reset(new GenHists(ctx, "02_GenHists_pT_Cut_unmatched", jet_label_gen));
  // h_GenHists3.reset(new GenHists(ctx, "03_GenHists_deltaR", jet_label_gen));
  // h_GenHists3_m.reset(new GenHists(ctx, "03_GenHists_deltaR_matched", jet_label_gen));
  // h_GenHists3_u.reset(new GenHists(ctx, "03_GenHists_deltaR_unmatched", jet_label_gen));
  // h_GenHists4.reset(new GenHists(ctx, "04_GenHists_masscut", jet_label_gen));
  // h_GenHists4_m.reset(new GenHists(ctx, "04_GenHists_masscut_matched", jet_label_gen));
  // h_GenHists4_m_XCone.reset(new GenHists(ctx, "04_GenHists_masscut_matched_XCone", jet_label_gen));
  // h_GenHists4_u.reset(new GenHists(ctx, "04_GenHists_masscut_unmatched", jet_label_gen));
  // h_GenHists4_u_XCone.reset(new GenHists(ctx, "04_GenHists_masscut_unmatched_XCone", jet_label_gen));
  // h_GenHists5_LowMass.reset(new GenHists(ctx, "05_GenHists_LowMass", jet_label_gen));
  // h_GenHists5_PeakMass.reset(new GenHists(ctx, "05_GenHists_PeakMass", jet_label_gen));
  // h_GenHists5_HighMass.reset(new GenHists(ctx, "05_GenHists_HighMass", jet_label_gen));

  // Topjet GenHists
  // h_GenHists0a_top.reset(new GenHists_topjet(ctx, "00a_GenHists_top_before_Cleaner"));
  // h_GenHists0b_top.reset(new GenHists_topjet(ctx, "00b_GenHists_top_after_Cleaner"));
  // h_GenHists1_top.reset(new GenHists_topjet(ctx, "01_GenHists_top_JetN"));
  // h_GenHists1_m_top.reset(new GenHists_topjet(ctx, "01_GenHists_top_JetN_matched"));
  // h_GenHists1_u_top.reset(new GenHists_topjet(ctx, "01_GenHists_top_JetN_unmatched"));
  // h_GenHists2_top.reset(new GenHists_topjet(ctx, "02_GenHists_top_pT_Cut"));
  // h_GenHists2_m_top.reset(new GenHists_topjet(ctx, "02_GenHists_top_pT_Cut_matched"));
  // h_GenHists2_u_top.reset(new GenHists_topjet(ctx, "02_GenHists_top_pT_Cut_unmatched"));
  // h_GenHists3_top.reset(new GenHists_topjet(ctx, "03_GenHists_top_deltaR"));
  // h_GenHists3_m_top.reset(new GenHists_topjet(ctx, "03_GenHists_top_deltaR_matched"));
  // h_GenHists3_u_top.reset(new GenHists_topjet(ctx, "03_GenHists_top_deltaR_unmatched"));
  // h_GenHists4_top.reset(new GenHists_topjet(ctx, "04_GenHists_top_masscut"));
  // h_GenHists4_m_top.reset(new GenHists_topjet(ctx, "04_GenHists_top_masscut_matched"));
  // h_GenHists4_u_top.reset(new GenHists_topjet(ctx, "04_GenHists_top_masscut_unmatched"));
  // h_GenHists4_xconeN5.reset(new GenHists_xconeN5(ctx, "04_GenHists_xconeN5", "xconeN5R4_gen"));
  // h_GenHists5_xconeN5_LowMass.reset(new GenHists_xconeN5(ctx, "05_GenHists_xconeSub_LowMass", "xcone23_gen_subjets"));
  // h_GenHists5_xconeN5_PeakMass.reset(new GenHists_xconeN5(ctx, "05_GenHists_xconeSub_PeakMass", "xcone23_gen_subjets"));
  // h_GenHists5_xconeN5_HighMass.reset(new GenHists_xconeN5(ctx, "05_GenHists_xconeSub_HighMass", "xcone23_gen_subjets"));

  //RECO
  // h_RecHists0a.reset(new RecoHists(ctx, "00a_RecHists_before_Cleaner", jet_label_rec));
  // h_RecHists0b.reset(new RecoHists(ctx, "00b_RecHists_after_Cleaner", jet_label_rec));
  // h_RecHists0c.reset(new RecoHists(ctx, "00c_RecHists_LepSel", jet_label_rec));
  // h_RecHists0d.reset(new RecoHists(ctx, "00d_RecHists_background_reduction", jet_label_rec));
  // h_RecHists1.reset(new RecoHists(ctx, "01_RecHists_JetN", jet_label_rec));
  // h_RecHists2.reset(new RecoHists(ctx, "02_RecHists_pT_Cut", jet_label_rec));
  // h_RecHists3.reset(new RecoHists(ctx, "03_RecHists_deltaR", jet_label_rec));
  // h_RecHists4.reset(new RecoHists(ctx, "04_RecHists_masscut", jet_label_rec));

  //Topjet Reco
  h_RecHists0a_top.reset(new RecoHists_topjet(ctx, "00a_RecHists_top_before_Cleaner", "topjets"));
  h_RecHists0b_top.reset(new RecoHists_topjet(ctx, "00b_RecHists_top_after_Cleaner", "topjets"));
  h_RecHists0c_top.reset(new RecoHists_topjet(ctx, "00c_RecHists_top_LepSel", "topjets"));
  h_RecHists0d_top.reset(new RecoHists_topjet(ctx, "00d_RecHists_top_background_reduction", "topjets"));
  h_RecHists1_top.reset(new RecoHists_topjet(ctx, "01_RecHists_top_JetN", "topjets"));
  h_RecHists2_top.reset(new RecoHists_topjet(ctx, "02_RecHists_top_pT_Cut", "topjets"));
  h_RecHists3_top.reset(new RecoHists_topjet(ctx, "03_RecHists_top_deltaR", "topjets"));
  h_RecHists4_top.reset(new RecoHists_topjet(ctx, "04_RecHists_top_masscut", "topjets"));

  // Reco Lepton
  h_Elec.reset(new ElectronHists(ctx, "ElecHist"));
  h_Muon.reset(new MuonHists(ctx, "MuonHist"));
  h_Elec2.reset(new ElectronHists(ctx, "ElecHist_after"));
  h_Muon2.reset(new MuonHists(ctx, "MuonHist_after"));

  //GEN+RECO (resolution)
  // h_RecGenHists0.reset(new RecoGenHists(ctx, "RecGenHists_all", jet_label_rec, jet_label_gen));
  // h_RecGenHists1.reset(new RecoGenHists(ctx, "RecGenHists_000to100", jet_label_rec, jet_label_gen));
  // h_RecGenHists2.reset(new RecoGenHists(ctx, "RecGenHists_100to150", jet_label_rec, jet_label_gen));
  // h_RecGenHists3.reset(new RecoGenHists(ctx, "RecGenHists_150to200", jet_label_rec, jet_label_gen));
  // h_RecGenHists4.reset(new RecoGenHists(ctx, "RecGenHists_200to250", jet_label_rec, jet_label_gen));
  // h_RecGenHists5.reset(new RecoGenHists(ctx, "RecGenHists_250to300", jet_label_rec, jet_label_gen));
  // h_RecGenHists6.reset(new RecoGenHists(ctx, "RecGenHists_300to500", jet_label_rec, jet_label_gen));
  // h_RecGenHists0_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_all"));
  // h_RecGenHists1_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_000to100"));
  // h_RecGenHists2_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_100to150"));
  // h_RecGenHists3_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_150to200"));
  // h_RecGenHists4_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_200to250"));
  // h_RecGenHists5_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_250to300"));
  // h_RecGenHists6_top.reset(new RecoGenHists_topjet(ctx, "RecGenHists_top_300to500"));

  // Clustering Hists
  Cluster01.reset(new ClusteringHists(ctx, "JetDisplay_event01", 0.0));
  Cluster02.reset(new ClusteringHists(ctx, "JetDisplay_event02", 0.0));
  Cluster03.reset(new ClusteringHists(ctx, "JetDisplay_event03", 0.0));
  Cluster04.reset(new ClusteringHists(ctx, "JetDisplay_event04", 0.0));
  Cluster05.reset(new ClusteringHists(ctx, "JetDisplay_event05", 0.0));
  Cluster06.reset(new ClusteringHists(ctx, "JetDisplay_event06", 0.0));
  Cluster07.reset(new ClusteringHists(ctx, "JetDisplay_event07", 0.0));
  Cluster08.reset(new ClusteringHists(ctx, "JetDisplay_event08", 0.4));
  Cluster09.reset(new ClusteringHists(ctx, "JetDisplay_event09", 0.0));
  Cluster10.reset(new ClusteringHists(ctx, "JetDisplay_event10", 0.0));
  Cluster11.reset(new ClusteringHists(ctx, "JetDisplay_event11", 0.0));
  Cluster12.reset(new ClusteringHists(ctx, "JetDisplay_event12", 0.0));
  Cluster13.reset(new ClusteringHists(ctx, "JetDisplay_event13", 0.0));
  Cluster14.reset(new ClusteringHists(ctx, "JetDisplay_event14", 0.0));
  Cluster15.reset(new ClusteringHists(ctx, "JetDisplay_event15", 0.0));
  Cluster16.reset(new ClusteringHists(ctx, "JetDisplay_event16", 0.0));
  Cluster17.reset(new ClusteringHists(ctx, "JetDisplay_event17", 0.0));
  Cluster18.reset(new ClusteringHists(ctx, "JetDisplay_event18", 0.0));
  Cluster19.reset(new ClusteringHists(ctx, "JetDisplay_event19", 0.0));
  Cluster20.reset(new ClusteringHists(ctx, "JetDisplay_event20", 0.0));

  //// EVENT SELECTION
  // GEN
  // n_genjets.reset(new NGenJets(ctx, jet_label_gen, 50, 2, 2));  // ==2 jets with pt > 150
  // topjetpt.reset(new LeadingJetPT(ctx, jet_label_gen, 400)); // leading jet pt > 400
  // deltaR.reset(new DeltaRCut(ctx, jet_label_gen, jet_radius_dR_Cut));
  // deltaPhi.reset(new DeltaPhiCut(ctx, jet_label_gen, 1.0));
  // masscut.reset(new MassCut(ctx, jet_label_gen));
  // matching.reset(new Matching(ctx, jet_label_gen, jet_radius));
  // matching_XCone.reset(new Matching_XCone23(ctx, "gen_xcone33subjets_1"));

  // GEN Top
  // n_genjets_top.reset(new NGenTopJets(ctx, 200, 2, 2));  // ==2 jets with pt > 200 (topjets are produced >~ 180)
  // topjetpt_top.reset(new LeadingTopJetPT(ctx, 400)); // leading jet pt > 400
  // deltaR_top.reset(new DeltaRCut_top(ctx, 0.8));
  // masscut_top.reset(new MassCut_top(ctx));
  // matching_top.reset(new Matching_top(ctx, 0.8));  ////

  // RECO
  // n_recjets.reset(new NRecoJets(ctx, jet_label_rec, 50, 2, 2));  // ==2 jets with pt > 150
  // topjetpt_rec.reset(new LeadingRecoJetPT(ctx, jet_label_rec, 400)); // leading jet pt > 400
  // met_sel.reset(new METCut  (20, uhh2::infinity));
  // htlep_sel.reset(new HTlepCut(100, uhh2::infinity));
  // twodcut_sel.reset(new TwoDCut1(.4, 40.));
  // deltaR_rec.reset(new DeltaRCutReco(ctx, jet_label_rec, jet_radius_dR_Cut));
  // masscut_rec.reset(new MassCutReco(ctx, jet_label_rec));
  ////

  // RECO Top
  // n_recjets_top.reset(new NRecoJets_topjet(ctx, 150, 2, 2));  // ==2 jets with pt > 150
  // topjetpt_rec_top.reset(new LeadingRecoJetPT_topjet(ctx, 400)); // leading jet pt > 400
  // deltaR_rec_top.reset(new DeltaRCutReco_topjet(ctx, 1.5));
  // masscut_rec_top.reset(new MassCutReco_topjet(ctx));
  ////


  // Selection for Mass binning
  // massbin1.reset(new MassCutGen1(ctx, jet_label_gen,   0, 100));
  // massbin2.reset(new MassCutGen1(ctx, jet_label_gen, 100, 150));
  // massbin3.reset(new MassCutGen1(ctx, jet_label_gen, 150, 200));
  // massbin4.reset(new MassCutGen1(ctx, jet_label_gen, 200, 250));
  // massbin5.reset(new MassCutGen1(ctx, jet_label_gen, 250, 300));
  // massbin6.reset(new MassCutGen1(ctx, jet_label_gen, 300, 500));
  // massbin1_top.reset(new MassCutGen1_top(ctx,  50, 100));
  // massbin2_top.reset(new MassCutGen1_top(ctx, 100, 150));
  // massbin3_top.reset(new MassCutGen1_top(ctx, 150, 200));
  // massbin4_top.reset(new MassCutGen1_top(ctx, 200, 250));
  // massbin5_top.reset(new MassCutGen1_top(ctx, 250, 300));
  // massbin6_top.reset(new MassCutGen1_top(ctx, 300, 500));



////

}

bool MTopJetMCDisplayPostSelectionModule::process(uhh2::Event& event){

  // cout<<event.event<<endl;
  //  COMMON MODULES
  // ================ set to true / false to run analysis =============================
  // bool do_gensel = true;
  // bool do_gensel_top = false;
  // bool do_recsel = false;
  // bool do_recsel_top = false;
  // ==================================================================================

  std::vector<int> selectedEvents;
  selectedEvents.push_back(38649894);
  selectedEvents.push_back(38650347);
  selectedEvents.push_back(38657926);
  selectedEvents.push_back(38681906);
  selectedEvents.push_back(38682050);
  selectedEvents.push_back(40780514);
  selectedEvents.push_back(40783315);
  selectedEvents.push_back(40783616);
  selectedEvents.push_back(40788695);
  selectedEvents.push_back(40795633);
  selectedEvents.push_back(41171186);
  selectedEvents.push_back(41177699);
  selectedEvents.push_back(41189702);
  selectedEvents.push_back(42788581);
  selectedEvents.push_back(42788823);
  selectedEvents.push_back(43569434);
  selectedEvents.push_back(51074160);
  selectedEvents.push_back(59483825);
  selectedEvents.push_back(60638255);
  selectedEvents.push_back(62705250);
  selectedEvents.push_back(99190083);

  bool eventisselected = false;
  for(auto evtnr: selectedEvents){
    if(event.event == evtnr) eventisselected = true;
  }
  if(!eventisselected) return false;

  if(foundallhists) return false;

  eventcounter ++;
  cout << eventcounter;

  ttgenprod->process(event);
  jetcluster->process(event);

  vector<Jet> jets = event.get(h_jets);

  if(jets[0].pt() > 400){
    if(jets[0].v4().M() > jets[1].v4().M()){
      cout << "  <-- passed; " << event.event << endl;
      histcounter++;
      if(histcounter ==  1) Cluster01->fill(event);
      if(histcounter ==  2) Cluster02->fill(event);
      if(histcounter ==  3) Cluster03->fill(event);
      if(histcounter ==  4) Cluster04->fill(event);
      if(histcounter ==  5) Cluster05->fill(event);
      if(histcounter ==  6) Cluster06->fill(event);
      if(histcounter ==  7) Cluster07->fill(event);
      if(histcounter ==  8) Cluster08->fill(event);
      if(histcounter ==  9) Cluster09->fill(event);
      if(histcounter == 10) Cluster10->fill(event);
      if(histcounter == 11) Cluster11->fill(event);
      if(histcounter == 12) Cluster12->fill(event);
      if(histcounter == 13) Cluster13->fill(event);
      if(histcounter == 14) Cluster14->fill(event);
      if(histcounter == 15) Cluster15->fill(event);
      if(histcounter == 16) Cluster16->fill(event);
      if(histcounter == 17) Cluster17->fill(event);
      if(histcounter == 18) Cluster18->fill(event);
      if(histcounter == 19) Cluster19->fill(event);
      if(histcounter == 20) Cluster20->fill(event);
      if(histcounter > 20) foundallhists = true;
    }
    else{
      cout << endl;
      return false;
    }
  }
  else{
    cout << endl;
    return false;
  }



  // muoSR_cleaner->process(event);
  // sort_by_pt<Muon>(*event.muons);
  //
  // eleSR_cleaner->process(event);
  // sort_by_pt<Electron>(*event.electrons);
  // cleaner_gen->process(event); //lepton cleaner needed for XCone merge

  // h_GenHists0a->fill(event);
  // h_GenHists0a_top->fill(event);
  // h_RecHists0a->fill(event);
  // h_RecHists0a_top->fill(event);
  // ---------------------------------------------------------------------------------------
  // ---------------------------- apply Gen Selection --------------------------------------
  // ---------------------------------------------------------------------------------------

  // Cleaner
  //cleaner_gen->process(event);
  // cleaner_topgen->process(event);

  // Cleaner for reco objects
  // muoSR_cleaner->process(event);
  // sort_by_pt<Muon>(*event.muons);
  //
  // eleSR_cleaner->process(event);
  // sort_by_pt<Electron>(*event.electrons);
  //
  // cleaner_rec->process(event);

  // topjet_corrector->process(event);
  // topjet_subjet_corrector->process(event);
  // // if(topjetER_smearer.get()) topjetER_smearer->process(event);
  // topjet_cleaner->process(event);
  // topjetlepton_cleaner->process(event); // cleaner by key matching
  // cleaner_toprec->process(event);
  // topjet_cleaner->process(event);
  // sort_by_pt<TopJet>(*event.topjets);

  // h_GenHists0b->fill(event);
  // h_GenHists0b_top->fill(event);
  // h_RecHists0b->fill(event);
  // h_RecHists0b_top->fill(event);

  // =============== ==2 Jets with pT > 150/200 =======================================
  // if(do_gensel && !(n_genjets->passes(event))) return false;
  // if(do_gensel_top && !(n_genjets_top->passes(event))) return false;
  // h_GenHists1->fill(event);
  // h_GenHists1_top->fill(event);
  // // apply matching
  // if(do_gensel  && matching_XCone->passes(event)){
  //   h_GenHists1_m->fill(event);
  //   cout << event.event << endl;
  //   // if() Cluster01->fill(event);
  // }
  // if(do_gensel_top  && matching_top->passes(event)){
  //   h_GenHists1_m_top->fill(event);
  // }
  // if(do_gensel  && !(matching_XCone->passes(event))){
  //   h_GenHists1_u->fill(event);
  // }
  // if(do_gensel_top  && !(matching_top->passes(event))){
  //   h_GenHists1_u_top->fill(event);
  // }
  // ===================================================================================

  // ================ pT cut on leading Jet ============================================
  // if(do_gensel  && !(topjetpt->passes(event))) return false;
  // if(do_gensel_top  && !(topjetpt_top->passes(event))) return false;
  // h_GenHists2->fill(event);
  // h_GenHists2_top->fill(event);
  // // apply matching
  // if(do_gensel  && matching_XCone->passes(event)){
  //   h_GenHists2_m->fill(event);
  // }
  // if(do_gensel_top  && matching_top->passes(event)){
  //   h_GenHists2_m_top->fill(event);
  // }
  // if(do_gensel  && !(matching_XCone->passes(event))){
  //   h_GenHists2_u->fill(event);
  // }
  // if(do_gensel_top  && !(matching_top->passes(event))){
  //   h_GenHists2_u_top->fill(event);
  // }
  // ===================================================================================


  // ================ deltaR(lepton, 2nd jet) < jet radius =============================
  // if(do_gensel  && !(deltaR->passes(event))) return false;
  // if(do_gensel  && !(deltaPhi->passes(event))) return false;
  // if(do_gensel_top  && !(deltaR_top->passes(event))) return false;
  // h_GenHists3->fill(event);
  // h_GenHists3_top->fill(event);
  // // apply matching
  // if(do_gensel  && matching_XCone->passes(event)){
  //   h_GenHists3_m->fill(event);
  // }
  // if(do_gensel_top  && matching_top->passes(event)){
  //   h_GenHists3_m_top->fill(event);
  // }
  // if(do_gensel  && !(matching_XCone->passes(event))){
  //   h_GenHists3_u->fill(event);
  // }
  // if(do_gensel_top  && !(matching_top->passes(event))){
  //   h_GenHists3_u_top->fill(event);
  // }

  // cout<<"processing event Number: "<<event.event<<endl;

  // ===================================================================================


  // ================ m(1st jet) > m(2nd jet + lepton)==================================
  // if(do_gensel  && !(masscut->passes(event))) return false;
  // if(do_gensel_top  && !(masscut_top->passes(event))) return false;
  // h_GenHists4->fill(event);
  // h_GenHists4_top->fill(event);
  // h_GenHists4_xconeN5->fill(event);

  // apply matching
  // if(do_gensel  && matching_XCone->passes(event)){
  //   h_GenHists4_m->fill(event);
  // }
  // if(do_gensel_top  && matching_top->passes(event)){
  //   h_GenHists4_m_top->fill(event);
  // }
  // if(do_gensel  && !(matching_XCone->passes(event))){
  //   h_GenHists4_u->fill(event);
  // }
  // if(do_gensel_top  && !(matching_top->passes(event))){
  //   h_GenHists4_u_top->fill(event);
  // }

  // ===================================================================================

  //
  // h_Elec->fill(event);
  // h_Muon->fill(event);

  // ================ mass bins for resolution studies =================================

  // if(do_gensel){
  //   h_RecGenHists0->fill(event);
  //   if(massbin1->passes(event)) h_RecGenHists1->fill(event);
  //   if(massbin2->passes(event)) h_RecGenHists2->fill(event);
  //   if(massbin3->passes(event)) h_RecGenHists3->fill(event);
  //   if(massbin4->passes(event)) h_RecGenHists4->fill(event);
  //   if(massbin5->passes(event)) h_RecGenHists5->fill(event);
  //   if(massbin6->passes(event)) h_RecGenHists6->fill(event);
  // }
  // if(do_gensel_top){
  //   h_RecGenHists0_top->fill(event);
  //   if(massbin1_top->passes(event)) h_RecGenHists1_top->fill(event);
  //   if(massbin2_top->passes(event)) h_RecGenHists2_top->fill(event);
  //   if(massbin3_top->passes(event)) h_RecGenHists3_top->fill(event);
  //   if(massbin4_top->passes(event)) h_RecGenHists4_top->fill(event);
  //   if(massbin5_top->passes(event)) h_RecGenHists5_top->fill(event);
  //   if(massbin6_top->passes(event)) h_RecGenHists6_top->fill(event);
  // }
  // ===================================================================================

  // ---------------------------------------------------------------------------------------
  // ---------------------------- apply Reco Selection -------------------------------------
  // ---------------------------------------------------------------------------------------

  //select one elec/muon
  // if(do_recsel){
  //   const bool pass_lep1 = (((event.muons->size() == 1)&&(event.electrons->size() == 0)) ||((event.muons->size() == 0)&&(event.electrons->size() == 1)));
  //   if(!pass_lep1) return false;
  // }
  //
  //
  // h_RecHists0c->fill(event);
  // h_RecHists0c_top->fill(event);
  // h_Elec2->fill(event);
  // h_Muon2->fill(event);

  // =====================================================================================================
  // =====================================================================================================
  // =============== FIX: Reco Selection should be based on AK4 jets! ====================================
  // =====================================================================================================
  // =====================================================================================================

  // reduction of background (MET, 2D, Ht_lep)
  // if(do_recsel && !(met_sel->passes(event))) return false;
  // if(do_recsel && !(htlep_sel->passes(event))) return false;




  /* lepton-2Dcut variables */
  // if(do_recsel){
  //   const bool pass_twodcut = twodcut_sel->passes(event);  {

  //     for(auto& muo : *event.muons){

  // 	float    dRmin, pTrel;
  // 	std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);

  // 	muo.set_tag(Muon::twodcut_dRmin, dRmin);
  // 	muo.set_tag(Muon::twodcut_pTrel, pTrel);
  //     }

  //     for(auto& ele : *event.electrons){

  // 	float    dRmin, pTrel;
  // 	std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

  // 	ele.set_tag(Electron::twodcut_dRmin, dRmin);
  // 	ele.set_tag(Electron::twodcut_pTrel, pTrel);
  //     }
  //   }

  // if(do_recsel && !(pass_twodcut)) return false;
  // }

  // h_RecHists0d->fill(event);
  // h_RecHists0d_top->fill(event);


  // // ==2 Jets with pT > 150
  // if(do_recsel && !(n_recjets->passes(event))) return false;
  // if(do_recsel_top && !(n_recjets_top->passes(event))) return false;
  // h_RecHists1->fill(event);
  // h_RecHists1_top->fill(event);

  // // pT cut on leading Jet
  // if(do_recsel && !(topjetpt_rec->passes(event))) return false;
  // if(do_recsel_top && !(topjetpt_rec_top->passes(event))) return false;
  // h_RecHists2->fill(event);
  // h_RecHists2_top->fill(event);

  // // deltaR(lepton, 2nd jet) < jet radius
  // // if(do_recsel && !(deltaR_rec->passes(event))) return false;
  // if(do_recsel_top && !(deltaR_rec_top->passes(event))) return false;
  // h_RecHists3->fill(event);
  // h_RecHists3_top->fill(event);

  // // m(1st jet) > m(2nd jet + lepton)
  // if(do_recsel && !(masscut_rec->passes(event))) return false;
  // if(do_recsel_top && !(masscut_rec_top->passes(event))) return false;
  // h_RecHists4->fill(event);
  // h_RecHists4_top->fill(event);
  // ---------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------

  return true; //false to delete all collections and keep hists
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetMCDisplayPostSelectionModule)
