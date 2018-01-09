#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/MCWeight.h>

#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/RecoSelections_topjet.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/RecoHists_xcone.h>
#include <UHH2/MTopJet/include/RecoHists_topjet.h>

#include <UHH2/MTopJet/include/GenHists_xcone.h>
#include <UHH2/MTopJet/include/GenHists_particles.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/SubjetHists_xcone.h>
#include "UHH2/MTopJet/include/RecoGenHists_subjets.h"
#include "UHH2/MTopJet/include/RecoGenHists_ak4.h"
#include "UHH2/MTopJet/include/CorrectionHists_subjets.h"
#include "UHH2/MTopJet/include/CorrectionFactor.h"
#include "UHH2/MTopJet/include/tt_width_reweight.h"

#include <vector>

/*
*******************************************************************
**************** TO DO ********************************************
*******************************************************************
- SF apllied to TTbar!
*******************************************************************
*******************************************************************
*/

class MTopJetPostSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners & Correctors
  // std::unique_ptr<uhh2::AnalysisModule> Correction;

  // selections
  std::unique_ptr<uhh2::Selection> njet_had;
  std::unique_ptr<uhh2::Selection> njet_lep;
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> pt350_gensel;
  std::unique_ptr<uhh2::Selection> pt350_sel;
  std::unique_ptr<uhh2::Selection> pt750_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;
  std::unique_ptr<uhh2::Selection> pt_gensel;
  std::unique_ptr<uhh2::Selection> mass_gensel;
  std::unique_ptr<uhh2::Selection> pt_gensel23;
  std::unique_ptr<uhh2::Selection> mass_gensel23;
  std::unique_ptr<uhh2::Selection> matched_sub;
  std::unique_ptr<uhh2::Selection> matched_fat;


  std::unique_ptr<uhh2::Selection> njet_ak8;
  std::unique_ptr<uhh2::Selection> deltaR_ak8;
  std::unique_ptr<uhh2::Selection> pt_sel_ak8;
  std::unique_ptr<uhh2::Selection> mass_ak8;


  // get weight (with all SF and weight applied in previous cycle)
  Event::Handle<double>h_weight;

  // use top pt reweight module instead of derived sf?
  std::unique_ptr<TopPtReweight> ttbar_reweight;

  // handles for output
  Event::Handle<bool>h_gensel_2;
  Event::Handle<bool>h_recsel_2;
  Event::Handle<bool>h_matched;
  Event::Handle<bool>h_measure_rec;
  Event::Handle<bool>h_gensel23;
  Event::Handle<bool>h_measure_gen;
  Event::Handle<bool>h_nomass_rec;
  Event::Handle<bool>h_nomass_gen;
  Event::Handle<bool>h_pt350_gen;
  Event::Handle<bool>h_pt350_rec;
  Event::Handle<bool>h_nobtag_rec;
  Event::Handle<bool>h_ttbar;
  Event::Handle<double>h_ttbar_SF;
  Event::Handle<double>h_mass_gen33;
  Event::Handle<double>h_mass_gen23;
  Event::Handle<double>h_mass_rec;
  Event::Handle<double>h_pt_gen33;
  Event::Handle<double>h_pt_gen23;
  Event::Handle<double>h_pt_rec;
  Event::Handle<double>h_genweight;
  Event::Handle<double>h_recweight;
  Event::Handle<std::vector<Jet>>h_recjets_had;
  Event::Handle<std::vector<Particle>>h_genjets23_had;
  Event::Handle<std::vector<Particle>>h_genjets33_had;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  //width reweight
  std::unique_ptr<AnalysisModule> width_reweight;

  //scale variation
  std::unique_ptr<AnalysisModule> scale_variation;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_Muon_PreSel, h_MTopJet_PreSel, h_Jets_PreSel, h_XCone_cor_PreSel;

  std::unique_ptr<Hists> h_Muon, h_MTopJet;
  std::unique_ptr<Hists> h_CorrectionHists, h_CorrectionHists_after;
  std::unique_ptr<Hists> h_RecGenHists_ak4, h_RecGenHists_ak4_noJEC;

  std::unique_ptr<Hists> h_XCone_cor, h_XCone_jec, h_XCone_raw;
  std::unique_ptr<Hists> h_XCone_cor_SF, h_XCone_jec_SF, h_XCone_raw_SF;
  std::unique_ptr<Hists> h_XCone_cor_pt350, h_XCone_cor_noptcut;
  std::unique_ptr<Hists> h_XCone_cor_subjets, h_XCone_jec_subjets, h_XCone_raw_subjets;
  std::unique_ptr<Hists> h_XCone_cor_subjets_SF, h_XCone_jec_subjets_SF, h_XCone_raw_subjets_SF;
  std::unique_ptr<Hists> h_XCone_cor_Sel_noSel, h_XCone_cor_Sel_noMass, h_XCone_cor_Sel_nobtag, h_XCone_cor_Sel_pt350;

  std::unique_ptr<Hists> h_XCone_cor_m, h_XCone_cor_u, h_XCone_cor_m_fat, h_XCone_cor_u_fat;

  std::unique_ptr<Hists> h_XCone_cor_PUlow, h_XCone_cor_PUmid,h_XCone_cor_PUhigh;
  std::unique_ptr<Hists> h_XCone_cor_PUlow_subjets, h_XCone_cor_PUmid_subjets, h_XCone_cor_PUhigh_subjets;


  std::unique_ptr<Hists> h_RecGenHists_lowPU, h_RecGenHists_medPU, h_RecGenHists_highPU, h_RecGenHists_lowPU_noJEC, h_RecGenHists_medPU_noJEC, h_RecGenHists_highPU_noJEC, h_RecGenHists_RecOnly_corr;
  std::unique_ptr<Hists> h_XCone_GEN_RecOnly, h_XCone_GEN_GenOnly, h_XCone_GEN_Both;
  std::unique_ptr<Hists> h_RecGenHists_GenOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly_noJEC;
  std::unique_ptr<Hists> h_RecGenHists_Both, h_RecGenHists_Both_corr;
  std::unique_ptr<Hists> h_RecGenHists_subjets, h_RecGenHists_subjets_noJEC, h_RecGenHists_subjets_corrected;
  std::unique_ptr<Hists> h_GenParticles_RecOnly, h_GenParticles_GenOnly, h_GenParticles_Both;
  std::unique_ptr<Hists> h_GenParticles_SM, h_GenParticles_newWidth;

  std::unique_ptr<Hists> h_XCone_cor_migration_pt, h_XCone_cor_migration_pt350, h_XCone_cor_migration_mass, h_XCone_cor_migration_btag;
  std::unique_ptr<Hists> h_750_ak8, h_750_xcone;

  std::unique_ptr<Hists> h_XCone_GEN_Sel_measurement, h_XCone_GEN_Sel_noMass, h_XCone_GEN_Sel_pt350;

  bool isMC;    //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part
  int counter;
};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  /*************************** CONFIGURATION **********************************************************************************/ 
  isMC = (ctx.get("dataset_type") == "MC");
  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700" || ctx.get("dataset_version") == "TTbar_Mtt0700to1000"  || ctx.get("dataset_version") == "TTbar_Mtt1000toInft") isTTbar = true;
  else  isTTbar = false;

  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  if(isTTbar) ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {
    std::string log("MTopJetPostSelectionModule::MTopJetPostSelectionModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";
    throw std::runtime_error(log);
  }

  counter = 0;

  // Top width reweight
  width_reweight.reset(new tt_width_reweight(ctx, 8.0));

  // Top PT reweight
  ttbar_reweight.reset(new TopPtReweight(ctx,0.0615,-0.0005,"","weight_ttbar",true)); // 13 TeV

  //scale variation
  scale_variation.reset(new MCScaleVariation(ctx));


  // get handle for weight
  h_weight=ctx.get_handle<double>("weight");

  // write output
  h_gensel_2 = ctx.get_handle<bool>("passed_gensel_2");
  h_recsel_2 = ctx.get_handle<bool>("passed_recsel_2");
  h_matched = ctx.declare_event_output<bool>("matched");
  h_measure_rec = ctx.declare_event_output<bool>("passed_measurement_rec");
  h_measure_gen = ctx.declare_event_output<bool>("passed_measurement_gen");
  h_nomass_rec = ctx.declare_event_output<bool>("passed_massmigration_rec");
  h_nomass_gen = ctx.declare_event_output<bool>("passed_massmigration_gen");
  h_pt350_rec = ctx.declare_event_output<bool>("passed_ptmigration_rec");
  h_pt350_gen = ctx.declare_event_output<bool>("passed_ptmigration_gen");
  h_nobtag_rec = ctx.declare_event_output<bool>("passed_btagmigration_rec");
  h_gensel23 = ctx.declare_event_output<bool>("passed_gensel23_full");
  h_ttbar = ctx.declare_event_output<bool>("is_TTbar");
  h_ttbar_SF = ctx.declare_event_output<double>("TTbar_SF");
  h_mass_gen23 = ctx.declare_event_output<double>("Mass_Gen23");
  h_mass_gen33 = ctx.declare_event_output<double>("Mass_Gen33");
  h_mass_rec = ctx.declare_event_output<double>("Mass_Rec");
  h_pt_gen23 = ctx.declare_event_output<double>("Pt_Gen23");
  h_pt_gen33 = ctx.declare_event_output<double>("Pt_Gen33");
  h_pt_rec = ctx.declare_event_output<double>("Pt_Rec");
  h_genweight = ctx.declare_event_output<double>("gen_weight");
  h_recweight = ctx.declare_event_output<double>("rec_weight");

  h_recjets_had = ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined_Corrected");
  if(isMC) h_genjets23_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone23_had_Combined");
  if(isMC) h_genjets33_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone33_had_Combined");

  /*************************** Setup Subjet Corrector **********************************************************************************/ 
  // Correction.reset(new CorrectionFactor(ctx, "XConeTopJets_Corrected"));
  /*************************** Setup Selections **********************************************************************************/ 
  // RECO Selection
  // chose: XCone33_had_Combined_Corrected,  XCone33_had_Combined_noJEC,  XCone33_had_Combined
  const std::string& jet_label_had = "XCone33_had_Combined_Corrected";
  const std::string& jet_label_lep = "XCone33_lep_Combined_Corrected";

  njet_had.reset(new NJetXCone(ctx, jet_label_had, 1));
  njet_lep.reset(new NJetXCone(ctx, jet_label_lep, 1));
 
  pt_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 400));
  pt350_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 350));
  mass_sel.reset(new MassCutXCone(ctx, jet_label_had, jet_label_lep));

  // GEN Selection
  pt_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 400));
  pt350_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 350));
  mass_gensel.reset(new MassCut_gen(ctx, "GEN_XCone33_had_Combined", "GEN_XCone33_lep_Combined"));
  pt_gensel23.reset(new LeadingJetPT_gen(ctx, "GEN_XCone23_had_Combined", 400));
  mass_gensel23.reset(new MassCut_gen(ctx, "GEN_XCone23_had_Combined", "GEN_XCone23_lep_Combined"));

  // Selection for matching reco jets to gen particles
  if(isTTbar) matched_sub.reset(new Matching_XCone33(ctx, true));
  if(isTTbar) matched_fat.reset(new Matching_XCone33(ctx, false));

  // AK8 Selection / 750 GeV
  njet_ak8.reset(new NRecoJets_topjet(ctx, 200, 2, 2));
  deltaR_ak8.reset(new DeltaRCutReco_topjet(ctx, 1.2));
  pt_sel_ak8.reset(new LeadingRecoJetPT_topjet(ctx, 750));
  mass_ak8.reset(new MassCutReco_topjet(ctx));
  pt750_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 750));

  /*************************** Set up Hists classes **********************************************************************************/

  //750GeV hists
  h_750_ak8.reset(new RecoHists_topjet(ctx, "750_ak8", "topjets"));
  h_750_xcone.reset(new RecoHists_xcone(ctx, "750_xcone", "cor"));

  // XCone Combined Jet
  h_XCone_raw.reset(new RecoHists_xcone(ctx, "XCone_raw", "raw"));
  h_XCone_cor.reset(new RecoHists_xcone(ctx, "XCone_cor", "cor"));
  h_XCone_jec.reset(new RecoHists_xcone(ctx, "XCone_jec", "jec"));
  h_XCone_raw_SF.reset(new RecoHists_xcone(ctx, "XCone_raw_SF", "raw"));
  h_XCone_cor_SF.reset(new RecoHists_xcone(ctx, "XCone_cor_SF", "cor"));
  h_XCone_jec_SF.reset(new RecoHists_xcone(ctx, "XCone_jec_SF", "jec"));

  h_XCone_cor_pt350.reset(new RecoHists_xcone(ctx, "XCone_cor_pt350", "cor"));
  h_XCone_cor_noptcut.reset(new RecoHists_xcone(ctx, "XCone_cor_noptcut", "cor"));

  // XCone Subjets
  h_XCone_jec_subjets.reset(new SubjetHists_xcone(ctx, "XCone_jec_subjets", "jec"));
  h_XCone_raw_subjets.reset(new SubjetHists_xcone(ctx, "XCone_raw_subjets", "raw"));
  h_XCone_cor_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_subjets", "cor"));
  h_XCone_jec_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_jec_subjets_SF", "jec"));
  h_XCone_raw_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_raw_subjets_SF", "raw"));
  h_XCone_cor_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_cor_subjets_SF", "cor"));

  // migrations from not passing rec cuts
  h_XCone_cor_migration_pt.reset(new RecoHists_xcone(ctx, "XCone_cor_migration_pt", "cor"));
  h_XCone_cor_migration_pt350.reset(new RecoHists_xcone(ctx, "XCone_cor_migration_pt350", "cor"));
  h_XCone_cor_migration_mass.reset(new RecoHists_xcone(ctx, "XCone_cor_migration_mass", "cor"));
  h_XCone_cor_migration_btag.reset(new RecoHists_xcone(ctx, "XCone_cor_migration_btag", "cor"));

  // Different REC Selection applied
  h_XCone_cor_Sel_noSel.reset(new RecoHists_xcone(ctx, "XCone_cor_Sel_noSel", "cor"));
  h_XCone_cor_Sel_noMass.reset(new RecoHists_xcone(ctx, "XCone_cor_Sel_noMass", "cor"));
  h_XCone_cor_Sel_pt350.reset(new RecoHists_xcone(ctx, "XCone_cor_Sel_pt350", "cor"));
  h_XCone_cor_Sel_nobtag.reset(new RecoHists_xcone(ctx, "XCone_cor_Sel_nobtag", "cor"));

  // PU dependence
  h_XCone_cor_PUlow.reset(new RecoHists_xcone(ctx, "XCone_cor_PUlow", "cor"));
  h_XCone_cor_PUmid.reset(new RecoHists_xcone(ctx, "XCone_cor_PUmid", "cor"));
  h_XCone_cor_PUhigh.reset(new RecoHists_xcone(ctx, "XCone_cor_PUhigh", "cor"));

  h_XCone_cor_PUlow_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUlow_subjets", "cor"));
  h_XCone_cor_PUmid_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUmid_subjets", "cor"));
  h_XCone_cor_PUhigh_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUhigh_subjets", "cor"));

  // Matching hists
  if(isTTbar) h_XCone_cor_m.reset(new RecoHists_xcone(ctx, "XCone_cor_matched", "cor"));
  if(isTTbar) h_XCone_cor_u.reset(new RecoHists_xcone(ctx, "XCone_cor_unmatched", "cor"));
  if(isTTbar) h_XCone_cor_m_fat.reset(new RecoHists_xcone(ctx, "XCone_cor_matched_fat", "cor"));
  if(isTTbar) h_XCone_cor_u_fat.reset(new RecoHists_xcone(ctx, "XCone_cor_unmatched_fat", "cor"));

  h_MTopJet.reset(new MTopJetHists(ctx, "EventHists"));
  h_Muon.reset(new MuonHists(ctx, "MuonHists"));

  h_MTopJet_PreSel.reset(new MTopJetHists(ctx, "EventHists_PreSel"));
  h_Muon_PreSel.reset(new MuonHists(ctx, "MuonHits_PreSel"));
  h_Jets_PreSel.reset(new JetHists(ctx, "JetHits_PreSel"));
  h_XCone_cor_PreSel.reset(new RecoHists_xcone(ctx, "XCone_cor_PreSel", "cor"));

  if(isMC){
    if(isTTbar) h_CorrectionHists.reset(new CorrectionHists_subjets(ctx, "CorrectionHists", "jec"));
    if(isTTbar) h_CorrectionHists_after.reset(new CorrectionHists_subjets(ctx, "CorrectionHists_after", "cor"));

    h_XCone_GEN_Sel_measurement.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_measurement"));
    h_XCone_GEN_Sel_noMass.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_noMass"));
    h_XCone_GEN_Sel_pt350.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_pt350"));

    h_XCone_GEN_GenOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly"));
    h_XCone_GEN_RecOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_RecOnly"));
    h_XCone_GEN_Both.reset(new GenHists_xcone(ctx, "XCone_GEN_Both"));

    h_GenParticles_GenOnly.reset(new GenHists_particles(ctx, "GenParticles_GenOnly"));
    h_GenParticles_RecOnly.reset(new GenHists_particles(ctx, "GenParticles_RecOnly"));
    h_GenParticles_Both.reset(new GenHists_particles(ctx, "GenParticles_Both"));

    h_GenParticles_SM.reset(new GenHists_particles(ctx, "GenParticles_SM"));
    h_GenParticles_newWidth.reset(new GenHists_particles(ctx, "GenParticles_newWidth"));

    h_RecGenHists_subjets.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets", "jec"));
    h_RecGenHists_subjets_noJEC.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC", "raw"));
    h_RecGenHists_subjets_corrected.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_corrected", "cor"));

    h_RecGenHists_ak4.reset(new RecoGenHists_ak4(ctx, "RecGenHists_ak4", true));
    h_RecGenHists_ak4_noJEC.reset(new RecoGenHists_ak4(ctx, "RecGenHists_ak4_noJEC", false));

    h_RecGenHists_lowPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_lowPU", "jec"));
    h_RecGenHists_medPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_medPU", "jec"));
    h_RecGenHists_highPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_highPU", "jec"));
    h_RecGenHists_lowPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_lowPU_noJEC", "raw"));
    h_RecGenHists_medPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_medPU_noJEC", "raw"));
    h_RecGenHists_highPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_highPU_noJEC", "raw"));

    // "true" for RecGenHists means to use jets with JEC applied
    h_RecGenHists_GenOnly.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly", "jec"));
    h_RecGenHists_RecOnly.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly", "jec"));
    h_RecGenHists_RecOnly_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_noJEC", "raw"));
    h_RecGenHists_RecOnly_corr.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_corrected", "cor"));
    h_RecGenHists_Both.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both", "jec"));
    h_RecGenHists_Both_corr.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_corrected", "cor"));

  }
  /*********************************************************************************************************************************/ 

}

bool MTopJetPostSelectionModule::process(uhh2::Event& event){

  // get bools for selections from root files
  bool passed_recsel;
  passed_recsel = event.get(h_recsel_2);
  bool passed_gensel33;
  passed_gensel33 = event.get(h_gensel_2);
  bool passed_gensel23;
  passed_gensel23 = event.get(h_gensel_2);
  ////


  // check if event has one had and one lep jet
  if( !(njet_had->passes(event)) ) return false;
  if( !(njet_lep->passes(event)) ) return false;
  ////

  // fill ttbargen class
  if(isTTbar) ttgenprod->process(event);
  ////

  /***************************  get jets to write mass *****************************************************************************************************/ 

  std::vector<Jet> rec_hadjets = event.get(h_recjets_had);
  double mass_rec = rec_hadjets.at(0).v4().M();
  double pt_rec = rec_hadjets.at(0).v4().Pt();
  event.set(h_mass_rec, mass_rec);
  event.set(h_pt_rec, pt_rec);

  if(isMC){
    std::vector<Particle> gen_hadjets23 = event.get(h_genjets23_had);
    std::vector<Particle> gen_hadjets33 = event.get(h_genjets33_had);
    double mass_gen23 = gen_hadjets23.at(0).v4().M();
    double mass_gen33 = gen_hadjets33.at(0).v4().M();
    double pt_gen23 = gen_hadjets23.at(0).v4().Pt();
    double pt_gen33 = gen_hadjets33.at(0).v4().Pt();
    event.set(h_mass_gen23, mass_gen23);
    event.set(h_mass_gen33, mass_gen33);
    event.set(h_pt_gen23, pt_gen23);
    event.set(h_pt_gen33, pt_gen33);
  }
  else{
    event.set(h_mass_gen23, 0.); // set gen mass to 0 for data
    event.set(h_mass_gen33, 0.); // set gen mass to 0 for data
    event.set(h_pt_gen23, 0.);   // set gen pt to 0 for data
    event.set(h_pt_gen33, 0.);   // set gen pt to 0 for data
  }

  /***************************  apply weight *****************************************************************************************************/

  bool do_scale = false;             // change scales mu_R and mu_F?
  bool reweight_ttbar = false;       // apply ttbar reweight?
  bool scale_ttbar = true;           // match MC and data cross-section (for plots only)?
  double SF_tt = 0.75;

  // get lumi weight = genweight (inkl scale variation)
  if(do_scale) scale_variation->process(event);
  double gen_weight = event.weight;

  // now get full weight from prev. Selection (weight = gen_weight * rec_weight)
  event.weight = event.get(h_weight);

  if(isTTbar)h_GenParticles_SM->fill(event);

  // choose if tt bar sample width should be reweighted
  bool do_width_reweight = false;
  if(do_width_reweight && isTTbar){
    width_reweight->process(event);
    gen_weight = event.weight; // store again as gen weight
    h_GenParticles_newWidth->fill(event);
  }

  if(do_scale) scale_variation->process(event);

  // rec weight is now:
  if(reweight_ttbar) ttbar_reweight->process(event);
  double rec_weight;
  if(gen_weight==0)rec_weight = 0;
  else rec_weight = (event.weight)/gen_weight;

  /** b-tagging *********************/
  int jetbtagN(0);
  bool passed_btag;
  for(const auto& j : *event.jets) if(CSVBTag(CSVBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN >= 1) passed_btag = true;
  else passed_btag = false;
  /**********************************/

  bool passed_presel_rec = (passed_recsel && passed_btag);

  /*************************** Events have to pass topjet pt > 400 & Mass_jet1 > Mass_jet2 *******************************************************/

  bool pass_measurement_rec;
  bool pass_pt350migration_rec;
  bool pass_ptmigration_rec;
  bool pass_massmigration_rec;
  bool pass_btagmigration_rec;

  if(passed_recsel && pt_sel->passes(event) && mass_sel->passes(event) && passed_btag) pass_measurement_rec = true;
  else pass_measurement_rec = false;

  if(passed_recsel && !pt_sel->passes(event) && pt350_sel->passes(event) && mass_sel->passes(event) && passed_btag ) pass_pt350migration_rec = true;
  else pass_pt350migration_rec = false;

  if(passed_recsel && !pt_sel->passes(event) && mass_sel->passes(event) && passed_btag ) pass_ptmigration_rec = true;
  else pass_ptmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && !mass_sel->passes(event) && passed_btag) pass_massmigration_rec = true;
  else pass_massmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && mass_sel->passes(event) && !passed_btag) pass_btagmigration_rec = true;
  else pass_btagmigration_rec = false;

  /*************************** Selection again on generator level (data events will not pass gen sel but will be stored if they pass rec sel)  ***/ 
  bool pass_measurement_gen;
  bool pass_pt350migration_gen;
  bool pass_massmigration_gen;

  if(isMC){
    if(passed_gensel33 && pt_gensel->passes(event) && mass_gensel->passes(event)) pass_measurement_gen = true;
    else pass_measurement_gen = false;

    if(passed_gensel33 && !pt_gensel->passes(event) && pt350_gensel->passes(event) && mass_gensel->passes(event)) pass_pt350migration_gen = true;
    else pass_pt350migration_gen = false;

    if(passed_gensel33 && pt_gensel->passes(event) && !mass_gensel->passes(event)) pass_massmigration_gen = true;
    else pass_massmigration_gen = false;

  }

  /***************************  750GeV SELECTION ***********************************************************************************/
  bool passes_750_ak8 = false;
  bool passes_750_xcone = false;

  if(passed_presel_rec && njet_ak8->passes(event)){
    if(deltaR_ak8->passes(event) && pt_sel_ak8->passes(event) && mass_ak8->passes(event)) passes_750_ak8 = true;
  }
  if(pass_measurement_rec && pt750_sel->passes(event)) passes_750_xcone = true;
  /******************************************************************************************************************************/
 
  /*************************** Pile Up bools  ***************************************************************************************************/
  bool lowPU = (event.pvs->size() <= 10);
  bool midPU = (event.pvs->size() > 10 && event.pvs->size() <= 20);
  bool highPU = (event.pvs->size() > 20);

  /*************************** fill hists with reco sel applied ***********************************************************************************/ 

  bool is_matched_sub = false;
  bool is_matched_fat = false;


  // hists after PreSel on REC Level
  if(passed_presel_rec){
    h_MTopJet_PreSel->fill(event);
    h_Muon_PreSel->fill(event);
    h_Jets_PreSel->fill(event);
    h_XCone_cor_PreSel->fill(event);
  }

  // Hists for 750GeV phase space
  if(passes_750_ak8) h_750_ak8->fill(event);
  if(passes_750_xcone) h_750_xcone->fill(event);

  // hists to see events that are generated in measurement phase-space, but reconstructed outside
  if(pass_measurement_gen){
    if(pass_pt350migration_rec) h_XCone_cor_migration_pt350->fill(event);
    if(pass_ptmigration_rec)    h_XCone_cor_migration_pt->fill(event);
    if(pass_massmigration_rec)  h_XCone_cor_migration_mass->fill(event);
    if(pass_btagmigration_rec)  h_XCone_cor_migration_btag->fill(event);
  }


  // see all events reconstructed outside measurement phase-space
  if(pass_ptmigration_rec) h_XCone_cor_noptcut->fill(event);
  if(pass_pt350migration_rec) h_XCone_cor_pt350->fill(event);
  if(pass_btagmigration_rec)  h_XCone_cor_Sel_nobtag->fill(event);
  if(pass_pt350migration_rec) h_XCone_cor_Sel_pt350->fill(event);
  if(pass_massmigration_rec)  h_XCone_cor_Sel_noMass->fill(event);

  if(pass_measurement_gen)    h_XCone_GEN_Sel_measurement->fill(event);
  if(pass_pt350migration_gen) h_XCone_GEN_Sel_pt350->fill(event);
  if(pass_massmigration_gen)  h_XCone_GEN_Sel_noMass->fill(event);

  if(pass_measurement_rec){
    h_XCone_raw->fill(event);
    h_XCone_jec->fill(event);
    h_XCone_cor->fill(event);

    h_XCone_raw_subjets->fill(event);
    h_XCone_jec_subjets->fill(event);
    h_XCone_cor_subjets->fill(event);

    if(isTTbar && scale_ttbar) event.weight *= SF_tt;

    h_XCone_raw_SF->fill(event);
    h_XCone_jec_SF->fill(event);
    h_XCone_cor_SF->fill(event);

    h_XCone_raw_subjets_SF->fill(event);
    h_XCone_jec_subjets_SF->fill(event);
    h_XCone_cor_subjets_SF->fill(event);

    h_MTopJet->fill(event);
    h_Muon->fill(event);

    if(lowPU){
      h_XCone_cor_PUlow_subjets->fill(event);
      h_XCone_cor_PUlow->fill(event);
      if(isMC)h_RecGenHists_lowPU->fill(event);
      if(isMC)h_RecGenHists_lowPU_noJEC->fill(event);
    }
    if(midPU){
      h_XCone_cor_PUmid_subjets->fill(event);
      h_XCone_cor_PUmid->fill(event);
      if(isMC)h_RecGenHists_medPU->fill(event);
      if(isMC)h_RecGenHists_medPU_noJEC->fill(event);
    }
    if(highPU){
      h_XCone_cor_PUhigh_subjets->fill(event);
      h_XCone_cor_PUhigh->fill(event);
      if(isMC)h_RecGenHists_highPU->fill(event);
      if(isMC)h_RecGenHists_highPU_noJEC->fill(event);
    }
    
    if(isTTbar){
      is_matched_sub = matched_sub->passes(event);
      is_matched_fat = matched_fat->passes(event);
      if(is_matched_sub) h_XCone_cor_m->fill(event);
      else h_XCone_cor_u->fill(event);
      if(is_matched_fat) h_XCone_cor_m_fat->fill(event);
      else h_XCone_cor_u_fat->fill(event);
    }

    if(isMC){
      h_XCone_GEN_RecOnly->fill(event);
      if(isTTbar) h_GenParticles_RecOnly->fill(event);
      h_RecGenHists_RecOnly->fill(event);
      h_RecGenHists_RecOnly_noJEC->fill(event);
      h_RecGenHists_RecOnly_corr->fill(event);
      h_RecGenHists_subjets->fill(event);
      h_RecGenHists_subjets_noJEC->fill(event);
      h_RecGenHists_subjets_corrected->fill(event);
      if(isTTbar) h_CorrectionHists->fill(event);
      if(isTTbar) h_CorrectionHists_after->fill(event);
      h_RecGenHists_ak4->fill(event);
      h_RecGenHists_ak4_noJEC->fill(event);
    }
  }

  /*************************** fill hists with gen sel applied *************************************************************************************/ 
  if(pass_measurement_gen){
    h_XCone_GEN_GenOnly->fill(event);
    if(isTTbar) h_GenParticles_GenOnly->fill(event);
    h_RecGenHists_GenOnly->fill(event);
  } 

  /*************************** fill hists with reco+gen selection applied **************************************************************************/ 
  if(pass_measurement_rec && pass_measurement_gen){
    h_XCone_GEN_Both->fill(event);
    if(isTTbar) h_GenParticles_Both->fill(event);
    h_RecGenHists_Both->fill(event);
    h_RecGenHists_Both_corr->fill(event);
  } 

  /*************************** write bools for passing selections **********************************************************************************/ 

  event.set(h_ttbar, isTTbar);
  event.set(h_matched, is_matched_sub);
  event.set(h_ttbar_SF, SF_tt);
  event.set(h_genweight, gen_weight);
  event.set(h_recweight, rec_weight);
  event.set(h_gensel23, passed_gensel23);

  event.set(h_measure_rec, pass_measurement_rec);
  event.set(h_measure_gen, pass_measurement_gen);
  event.set(h_pt350_rec, pass_pt350migration_rec);
  event.set(h_pt350_gen, pass_pt350migration_gen);
  event.set(h_nomass_rec, pass_massmigration_rec);
  event.set(h_nomass_gen, pass_massmigration_gen);
  event.set(h_nobtag_rec, pass_btagmigration_rec);

  /*************************** only store events that survive one of the selections (use looser pt cut) ****************************************************************/
  bool in_migrationmatrix = (pass_measurement_rec || pass_measurement_gen || pass_pt350migration_rec || pass_pt350migration_gen || pass_massmigration_rec || pass_massmigration_gen || pass_btagmigration_rec);

  if(!in_migrationmatrix) return false;
  else return true;
  
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
