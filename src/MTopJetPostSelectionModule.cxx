#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/JetCorrectionSets.h>
#include "UHH2/common/include/YearRunSwitchers.h"

#include <UHH2/MTopJet/include/CorrectionHists_subjets.h>
#include <UHH2/MTopJet/include/ControlHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/GenHists_xcone.h>
#include <UHH2/MTopJet/include/GenHists_particles.h>
#include <UHH2/MTopJet/include/GenHists_process.h>
#include <UHH2/MTopJet/include/LeptonicTop_Hists.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/PDFHists.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone.h>
#include <UHH2/MTopJet/include/RecoGenHists_subjets.h>
#include <UHH2/MTopJet/include/RecoGenHists_ak4.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone_topjet.h>
#include <UHH2/MTopJet/include/RecoHists_xcone.h>
#include <UHH2/MTopJet/include/RecoHists_puppi.h>
#include <UHH2/MTopJet/include/SubjetHists_xcone.h>
#include <UHH2/MTopJet/include/StoreBJet.h>

#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/CorrectionFactor.h>
#include <UHH2/MTopJet/include/ElecTriggerSF.h>
#include <UHH2/MTopJet/include/CorrectionFactor_JMS.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include "UHH2/MTopJet/include/PartonShowerWeight.h"
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/tt_width_reweight.h>

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"


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
  JetId DeepjetMedium, DeepjetLoose, DeepjetTight;
  // cleaners & Correctors
  // std::unique_ptr<uhh2::AnalysisModule> Correction;

  // selections
  std::unique_ptr<uhh2::Selection> njet_had;
  std::unique_ptr<uhh2::Selection> njet_lep;
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> pt450_sel;
  std::unique_ptr<uhh2::Selection> pt530_sel;
  std::unique_ptr<uhh2::Selection> pt2_sel;
  std::unique_ptr<uhh2::Selection> eta_sel;
  std::unique_ptr<uhh2::Selection> pt350_gensel;
  std::unique_ptr<uhh2::Selection> pt350_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;
  std::unique_ptr<uhh2::Selection> pt_gensel;
  std::unique_ptr<uhh2::Selection> pt2_gensel;
  std::unique_ptr<uhh2::Selection> mass_gensel;
  std::unique_ptr<uhh2::Selection> matched_sub_GEN;
  std::unique_ptr<uhh2::Selection> matched_fat_GEN;
  std::unique_ptr<uhh2::Selection> matched_sub;
  std::unique_ptr<uhh2::Selection> matched_fat;
  std::unique_ptr<uhh2::Selection> subjet_quality;
  std::unique_ptr<uhh2::Selection> subjet_quality_gen;
  std::unique_ptr<uhh2::Selection> lepton_sel;
  std::unique_ptr<uhh2::Selection> muon_highpt_sel;
  std::unique_ptr<uhh2::Selection> lepton_sel_gen;
  std::unique_ptr<uhh2::Selection> lepton_Nsel_gen;


  // get weight (with all SF and weight applied in previous cycle)
  Event::Handle<double>h_weight;


  // use top pt reweight module instead of derived sf?
  std::unique_ptr<TopPtReweight> ttbar_reweight;

  // handles for output
  Event::Handle<double>h_musf_central;
  Event::Handle<TTbarGen>h_ttbargen;
  Event::Handle<std::vector<TopJet>>h_hadjets_raw;
  Event::Handle<std::vector<TopJet>>h_hadjets;
  Event::Handle<std::vector<TopJet>>h_lepjets;
  Event::Handle<std::vector<TopJet>>h_fatjets;
  Event::Handle<std::vector<TopJet>>h_fatjets_raw;
  Event::Handle<bool>h_gensel_2;
  Event::Handle<bool>h_recsel_2;

  // store weights (or factor to get to these weights)
  Event::Handle<float> h_muid_up;
  Event::Handle<float> h_muid_down;
  Event::Handle<float> h_mutr_up;
  Event::Handle<float> h_mutr_down;
  Event::Handle<float> h_elid_up;
  Event::Handle<float> h_elid_down;
  Event::Handle<float> h_eltr_up;
  Event::Handle<float> h_eltr_down;
  Event::Handle<float> h_elreco_up;
  Event::Handle<float> h_elreco_down;
  Event::Handle<float> h_pu_up;
  Event::Handle<float> h_pu_down;
  Event::Handle<float> h_btag_up;
  Event::Handle<float> h_btag_down;
  Event::Handle<float> h_scale_upup;
  Event::Handle<float> h_scale_upnone;
  Event::Handle<float> h_scale_noneup;
  Event::Handle<float> h_scale_nonedown;
  Event::Handle<float> h_scale_downnone;
  Event::Handle<float> h_scale_downdown;
  Event::Handle<float> h_fsr_upsqrt2;
  Event::Handle<float> h_fsr_up2;
  Event::Handle<float> h_fsr_up4;
  Event::Handle<float> h_fsr_downsqrt2;
  Event::Handle<float> h_fsr_down2;
  Event::Handle<float> h_fsr_down4;
  Event::Handle<float> h_isr_up2;
  Event::Handle<float> h_isr_down2;
  ////

  Event::Handle<bool>h_matched;
  Event::Handle<bool>h_measure_rec;
  Event::Handle<bool>h_measure_gen;
  Event::Handle<bool>h_nomass_rec;
  Event::Handle<bool>h_nomass_gen;
  Event::Handle<bool>h_ptsub_rec;
  Event::Handle<bool>h_ptsub_gen;
  Event::Handle<bool>h_ptlepton_rec;
  Event::Handle<bool>h_ptlepton_gen;
  Event::Handle<bool>h_pt350_gen;
  Event::Handle<bool>h_pt350_rec;
  Event::Handle<bool>h_nobtag_rec;
  Event::Handle<bool>h_ttbar;
  Event::Handle<double>h_ttbar_SF;
  Event::Handle<double>h_mass_gen33;
  Event::Handle<double>h_mass_rec;
  Event::Handle<double>h_ptsub1_rec;
  Event::Handle<double>h_ptsub2_rec;
  Event::Handle<double>h_ptsub3_rec;
  Event::Handle<double>h_mW_rec;
  Event::Handle<double>h_mass_jms;
  Event::Handle<double>h_softdropmass_rec;
  Event::Handle<double>h_pt_gen33;
  Event::Handle<double>h_pt_rec;
  Event::Handle<double>h_genweight;
  Event::Handle<double>h_recweight;
  Event::Handle<double>h_genweight_ttfactor;
  Event::Handle<double>h_factor_2width;
  Event::Handle<double>h_factor_4width;
  Event::Handle<double>h_factor_8width;
  Event::Handle<std::vector<double>>h_bquark_pt;
  Event::Handle<std::vector<double>>h_pdf_weights;
  Event::Handle<std::vector<TopJet>>h_recjets_had;
  Event::Handle<std::vector<GenTopJet>>h_genjets33_had;
  Event::Handle<std::vector<GenTopJet>>h_genfatjets;

  std::unique_ptr<PartonShowerWeight> ps_weights;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // Weights
  Event::Handle<float>h_weight_muid;
  Event::Handle<float>h_weight_muid_up;
  Event::Handle<float>h_weight_muid_down;
  Event::Handle<float>h_weight_mutr;
  Event::Handle<float>h_weight_mutr_up;
  Event::Handle<float>h_weight_mutr_down;
  Event::Handle<float>h_weight_elid;
  Event::Handle<float>h_weight_elid_up;
  Event::Handle<float>h_weight_elid_down;
  Event::Handle<float>h_weight_eltr;
  Event::Handle<float>h_weight_eltr_up;
  Event::Handle<float>h_weight_eltr_down;
  Event::Handle<float>h_weight_elreco;
  Event::Handle<float>h_weight_elreco_up;
  Event::Handle<float>h_weight_elreco_down;
  Event::Handle<float>h_weight_pu;
  Event::Handle<float>h_weight_pu_up;
  Event::Handle<float>h_weight_pu_down;
  Event::Handle<float> h_weight_scale_upup;
  Event::Handle<float> h_weight_scale_upnone;
  Event::Handle<float> h_weight_scale_noneup;
  Event::Handle<float> h_weight_scale_nonedown;
  Event::Handle<float> h_weight_scale_downnone;
  Event::Handle<float> h_weight_scale_downdown;
  Event::Handle<float> h_weight_btag;
  // Event::Handle<float> h_weight_btag_jesup;
  // Event::Handle<float> h_weight_btag_jesdown;
  // Event::Handle<float> h_weight_btag_lfup;
  // Event::Handle<float> h_weight_btag_lfdown;
  // Event::Handle<float> h_weight_btag_hfup;
  // Event::Handle<float> h_weight_btag_hfdown;
  // Event::Handle<float> h_weight_btag_hfstats1up;
  // Event::Handle<float> h_weight_btag_hfstats1down;
  // Event::Handle<float> h_weight_btag_hfstats2up;
  // Event::Handle<float> h_weight_btag_hfstats2down;
  // Event::Handle<float> h_weight_btag_lfstats1up;
  // Event::Handle<float> h_weight_btag_lfstats1down;
  // Event::Handle<float> h_weight_btag_lfstats2up;
  // Event::Handle<float> h_weight_btag_lfstats2down;
  // Event::Handle<float> h_weight_btag_cferr1up;
  // Event::Handle<float> h_weight_btag_cferr1down;
  // Event::Handle<float> h_weight_btag_cferr2up;
  // Event::Handle<float> h_weight_btag_cferr2down;
  std::vector<Event::Handle<float>> h_weight_btag_upvars;
  std::vector<Event::Handle<float>> h_weight_btag_downvars;


  Event::Handle<float> h_ak8tau, h_ak8mass;

  //width reweight
  std::unique_ptr<tt_width_reweight> width2_reweight;
  std::unique_ptr<tt_width_reweight> width4_reweight;
  std::unique_ptr<tt_width_reweight> width8_reweight;

  //scale variation
  std::unique_ptr<AnalysisModule> scale_variation;
  std::unique_ptr<AnalysisModule> PUreweight;

  // Best fit
  std::unique_ptr<CorrectionFactor_JMS> BestFit;
  std::unique_ptr<JetMassScaleHists> h_XCone_JMS;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_Muon_PreSel01, h_Elec_PreSel01, h_MTopJet_PreSel01, h_Jets_PreSel01, h_XCone_cor_PreSel01;
  std::unique_ptr<Hists> h_Muon_PreSel02, h_Elec_PreSel02, h_MTopJet_PreSel02, h_Jets_PreSel02, h_XCone_cor_PreSel02;
  std::unique_ptr<Hists> h_Muon_PreSel03, h_Elec_PreSel03, h_MTopJet_PreSel03, h_Jets_PreSel03, h_XCone_cor_PreSel03;
  std::unique_ptr<Hists> h_Muon_PreSel03b, h_Elec_PreSel03b, h_MTopJet_PreSel03b, h_Jets_PreSel03b, h_XCone_cor_PreSel03b;
  std::unique_ptr<Hists> h_Muon_PreSel04, h_Elec_PreSel04, h_MTopJet_PreSel04, h_Jets_PreSel04, h_XCone_cor_PreSel04;

  std::unique_ptr<Hists> h_Muon, h_Elec, h_MTopJet, h_Jets;
  std::unique_ptr<Hists> h_CorrectionHists, h_CorrectionHists_after, h_CorrectionHists_WJets;
  std::unique_ptr<Hists> h_RecGenHists_ak4, h_RecGenHists_ak4_noJEC;
  std::unique_ptr<Hists> h_weights01, h_weights02, h_weights03, h_weights04;

  std::unique_ptr<Hists> h_XCone_cor, h_XCone_jec, h_XCone_raw;
  std::unique_ptr<Hists> h_XCone_puppi;
  std::unique_ptr<Hists> h_XCone_cor_SF, h_XCone_jec_SF, h_XCone_raw_SF;
  std::unique_ptr<Hists> h_XCone_cor_SF_pt400, h_XCone_cor_SF_pt450, h_XCone_cor_SF_pt530;
  std::unique_ptr<Hists> h_XCone_cor_pt350, h_XCone_cor_noptcut;
  std::unique_ptr<Hists> h_XCone_cor_subjets, h_XCone_jec_subjets, h_XCone_raw_subjets;
  std::unique_ptr<Hists> h_XCone_cor_subjets_SF, h_XCone_jec_subjets_SF, h_XCone_raw_subjets_SF;
  std::unique_ptr<Hists> h_XCone_cor_Sel_noSel, h_XCone_cor_Sel_noMass, h_XCone_cor_Sel_nobtag, h_XCone_cor_Sel_pt350, h_XCone_cor_Sel_ptsub, h_XCone_cor_subjets_Sel_ptsub;

  std::unique_ptr<Hists> h_XCone_cor_m, h_XCone_cor_u, h_XCone_cor_m_fat, h_XCone_cor_u_fat;

  std::unique_ptr<Hists> h_XCone_cor_PUlow, h_XCone_cor_PUmid,h_XCone_cor_PUhigh;
  std::unique_ptr<Hists> h_XCone_cor_PUlow_subjets, h_XCone_cor_PUmid_subjets, h_XCone_cor_PUhigh_subjets;

  std::unique_ptr<Hists> h_RecGenHists_lowPU, h_RecGenHists_medPU, h_RecGenHists_highPU, h_RecGenHists_lowPU_noJEC, h_RecGenHists_medPU_noJEC, h_RecGenHists_highPU_noJEC, h_RecGenHists_RecOnly_corr;
  std::unique_ptr<Hists> h_XCone_GEN_RecOnly, h_XCone_GEN_Both;
  std::unique_ptr<Hists> h_XCone_GEN_ST;
  std::unique_ptr<Hists> h_XCone_GEN_GenOnly, h_XCone_GEN_GenOnly_matched, h_XCone_GEN_GenOnly_unmatched, h_XCone_GEN_GenOnly_matched_fat, h_XCone_GEN_GenOnly_unmatched_fat;
  std::unique_ptr<Hists> h_RecGenHists_GenSelRecInfo, h_RecGenHists_GenSelRecInfo_lowPU, h_RecGenHists_GenSelRecInfo_midPU, h_RecGenHists_GenSelRecInfo_highPU;
  std::unique_ptr<Hists> h_RecGenHists_GenSelRecInfo_matched, h_RecGenHists_GenSelRecInfo_matched_lowPU, h_RecGenHists_GenSelRecInfo_matched_midPU, h_RecGenHists_GenSelRecInfo_matched_highPU;
  std::unique_ptr<Hists> h_RecGenHists_GenOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly_noJEC;
  std::unique_ptr<Hists> h_RecGenHists_Both, h_RecGenHists_Both_corr;
  std::unique_ptr<Hists> h_RecGenHists_subjets, h_RecGenHists_subjets_noJEC, h_RecGenHists_subjets_corrected, h_RecGenHists_subjets_matched;
  std::unique_ptr<Hists> h_RecGenHistsST_subjets, h_RecGenHistsST_subjets_noJEC, h_RecGenHistsST_subjets_corrected;
  std::unique_ptr<Hists> h_RecGenHists_subjets_lowPU, h_RecGenHists_subjets_medPU, h_RecGenHists_subjets_highPU;
  std::unique_ptr<Hists> h_RecGenHists_subjets_noJEC_lowPU, h_RecGenHists_subjets_noJEC_medPU, h_RecGenHists_subjets_noJEC_highPU;
  std::unique_ptr<Hists> h_RecGenHists_subjets_WJets, h_RecGenHists_subjets_noJEC_WJets;
  std::unique_ptr<Hists> h_GenParticles_RecOnly, h_GenParticles_GenOnly, h_GenParticles_Both;
  std::unique_ptr<Hists> h_GenParticles_SM, h_GenParticles_newWidth;
  std::unique_ptr<Hists> h_GenProcess;

  std::unique_ptr<Hists> h_XCone_cor_migration_pt, h_XCone_cor_migration_pt350, h_XCone_cor_migration_mass, h_XCone_cor_migration_btag;
  std::unique_ptr<Hists> h_750_xcone;

  std::unique_ptr<Hists> h_XCone_GEN_Sel_measurement, h_XCone_GEN_Sel_noMass, h_XCone_GEN_Sel_pt350, h_XCone_GEN_Sel_ptsub;
  std::unique_ptr<Hists> h_PDFHists;

  std::unique_ptr<RecoGenHists_xcone_topjet> h_comparison_topjet_xcone_pass_gen, h_comparison_topjet_xcone_pass_rec, h_comparison_topjet_xcone_pass_genrec, h_comparison_topjet_xcone;
  std::unique_ptr<RecoGenHists_xcone_topjet> h_comparison_topjet_xcone_pass_rec_masscut_120, h_comparison_topjet_xcone_pass_rec_masscut_130, h_comparison_topjet_xcone_pass_rec_masscut_140, h_comparison_topjet_xcone_pass_rec_masscut_150;

  std::unique_ptr<Hists> h_gen_weights, h_gen_weights_pass_gen, h_gen_weights_massbin_145, h_gen_weights_massbin_275;
  std::unique_ptr<Hists> h_gen_weight_300, h_gen_weight_gensel_300;
  std::unique_ptr<Hists> h_leptonictop, h_leptonictop_SF;

  std::unique_ptr<Hists> h_wrong_events_mtop_17;

  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  std::unique_ptr<uhh2::AnalysisModule> BTagScaleFactors, BTagReshape, CreateBTagJets;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;
  std::unique_ptr<uhh2::AnalysisModule> ele_id_SF, ele_trigger_SF, ele_reco_SF;

  BTag::algo btag_algo;
  BTag::wp wp_loose, wp_medium, wp_tight;

  // Jet Correction Stuff
  std::unique_ptr<JetCorrections_xcone> JetCorrections;
  std::unique_ptr<JER_Smearer_xcone> JERSmearing;
  std::unique_ptr<uhh2::AnalysisModule> Correction;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_corrected;
  std::unique_ptr<CombineXCone33> wmass_indices;
  ////

  string BTag_variation ="central";
  string MuScale_variation ="nominal";
  string MuTrigger_variation ="nominal";
  string ElID_variation ="nominal";
  string ElReco_variation ="nominal";
  string ElTrigger_variation ="nominal";
  string PU_variation ="nominal";
  TString jms_direction ="nominal";
  string corvar ="nominal";
  Year year;

  bool isMC;    //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part
  bool debug;
  bool isMTop;

  bool year_16;
  bool year_17;
  bool year_18;

  bool do_ps = false;
};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  debug = false;
  /*
  .██████ ████████ ██   ██
  ██         ██     ██ ██
  ██         ██      ███
  ██         ██     ██ ██
  .██████    ██    ██   ██
  */
  if(debug) cout << "\nStart Module - CTX\n";

  year_16 = false;
  year_17 = false;
  year_18 = false;
  year = extract_year(ctx);

  if(year == Year::is2016v3) year_16 = true;
  else if(year == Year::is2017v2) year_17 = true;
  else if(year == Year::is2018) year_18 = true;
  else throw runtime_error("In PostSelectionModule: This Event is not from 2016v3, 2017v2 or 2018!");

  /*************************** CONFIGURATION **********************************************************************************/
  isMC = (ctx.get("dataset_type") == "MC");

  TString dataset_version = (TString) ctx.get("dataset_version");
  if(dataset_version.Contains("TTbar") || dataset_version.Contains("TTTo")) isTTbar = true;
  else  isTTbar = false;

  // Avoid empty fatjets for JMS
  if(dataset_version.Contains("_mtop")&&year_17) isMTop = true;
  else                                           isMTop = false;

  if(dataset_version.Contains("TTbar") && !dataset_version.Contains("_mtop") && !year_16) do_ps = true;

  if(isMTop) h_wrong_events_mtop_17.reset(new CountingEventHists(ctx, "Wrong_17_mTop_events"));

  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  if(isTTbar) ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {
    std::string log("MTopJetPostSelectionModule::MTopJetPostSelectionModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";
    throw std::runtime_error(log);
  }

  if(debug) cout << "\nSome reweight\n";

  // Top width reweigh
  width2_reweight.reset(new tt_width_reweight(ctx, 2.0));
  width4_reweight.reset(new tt_width_reweight(ctx, 4.0));
  width8_reweight.reset(new tt_width_reweight(ctx, 8.0));

  // Top PT reweight
  ttbar_reweight.reset(new TopPtReweight(ctx,0.0615,-0.0005,"","weight_ttbar",true)); // 13 TeV

  //scale variation
  scale_variation.reset(new MCScaleVariation(ctx));

  // PU reweighting
  PU_variation = ctx.get("PU_variation","central");
  PUreweight.reset(new MCPileupReweight(ctx, PU_variation));

  // PS reweighting
  string PS_variation = "central";
  PS_variation = ctx.get("PS_variation", "central");
  ps_weights.reset(new PartonShowerWeight(ctx, PS_variation));

  h_recjets_had = ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  if(isTTbar) h_genjets33_had = ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");
  h_genfatjets = ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  // weight handles
  h_weight_muid = ctx.get_handle<float>("weight_sfmu_tightID");
  h_weight_muid_up = ctx.get_handle<float>("weight_sfmu_tightID_up");
  h_weight_muid_down = ctx.get_handle<float>("weight_sfmu_tightID_down");
  h_weight_mutr = ctx.get_handle<float>("weight_sfmu_trigger");
  h_weight_mutr_up = ctx.get_handle<float>("weight_sfmu_trigger_up");
  h_weight_mutr_down = ctx.get_handle<float>("weight_sfmu_trigger_down");
  h_weight_elid = ctx.get_handle<float>("weight_sfelec_tightID");
  h_weight_elid_up = ctx.get_handle<float>("weight_sfelec_tightID_up");
  h_weight_elid_down = ctx.get_handle<float>("weight_sfelec_tightID_down");
  h_weight_elreco = ctx.get_handle<float>("weight_sfelec_reco");
  h_weight_elreco_up = ctx.get_handle<float>("weight_sfelec_reco_up");
  h_weight_elreco_down = ctx.get_handle<float>("weight_sfelec_reco_down");
  h_weight_eltr = ctx.get_handle<float>("weight_sfelec_trigger");
  h_weight_eltr_up = ctx.get_handle<float>("weight_sfelec_trigger_up");
  h_weight_eltr_down = ctx.get_handle<float>("weight_sfelec_trigger_down");
  h_weight_pu = ctx.get_handle<float>("weight_pu");
  h_weight_pu_up = ctx.get_handle<float>("weight_pu_up");
  h_weight_pu_down = ctx.get_handle<float>("weight_pu_down");
  h_weight_scale_upup = ctx.get_handle<float>("weight_murmuf_upup");
  h_weight_scale_upnone = ctx.get_handle<float>("weight_murmuf_upnone");
  h_weight_scale_noneup = ctx.get_handle<float>("weight_murmuf_noneup");
  h_weight_scale_nonedown = ctx.get_handle<float>("weight_murmuf_nonedown");
  h_weight_scale_downnone = ctx.get_handle<float>("weight_murmuf_downnone");
  h_weight_scale_downdown = ctx.get_handle<float>("weight_murmuf_downdown");

  h_weight_btag = ctx.get_handle<float>("weight_btagdisc_central");

  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_jesup"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfup"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfup"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfstats1up"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfstats2up"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfstats1up"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfstats2up"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_cferr1up"));
  h_weight_btag_upvars.push_back(ctx.get_handle<float>("weight_btagdisc_cferr2up"));

  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_jesdown"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfdown"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfdown"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfstats1down"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_hfstats2down"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfstats1down"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_lfstats2down"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_cferr1down"));
  h_weight_btag_downvars.push_back(ctx.get_handle<float>("weight_btagdisc_cferr2down"));



  /*************************** Setup Subjet Corrector **********************************************************************************/
  // Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected"));
  /*************************** Setup Selections ****************************************************************************************/
  if(debug) cout << "\nSetup Selection (RECO)\n";

  // RECO Selection
  // chose: XCone33_had_Combined_Corrected,  XCone33_had_Combined_noJEC,  XCone33_had_Combined
  const std::string& jet_label_had = "XCone33_had_Combined_Corrected";
  const std::string& jet_label_lep = "XCone33_lep_Combined_Corrected";
  njet_had.reset(new NJetXCone(ctx, jet_label_had, 1));
  njet_lep.reset(new NJetXCone(ctx, jet_label_lep, 1));


  subjet_quality.reset(new SubjetQuality(ctx, jet_label_had, 30, 2.5));
  pt_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 400));
  pt450_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 450));
  pt530_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 530));
  pt2_sel.reset(new LeadingRecoJetPT(ctx, jet_label_lep, 10));
  eta_sel.reset(new LeadingRecoJetETA(ctx, jet_label_had, 2.5));
  pt350_sel.reset(new LeadingRecoJetPT(ctx, jet_label_had, 350));
  mass_sel.reset(new MassCutXCone(ctx, jet_label_had, jet_label_lep));

  MuonId muid = AndId<Muon>(MuonID(Muon::Tight), PtEtaCut(60., 2.4));
  ElectronId eleid;
  if(year_16) eleid = AndId<Electron>(PtEtaSCCut(60., 2.4), ElectronID_Summer16_tight_noIso);
  else        eleid = AndId<Electron>(PtEtaSCCut(60., 2.4), ElectronID_Fall17_tight_noIso);

  if(channel_ == muon) lepton_sel.reset(new NMuonSelection(1, -1, muid));
  else                 lepton_sel.reset(new NElectronSelection(1, -1, eleid));

  if(debug) cout << "\nSetup Selection (GEN)\n";
  // GEN Selection
  if(channel_ == muon){
    lepton_sel_gen.reset(new GenMuonSel(ctx, 60.));
    lepton_Nsel_gen.reset(new GenMuonCount(ctx)); // Count is confusing here, since it only chooses the one lepton from the SemiLepDecay
  }
  else{
    lepton_sel_gen.reset(new GenElecSel(ctx, 60.));
    lepton_Nsel_gen.reset(new GenElecCount(ctx));
  }
  subjet_quality_gen.reset(new SubjetQuality_gen(ctx, "GEN_XCone33_had_Combined", 30, 2.5));
  pt_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 400));
  pt2_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_lep_Combined", 10));
  pt350_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 350));
  mass_gensel.reset(new MassCut_gen(ctx, "GEN_XCone33_had_Combined", "GEN_XCone33_lep_Combined"));
  matched_sub_GEN.reset(new Matching_XCone33GEN(ctx, "GEN_XCone33_had_Combined", true));
  matched_fat_GEN.reset(new Matching_XCone33GEN(ctx, "GEN_XCone33_had_Combined", false));

  // Selection for matching reco jets to gen particles
  if(isTTbar) matched_sub.reset(new Matching_XCone33(ctx, true));
  if(isTTbar) matched_fat.reset(new Matching_XCone33(ctx, false));

  if(debug) cout << "\nBTAG\n";
  //B-Tagging
  btag_algo = BTag::DEEPJET;
  wp_loose = BTag::WP_LOOSE;
  wp_medium = BTag::WP_MEDIUM;
  wp_tight = BTag::WP_TIGHT;

  DeepjetLoose = BTag(btag_algo, wp_loose);
  DeepjetMedium = BTag(btag_algo, wp_medium);
  DeepjetTight = BTag(btag_algo, wp_tight);

  if(debug) cout << "\nScale factors\n";
  // Scale factors
  BTag_variation = ctx.get("BTag_variation","central");
  MuScale_variation = ctx.get("MuScale_variation","nominal");
  MuTrigger_variation = ctx.get("MuTrigger_variation","nominal");
  ElID_variation = ctx.get("ElID_variation","nominal");
  ElReco_variation = ctx.get("ElReco_variation","nominal");
  ElTrigger_variation = ctx.get("ElTrigger_variation","nominal");

  string bjet_collection = "highest_deepjet_jet";
  CreateBTagJets.reset(new StoreBJet(ctx, bjet_collection));
  BTagReshape.reset(new MCBTagDiscriminantReweighting(ctx, btag_algo, bjet_collection, BTag_variation, "iterativefit", "", "BTagCalibration"));
  // BTagScaleFactors.reset(new MCBTagScaleFactor(ctx, btag_algo, wp_tight,"jets",BTag_variation,"comb"));

  if(year_16){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root","NUM_TightID_DEN_genTracks_eta_pt",1, "tightID", false, MuScale_variation));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "trigger", false, MuTrigger_variation));
    ele_trigger_SF.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/electrSF", ElTrigger_variation));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "reco", ElReco_variation));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "tightID", ElID_variation));
  }
  else if(year_17){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonID_94X_RunBCDEF_SF_ID.root","NUM_TightID_DEN_genTracks_pt_abseta",1, "tightID", true, MuScale_variation));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","Mu50_PtEtaBins/pt_abseta_ratio",1, "trigger", true, MuTrigger_variation));
    ele_trigger_SF.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/electrSF", ElTrigger_variation));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2017/2017_ElectronTight.root", 1.0, "tightID", ElID_variation));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", 1.0, "reco", ElReco_variation));
  }
  else if(year_18){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_ID_SF_RunABCD.root","NUM_TightID_DEN_TrackerMuons_pt_abseta",1, "tightID", true, MuScale_variation));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_Trigger_Eff_SF_AfterMuonHLTUpdate.root","Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/pt_abseta_ratio",1, "trigger", true, MuTrigger_variation));
    ele_trigger_SF.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/electrSF", ElTrigger_variation));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2018/2018_ElectronTight.root", 1.0, "tightID", ElID_variation));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/common/data/2018/egammaEffi.txt_EGM2D_updatedAll.root", 1.0, "reco", ElReco_variation));
  }
  else throw runtime_error("In PostSelectionModule: There is no Event from 2016_v3, 2017_v2 or 2018!");

  // UPDATE THIS FOR 2017 + 2018!!

  /*************************** Set up Hists classes **********************************************************************************/
  if(debug) cout << "\nSetup Hists \n";
  //BTagMCEfficiencyHists
  BTagEffHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",DeepjetTight));
  //LeptonicTop_Hists
  h_leptonictop.reset(new LeptonicTop_Hists(ctx, "LeptonicTop"));
  h_leptonictop_SF.reset(new LeptonicTop_Hists(ctx, "LeptonicTop_SF"));

  //750GeV hists
  h_750_xcone.reset(new RecoHists_xcone(ctx, "750_xcone", "cor"));

  // XCone Combined Jet
  h_XCone_raw.reset(new RecoHists_xcone(ctx, "XCone_raw", "raw"));
  h_XCone_cor.reset(new RecoHists_xcone(ctx, "XCone_cor", "cor"));
  h_XCone_jec.reset(new RecoHists_xcone(ctx, "XCone_jec", "jec"));
  h_XCone_raw_SF.reset(new RecoHists_xcone(ctx, "XCone_raw_SF", "raw"));
  h_XCone_cor_SF.reset(new RecoHists_xcone(ctx, "XCone_cor_SF", "cor"));
  h_XCone_jec_SF.reset(new RecoHists_xcone(ctx, "XCone_jec_SF", "jec"));

  h_XCone_cor_SF_pt400.reset(new RecoHists_xcone(ctx, "XCone_cor_SF_pt400", "cor"));
  h_XCone_cor_SF_pt450.reset(new RecoHists_xcone(ctx, "XCone_cor_SF_pt450", "cor"));
  h_XCone_cor_SF_pt530.reset(new RecoHists_xcone(ctx, "XCone_cor_SF_pt530", "cor"));

  h_XCone_cor_pt350.reset(new RecoHists_xcone(ctx, "XCone_cor_pt350", "cor"));
  h_XCone_cor_noptcut.reset(new RecoHists_xcone(ctx, "XCone_cor_noptcut", "cor"));

  h_XCone_puppi.reset(new RecoHists_puppi(ctx, "XCone_puppi"));

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
  h_XCone_cor_Sel_ptsub.reset(new RecoHists_xcone(ctx, "XCone_cor_Sel_ptsub", "cor"));
  h_XCone_cor_subjets_Sel_ptsub.reset(new SubjetHists_xcone(ctx, "h_XCone_cor_subjets_Sel_ptsub", "cor"));

  // PU dependence
  h_XCone_cor_PUlow.reset(new RecoHists_xcone(ctx, "XCone_cor_PUlow", "cor"));
  h_XCone_cor_PUmid.reset(new RecoHists_xcone(ctx, "XCone_cor_PUmid", "cor"));
  h_XCone_cor_PUhigh.reset(new RecoHists_xcone(ctx, "XCone_cor_PUhigh", "cor"));

  h_XCone_cor_PUlow_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUlow_subjets", "cor"));
  h_XCone_cor_PUmid_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUmid_subjets", "cor"));
  h_XCone_cor_PUhigh_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_PUhigh_subjets", "cor"));

  // Matching hists
  h_XCone_cor_m.reset(new RecoHists_xcone(ctx, "XCone_cor_matched", "cor"));
  h_XCone_cor_u.reset(new RecoHists_xcone(ctx, "XCone_cor_unmatched", "cor"));
  h_XCone_cor_m_fat.reset(new RecoHists_xcone(ctx, "XCone_cor_matched_fat", "cor"));
  h_XCone_cor_u_fat.reset(new RecoHists_xcone(ctx, "XCone_cor_unmatched_fat", "cor"));

  // Weight Hists
  h_weights01.reset(new WeightHists(ctx, "WeightHists01_noSF"));
  h_weights02.reset(new WeightHists(ctx, "WeightHists02_PU"));
  h_weights03.reset(new WeightHists(ctx, "WeightHists03_Lepton"));
  h_weights04.reset(new WeightHists(ctx, "WeightHists04_BTag"));

  // AK8 N-subjetiness
  h_comparison_topjet_xcone.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone", isTTbar, 0));
  h_comparison_topjet_xcone_pass_gen.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_gen", isTTbar, 0));
  h_comparison_topjet_xcone_pass_rec.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_rec", isTTbar, 0));
  h_comparison_topjet_xcone_pass_genrec.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_genrec", isTTbar, 0));

  h_comparison_topjet_xcone_pass_rec_masscut_120.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_rec_masscut_120", isTTbar, 120));
  h_comparison_topjet_xcone_pass_rec_masscut_130.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_rec_masscut_130", isTTbar, 130));
  h_comparison_topjet_xcone_pass_rec_masscut_140.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_rec_masscut_140", isTTbar, 140));
  h_comparison_topjet_xcone_pass_rec_masscut_150.reset(new RecoGenHists_xcone_topjet(ctx, "comparison_topjet_xcone_pass_rec_masscut_150", isTTbar, 150));

  // Control Hists for FSRx_4
  if(isTTbar){
    h_gen_weights.reset(new GenWeightRangeHists(ctx, "Gen_Weights"));
    h_gen_weights_massbin_145.reset(new GenWeightRangeHists(ctx, "Gen_Weights_Massbin_145"));
    h_gen_weights_pass_gen.reset(new GenWeightRangeHists(ctx, "Gen_Weights_pass_gen"));
    h_gen_weights_massbin_275.reset(new GenWeightRangeHists(ctx, "Gen_Weights_Massbin_275"));

    h_gen_weight_300.reset(new CountingEventHists(ctx, "Gen_Weight_300"));
    h_gen_weight_gensel_300.reset(new CountingEventHists(ctx, "Gen_Weight_gensel_300"));
  }

  h_MTopJet.reset(new MTopJetHists(ctx, "EventHists"));
  h_Muon.reset(new MuonHists(ctx, "MuonHists"));
  h_Elec.reset(new ElectronHists(ctx, "ElecHists"));
  h_Jets.reset(new JetHists(ctx, "JetHits"));

  h_MTopJet_PreSel01.reset(new MTopJetHists(ctx, "PreSel01_Event"));
  h_Muon_PreSel01.reset(new MuonHists(ctx, "PreSel01_Muon"));
  h_Elec_PreSel01.reset(new ElectronHists(ctx, "PreSel01_Elec"));
  h_Jets_PreSel01.reset(new JetHists(ctx, "PreSel01_Jet"));
  h_XCone_cor_PreSel01.reset(new RecoHists_xcone(ctx, "PreSel01_XCone", "cor"));

  h_MTopJet_PreSel02.reset(new MTopJetHists(ctx, "PreSel02_Event"));
  h_Muon_PreSel02.reset(new MuonHists(ctx, "PreSel02_Muon"));
  h_Elec_PreSel02.reset(new ElectronHists(ctx, "PreSel02_Elec"));
  h_Jets_PreSel02.reset(new JetHists(ctx, "PreSel02_Jet"));
  h_XCone_cor_PreSel02.reset(new RecoHists_xcone(ctx, "PreSel02_XCone", "cor"));

  h_MTopJet_PreSel03.reset(new MTopJetHists(ctx, "PreSel03_Event"));
  h_Muon_PreSel03.reset(new MuonHists(ctx, "PreSel03_Muon"));
  h_Elec_PreSel03.reset(new ElectronHists(ctx, "PreSel03_Elec"));
  h_Jets_PreSel03.reset(new JetHists(ctx, "PreSel03_Jet"));
  h_XCone_cor_PreSel03.reset(new RecoHists_xcone(ctx, "PreSel03_XCone", "cor"));

  h_MTopJet_PreSel03b.reset(new MTopJetHists(ctx, "PreSel03b_Event"));
  h_Muon_PreSel03b.reset(new MuonHists(ctx, "PreSel03b_Muon"));
  h_Elec_PreSel03b.reset(new ElectronHists(ctx, "PreSel03b_Elec"));
  h_Jets_PreSel03b.reset(new JetHists(ctx, "PreSel03b_Jet"));
  h_XCone_cor_PreSel03b.reset(new RecoHists_xcone(ctx, "PreSel03b_XCone", "cor"));

  h_MTopJet_PreSel04.reset(new MTopJetHists(ctx, "PreSel04_Event"));
  h_Muon_PreSel04.reset(new MuonHists(ctx, "PreSel04_Muon"));
  h_Elec_PreSel04.reset(new ElectronHists(ctx, "PreSel04_Elec"));
  h_Jets_PreSel04.reset(new JetHists(ctx, "PreSel04_Jet"));
  h_XCone_cor_PreSel04.reset(new RecoHists_xcone(ctx, "PreSel04_XCone", "cor"));

  h_XCone_JMS.reset(new JetMassScaleHists(ctx, "JetMassScaleHists"));

  if(debug) cout << "\nMC Hists\n";
  if(isMC){

    h_CorrectionHists.reset(new CorrectionHists_subjets(ctx, "CorrectionHists", "jec"));
    h_CorrectionHists_after.reset(new CorrectionHists_subjets(ctx, "CorrectionHists_after", "cor"));
    h_CorrectionHists_WJets.reset(new CorrectionHists_subjets(ctx, "CorrectionHists_WJets", "jec"));


    h_XCone_GEN_Sel_measurement.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_measurement"));
    h_XCone_GEN_Sel_noMass.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_noMass"));
    h_XCone_GEN_Sel_pt350.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_pt350"));
    h_XCone_GEN_Sel_ptsub.reset(new GenHists_xcone(ctx, "XCone_GEN_Sel_ptsub"));

    h_XCone_GEN_GenOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly"));
    h_XCone_GEN_GenOnly_matched.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly_matched"));
    h_XCone_GEN_GenOnly_unmatched.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly_unmatched"));
    h_XCone_GEN_GenOnly_matched_fat.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly_matched_fat"));
    h_XCone_GEN_GenOnly_unmatched_fat.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly_unmatched_fat"));
    h_XCone_GEN_RecOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_RecOnly"));
    h_XCone_GEN_Both.reset(new GenHists_xcone(ctx, "XCone_GEN_Both"));

    h_XCone_GEN_ST.reset(new GenHists_xcone(ctx, "XCone_GEN_ST"));

    h_GenParticles_GenOnly.reset(new GenHists_particles(ctx, "GenParticles_GenOnly"));
    h_GenParticles_RecOnly.reset(new GenHists_particles(ctx, "GenParticles_RecOnly"));
    h_GenParticles_Both.reset(new GenHists_particles(ctx, "GenParticles_Both"));

    h_GenParticles_SM.reset(new GenHists_particles(ctx, "GenParticles_SM"));
    h_GenParticles_newWidth.reset(new GenHists_particles(ctx, "GenParticles_newWidth"));

    h_GenProcess.reset(new GenHists_process(ctx, "GenProcess"));

    h_RecGenHists_subjets.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets", "jec"));
    h_RecGenHists_subjets_matched.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_matched", "jec"));
    h_RecGenHists_subjets_lowPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_lowPU", "jec"));
    h_RecGenHists_subjets_medPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_medPU", "jec"));
    h_RecGenHists_subjets_highPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_highPU", "jec"));
    h_RecGenHists_subjets_noJEC.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC", "raw"));
    h_RecGenHists_subjets_noJEC_lowPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC_lowPU", "raw"));
    h_RecGenHists_subjets_noJEC_medPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC_medPU", "raw"));
    h_RecGenHists_subjets_noJEC_highPU.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC_highPU", "raw"));
    h_RecGenHists_subjets_corrected.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_corrected", "cor"));

    h_RecGenHists_subjets_WJets.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_WJets", "jec"));
    h_RecGenHists_subjets_noJEC_WJets.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC_WJets", "raw"));

    h_RecGenHistsST_subjets.reset(new RecoGenHists_subjets(ctx, "RecGenHistsST_subjets", "jec"));
    h_RecGenHistsST_subjets_noJEC.reset(new RecoGenHists_subjets(ctx, "RecGenHistsST_subjets_noJEC", "raw"));
    h_RecGenHistsST_subjets_corrected.reset(new RecoGenHists_subjets(ctx, "RecGenHistsST_subjets_corrected", "cor"));


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
    h_RecGenHists_GenSelRecInfo.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo", "cor"));
    h_RecGenHists_GenSelRecInfo_lowPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_lowPU", "cor"));
    h_RecGenHists_GenSelRecInfo_midPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_midPU", "cor"));
    h_RecGenHists_GenSelRecInfo_highPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_highPU", "cor"));
    h_RecGenHists_GenSelRecInfo_matched.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_matched", "cor"));
    h_RecGenHists_GenSelRecInfo_matched_lowPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_matched_lowPU", "cor"));
    h_RecGenHists_GenSelRecInfo_matched_midPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_matched_midPU", "cor"));
    h_RecGenHists_GenSelRecInfo_matched_highPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenSelRecInfo_matched_highPU", "cor"));
  }

  // PDF hists
  h_PDFHists.reset(new PDFHists(ctx, "PDFHists"));

  // Best Fit
  if(isMC){
    if(debug) cout << "\nBestFit\n";
    jms_direction =  ctx.get("JetMassScale_direction","nominal");
    BestFit.reset(new CorrectionFactor_JMS(ctx, "XCone33_had_Combined_noJEC", "genXCone33TopJets", year));
  }

  // Jet Correction Stuff
  corvar = ctx.get("JetCorrection_direction","nominal"); // "nominal" nur wenn nicht up oder down
  // correct subjets (JEC + additional correction)
  JetCorrections.reset(new JetCorrections_xcone());
  JetCorrections->init(ctx, "xconeCHS");
  // smear jets after Correction
  JERSmearing.reset(new JER_Smearer_xcone());
  JERSmearing->init(ctx, "xconeCHS", "genXCone33TopJets", "sub");

  // Xcone correction
  if(year_16) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2016"));
  else if(year_17) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2017"));
  else if(year_18) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2018"));
  else throw runtime_error("In PostSelectionModule: There is no Event from 2016_v2, 2017_v2 or 2018!");

  // combine jets after correction
  jetprod_reco.reset(new CombineXCone33(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "xconeCHS"));
  jetprod_reco_corrected.reset(new CombineXCone33(ctx, "XCone33_had_Combined_Corrected", "XCone33_lep_Combined_Corrected", "xconeCHS_Corrected"));

  wmass_indices.reset(new CombineXCone33(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "xconeCHS"));

  // get handles
  if(debug) cout << "\nHandels\n";
  h_weight = ctx.get_handle<double>("weight");
  h_gensel_2 = ctx.get_handle<bool>("passed_gensel_2");
  h_recsel_2 = ctx.get_handle<bool>("passed_recsel_2");
  h_musf_central = ctx.get_handle<double>("passed_recsel_2");
  // h_weight_btag = ctx.get_handle<float>("weight_btag");

  h_hadjets_raw=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_noJEC");

  h_hadjets    =ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  h_lepjets    =ctx.get_handle<std::vector<TopJet>>("XCone33_lep_Combined_Corrected");
  h_fatjets    =ctx.get_handle<std::vector<TopJet>>("xconeCHS_Corrected");
  h_fatjets_raw=ctx.get_handle<std::vector<TopJet>>("xconeCHS");


  // undeclare event output (jet collections etc) to get small root files
  ctx.undeclare_all_event_output();
  // declare event output used for unfolding
  h_matched = ctx.declare_event_output<bool>("matched");
  h_measure_rec = ctx.declare_event_output<bool>("passed_measurement_rec");
  h_measure_gen = ctx.declare_event_output<bool>("passed_measurement_gen");
  h_nomass_rec = ctx.declare_event_output<bool>("passed_massmigration_rec");
  h_nomass_gen = ctx.declare_event_output<bool>("passed_massmigration_gen");
  h_ptsub_rec = ctx.declare_event_output<bool>("passed_subptmigration_rec");
  h_ptsub_gen = ctx.declare_event_output<bool>("passed_subptmigration_gen");
  h_pt350_rec = ctx.declare_event_output<bool>("passed_ptmigration_rec");
  h_pt350_gen = ctx.declare_event_output<bool>("passed_ptmigration_gen");
  h_nobtag_rec = ctx.declare_event_output<bool>("passed_btagmigration_rec");
  h_ptlepton_rec = ctx.declare_event_output<bool>("passed_leptonptmigration_rec");
  h_ptlepton_gen = ctx.declare_event_output<bool>("passed_leptonptmigration_gen");;
  h_ttbar = ctx.declare_event_output<bool>("is_TTbar");
  h_ttbar_SF = ctx.declare_event_output<double>("TTbar_SF");
  h_mass_gen33 = ctx.declare_event_output<double>("Mass_Gen33");
  h_mass_rec = ctx.declare_event_output<double>("Mass_Rec_old");
  h_mass_jms = ctx.declare_event_output<double>("Mass_Rec");
  h_ptsub1_rec = ctx.declare_event_output<double>("ptsub1");
  h_ptsub2_rec = ctx.declare_event_output<double>("ptsub2");
  h_ptsub3_rec = ctx.declare_event_output<double>("ptsub3");
  h_mW_rec = ctx.declare_event_output<double>("mW");

  h_pt_gen33 = ctx.declare_event_output<double>("Pt_Gen33");
  h_pt_rec = ctx.declare_event_output<double>("Pt_Rec");
  h_genweight = ctx.declare_event_output<double>("gen_weight");
  h_recweight = ctx.declare_event_output<double>("rec_weight");
  h_genweight_ttfactor = ctx.declare_event_output<double>("gen_weight_ttfactor");
  h_factor_2width = ctx.declare_event_output<double>("factor_2width");
  h_factor_4width = ctx.declare_event_output<double>("factor_4width");
  h_factor_8width = ctx.declare_event_output<double>("factor_8width");
  h_pdf_weights = ctx.declare_event_output<vector<double>>("pdf_weights");
  h_bquark_pt = ctx.declare_event_output<vector<double>>("bquark_pt");
  h_softdropmass_rec = ctx.declare_event_output<double>("SoftDropMass_Rec");
  h_muid_up = ctx.declare_event_output<float>("sf_muid_up");
  h_muid_down = ctx.declare_event_output<float>("sf_muid_down");
  h_mutr_up = ctx.declare_event_output<float>("sf_mutr_up");
  h_mutr_down = ctx.declare_event_output<float>("sf_mutr_down");
  h_elid_up = ctx.declare_event_output<float>("sf_elid_up");
  h_elid_down = ctx.declare_event_output<float>("sf_elid_down");
  h_eltr_up = ctx.declare_event_output<float>("sf_eltr_up");
  h_eltr_down = ctx.declare_event_output<float>("sf_eltr_down");
  h_elreco_up = ctx.declare_event_output<float>("sf_elreco_up");
  h_elreco_down = ctx.declare_event_output<float>("sf_elreco_down");
  h_pu_up = ctx.declare_event_output<float>("sf_pu_up");
  h_pu_down = ctx.declare_event_output<float>("sf_pu_down");
  h_btag_up = ctx.declare_event_output<float>("sf_btag_up");
  h_btag_down = ctx.declare_event_output<float>("sf_btag_down");
  h_scale_upup = ctx.declare_event_output<float>("sf_scale_upup");
  h_scale_upnone = ctx.declare_event_output<float>("sf_scale_upnone");
  h_scale_noneup = ctx.declare_event_output<float>("sf_scale_noneup");
  h_scale_nonedown = ctx.declare_event_output<float>("sf_scale_nonedown");
  h_scale_downnone = ctx.declare_event_output<float>("sf_scale_downnone");
  h_scale_downdown = ctx.declare_event_output<float>("sf_scale_downdown");
  h_fsr_upsqrt2 = ctx.declare_event_output<float>("sf_fsr_upsqrt2");
  h_fsr_up2 = ctx.declare_event_output<float>("sf_fsr_up2");
  h_fsr_up4 = ctx.declare_event_output<float>("sf_fsr_up4");
  h_fsr_downsqrt2 = ctx.declare_event_output<float>("sf_fsr_downsqrt2");
  h_fsr_down2 = ctx.declare_event_output<float>("sf_fsr_down2");
  h_fsr_down4 = ctx.declare_event_output<float>("sf_fsr_down4");
  h_isr_up2 = ctx.declare_event_output<float>("sf_isr_up2");
  h_isr_down2 = ctx.declare_event_output<float>("sf_isr_down2");
  h_ak8tau = ctx.declare_event_output<float>("tau32_ak8_had");
  h_ak8mass = ctx.declare_event_output<float>("mass_ak8_had");


  /*********************************************************************************************************************************/

}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool MTopJetPostSelectionModule::process(uhh2::Event& event){

  /***************************  First do JEC/JER/XCONE *******************************************************************************/

  JetCorrections->process(event);          // apply AK4 JEC to subjets of the original Fatjet Collection
  JERSmearing->process(event);             // apply JER smearing to subjets
  jetprod_reco->process(event);            // now store sum of 'jec' subjets
  Correction->process(event);              // apply additional correction (a new 'cor' TopJet Collection is generated)
  jetprod_reco_corrected->process(event);  // finally store sum of 'cor' subjets


  /************************************************************************************************************************************/

  if(debug) cout << "\nStart process --------------------------------------------------------------\n";
  if(isMTop){
    std::vector<TopJet> test_fatjets = event.get(h_fatjets);
    for(auto & jet: test_fatjets){
      if (jet.subjets().size() ==0){
        h_wrong_events_mtop_17->fill(event);
        return false;
      }
    }
  }

  if(isMC) h_GenProcess->fill(event);

  // get bools for selections from root files
  bool passed_recsel;
  bool passed_gensel33;
  passed_recsel = event.get(h_recsel_2);
  passed_gensel33 = event.get(h_gensel_2);

  // check if event has one had and one lep jet
  if( !(njet_had->passes(event)) ) return false;
  if( !(njet_lep->passes(event)) ) return false;

  ////
  // fill ttbargen class
  if(isTTbar) ttgenprod->process(event);
  ////
  /***************************  get jets to write mass *****************************************************************************************************/
  std::vector<TopJet> rec_hadjets = event.get(h_recjets_had);
  double mass_rec = rec_hadjets.at(0).v4().M();
  double pt_rec = rec_hadjets.at(0).v4().Pt();
  event.set(h_mass_rec, mass_rec);
  event.set(h_pt_rec, pt_rec);
  vector<Jet> subjets = rec_hadjets.at(0).subjets();
  event.set(h_ptsub1_rec, subjets.at(0).v4().Pt());
  event.set(h_ptsub2_rec, subjets.at(1).v4().Pt());
  event.set(h_ptsub3_rec, subjets.at(2).v4().Pt());

  // JetMassScale --------------------------------------------------------------
  if(debug) cout << "JMS\n";
  vector<double> points;
  double mass_jms =0;
  double mass_wjms = 0;
  vector<int> WSubjetIndex = wmass_indices->GetWSubjetsIndices(event);
  TopJet wjet = wmass_indices->GetHadronicWJet(event, WSubjetIndex);
  double mass_wjet_rec = wjet.v4().M();

  if(isMC){
    // Data+MC
    if     (jms_direction == "nominal")  points = {0.7495, -0.3392}; // BestFit point
    else if(jms_direction == "upup")     points = {0.887884, -0.291257}; // clostest point, right
    else if(jms_direction == "updown")   points = {0.953384, -0.923157}; // furthest point, down
    else if(jms_direction == "downup")   points = {0.533584,  0.256443}; // furthest point, up
    else if(jms_direction == "downdown") points = {0.610384, -0.389157}; // clostest point, left
    else throw runtime_error("Your JetMassScale in the Config-File PostSel is not set correctly (nominal, upup, updown, downup, downdown)");
    if(debug) cout << "JMS - get_mass_BestFit\n";
    mass_jms = BestFit->get_mass_BestFit(event, points);
    event.set(h_mass_jms, mass_jms);
    mass_wjms = BestFit->get_wmass_BestFit(event, points, WSubjetIndex);
    event.set(h_mW_rec, mass_wjms);
  }
  else{
    event.set(h_mW_rec, mass_wjet_rec);
    event.set(h_mass_jms, mass_rec);
  }

  if(debug) cout << "Gen Mass\n";
  double mass_gen33 = 0;
  double pt_gen33 = 0;
  if(isTTbar){
    std::vector<GenTopJet> gen_hadjets33 = event.get(h_genjets33_had);
    if(gen_hadjets33.size()>0){
      mass_gen33 = gen_hadjets33.at(0).v4().M();
      pt_gen33 = gen_hadjets33.at(0).v4().Pt();
    }
    event.set(h_mass_gen33, mass_gen33);
    event.set(h_pt_gen33, pt_gen33);
  }
  else{
    event.set(h_mass_gen33, 0.); // set gen mass to 0 for data
    event.set(h_pt_gen33, 0.);   // set gen pt to 0 for data
  }

  if(debug) cout << "Write SD mass\n";
  // SoftDropMass for HadJet
  std::vector<TopJet> hadjets_raw = event.get(h_hadjets_raw);
  std::vector<TopJet> hadjets     = event.get(h_hadjets);
  std::vector<TopJet> lepjets     = event.get(h_lepjets);
  std::vector<TopJet> fatjets     = event.get(h_fatjets);
  TopJet hadjet;
  if(deltaR(lepjets.at(0), fatjets.at(0)) < deltaR(hadjets.at(0), fatjets.at(0))) hadjet = fatjets.at(1);
  else hadjet = fatjets.at(0);
  event.set(h_softdropmass_rec, hadjet.softdropmass());

  // also set bquark pt
  vector<double> bquark_pt = {0.0, 0.0};
  if(isTTbar){
    if(passed_gensel33){
      const auto & ttbargen = event.get(h_ttbargen);
      GenParticle bot1 = ttbargen.bTop();
      GenParticle bot2 = ttbargen.bAntitop();
      bquark_pt[0] = bot1.pt();
      bquark_pt[1] = bot2.pt();
    }
  }
  event.set(h_bquark_pt, bquark_pt);

  /***************************  apply weight *****************************************************************************************************/
  // bool reweight_ttbar = false;       // apply ttbar reweight?
  bool scale_ttbar = true;           // match MC and data cross-section (for plots only)?
  double SF_tt=0;
  if(year_16)      SF_tt = 0.728693;
  else if(year_17) SF_tt = 0.896576;
  else if(year_18) SF_tt = 0.847446;
  else throw runtime_error{"Year is wrong in PostSel: ttbar scale factor"}; // dummy
  // get lumi weight = genweight (inkl scale variation)
  // now get full weight from prev. Selection (weight = gen_weight * rec_weight)
  event.weight = event.get(h_weight);
  if(debug) cout << "Vary scale\n";
  scale_variation->process(event);
  if(isMC){
    // For the scale variations, the SF is directly stored, so just put them in own handles
    event.set(h_scale_upup, event.get(h_weight_scale_upup));
    event.set(h_scale_upnone, event.get(h_weight_scale_upnone));
    event.set(h_scale_noneup, event.get(h_weight_scale_noneup));
    event.set(h_scale_nonedown, event.get(h_weight_scale_nonedown));
    event.set(h_scale_downnone, event.get(h_weight_scale_downnone));
    event.set(h_scale_downdown, event.get(h_weight_scale_downdown));
  }
  else{
    event.set(h_scale_upup, 1.);
    event.set(h_scale_upnone, 1.);
    event.set(h_scale_noneup, 1.);
    event.set(h_scale_nonedown, 1.);
    event.set(h_scale_downnone, 1.);
    event.set(h_scale_downdown, 1.);
  }

  if(debug) cout << "Vary PS\n";
  if(isTTbar && do_ps){
    ps_weights->process(event);
    vector<float> ps_factors = ps_weights->get_factors(event);
    if(ps_factors.size()>0){
      event.set(h_fsr_upsqrt2, ps_factors[0]);
      event.set(h_fsr_up2, ps_factors[1]);
      event.set(h_fsr_up4, ps_factors[2]);
      event.set(h_fsr_downsqrt2, ps_factors[3]);
      event.set(h_fsr_down2, ps_factors[4]);
      event.set(h_fsr_down4, ps_factors[5]);
      event.set(h_isr_up2, ps_factors[6]);
      event.set(h_isr_down2, ps_factors[7]);
    }
    else{
      event.set(h_fsr_upsqrt2, 1.);
      event.set(h_fsr_up2, 1.);
      event.set(h_fsr_up4, 1.);
      event.set(h_fsr_downsqrt2, 1.);
      event.set(h_fsr_down2, 1.);
      event.set(h_fsr_down4, 1.);
      event.set(h_isr_up2, 1.);
      event.set(h_isr_down2, 1.);
    }
  }
  else{
    event.set(h_fsr_upsqrt2, 1.);
    event.set(h_fsr_up2, 1.);
    event.set(h_fsr_up4, 1.);
    event.set(h_fsr_downsqrt2, 1.);
    event.set(h_fsr_down2, 1.);
    event.set(h_fsr_down4, 1.);
    event.set(h_isr_up2, 1.);
    event.set(h_isr_down2, 1.);
  }
  double gen_weight = event.weight;

  // select events by gen weight - Control stage for FSRup_4
  bool GenWeight300 = false;
  if(isTTbar && gen_weight>300){
    GenWeight300 = true;
    h_gen_weight_300->fill(event);
  }

  // FILL CONTROL PLOTS
  if(debug) cout << "\nFill PreSel01\n";
  if(passed_recsel){
    h_MTopJet_PreSel01->fill(event);
    h_Muon_PreSel01->fill(event);
    h_Elec_PreSel01->fill(event);
    h_Jets_PreSel01->fill(event);
    h_XCone_cor_PreSel01->fill(event);
  }
  if(isTTbar)h_GenParticles_SM->fill(event);

  // choose if tt bar sample width should be reweighted
  double factor_2w, factor_4w, factor_8w;
  if(isTTbar){
    factor_2w = width2_reweight->get_factor(event);
    factor_4w = width4_reweight->get_factor(event);
    factor_8w = width8_reweight->get_factor(event);
  }
  else{
    factor_2w = 1.0;
    factor_4w = 1.0;
    factor_8w = 1.0;
  }

  // if(reweight_ttbar) ttbar_reweight->process(event);
  h_weights01->fill(event);

  /** store weights for pdf variation *********************/
  if(debug) cout << "\nVary PDF\n";
  std::vector<double> pdf_weights;
  if(isTTbar){
    if(event.genInfo->systweights().size()){
      for(int i=0; i<100; i++){
        double pdf_weight = event.genInfo->systweights().at(i+9);
        pdf_weights.push_back(pdf_weight);
      }
    }
  }
  /** PU Reweighting *********************/
  if(debug) cout << "\nVary PU\n";
  PUreweight->process(event);
  event.set(h_pu_up, event.get(h_weight_pu_up)/event.get(h_weight_pu));
  event.set(h_pu_down, event.get(h_weight_pu_down)/event.get(h_weight_pu));
  h_weights02->fill(event);
  // FILL CONTROL PLOTS
  if(debug) cout << "\nFill PreSel02\n";
  if(passed_recsel){
    h_MTopJet_PreSel02->fill(event);
    h_Muon_PreSel02->fill(event);
    h_Elec_PreSel02->fill(event);
    h_Jets_PreSel02->fill(event);
    h_XCone_cor_PreSel02->fill(event);
  }

  /** muon SF *********************/
  if(channel_ == muon){
    if(debug) cout << "\nMuon SF\n";
    muo_tight_noniso_SF->process(event);
    muo_trigger_SF->process(event);
    // now store factors to get all the weights in unfolding
    event.set(h_muid_up, event.get(h_weight_muid_up)/event.get(h_weight_muid));
    event.set(h_muid_down, event.get(h_weight_muid_down)/event.get(h_weight_muid));
    event.set(h_mutr_up, event.get(h_weight_mutr_up)/event.get(h_weight_mutr));
    event.set(h_mutr_down, event.get(h_weight_mutr_down)/event.get(h_weight_mutr));
    event.set(h_elid_up, 1.);
    event.set(h_elid_down, 1.);
    event.set(h_elreco_up, 1.);
    event.set(h_elreco_down, 1.);
    event.set(h_eltr_up, 1.);
    event.set(h_eltr_down, 1.);
  }
  else{
    if(debug) cout << "\nElectron SF\n";
    ele_id_SF->process(event);
    ele_reco_SF->process(event);
    ele_trigger_SF->process(event);
    event.set(h_elreco_up, event.get(h_weight_elreco_up)/event.get(h_weight_elreco));
    event.set(h_elreco_down, event.get(h_weight_elreco_down)/event.get(h_weight_elreco));
    event.set(h_elid_up, event.get(h_weight_elid_up)/event.get(h_weight_elid));
    event.set(h_elid_down, event.get(h_weight_elid_down)/event.get(h_weight_elid));
    event.set(h_eltr_up, event.get(h_weight_eltr_up)/event.get(h_weight_eltr));
    event.set(h_eltr_down, event.get(h_weight_eltr_down)/event.get(h_weight_eltr));
    event.set(h_muid_up, 1.);
    event.set(h_muid_down, 1.);
    event.set(h_mutr_up, 1.);
    event.set(h_mutr_down, 1.);
  }


  h_weights03->fill(event);

  // FILL CONTROL PLOTS
  if(debug) cout << "\nFill PreSel03\n";
  if(passed_recsel){
    h_MTopJet_PreSel03->fill(event);
    h_Muon_PreSel03->fill(event);
    h_Elec_PreSel03->fill(event);
    h_Jets_PreSel03->fill(event);
    h_XCone_cor_PreSel03->fill(event);
    if(!event.isRealData) BTagEffHists->fill(event);
  }

  /** b-tagging *********************/
  // first do cvs reshape (instead of SF)
  if(debug) cout << "\nb tag reshape\n";
  CreateBTagJets->process(event); // First store a jet vector with only the highest btag
  BTagReshape->process(event);
  float btag_shift_up2 = 0;
  float btag_shift_down2 = 0;
  float btag_central = event.get(h_weight_btag);

  for(auto h: h_weight_btag_upvars){
    double sf_sys = event.get(h);
    if(sf_sys == 1.0) sf_sys = btag_central;
    btag_shift_up2 += pow(sf_sys-btag_central,2);
  }
  for(auto h: h_weight_btag_downvars){
    double sf_sys = event.get(h);
    if(sf_sys == 1.0) sf_sys = btag_central;
    btag_shift_down2 += pow(sf_sys-btag_central,2);
  }

  float btag_factor_up = (btag_central+sqrt(btag_shift_up2))/btag_central;
  float btag_factor_down = (btag_central-sqrt(btag_shift_up2))/btag_central;
  event.set(h_btag_up, btag_factor_up);
  event.set(h_btag_down, btag_factor_down);


  // cout << "=========== " << endl;
  // cout << "central = " <<  btag_central << endl;
  // cout << "shift up = " <<  sqrt(btag_shift_up2) << endl;
  // cout << "shift down = " <<  sqrt(btag_shift_down2) << endl;
  // cout << "factor up = " <<  btag_factor_up << endl;
  // cout << "factor down = " <<  btag_factor_down << endl;

  if(passed_recsel){
    h_MTopJet_PreSel03b->fill(event);
    h_Muon_PreSel03b->fill(event);
    h_Elec_PreSel03b->fill(event);
    h_Jets_PreSel03b->fill(event);
    h_XCone_cor_PreSel03b->fill(event);
  }
  bool passed_btag;
  bool passed_btag_medium;
  bool passed_btag_loose;
  int jetbtagN(0);
  int jetbtagN_medium(0);
  int jetbtagN_loose(0);
  for(const auto& j : *event.jets){
    if(DeepjetTight(j, event)) ++jetbtagN;
    if(DeepjetMedium(j, event)) ++jetbtagN_medium;
    if(DeepjetLoose(j, event)) ++jetbtagN_loose;
  }
  if(jetbtagN >= 1) passed_btag = true;
  else passed_btag = false;
  if(jetbtagN_medium >= 1) passed_btag_medium = true;
  else passed_btag_medium = false;
  if(jetbtagN_loose >= 1) passed_btag_loose = true;
  else passed_btag_loose = false;

  h_weights04->fill(event);

  /** calculate recweight now *********************/
  double rec_weight;
  if(gen_weight==0)rec_weight = 0;
  else rec_weight = (event.weight)/gen_weight;

  /**********************************/
  bool passed_presel_rec = (passed_recsel && passed_btag);
  // FILL CONTROL PLOTS
  if(passed_presel_rec){
    h_MTopJet_PreSel04->fill(event);
    h_Muon_PreSel04->fill(event);
    h_Elec_PreSel04->fill(event);
    h_Jets_PreSel04->fill(event);
    h_XCone_cor_PreSel04->fill(event);
  }
  /*************************** Events have to pass topjet pt > 400 & Mass_jet1 > Mass_jet2 *******************************************************/
  bool pass_measurement_rec;
  bool pass_pt350migration_rec;
  bool pass_ptmigration_rec;
  bool pass_massmigration_rec;
  bool pass_btagmigration_rec;
  bool pass_WJets_sel;
  bool pass_subptmigration_rec;
  bool pass_leptonptmigration_rec;
  if(debug) cout << "\nStart Main Selection\n";

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && passed_btag && subjet_quality->passes(event) && lepton_sel->passes(event)) pass_measurement_rec = true;
  else pass_measurement_rec = false;

  if(passed_recsel && !pt_sel->passes(event) && pt2_sel->passes(event) && pt350_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && passed_btag && subjet_quality->passes(event) && lepton_sel->passes(event) ) pass_pt350migration_rec = true;
  else pass_pt350migration_rec = false;

  if(passed_recsel && !pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && passed_btag && subjet_quality->passes(event) && lepton_sel->passes(event)) pass_ptmigration_rec = true;
  else pass_ptmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && !mass_sel->passes(event) && passed_btag && subjet_quality->passes(event) && lepton_sel->passes(event)) pass_massmigration_rec = true;
  else pass_massmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && !passed_btag && passed_btag_medium && subjet_quality->passes(event) && lepton_sel->passes(event)) pass_btagmigration_rec = true;
  else pass_btagmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event)  && eta_sel->passes(event) && !passed_btag_loose && subjet_quality->passes(event) && lepton_sel->passes(event)) pass_WJets_sel = true;
  else pass_WJets_sel = false;

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && passed_btag && !subjet_quality->passes(event) && lepton_sel->passes(event)) pass_subptmigration_rec = true;
  else pass_subptmigration_rec = false;

  if(passed_recsel && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && passed_btag && subjet_quality->passes(event) && !lepton_sel->passes(event)) pass_leptonptmigration_rec = true;
  else pass_leptonptmigration_rec = false;
  /*************************** Selection again on generator level (data events will not pass gen sel but will be stored if they pass rec sel)  *****/
  bool pass_measurement_gen;
  bool pass_pt350migration_gen;
  bool pass_massmigration_gen;
  bool pass_subptmigration_gen;
  bool pass_leptonptmigration_gen;
  if(isTTbar){
    if(!lepton_Nsel_gen->passes(event)) passed_gensel33 = false;

    if(passed_gensel33 && pt_gensel->passes(event) && pt2_gensel->passes(event) && mass_gensel->passes(event) && subjet_quality_gen->passes(event) && lepton_sel_gen->passes(event)) pass_measurement_gen = true;
    else pass_measurement_gen = false;

    if(passed_gensel33 && !pt_gensel->passes(event) && pt2_gensel->passes(event) && pt350_gensel->passes(event) && mass_gensel->passes(event) && subjet_quality_gen->passes(event) && lepton_sel_gen->passes(event)) pass_pt350migration_gen = true;
    else pass_pt350migration_gen = false;

    if(passed_gensel33 && pt_gensel->passes(event) && pt2_gensel->passes(event) && !mass_gensel->passes(event) && subjet_quality_gen->passes(event) && lepton_sel_gen->passes(event)) pass_massmigration_gen = true;
    else pass_massmigration_gen = false;

    if(passed_gensel33 && pt_gensel->passes(event) && pt2_gensel->passes(event) && mass_gensel->passes(event) && !subjet_quality_gen->passes(event) && lepton_sel_gen->passes(event)) pass_subptmigration_gen = true;
    else pass_subptmigration_gen = false;

    if(passed_gensel33 && pt_gensel->passes(event) && pt2_gensel->passes(event) && mass_gensel->passes(event) && subjet_quality_gen->passes(event) && !lepton_sel_gen->passes(event)) pass_leptonptmigration_gen = true;
    else pass_leptonptmigration_gen = false;
  }
  else if(!isTTbar){
    pass_measurement_gen = false;
    pass_pt350migration_gen = false;
    pass_massmigration_gen = false;
    pass_subptmigration_gen = false;
    pass_leptonptmigration_gen = false;
  }

  /*************************************************************************************************************************************************/
  /* In two steps because Event with gen_weight > 300 should be deleted completly but still be shown if it passes the GenSelection *****************/
  if(GenWeight300 && pass_measurement_gen) h_gen_weight_gensel_300->fill(event);
  if(GenWeight300) return false;

  /*************************** Pile Up bools  ******************************************************************************************************/
  bool lowPU = (event.pvs->size() <= 10);
  bool midPU = (event.pvs->size() > 10 && event.pvs->size() <= 20);
  bool highPU = (event.pvs->size() > 20);
  bool is_matched_sub = false;
  bool is_matched_fat = false;

  /*************************** fill hists with reco sel applied ************************************************************************************/
  if(debug) cout << "\nFill Main Hists\n";
  h_comparison_topjet_xcone->fill(event);

  // hists to see events that are generated in measurement phase-space, but reconstructed outside
  if(debug) cout << "\npass_measurement_gen\n";
  if(pass_measurement_gen){
    if(pass_pt350migration_rec) h_XCone_cor_migration_pt350->fill(event);
    if(pass_ptmigration_rec)    h_XCone_cor_migration_pt->fill(event);
    if(pass_massmigration_rec)  h_XCone_cor_migration_mass->fill(event);
    if(pass_btagmigration_rec)  h_XCone_cor_migration_btag->fill(event);
  }
  // see all events reconstructed outside measurement phase-space
  if(debug) cout << "\nSingle passes\n";
  if(pass_ptmigration_rec) h_XCone_cor_noptcut->fill(event);
  if(pass_pt350migration_rec) h_XCone_cor_pt350->fill(event);
  if(pass_btagmigration_rec)  h_XCone_cor_Sel_nobtag->fill(event);
  if(pass_pt350migration_rec) h_XCone_cor_Sel_pt350->fill(event);
  if(pass_massmigration_rec)  h_XCone_cor_Sel_noMass->fill(event);
  if(pass_subptmigration_rec){
    h_XCone_cor_Sel_ptsub->fill(event);
    h_XCone_cor_subjets_Sel_ptsub->fill(event);
  }
  if(debug) cout << "\nSinge passes - substep\n";
  if(pass_measurement_gen)    h_XCone_GEN_Sel_measurement->fill(event);
  if(pass_pt350migration_gen) h_XCone_GEN_Sel_pt350->fill(event);
  if(pass_massmigration_gen)  h_XCone_GEN_Sel_noMass->fill(event);
  if(pass_subptmigration_gen) h_XCone_GEN_Sel_ptsub->fill(event);
  if(pass_WJets_sel && isTTbar){
    h_CorrectionHists_WJets->fill(event);
    h_RecGenHists_subjets_WJets->fill(event);
    h_RecGenHists_subjets_noJEC_WJets->fill(event);
  }
  // fill resolution hists here without the subjet pt cut
  if(debug) cout << "\npass_measurement_gen or pass_subptmigration_gen\n";
  if(pass_measurement_gen || pass_subptmigration_gen){
    if(lowPU){
      if(isTTbar)h_RecGenHists_lowPU->fill(event);
      if(isTTbar)h_RecGenHists_lowPU_noJEC->fill(event);
      if(isTTbar)h_RecGenHists_subjets_lowPU->fill(event);
      if(isTTbar)h_RecGenHists_subjets_noJEC_lowPU->fill(event);
    }
    if(midPU){
      if(isTTbar)h_RecGenHists_medPU->fill(event);
      if(isTTbar)h_RecGenHists_medPU_noJEC->fill(event);
      if(isTTbar)h_RecGenHists_subjets_medPU->fill(event);
      if(isTTbar)h_RecGenHists_subjets_noJEC_medPU->fill(event);
    }
    if(highPU){
      if(isTTbar)h_RecGenHists_highPU->fill(event);
      if(isTTbar)h_RecGenHists_highPU_noJEC->fill(event);
      if(isTTbar)h_RecGenHists_subjets_highPU->fill(event);
      if(isTTbar)h_RecGenHists_subjets_noJEC_highPU->fill(event);
    }
    h_RecGenHists_RecOnly->fill(event);
    h_RecGenHists_RecOnly_noJEC->fill(event);
    h_RecGenHists_RecOnly_corr->fill(event);
    if(matched_fat->passes(event)) h_RecGenHists_subjets_matched->fill(event);
    h_RecGenHists_subjets->fill(event);
    h_RecGenHists_subjets_noJEC->fill(event);
    h_RecGenHists_subjets_corrected->fill(event);
    if(isTTbar) h_CorrectionHists->fill(event);
    if(isTTbar) h_CorrectionHists_after->fill(event);
    h_RecGenHists_ak4->fill(event);
    h_RecGenHists_ak4_noJEC->fill(event);
  }
  if(isMC){
    std::vector<GenTopJet> gen_fatjets = event.get(h_genfatjets);
    double pt = gen_fatjets.at(0).v4().Pt();
    double mjet = gen_fatjets.at(0).v4().M();
    if(pt > 400 && mjet > 120){
      // h_XCone_GEN_ST->fill(event);
      h_RecGenHistsST_subjets->fill(event);
      h_RecGenHistsST_subjets_noJEC->fill(event);
      h_RecGenHistsST_subjets_corrected->fill(event);
    }
  }

  event.set(h_ak8tau, -1.); // set default value
  event.set(h_ak8mass, -1.); // set default value

  if(debug) cout << "pass_measurement_rec\n";
  if(pass_measurement_rec){
    h_PDFHists->fill(event);

    h_XCone_raw->fill(event);
    h_XCone_jec->fill(event);
    h_XCone_cor->fill(event);
    h_XCone_puppi->fill(event);
    if(isMC){
      // cout << "isMC\n";
      h_XCone_JMS->fill_mass(event, mass_jms);
    }
    else{
      // cout << "isNotMC\n";
      h_XCone_JMS->fill_mass(event, mass_rec);
    }
    h_leptonictop->fill(event);
    h_XCone_raw_subjets->fill(event);
    h_XCone_jec_subjets->fill(event);
    h_XCone_cor_subjets->fill(event);

    h_comparison_topjet_xcone_pass_rec->fill(event);
    h_comparison_topjet_xcone_pass_rec_masscut_120->fill(event);
    h_comparison_topjet_xcone_pass_rec_masscut_130->fill(event);
    h_comparison_topjet_xcone_pass_rec_masscut_140->fill(event);
    h_comparison_topjet_xcone_pass_rec_masscut_150->fill(event);
    double tau32 = h_comparison_topjet_xcone_pass_rec_masscut_140->get_tau32();
    double mass = h_comparison_topjet_xcone_pass_rec_masscut_140->get_mass();
    event.set(h_ak8tau, tau32);
    event.set(h_ak8mass, mass);

    if(isTTbar && scale_ttbar) event.weight *= SF_tt;

    h_XCone_raw_SF->fill(event);
    h_XCone_jec_SF->fill(event);
    h_XCone_cor_SF->fill(event);

    h_leptonictop_SF->fill(event);

    if(!pt450_sel->passes(event)) h_XCone_cor_SF_pt400->fill(event);
    if(pt450_sel->passes(event) && !pt530_sel->passes(event)) h_XCone_cor_SF_pt450->fill(event);
    if(pt530_sel->passes(event)) h_XCone_cor_SF_pt530->fill(event);

    h_XCone_raw_subjets_SF->fill(event);
    h_XCone_jec_subjets_SF->fill(event);
    h_XCone_cor_subjets_SF->fill(event);

    h_MTopJet->fill(event);
    h_Muon->fill(event);
    h_Elec->fill(event);
    h_Jets->fill(event);

    if(lowPU){
      h_XCone_cor_PUlow_subjets->fill(event);
      h_XCone_cor_PUlow->fill(event);
    }
    if(midPU){
      h_XCone_cor_PUmid_subjets->fill(event);
      h_XCone_cor_PUmid->fill(event);
    }
    if(highPU){
      h_XCone_cor_PUhigh_subjets->fill(event);
      h_XCone_cor_PUhigh->fill(event);
    }

    if(isTTbar){
      is_matched_sub = matched_sub->passes(event);
      is_matched_fat = matched_fat->passes(event);
      if(is_matched_sub) h_XCone_cor_m->fill(event);
      else h_XCone_cor_u->fill(event);
      if(is_matched_fat) h_XCone_cor_m_fat->fill(event);
      else h_XCone_cor_u_fat->fill(event);
    }

    if(isTTbar){
      h_XCone_GEN_RecOnly->fill(event);
      if(isTTbar) h_GenParticles_RecOnly->fill(event);
    }
  }

  /*************************** fill hists with gen sel applied *************************************************************************************/
  if(debug) cout << "pass_measurement_gen\n";
  if(pass_measurement_gen){

    // The next 3 Histograms are using the gen_weight of the event.
    // To pass the weight properly through the event.weight, this method is used.
    // I am open for a more elegant way !!
    // double original_weight = event.weight;
    // event.weight = gen_weight;
    // h_gen_weights_pass_gen->fill(event);
    // if(mass_gen33 < 150 && mass_gen33 > 140) h_gen_weights_massbin_145->fill(event); // 2017
    // if(mass_gen33 < 280 && mass_gen33 > 270) h_gen_weights_massbin_275->fill(event); // 2018
    // event.weight = original_weight;

    if(isTTbar) h_comparison_topjet_xcone_pass_gen->fill(event);
    if(isTTbar) h_GenParticles_GenOnly->fill(event);

    h_XCone_GEN_GenOnly->fill(event);
    if(matched_sub_GEN->passes(event)) h_XCone_GEN_GenOnly_matched->fill(event);
    else                               h_XCone_GEN_GenOnly_unmatched->fill(event);
    if(matched_fat_GEN->passes(event)) h_XCone_GEN_GenOnly_matched_fat->fill(event);
    else                               h_XCone_GEN_GenOnly_unmatched_fat->fill(event);
    h_RecGenHists_GenOnly->fill(event);
    // fill this only if event passes at least one reco sideband or measurement phase space
    if(pass_measurement_rec || pass_pt350migration_rec || pass_massmigration_rec || pass_btagmigration_rec || pass_subptmigration_rec || pass_leptonptmigration_rec){
      h_RecGenHists_GenSelRecInfo->fill(event);
      if(lowPU)       h_RecGenHists_GenSelRecInfo_lowPU->fill(event);
      else if(midPU)  h_RecGenHists_GenSelRecInfo_midPU->fill(event);
      else if(highPU) h_RecGenHists_GenSelRecInfo_highPU->fill(event);
      if(matched_sub_GEN->passes(event)){
        h_RecGenHists_GenSelRecInfo_matched->fill(event);
        if(lowPU)       h_RecGenHists_GenSelRecInfo_matched_lowPU->fill(event);
        else if(midPU)  h_RecGenHists_GenSelRecInfo_matched_midPU->fill(event);
        else if(highPU) h_RecGenHists_GenSelRecInfo_matched_highPU->fill(event);
      }
    }
  }


  /*************************** fill hists with reco+gen selection applied **************************************************************************/
  if(debug) cout << "pass_measurement_rec && pass_measurement_gen\n";
  if(pass_measurement_rec && pass_measurement_gen){
    if(isTTbar) h_comparison_topjet_xcone_pass_genrec->fill(event);
    if(isTTbar) h_GenParticles_Both->fill(event);
    h_XCone_GEN_Both->fill(event);
    h_RecGenHists_Both->fill(event);
    h_RecGenHists_Both_corr->fill(event);
  }

  /*************************** also store factor from ttbar reweighting ****************************************************************************/
  double ttfactor = 1.0;
  if(isTTbar){
    double weight_before = event.weight;
    ttbar_reweight->process(event);
    double weight_after = event.weight;
    ttfactor = weight_after / weight_before;
  }

  /*************************** write bools for passing selections **********************************************************************************/
  if(debug) cout << "Event Set\n";

  event.set(h_ttbar, isTTbar);
  event.set(h_matched, is_matched_sub);
  event.set(h_ttbar_SF, SF_tt);
  event.set(h_genweight, gen_weight);
  event.set(h_recweight, rec_weight);
  event.set(h_genweight_ttfactor, ttfactor);

  event.set(h_measure_rec, pass_measurement_rec);
  event.set(h_measure_gen, pass_measurement_gen);
  event.set(h_pt350_rec, pass_pt350migration_rec);
  event.set(h_pt350_gen, pass_pt350migration_gen);
  event.set(h_nomass_rec, pass_massmigration_rec);
  event.set(h_nomass_gen, pass_massmigration_gen);
  event.set(h_nobtag_rec, pass_btagmigration_rec);
  event.set(h_ptsub_rec, pass_subptmigration_rec);
  event.set(h_ptsub_gen, pass_subptmigration_gen);
  event.set(h_ptlepton_rec, pass_leptonptmigration_rec);
  event.set(h_ptlepton_gen, pass_leptonptmigration_gen);

  event.set(h_factor_2width, factor_2w);
  event.set(h_factor_4width, factor_4w);
  event.set(h_factor_8width, factor_8w);

  event.set(h_pdf_weights, pdf_weights);
  /*************************** only store events that survive one of the selections (use looser pt cut) ****************************************************************/
  if(debug) cout << "Event Set\n";
  bool in_migrationmatrix = ( pass_measurement_rec ||
    pass_measurement_gen ||
    pass_pt350migration_rec ||
    pass_pt350migration_gen ||
    pass_massmigration_rec ||
    pass_massmigration_gen ||
    pass_btagmigration_rec ||
    pass_subptmigration_rec ||
    pass_subptmigration_gen ||
    pass_leptonptmigration_rec ||
    pass_leptonptmigration_gen   );
    if(debug) cout << "Event Set\n";

    if(!in_migrationmatrix) return false;







    else return true;

  }

  UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
