#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/TopPtReweight.h>
#include "UHH2/common/include/YearRunSwitchers.h"

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CutHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/NonClosureUncertainty.h>
#include <UHH2/MTopJet/include/RemoveLepton.h>

#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/CorrectionFactor.h>


/*
*******************************************************************
**************** TO DO ********************************************
*******************************************************************
- b tagging SF
*******************************************************************
*******************************************************************
*/

class MTopJetSelectionModule : public ModuleBASE {

public:
  explicit MTopJetSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners & Correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;


  std::unique_ptr<JetCleaner> jet_cleaner1;
  std::unique_ptr<JetCleaner> jet_cleaner2;
  std::unique_ptr<JetCorrections_xcone> JetCorrections;
  std::unique_ptr<JER_Smearer_xcone> JERSmearing;
  std::unique_ptr<uhh2::AnalysisModule> Correction;
  std::unique_ptr<uhh2::AnalysisModule> NonClosureSYS;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> remove_lepton_rec;
  std::unique_ptr<uhh2::AnalysisModule> remove_lepton_gen33;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_pupppi;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_noJEC;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_corrected;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen33;
  std::unique_ptr<uhh2::AnalysisModule> copy_jet;
  std::unique_ptr<uhh2::Selection> trigger_mu_A;
  std::unique_ptr<uhh2::Selection> trigger_mu_B;
  std::unique_ptr<uhh2::Selection> trigger_el_A;
  std::unique_ptr<uhh2::Selection> trigger_el_B;
  std::unique_ptr<uhh2::Selection> trigger_el_C;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> elec_sel_120;
  std::unique_ptr<uhh2::Selection> elec_sel_triggerA;
  std::unique_ptr<uhh2::Selection> elec_sel_triggerB;
  std::unique_ptr<uhh2::Selection> elec_etaveto;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;

  Event::Handle<bool>h_gensel;
  Event::Handle<bool>h_recsel;
  Event::Handle<bool>h_gensel_2;
  Event::Handle<bool>h_recsel_2;



  Event::Handle<std::vector<TopJet>>h_fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen33fatjets;

  // just for testing
  std::unique_ptr<TopPtReweight> ttbar_reweight;


  // write output
  std::unique_ptr<uhh2::AnalysisModule> output;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_PreSel_event,  h_PreSel_elec, h_PreSel_muon, h_PreSel_jets, h_PreSel_lumi;
  std::unique_ptr<Hists> h_Trigger_event,  h_Trigger_elec, h_Trigger_muon, h_Trigger_jets, h_Trigger_lumi;
  std::unique_ptr<Hists> h_Cleaner_event,  h_Cleaner_elec, h_Cleaner_muon, h_Cleaner_jets, h_Cleaner_lumi;
  std::unique_ptr<Hists> h_Lepton_event,  h_Lepton_elec, h_Lepton_muon, h_Lepton_jets, h_Lepton_lumi;
  std::unique_ptr<Hists> h_TwoD_event,  h_TwoD_elec, h_TwoD_muon, h_TwoD_jets, h_TwoD_lumi;
  std::unique_ptr<Hists> h_Jet_event,  h_Jet_elec, h_Jet_muon, h_Jet_jets, h_Jet_lumi;
  std::unique_ptr<Hists> h_bTag_event,  h_bTag_elec, h_bTag_muon, h_bTag_jets, h_bTag_lumi;
  std::unique_ptr<Hists> h_Side_event,  h_Side_elec, h_Side_muon, h_Side_jets, h_Side_lumi;
  std::unique_ptr<Hists> h_HTlep_event,  h_HTlep_elec, h_HTlep_muon, h_HTlep_jets, h_HTlep_lumi;
  std::unique_ptr<Hists> h_MET_event,  h_MET_elec, h_MET_muon, h_MET_jets, h_MET_lumi;
  std::unique_ptr<Hists> h_ttbar_reweight_event, h_ttbar_reweight_elec, h_ttbar_reweight_muon, h_ttbar_reweight_jets, h_ttbar_reweight_lumi;
  std::unique_ptr<Hists> h_cuts_all, h_cuts_all_but_2D, h_cuts_all_but_met, h_cuts_all_but_btag;


  bool isMC; //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part
  bool recsel_only;
  bool isElectronStream;
  bool isPhotonStream;

  string corvar ="nominal";
  string nonclosureSYS = "false";
  bool do_nonClosure = false;

  JetId btag_tight;

  Year year;

};

/*
███    ███  ██████  ██████  ██    ██ ██      ███████
████  ████ ██    ██ ██   ██ ██    ██ ██      ██
██ ████ ██ ██    ██ ██   ██ ██    ██ ██      █████
██  ██  ██ ██    ██ ██   ██ ██    ██ ██      ██
██      ██  ██████  ██████   ██████  ███████ ███████
*/



MTopJetSelectionModule::MTopJetSelectionModule(uhh2::Context& ctx){

  recsel_only = false;
  year = extract_year(ctx); // Ask for the year of Event

  //// CONFIGURATION
  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700_2016v3"  ||
  ctx.get("dataset_version") == "TTbar_Mtt0700to1000_2016v3"  ||
  ctx.get("dataset_version") == "TTbar_Mtt1000toInft_2016v3"  ||
  ctx.get("dataset_version") == "TTbar_mtop1695_2016v3"       ||
  ctx.get("dataset_version") == "TTbar_mtop1715_2016v3"       ||
  ctx.get("dataset_version") == "TTbar_mtop1735_2016v3"       ||
  ctx.get("dataset_version") == "TTbar_mtop1755_2016v3"       ||
  ctx.get("dataset_version") == "TTbar_fsrup_2016v3"          ||
  ctx.get("dataset_version") == "TTbar_fsrdown_2016v3"        ||
  ctx.get("dataset_version") == "TTbar_isrup_2016v3"          ||
  ctx.get("dataset_version") == "TTbar_isrdown_2016v3"        ||
  ctx.get("dataset_version") == "TTbar_hdampup_2016v3"        ||
  ctx.get("dataset_version") == "TTbar_Mtt0000to0700_2L2Nu_2017v2"          ||
  ctx.get("dataset_version") == "TTbar_Mtt0000to0700_SemiLep_2017v2"        ||
  ctx.get("dataset_version") == "TTbar_Mtt0000to0700_Hadronic_2017v2"        ||
  ctx.get("dataset_version") == "TTbar_hdampdown_2016v3"      ) isTTbar = true;
  else  isTTbar = false;

  if(ctx.get("dataset_version") == "SingleElecB" ||
  ctx.get("dataset_version") == "SingleElecC" ||
  ctx.get("dataset_version") == "SingleElecD" ||
  ctx.get("dataset_version") == "SingleElecE" ||
  ctx.get("dataset_version") == "SingleElecF" ||
  ctx.get("dataset_version") == "SingleElecG" ||
  ctx.get("dataset_version") == "SingleElecHver2" ||
  ctx.get("dataset_version") == "SingleElecHver3") isElectronStream = true;
  else isElectronStream = false;

  if(ctx.get("dataset_version") == "SinglePhotonB" ||
  ctx.get("dataset_version") == "SinglePhotonC" ||
  ctx.get("dataset_version") == "SinglePhotonD" ||
  ctx.get("dataset_version") == "SinglePhotonE" ||
  ctx.get("dataset_version") == "SinglePhotonF" ||
  ctx.get("dataset_version") == "SinglePhotonG" ||
  ctx.get("dataset_version") == "SinglePhotonHver2" ||
  ctx.get("dataset_version") == "SinglePhotonHver3") isPhotonStream = true;
  else isPhotonStream = false;

  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  if(isTTbar) ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  if(isTTbar) h_gensel = ctx.get_handle<bool>("passed_gensel");
  h_recsel = ctx.get_handle<bool>("passed_recsel");
  h_gensel_2 = ctx.declare_event_output<bool>("passed_gensel_2");
  h_recsel_2 = ctx.declare_event_output<bool>("passed_recsel_2");
  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  h_gen33fatjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  isMC = (ctx.get("dataset_type") == "MC");

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {

    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";

    throw std::runtime_error(log);
  }

  // remove lepton
  remove_lepton_rec.reset(new RemoveLepton(ctx, "xconeCHS"));
  if(isMC){
    remove_lepton_gen33.reset(new RemoveLeptonGen(ctx, "genXCone33TopJets"));
  }

  // combine XCone
  jetprod_reco.reset(new CombineXCone33(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "xconeCHS"));
  jetprod_reco_noJEC.reset(new CombineXCone33(ctx, "XCone33_had_Combined_noJEC", "XCone33_lep_Combined_noJEC", "xconeCHS" ));
  jetprod_reco_corrected.reset(new CombineXCone33(ctx, "XCone33_had_Combined_Corrected", "XCone33_lep_Combined_Corrected", "xconeCHS_Corrected"));
  copy_jet.reset(new CopyJets(ctx, "xconeCHS", "xconeCHS_noJEC"));

  if(isMC){
    jetprod_gen33.reset(new CombineXCone33_gen(ctx));
  }
  ////

  // write output
  output.reset(new WriteOutput(ctx));
  ////
  // just for testing
  // ttbar_reweight.reset(new TopPtReweight(ctx,0.159,-0.00141,"","weight_ttbar",true,0.9910819)); // 8 TeV
  ttbar_reweight.reset(new TopPtReweight(ctx,0.0615,-0.0005,"","weight_ttbar",true)); // 13 TeV
  // double jecsysfactor = 1.0;

  corvar = ctx.get("JetCorrection_direction","nominal");
  nonclosureSYS = ctx.get("NonClosureUncertainty","false");
  if(nonclosureSYS == "true") do_nonClosure = true;
  // correct subjets (JEC + additional correction)
  JetCorrections.reset(new JetCorrections_xcone());
  JetCorrections->init(ctx, "xconeCHS");
  // smear jets after Correction
  JERSmearing.reset(new JER_Smearer_xcone());
  JERSmearing->init(ctx, "xconeCHS", "genXCone33TopJets", "sub");
  Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, false));
  NonClosureSYS.reset(new NonClosureUncertainty(ctx));
  //// EVENT SELECTION

  // define IDs
  MuonId muid = AndId<Muon>(MuonID(Muon::Tight), PtEtaCut(55., 2.4));
  // this is only used for cleaner and electron veto
  ElectronId eleid_noiso55  = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Summer16_tight_noIso);
  // this is used to decide which ele trigger is used
  ElectronId eleid_noiso120 = AndId<Electron>(PtEtaSCCut(120., 2.4), ElectronID_Summer16_tight_noIso);
  // this is used in combination with iso trigger
  ElectronId eleid_iso55    = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Summer16_tight);
  // jet ids
  JetId jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));
  ////

  // define Trigger
  trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
  trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");
  trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele27_WPTight_Gsf_v*");
  trigger_el_B = uhh2::make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
  trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon175_v*");

  /*Only select event with exacly 1 muon or electron */
  if(channel_ == elec){
    muon_sel.reset(new NMuonSelection(0, 0, muid));
    elec_sel.reset(new NElectronSelection(1, 1, eleid_noiso55));
  }
  else if (channel_ == muon){
    muon_sel.reset(new NMuonSelection(1, 1, muid));
    elec_sel.reset(new NElectronSelection(0, 0, eleid_noiso55));
  }
  elec_etaveto.reset(new ElectronEtaVeto(1.44, 1.57));
  elec_sel_triggerA.reset(new NElectronSelection(1, 1, eleid_iso55));
  elec_sel_120.reset(new NElectronSelection(1, 1, eleid_noiso120));

  double metcut = 0;
  if(channel_ == muon)      metcut = 50;
  else if(channel_ == elec) metcut = 50;
  met_sel.reset(new METCut  (metcut , uhh2::infinity));
  twodcut_sel.reset(new TwoDCut1(0.4, 40));
  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));

  //// Obj Cleaning

  common.reset(new CommonModules());
  common->set_HTjetid(jetid_cleaner);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  common->disable_mcpileupreweight(); // do this in PostSel
  // common->disable_jersmear();
  // common->disable_jec();
  common->init(ctx);

  muoSR_cleaner.reset(new     MuonCleaner(muid));
  eleSR_cleaner.reset(new ElectronCleaner(eleid_noiso55));

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //
  btag_tight = DeepJetBTag(DeepJetBTag::WP_TIGHT);
  BTagEffHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",btag_tight));

  //// set up Hists classes:
  h_cuts_all.reset(new CutHists(ctx, "h_cuts_all"));
  h_cuts_all_but_2D.reset(new CutHists(ctx, "h_cuts_all_but_2D"));
  h_cuts_all_but_met.reset(new CutHists(ctx, "h_cuts_all_but_met"));
  h_cuts_all_but_btag.reset(new CutHists(ctx, "h_cuts_all_but_btag"));

  h_PreSel_event.reset(new MTopJetHists(ctx, "00_PreSel_Event"));
  h_PreSel_elec.reset(new ElectronHists(ctx, "00_PreSel_Elec"));
  h_PreSel_muon.reset(new MuonHists(ctx, "00_PreSel_Muon"));
  h_PreSel_jets.reset(new JetHists(ctx, "00_PreSel_Jets"));
  h_PreSel_lumi.reset(new LuminosityHists(ctx, "00_PreSel_lumi"));

  h_Cleaner_event.reset(new MTopJetHists(ctx, "01_Cleaner_Event"));
  h_Cleaner_elec.reset(new ElectronHists(ctx, "01_Cleaner_Elec"));
  h_Cleaner_muon.reset(new MuonHists(ctx, "01_Cleaner_Muon"));
  h_Cleaner_jets.reset(new JetHists(ctx, "01_Cleaner_Jets"));
  h_Cleaner_lumi.reset(new LuminosityHists(ctx, "01_Cleaner_lumi"));

  h_Trigger_event.reset(new MTopJetHists(ctx, "02_Trigger_Event"));
  h_Trigger_elec.reset(new ElectronHists(ctx, "02_Trigger_Elec"));
  h_Trigger_muon.reset(new MuonHists(ctx, "02_Trigger_Muon"));
  h_Trigger_jets.reset(new JetHists(ctx, "02_Trigger_Jets"));
  h_Trigger_lumi.reset(new LuminosityHists(ctx, "02_Trigger_lumi"));

  h_Lepton_event.reset(new MTopJetHists(ctx, "03_Lepton_Event"));
  h_Lepton_elec.reset(new ElectronHists(ctx, "03_Lepton_Elec"));
  h_Lepton_muon.reset(new MuonHists(ctx, "03_Lepton_Muon"));
  h_Lepton_jets.reset(new JetHists(ctx, "03_Lepton_Jets"));
  h_Lepton_lumi.reset(new LuminosityHists(ctx, "03_Lepton_lumi"));

  h_Jet_event.reset(new MTopJetHists(ctx, "04_Jet_Event"));
  h_Jet_elec.reset(new ElectronHists(ctx, "04_Jet_Elec"));
  h_Jet_muon.reset(new MuonHists(ctx, "04_Jet_Muon"));
  h_Jet_jets.reset(new JetHists(ctx, "04_Jet_Jets"));
  h_Jet_lumi.reset(new LuminosityHists(ctx, "04_Jet_lumi"));

  h_TwoD_event.reset(new MTopJetHists(ctx, "05_TwoD_Event"));
  h_TwoD_elec.reset(new ElectronHists(ctx, "05_TwoD_Elec"));
  h_TwoD_muon.reset(new MuonHists(ctx, "05_TwoD_Muon"));
  h_TwoD_jets.reset(new JetHists(ctx, "05_TwoD_Jets"));
  h_TwoD_lumi.reset(new LuminosityHists(ctx, "05_TwoD_lumi"));

  h_MET_event.reset(new MTopJetHists(ctx, "06_MET_Event"));
  h_MET_elec.reset(new ElectronHists(ctx, "06_MET_Elec"));
  h_MET_muon.reset(new MuonHists(ctx, "06_MET_Muon"));
  h_MET_jets.reset(new JetHists(ctx, "06_MET_Jets"));
  h_MET_lumi.reset(new LuminosityHists(ctx, "06_MET_lumi"));

  h_HTlep_event.reset(new MTopJetHists(ctx, "07_HTlep_Event"));
  h_HTlep_elec.reset(new ElectronHists(ctx, "07_HTlep_Elec"));
  h_HTlep_muon.reset(new MuonHists(ctx, "07_HTlep_Muon"));
  h_HTlep_jets.reset(new JetHists(ctx, "07_HTlep_Jets"));
  h_HTlep_lumi.reset(new LuminosityHists(ctx, "07_HTlep_lumi"));

  h_bTag_event.reset(new MTopJetHists(ctx, "08_bTag_Event"));
  h_bTag_elec.reset(new ElectronHists(ctx, "08_bTag_Elec"));
  h_bTag_muon.reset(new MuonHists(ctx, "08_bTag_Muon"));
  h_bTag_jets.reset(new JetHists(ctx, "08_bTag_jets"));
  h_bTag_lumi.reset(new LuminosityHists(ctx, "08_bTag_lumi"));

  h_Side_event.reset(new MTopJetHists(ctx, "Side_Event"));
  h_Side_elec.reset(new ElectronHists(ctx, "Side_Elec"));
  h_Side_muon.reset(new MuonHists(ctx, "Side_Muon"));
  h_Side_jets.reset(new JetHists(ctx, "Side_Jets"));
  h_Side_lumi.reset(new LuminosityHists(ctx, "Side_lumi"));

  h_ttbar_reweight_event.reset(new MTopJetHists(ctx, "ttbar_reweight_Event"));
  h_ttbar_reweight_elec.reset(new ElectronHists(ctx, "ttbar_reweight_Elec"));
  h_ttbar_reweight_muon.reset(new MuonHists(ctx, "ttbar_reweight_Muon"));
  h_ttbar_reweight_jets.reset(new JetHists(ctx, "ttbar_reweight_Jets"));
  h_ttbar_reweight_lumi.reset(new LuminosityHists(ctx, "ttbar_reweight_lumi"));
  //

}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool MTopJetSelectionModule::process(uhh2::Event& event){

  bool passed_gensel;
  if(isTTbar) passed_gensel = event.get(h_gensel);
  else passed_gensel = false;

  bool passed_recsel;
  if(isTTbar) passed_recsel = event.get(h_recsel);
  else passed_recsel = true;


  // fill ttbargen class
  if(isTTbar) ttgenprod->process(event);
  ////

  /* *********** Lepton Cleaner and Selection *********** */
  if(passed_recsel){
    muoSR_cleaner->process(event);
    if(event.muons->size() > 0) sort_by_pt<Muon>(*event.muons);

    eleSR_cleaner->process(event);
    if(event.electrons->size() > 0) sort_by_pt<Electron>(*event.electrons);
    const bool pass_lep1 = ((event.muons->size() > 0) || (event.electrons->size() > 0));
    if(!pass_lep1) passed_recsel = false;
  }
  ////

  /* *********** Cleaner ********** */
  if(!common->process(event)) return false;
  if(passed_recsel){
    jet_cleaner1->process(event);
    sort_by_pt<Jet>(*event.jets);
  }

  /* *********** at least 1 good primary vertex *********** */
  if(passed_recsel){
    if(!pv_sel->passes(event)) passed_recsel = false;
  }

  /* *********** scale factor muon *********** */
  // muo_tight_noniso_SF->process(event);
  if(passed_recsel){
    h_Cleaner_event->fill(event);
    h_Cleaner_elec->fill(event);
    h_Cleaner_muon->fill(event);
    h_Cleaner_jets->fill(event);
    h_Cleaner_lumi->fill(event);
  }

  /* *********** Trigger *********** */
  // for DATA until run 274954 -> use only Trigger A
  // for MC and DATA from 274954 -> use "A || B"
  // HLT_TkMu50_v is not availible for 2017&18. Alternative still needs to be inlcuded !!!!!!!!!!!!!!!!
  bool elec_is_isolated = false;
  if(channel_ == muon){
    if(passed_recsel){
      if(!isMC) {
        if(year == Year::is2016v3){
          if(event.run < 274954){
            if(!trigger_mu_A->passes(event)) passed_recsel = false;
          }
          else{
            if( !(trigger_mu_A->passes(event) || trigger_mu_B->passes(event)) ) passed_recsel = false;
          }
        }
        if(year == Year::is2017v2){
          if(!trigger_mu_A->passes(event)) passed_recsel = false;
        }
        if(year == Year::is2018){
          if(!trigger_mu_A->passes(event)) passed_recsel = false;
        }
      }
    }
  }

  if(year == Year::is2017v2){
    if(channel_ == muon){
      if(passed_recsel){
        if( !isMC)
        if(!trigger_mu_A->passes(event)) passed_recsel = false;
      }
    }
  }
  // only use triggerA and isolation if elec pt < 120
  // for pt > 120 use triggerB || triggerC
  else if(channel_ == elec){
    if(passed_recsel){
      if(!elec_sel_120->passes(event)){
        if(isPhotonStream) passed_recsel = false;
        if(!trigger_el_A->passes(event))      passed_recsel = false;
        if(!elec_sel_triggerA->passes(event)) passed_recsel = false;
        if(passed_recsel) elec_is_isolated = false;
      }
      else{
        if(isMC){
          if( !(trigger_el_B->passes(event) || trigger_el_C->passes(event)) ) passed_recsel = false;
        }
        else if(isElectronStream){
          if(!trigger_el_B->passes(event))  passed_recsel = false;
        }
        else if(isPhotonStream){
          if(trigger_el_B->passes(event))  passed_recsel = false;
          if(!trigger_el_C->passes(event))  passed_recsel = false;
        }
      }
    }
  }

  if(passed_recsel){
    h_Trigger_event->fill(event);
    h_Trigger_elec->fill(event);
    h_Trigger_muon->fill(event);
    h_Trigger_jets->fill(event);
    h_Trigger_lumi->fill(event);
  }
  /* *********** lEPTON Selection *********** */
  if(passed_recsel){
    bool pass_lepsel = (muon_sel->passes(event) && elec_sel->passes(event));
    // in elec channel use additional eta cut
    if(channel_ == elec){
      if(!elec_etaveto->passes(event)) pass_lepsel = false;
    }
    ////
    if(!pass_lepsel) passed_recsel = false;
    else{
      h_Lepton_event->fill(event);
      h_Lepton_elec->fill(event);
      h_Lepton_muon->fill(event);
      h_Lepton_jets->fill(event);
      h_Lepton_lumi->fill(event);
    }
  }
  ////



  bool presel = passed_recsel;

  /* *********** BTag Effi Hist *********** */
  if(passed_recsel){
    if(!event.isRealData) BTagEffHists->fill(event);
  }

  /* *********** lepton-2Dcut variables ***********  */
  bool pass_twodcut  = twodcut_sel->passes(event); {

    for(auto& muo : *event.muons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);

      muo.set_tag(Muon::twodcut_dRmin, dRmin);
      muo.set_tag(Muon::twodcut_pTrel, pTrel);
    }

    for(auto& ele : *event.electrons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

      ele.set_tag(Electron::twodcut_dRmin, dRmin);
      ele.set_tag(Electron::twodcut_pTrel, pTrel);
    }
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ------------------- NOCHMAL ANGUCKEN WENN E-CHANNEL ------------------------
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(elec_is_isolated) pass_twodcut = true; // do not do 2D cut for isolated electrons
  if(!pass_twodcut) passed_recsel = false;
  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);

  if(passed_recsel){
    h_TwoD_event->fill(event);
    h_TwoD_elec->fill(event);
    h_TwoD_muon->fill(event);
    h_TwoD_jets->fill(event);
    h_TwoD_lumi->fill(event);
  }
  ////

  /* *********** MET selection *********** */
  bool pass_met = met_sel->passes(event);
  if(!pass_met) passed_recsel = false;

  if(passed_recsel){
    h_MET_event->fill(event);
    h_MET_elec->fill(event);
    h_MET_muon->fill(event);
    h_MET_jets->fill(event);
    h_MET_lumi->fill(event);
  }

  /* *********** b-tag counter *********** */
  bool passed_btag = false;
  int jetbtagN(0);
  for(const auto& j : *event.jets) if(CSVBTag(CSVBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN == 0) passed_btag = false;
  else passed_btag = true;

  // BTagScaleFactors->process(event);

  if(passed_recsel && !passed_btag){
    h_Side_event->fill(event);
    h_Side_elec->fill(event);
    h_Side_muon->fill(event);
    h_Side_jets->fill(event);
    h_Side_lumi->fill(event);
  }

  // dont put bool here, set it in PostSel
  // if(!passed_btag) passed_recsel = false;

  if(passed_recsel){
    h_bTag_event->fill(event);
    h_bTag_elec->fill(event);
    h_bTag_muon->fill(event);
    h_bTag_jets->fill(event);
    h_bTag_lumi->fill(event);
  }

  if(presel && pass_twodcut && pass_met && passed_btag) h_cuts_all->fill(event);
  if(presel && !pass_twodcut && pass_met && passed_btag) h_cuts_all_but_2D->fill(event);
  if(presel && pass_twodcut && !pass_met && passed_btag) h_cuts_all_but_met->fill(event);
  if(presel && pass_twodcut && pass_met && !passed_btag) h_cuts_all_but_btag->fill(event);

  // only keep events that passed rec or gen solution
  if(recsel_only){
    if(!passed_recsel) return false;
  }
  if(!passed_gensel && !passed_recsel) return false;

  /* *********** now produce final XCone Jets and write output (especially weight) *********** */
  // store reco jets with and without JEC applied, and also copy uncorrected subjets
  std::vector<TopJet> jets = event.get(h_fatjets);
  if(jets.size() < 2) return false;
  if(!event.isRealData){
    std::vector<GenTopJet> gen33jets = event.get(h_gen33fatjets);
    if(gen33jets.size() < 2) return false;
  }

  // remove leptons from jets
  // cout << "before removing" << endl;
  remove_lepton_rec->process(event);
  if(!event.isRealData){
    remove_lepton_gen33->process(event);
  }
  // cout << "after removing" << endl;

  if(!event.isRealData){
    jetprod_gen33->process(event);
  }

  // Here all the Correction is happening
  jetprod_reco_noJEC->process(event);      // first store sum of 'raw' subjets
  copy_jet->process(event);                // copy 'raw' Fatets (with subjets) and name one copy 'noJEC'
  JetCorrections->process(event);          // apply AK4 JEC to subjets of the original Fatjet Collection
  JERSmearing->process(event);             // apply JER smearing to subjets
  if(do_nonClosure) NonClosureSYS->process(event);           // do nonClosure variation here
  jetprod_reco->process(event);            // now store sum of 'jec' subjets
  Correction->process(event);              // apply additional correction (a new 'cor' TopJet Collection is generated)
  jetprod_reco_corrected->process(event);  // finally store sum of 'cor' subjets
  ////

  output->process(event);

  /* *********** just a check ****************************************** */

  if(passed_recsel){
    ttbar_reweight->process(event);
    h_ttbar_reweight_event->fill(event);
    h_ttbar_reweight_elec->fill(event);
    h_ttbar_reweight_muon->fill(event);
    h_ttbar_reweight_jets->fill(event);
    h_ttbar_reweight_lumi->fill(event);
  }

  event.set(h_recsel_2, passed_recsel);
  event.set(h_gensel_2, passed_gensel);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetSelectionModule)
