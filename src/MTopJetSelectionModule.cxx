#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
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
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/TopPtReweight.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CutHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>

#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include "UHH2/MTopJet/include/CorrectionFactor.h"


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
  std::unique_ptr<uhh2::AnalysisModule> Correction;

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_noJEC;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_corrected;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen23;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen33;
  std::unique_ptr<uhh2::AnalysisModule> copy_jet;
  std::unique_ptr<uhh2::Selection> trigger_sel_A;
  std::unique_ptr<uhh2::Selection> trigger_sel_B;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> jet_sel;

  Event::Handle<bool>h_gensel;
  Event::Handle<bool>h_recsel;
  Event::Handle<bool>h_gensel_2;
  Event::Handle<bool>h_recsel_2;
 


  Event::Handle<std::vector<TopJet>>h_fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen33fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen23fatjets;

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
  bool jecsys;

  string corvar ="nominal";

};

MTopJetSelectionModule::MTopJetSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700"  || 
     ctx.get("dataset_version") == "TTbar_Mtt0700to1000"  || 
     ctx.get("dataset_version") == "TTbar_Mtt1000toInft"  ||
     ctx.get("dataset_version") == "TTbar_mtop1665"       ||
     ctx.get("dataset_version") == "TTbar_mtop1695_ext1"  ||
     ctx.get("dataset_version") == "TTbar_mtop1695_ext2"  ||
     ctx.get("dataset_version") == "TTbar_mtop1715"       ||
     ctx.get("dataset_version") == "TTbar_mtop1735"       ||
     ctx.get("dataset_version") == "TTbar_mtop1755"       ||
     ctx.get("dataset_version") == "TTbar_mtop1785"       ||
     ctx.get("dataset_version") == "TTbar_amcatnlo-pythia"||
     ctx.get("dataset_version") == "TTbar_powheg-herwig") isTTbar = true;
  else  isTTbar = false;

  if(isTTbar) h_gensel = ctx.get_handle<bool>("passed_gensel");
  h_recsel = ctx.get_handle<bool>("passed_recsel");
  h_gensel_2 = ctx.declare_event_output<bool>("passed_gensel_2");
  h_recsel_2 = ctx.declare_event_output<bool>("passed_recsel_2");
  h_fatjets=ctx.get_handle<std::vector<TopJet>>("XConeTopJets");
  h_gen23fatjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone23TopJets");
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

  // combine XCone
  jetprod_reco.reset(new CombineXCone33(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "XConeTopJets")); 
  jetprod_reco_noJEC.reset(new CombineXCone33(ctx, "XCone33_had_Combined_noJEC", "XCone33_lep_Combined_noJEC", "XConeTopJets" )); 
  jetprod_reco_corrected.reset(new CombineXCone33(ctx, "XCone33_had_Combined_Corrected", "XCone33_lep_Combined_Corrected", "XConeTopJets_Corrected")); 
  copy_jet.reset(new CopyJets(ctx, "XConeTopJets", "XConeTopJets_noJEC")); 

  if(isMC){
    jetprod_gen23.reset(new CombineXCone23_gen(ctx));
    jetprod_gen33.reset(new CombineXCone33_gen(ctx));
  }	
  ////

  // write output
  output.reset(new WriteOutput(ctx));
  ////

  // just for testing
  // ttbar_reweight.reset(new TopPtReweight(ctx,0.159,-0.00141,"","weight_ttbar",true,0.9910819)); // 8 TeV
  ttbar_reweight.reset(new TopPtReweight(ctx,0.0615,-0.0005,"","weight_ttbar",true)); // 13 TeV

  jecsys = false;
  // double jecsysfactor = 1.0;

  corvar = ctx.get("JetCorrection_direction","nominal");

  // correct subjets (JEC + additional correction)
  JetCorrections.reset(new JetCorrections_xcone());
  JetCorrections->init(ctx, "XConeTopJets");
  Correction.reset(new CorrectionFactor(ctx, "XConeTopJets_Corrected", corvar));

  //// EVENT SELECTION


  // define IDs
  MuonId muid = AndId<Muon>(MuonIDTight(), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Spring16_medium_noIso, PtEtaCut(55., 2.4));
  JetId jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  JetId jetid_selection = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4));
  ////

  // define Trigger
  if(channel_ == muon){
    trigger_sel_A = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
    trigger_sel_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");
  }
  ////

  /*Only select event with exacly 1 muon or electron */
  if(channel_ == elec){
    muon_sel.reset(new NMuonSelection(0, 0, muid));
    elec_sel.reset(new NElectronSelection(1, 1, eleid));
  }
  else if (channel_ == muon){
    muon_sel.reset(new NMuonSelection(1, 1, muid));
    elec_sel.reset(new NElectronSelection(0, 0, eleid));
  }

  jet_sel.reset(new NJetSelection(2, -1, jetid_selection));
  met_sel  .reset(new METCut  (50 , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(100, uhh2::infinity));
  twodcut_sel.reset(new TwoDCut1(0.4, 40));
  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));

  // if     (channel_ == elec) triangc_sel.reset(new TriangularCuts(M_PI/2., (M_PI/2.)/75.));
  // else if(channel_ == muon) triangc_sel.reset(new uhh2::AndSelection(ctx));
  ////


   //// Obj Cleaning

  common.reset(new CommonModules());
  common->set_HTjetid(jetid_cleaner);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  common->disable_mcpileupreweight(); // do this in PostSel
  common->init(ctx);
  

  muoSR_cleaner.reset(new     MuonCleaner(muid));
  eleSR_cleaner.reset(new ElectronCleaner(eleid));

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //
  BTagEffHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",CSVBTag::WP_TIGHT));

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

bool MTopJetSelectionModule::process(uhh2::Event& event){

  bool passed_gensel;
  if(isTTbar) passed_gensel = event.get(h_gensel);
  else passed_gensel = false;

  bool passed_recsel;
  if(isTTbar) passed_recsel = event.get(h_recsel);
  passed_recsel = true;

  // cout << "*******************************" << endl;
  // cout << "processing event nr. " << event.event<< endl;

  // if(passed_recsel){
  //   h_PreSel_event->fill(event);
  //   h_PreSel_elec->fill(event);
  //   h_PreSel_muon->fill(event);
  //   h_PreSel_jets->fill(event);
  //   h_PreSel_lumi->fill(event);
  // }
  ////

  /* *********** Lepton Cleaner and Selection *********** */
  if(passed_recsel){
    muoSR_cleaner->process(event);
    sort_by_pt<Muon>(*event.muons);

    eleSR_cleaner->process(event);
    sort_by_pt<Electron>(*event.electrons);
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
  if(channel_ == muon){
    if(passed_recsel){
      if( !isMC && event.run < 274954) {
	if(!trigger_sel_A->passes(event)) passed_recsel = false;
      }else{
	if( !(trigger_sel_A->passes(event) || trigger_sel_B->passes(event)) ) passed_recsel = false;
      }
    }
  }
  // if(channel_ == muon) muo_trigger_SF->process(event);

  if(passed_recsel){
    h_Trigger_event->fill(event);
    h_Trigger_elec->fill(event);
    h_Trigger_muon->fill(event);
    h_Trigger_jets->fill(event);
    h_Trigger_lumi->fill(event);
  }
  /* *********** lEPTON Selection *********** */
  if(passed_recsel){
    const bool pass_lepsel = (muon_sel->passes(event) && elec_sel->passes(event));
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

  /* *********** Jet Selection *********** */
  bool pass_jetsel = false;
  if(passed_recsel){
    pass_jetsel = (jet_sel->passes(event));
    if(!pass_jetsel) passed_recsel = false;
    else{
      h_Jet_elec->fill(event);
      h_Jet_muon->fill(event);
      h_Jet_jets->fill(event);
      h_Jet_lumi->fill(event);
    }
  }

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
  if(jecsys){
    if(!passed_recsel) return false;
  }
  if(!passed_gensel && !passed_recsel) return false;


  /* *********** now produce final XCone Jets and write output (especially weight) *********** */
  // store reco jets with and without JEC applied, and also copy uncorrected subjets
  std::vector<TopJet> jets = event.get(h_fatjets);
  if(jets.size() < 2) return false;
  if(!event.isRealData){
    std::vector<GenTopJet> gen23jets = event.get(h_gen23fatjets);
    if(gen23jets.size() < 2) return false;
    std::vector<GenTopJet> gen33jets = event.get(h_gen33fatjets);
    if(gen33jets.size() < 2) return false;
  }

  jetprod_reco_noJEC->process(event);
  copy_jet->process(event);
  JetCorrections->process(event);
  jetprod_reco->process(event);
  Correction->process(event);
  jetprod_reco_corrected->process(event);


  if(!event.isRealData){
    jetprod_gen23->process(event);
    jetprod_gen33->process(event);
  } 

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
