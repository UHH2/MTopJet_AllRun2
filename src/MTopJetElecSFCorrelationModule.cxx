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
#include <UHH2/common/include/TTbarGen.h>
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
#include <UHH2/MTopJet/include/ControlHists.h>
#include <UHH2/MTopJet/include/CutHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/NonClosureUncertainty.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/RemoveLepton.h>
#include <UHH2/MTopJet/include/ElecTriggerSF.h>

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

class MTopJetElecSFCorrelationModule : public ModuleBASE {

public:
  explicit MTopJetElecSFCorrelationModule(uhh2::Context&);
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

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections

  std::unique_ptr<uhh2::Selection> trigger_mu_A;
  std::unique_ptr<uhh2::Selection> trigger_mu_B;
  std::unique_ptr<uhh2::Selection> trigger_el_A;
  std::unique_ptr<uhh2::Selection> trigger_el_B;
  std::unique_ptr<uhh2::Selection> trigger_el_C;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel1;
  std::unique_ptr<uhh2::Selection> elec_sel2;
  std::unique_ptr<uhh2::Selection> elec_etaveto;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> sel_badhcal;

  Event::Handle<bool>h_recsel;
  std::unique_ptr<uhh2::Hists> h_all, h_muon, h_elec, h_both;

  std::unique_ptr<uhh2::AnalysisModule> ele_id_SF, ele_trigger_SFpt, ele_trigger_SFeta, ele_reco_SF;
  std::unique_ptr<uhh2::AnalysisModule> ele_trigger_SFetapt, ele_trigger_SFetaptUP, ele_trigger_SFetaptDOWN;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF, muo_trigger_SF_B;

  bool debug = false;
  bool isMC; //define here to use it in "process" part
  bool year_16;
  bool year_17;
  bool year_18;
  Year year;
};

MTopJetElecSFCorrelationModule::MTopJetElecSFCorrelationModule(uhh2::Context& ctx){

  if(debug) cout << "Get Year ... " << endl;
  year_16 = false;
  year_17 = false;
  year_18 = false;
  year = extract_year(ctx);

  if(year == Year::is2016v3) year_16 = true;
  else if(year == Year::is2017v2) year_17 = true;
  else if(year == Year::is2018) year_18 = true;
  else throw runtime_error("In PostSelectionModule: This Event is not from 2016v3, 2017v2 or 2018!");

  isMC = (ctx.get("dataset_type") == "MC");

  if(debug) cout << "Get ResSel ... " << endl;
  h_recsel = ctx.get_handle<bool>("passed_recsel");

  if(debug) cout << "Define IDs ... " << endl;
  // define IDs
  MuonId muid = AndId<Muon>(MuonID(Muon::Tight), PtEtaCut(55., 2.4));
  // this is only used for cleaner and electron veto
  ElectronId eleid_noiso55;
  if(year_16) eleid_noiso55 = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Summer16_tight_noIso);
  else        eleid_noiso55 = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Fall17_tight_noIso);
  // this is used to decide which ele trigger is used
  ElectronId eleid_noiso120;
  if(year_16) eleid_noiso120 = AndId<Electron>(PtEtaSCCut(120., 2.4), ElectronID_Summer16_tight_noIso);
  else        eleid_noiso120 = AndId<Electron>(PtEtaSCCut(120., 2.4), ElectronID_Fall17_tight_noIso);
  // this is used in combination with iso trigger
  ElectronId eleid_iso55;
  if(year_16)   eleid_iso55  = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Summer16_tight);
  else          eleid_iso55  = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Fall17_tight);
  // jet ids
  if(debug) cout << "Define Jet ID ... " << endl;

  JetId jetid_cleaner;
  // if(year_16) jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE_CHS), PtEtaCut(30.0, 2.4));
  jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));
  ////

  if(debug) cout << "Define Trigger ... " << endl;
  // define Trigger
  // trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu24_v*");
  // if(year == Year::is2017v2) trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu27_v*");
  // trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu24_v*");
  trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
  trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");
  if(year_16)      trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele27_WPTight_Gsf_v*");
  else if(year_17) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele35_WPTight_Gsf_v*");
  else if(year_18) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele32_WPTight_Gsf_v*");
  trigger_el_B = uhh2::make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
  if(year_16) trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon175_v*");
  else        trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon200_v*");

  // Scale Factors
  TString syear;
  if(year_16) syear = "2016";
  if(year_17) syear = "2017";
  if(year_18) syear = "2018";
  ele_trigger_SFpt.reset(new ElecTriggerSF(ctx, "nominal", "pt", syear));
  ele_trigger_SFeta.reset(new ElecTriggerSF(ctx, "nominal", "eta", syear));
  ele_trigger_SFetapt.reset(new ElecTriggerSF(ctx, "nominal", "eta_ptbins", syear));
  ele_trigger_SFetaptUP.reset(new ElecTriggerSF(ctx, "up", "eta_ptbins", syear));
  ele_trigger_SFetaptDOWN.reset(new ElecTriggerSF(ctx, "down", "eta_ptbins", syear));

  if(debug) cout << "Define SF ... " << endl;
  if(year_16){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root","NUM_TightID_DEN_genTracks_eta_pt",1, "tightID", false, "nominal"));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "trigger", false, "nominal"));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "reco", "nominal"));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2016/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "tightID", "nominal"));
  }
  else if(year_17){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonID_94X_RunBCDEF_SF_ID.root","NUM_TightID_DEN_genTracks_pt_abseta",1, "tightID", true, "nominal"));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","Mu50_PtEtaBins/pt_abseta_ratio",1, "trigger", true, "nominal"));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", 1.0, "reco", "nominal"));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/2017_ElectronTight.root", 1.0, "tightID", "nominal"));
  }
  else if(year_18){
    muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_ID_SF_RunABCD.root","NUM_TightID_DEN_TrackerMuons_pt_abseta",1, "tightID", true, "nominal"));
    muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/ScaleFactors/Muons/Muon_Trigger_SF_2018.root","Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/pt_abseta_ratio",1, "trigger", true, "nominal"));
    ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/egammaEffi.txt_EGM2D_updatedAll.root", 1.0, "reco", "nominal"));
    ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/2018_ElectronTight.root", 1.0, "tightID", "nominal"));
  }
  else throw runtime_error("In ElecSFModule: There is no Event from 2016_v3, 2017_v2 or 2018!");

  if(debug) cout << "Define Selection ... " << endl;
  /*Only select event with exacly 1 muon or electron */
  muon_sel.reset(new NMuonSelection(1, 1, muid));
  elec_sel1.reset(new NElectronSelection(1, 1, eleid_iso55));
  elec_sel2.reset(new NElectronSelection(1, 1, eleid_noiso120));

  elec_etaveto.reset(new ElectronEtaVeto(1.44, 1.57));

  met_sel  .reset(new METCut  (50 , uhh2::infinity));
  twodcut_sel.reset(new TwoDCut1(0.4, 40));
  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));
  sel_badhcal.reset(new BadHCALSelection(ctx));

  if(debug) cout << "Common Modules ... " << endl;
  //// Obj Cleaning
  common.reset(new CommonModules());
  common->set_HTjetid(jetid_cleaner);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  // common->disable_mcpileupreweight(); // do this in PostSel
  common->init(ctx);

  muoSR_cleaner.reset(new MuonCleaner(muid));
  eleSR_cleaner.reset(new ElectronCleaner(eleid_noiso55));

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  h_all.reset(new CountingEventHists(ctx, "no_trigger"));
  h_muon.reset(new CountingEventHists(ctx, "muon_trigger"));
  h_elec.reset(new CountingEventHists(ctx, "elec_trigger"));
  h_both.reset(new CountingEventHists(ctx, "both_trigger"));

}


bool MTopJetElecSFCorrelationModule::process(uhh2::Event& event){
  if(debug) cout << " ------------------------------------------- " << endl;
  if(debug) cout << " ----------------- NewEvent ---------------- " << endl;
  if(debug) cout << " ------------------------------------------- " << endl;

  bool passed_recsel = event.get(h_recsel);
  if(!passed_recsel) return false;

  /* *********** Apply Prefireing Weights *********** */
  if(debug) cout << "No Pre: " << event.weight << endl;
  if(!event.isRealData) event.weight *= event.prefiringWeight;
  if(debug) cout << "-- Pre: " << event.weight << endl;

  /* *********** Lepton Cleaner and Selection *********** */
  if(debug) cout << "Lepton Cleaner ... " << endl;
  muoSR_cleaner->process(event);
  if(event.muons->size() > 0) sort_by_pt<Muon>(*event.muons);

  eleSR_cleaner->process(event);
  if(event.electrons->size() > 0) sort_by_pt<Electron>(*event.electrons);


  const bool pass_lep1 = ((event.muons->size() > 0) || (event.electrons->size() > 0));
  if(!pass_lep1) return false;
  ////

  /* *********** Cleaner ********** */
  if(debug) cout << "Cleaner ... " << endl;
  if(!common->process(event)) return false;
  if(debug) cout << "Jet Cleaner ... " << endl;
  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);

  /* *********** at least 1 good primary vertex *********** */
  if(!pv_sel->passes(event)) return false;

  /* *********** LEPTON Selection *********** */
  if(debug) cout << "Lepton Selection ... " << endl;
  if(!muon_sel->passes(event)) return false;
  bool pass_elec1 = elec_sel1->passes(event);
  bool pass_elec2 = elec_sel2->passes(event);
  if(!pass_elec1 && !pass_elec2) return false;

  if(!elec_etaveto->passes(event)) return false;
  ////

  /* *********** lepton-2Dcut variables ***********  */
  if(debug) cout << "2D Cut ... " << endl;
  bool pass_twodcut  = twodcut_sel->passes(event); {

    for(auto& muo : *event.muons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);

      muo.set_tag(Muon::twodcut_dRmin, dRmin);
      muo.set_tag(Muon::twodcut_pTrel, pTrel);
    }
  }
  if(!pass_twodcut) return false;

  /* *********** MET selection *********** */
  if(debug) cout << "Apply SFs ... " << endl;
  // bool pass_met = met_sel->passes(event);
  // if(!pass_met) return false;

  /* *********** bad HCAL selection for 2018 *********** */
  if(!sel_badhcal->passes(event)) return false;

  // Muon Scale Factors
  muo_tight_noniso_SF->process(event);
  muo_trigger_SF->process(event);

  // Elec Scale Factors
  ele_id_SF->process(event);
  ele_reco_SF->process(event);

  // ===========================================================================
  //          Trigger
  // ===========================================================================
  /* *********** Trigger *********** */
  if(debug) cout << "Trigger... " << endl;
  // for DATA until run 274954 -> use only Trigger A
  // for MC and DATA from 274954 -> use "A || B"
  // HLT_TkMu50_v is not availible for 2017&18. Alternative still needs to be inlcuded !!!!!!!!!!!!!!!!

  // ===========================================================================
  //          Muon

  bool pass_muon=true;
  if(year == Year::is2016v3){
    if(!isMC && event.run < 274954) pass_muon = trigger_mu_A->passes(event);
    else pass_muon = trigger_mu_A->passes(event) || trigger_mu_B->passes(event);
  }
  if(year == Year::is2017v2) pass_muon = trigger_mu_A->passes(event);
  if(year == Year::is2018) pass_muon = trigger_mu_A->passes(event);

  // ===========================================================================
  //          Elec

  if(debug) cout << "Start Fill ... " << endl;
  bool pass_elec = true;
  // note that at least one of pass_elec1 and pass_elec2 has to be true here!
  if(pass_elec1 && !pass_elec2) pass_elec = trigger_el_A->passes(event);
  else{
    if(year == Year::is2016v3)  pass_elec = (trigger_el_B->passes(event) || trigger_el_C->passes(event));
    if(year == Year::is2017v2){
      // for MC event.run=1
      if(!isMC && event.run <= 299329) pass_elec = (trigger_el_A->passes(event) || trigger_el_C->passes(event));
      else                             pass_elec = (trigger_el_B->passes(event) || trigger_el_C->passes(event));
    }
    if(year == Year::is2018)  pass_elec = (trigger_el_B->passes(event) || trigger_el_C->passes(event));
  }


  // bool pass_elec = false;
  // if(pass_elec1) pass_elec = (trigger_el_A->passes(event) || trigger_el_C->passes(event)); // noHighTrigger

  if(debug) cout << "Fill hists ... " << endl;

  h_all->fill(event);
  if(pass_muon) h_muon->fill(event);
  if(pass_elec) h_elec->fill(event);
  if(pass_muon && pass_elec) h_both->fill(event);

  // if(debug) cout << "Event done ... " << endl;
  if(debug) cout << "Event done ... " << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetElecSFCorrelationModule)
