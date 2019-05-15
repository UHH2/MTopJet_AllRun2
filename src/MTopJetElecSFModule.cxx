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
#include <UHH2/MTopJet/include/CutHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/NonClosureUncertainty.h>
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

class MTopJetElecSFModule : public ModuleBASE {

public:
  explicit MTopJetElecSFModule(uhh2::Context&);
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

  Event::Handle<bool>h_recsel;
  Event::Handle<double>h_pt;
  Event::Handle<double>h_eta;
  Event::Handle<double>h_weight;
  Event::Handle<double>h_weight_SFpt;
  Event::Handle<double>h_weight_SFeta;
  Event::Handle<double>h_weight_SFetapt;
  Event::Handle<double>h_weight_SFetaptUP;
  Event::Handle<double>h_weight_SFetaptDOWN;
  Event::Handle<bool>h_passed;

  std::unique_ptr<uhh2::AnalysisModule> ele_id_SF, ele_trigger_SFpt, ele_trigger_SFeta, ele_reco_SF;
  std::unique_ptr<uhh2::AnalysisModule> ele_trigger_SFetapt, ele_trigger_SFetaptUP, ele_trigger_SFetaptDOWN;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;

  bool isMC; //define here to use it in "process" part


  std::unique_ptr<Hists> h_pass, h_all;
};

MTopJetElecSFModule::MTopJetElecSFModule(uhh2::Context& ctx){

  isMC = (ctx.get("dataset_type") == "MC");


  h_recsel = ctx.get_handle<bool>("passed_recsel");

  ctx.undeclare_all_event_output();

  h_eta = ctx.declare_event_output<double>("eta");
  h_pt = ctx.declare_event_output<double>("pt");
  h_weight = ctx.declare_event_output<double>("weight");
  h_weight_SFpt = ctx.declare_event_output<double>("weight_sfpt");
  h_weight_SFeta = ctx.declare_event_output<double>("weight_sfeta");
  h_weight_SFetapt = ctx.declare_event_output<double>("weight_sfetapt");
  h_weight_SFetaptUP = ctx.declare_event_output<double>("weight_sfetaptUP");
  h_weight_SFetaptDOWN = ctx.declare_event_output<double>("weight_sfetaptDOWN");
  h_passed = ctx.declare_event_output<bool>("passed");

  // define IDs
  MuonId muid = AndId<Muon>(MuonIDTight(), PtEtaCut(55., 2.4));
  // this is only used for cleaner and electron veto
  ElectronId eleid_noiso55  = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Spring16_tight_noIso);
  // this is used to decide which ele trigger is used
  ElectronId eleid_noiso120 = AndId<Electron>(PtEtaSCCut(120., 2.4), ElectronID_Spring16_tight_noIso);
  // this is used in combination with iso trigger
  ElectronId eleid_iso55    = AndId<Electron>(PtEtaSCCut(55., 2.4), ElectronID_Spring16_tight);
  // jet ids
  JetId jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  ////

  // define Trigger
  trigger_mu_A = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
  trigger_mu_B = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");
  trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele27_WPTight_Gsf_v*");
  trigger_el_B = uhh2::make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
  trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon175_v*");

  // Scale Factors
  ele_trigger_SFpt.reset(new ElecTriggerSF(ctx, "nominal", "pt"));
  ele_trigger_SFeta.reset(new ElecTriggerSF(ctx, "nominal", "eta"));
  ele_trigger_SFetapt.reset(new ElecTriggerSF(ctx, "nominal", "eta_ptbins"));
  ele_trigger_SFetaptUP.reset(new ElecTriggerSF(ctx, "up", "eta_ptbins"));
  ele_trigger_SFetaptDOWN.reset(new ElecTriggerSF(ctx, "down", "eta_ptbins"));
  ele_reco_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", "nominal"));
  ele_id_SF.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "", "nominal"));
  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID", true, "nominal"));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger", true, "nominal"));


  /*Only select event with exacly 1 muon or electron */
  muon_sel.reset(new NMuonSelection(1, 1, muid));
  elec_sel1.reset(new NElectronSelection(1, 1, eleid_iso55));
  elec_sel2.reset(new NElectronSelection(1, 1, eleid_noiso120));

  elec_etaveto.reset(new ElectronEtaVeto(1.44, 1.57));

  met_sel  .reset(new METCut  (50 , uhh2::infinity));
  twodcut_sel.reset(new TwoDCut1(0.4, 40));
  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));

  //// Obj Cleaning

  common.reset(new CommonModules());
  common->set_HTjetid(jetid_cleaner);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  // common->disable_mcpileupreweight(); // do this in PostSel
  common->init(ctx);


  muoSR_cleaner.reset(new     MuonCleaner(muid));
  eleSR_cleaner.reset(new ElectronCleaner(eleid_noiso55));

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));


  h_pass.reset(new ElectronHists(ctx, "pass_Elec"));
  h_all.reset(new ElectronHists(ctx, "all_Elec"));


}


bool MTopJetElecSFModule::process(uhh2::Event& event){

  bool passed_recsel = event.get(h_recsel);
  if(!passed_recsel) return false;

  /* *********** Lepton Cleaner and Selection *********** */
  muoSR_cleaner->process(event);
  if(event.muons->size() > 0) sort_by_pt<Muon>(*event.muons);

  eleSR_cleaner->process(event);
  if(event.electrons->size() > 0) sort_by_pt<Electron>(*event.electrons);


  const bool pass_lep1 = ((event.muons->size() > 0) || (event.electrons->size() > 0));
  if(!pass_lep1) return false;
  ////

  /* *********** Cleaner ********** */
  if(!common->process(event)) return false;

  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);

  /* *********** at least 1 good primary vertex *********** */
  if(!pv_sel->passes(event)) return false;


  /* *********** Trigger *********** */
  // for DATA until run 274954 -> use only Trigger A
  // for MC and DATA from 274954 -> use "A || B"
  if( !isMC && event.run < 274954) {
    if(!trigger_mu_A->passes(event)) return false;
  }
  else{
    if( !(trigger_mu_A->passes(event) || trigger_mu_B->passes(event)) ) return false;
  }

  /* *********** lEPTON Selection *********** */
  if(!muon_sel->passes(event)) return false;
  bool pass_elec1 = elec_sel1->passes(event);
  bool pass_elec2 = elec_sel2->passes(event);
  if(!pass_elec1 && !pass_elec2) return false;

  if(!elec_etaveto->passes(event)) return false;
  ////

  /* *********** lepton-2Dcut variables ***********  */
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
  // bool pass_met = met_sel->passes(event);
  // if(!pass_met) return false;

  // Muon Scale Factors
  muo_tight_noniso_SF->process(event);
  muo_trigger_SF->process(event);

  // Elec Scale Factors
  ele_id_SF->process(event);
  ele_reco_SF->process(event);

  // fill hists with all events
  h_all->fill(event);


  // HERE FILL PT AND ETA HISTS FOR PASSING AND NOT PASSING ELEC TRIGGER
  bool passed_elec_trigger = true;
  // note that at least one of pass_elec1 and pass_elec2 has to be true here!
  if(!pass_elec2){
    if(!trigger_el_A->passes(event)) passed_elec_trigger = false;
  }
  else{
    if( !(trigger_el_B->passes(event) || trigger_el_C->passes(event)) ) passed_elec_trigger = false;
  }

  event.set(h_pt, event.electrons->at(0).pt());
  event.set(h_eta, event.electrons->at(0).eta());
  event.set(h_weight, event.weight);
  event.set(h_passed, passed_elec_trigger);

  // now apply SF and store weight after SF
  // first store raw weight
  double wraw = event.weight;
  // use pt SF and store new weight
  ele_trigger_SFpt->process(event);
  event.set(h_weight_SFpt, event.weight);
  event.weight = wraw;
  // reset weight and do eta SF
  ele_trigger_SFeta->process(event);
  event.set(h_weight_SFeta, event.weight);
  event.weight = wraw;
  // reset weight and do eta_ptbin SF
  ele_trigger_SFetapt->process(event);
  event.set(h_weight_SFetapt, event.weight);
  event.weight = wraw;
  // up variation
  ele_trigger_SFetaptUP->process(event);
  event.set(h_weight_SFetaptUP, event.weight);
  event.weight = wraw;
  // down variation
  ele_trigger_SFetaptDOWN->process(event);
  event.set(h_weight_SFetaptDOWN, event.weight);

  if(passed_elec_trigger){
    // fill pass histograms
    h_pass->fill(event);
  }



  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetElecSFModule)
