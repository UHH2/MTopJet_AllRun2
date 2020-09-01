#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/CleaningModules.h>
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
#include <UHH2/MTopJet/include/ControlHists.h>
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

class MTopJetAllHadronicSelectionModule : public ModuleBASE {

public:
  explicit MTopJetAllHadronicSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  std::unique_ptr<CommonModules> common;
  // cleaners & Correctors
  std::unique_ptr<JetCleaner> jet_cleaner1;
  std::unique_ptr<JetCleaner> jet_cleaner2;
  std::unique_ptr<JetCorrections_xcone> JetCorrections;
  std::unique_ptr<JER_Smearer_xcone> JERSmearing;
  std::unique_ptr<uhh2::AnalysisModule> Correction;

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_pupppi;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_noJEC;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_corrected;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen33;
  std::unique_ptr<uhh2::AnalysisModule> copy_jet;
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> jet_sel;

  //HISTS

  Event::Handle<std::vector<TopJet>>h_fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen33fatjets;

  // just for testing
  std::unique_ptr<TopPtReweight> ttbar_reweight;

  // store Hist collection as member variables

  bool isMC; //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part

  string corvar = "nominal";

  JetId btag_tight;
  Year year;

};

MTopJetAllHadronicSelectionModule::MTopJetAllHadronicSelectionModule(uhh2::Context& ctx){

  //======================= YearSwitcher =======================================
  bool year_16 = false;
  bool year_17 = false;
  bool year_18 = false;
  year = extract_year(ctx);

  if(year == Year::is2016v3) year_16 = true;
  else if(year == Year::is2017v2) year_17 = true;
  else if(year == Year::is2018) year_18 = true;
  else throw runtime_error("In AllHadronicSelectionModule: This Event is not from 2016v3, 2017v2 or 2018!");

  /*************************** CONFIGURATION **********************************************************************************/
  isMC = (ctx.get("dataset_type") == "MC");
  TString dataset_version = (TString) ctx.get("dataset_version");
  if(dataset_version.Contains("TTbar")) isTTbar = true;
  else  isTTbar = false;


  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  h_gen33fatjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  // combine XCone
  jetprod_reco.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "xconeCHS"));
  jetprod_reco_noJEC.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined_noJEC", "XCone33_lep_Combined_noJEC", "xconeCHS" ));
  jetprod_reco_corrected.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined_Corrected", "XCone33_lep_Combined_Corrected", "xconeCHS_Corrected"));
  copy_jet.reset(new CopyJets(ctx, "xconeCHS", "xconeCHS_noJEC"));

  if(isMC) jetprod_gen33.reset(new CombineXConeAllHad_gen(ctx));
  ////


  // double jecsysfactor = 1.0;
  corvar = ctx.get("JetCorrection_direction","nominal");

  // correct subjets (JEC + additional correction)
  JetCorrections.reset(new JetCorrections_xcone());
  JetCorrections->init(ctx, "xconeCHS");

  if(year_16) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2016"));
  else if(year_17) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2017"));
  else if(year_18) Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true, "2018"));
  else throw runtime_error("In PostSelectionModule: There is no Event from 2016_v2, 2017_v2 or 2018!");

  JERSmearing.reset(new JER_Smearer_xcone());
  JERSmearing->init(ctx, "xconeCHS", "genXCone33TopJets", "sub");
  //// EVENT SELECTION


  // define IDs
  JetId jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30.0, 2.4));
  JetId jetid_selection = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(50.0, 2.4));
  ////

  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));


  //// Obj Cleaning

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //
  btag_tight = DeepJetBTag(DeepJetBTag::WP_TIGHT);
  BTagEffHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",btag_tight));

  //// set up Hists classes:

  common.reset(new CommonModules());
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  common->disable_mcpileupreweight(); // do this in PostSel
  common->init(ctx);

}

bool MTopJetAllHadronicSelectionModule::process(uhh2::Event& event){
  if(!common->process(event)) return false;
  /* *********** at least 1 good primary vertex *********** */
  if(!pv_sel->passes(event)) return false;

  /* *********** BTag Effi Hist *********** */
  if(!event.isRealData) BTagEffHists->fill(event);

  /* *********** b-tag counter *********** */
  // bool passed_btag = false;
  int jetbtagN(0);
  for(const auto& j : *event.jets) if(DeepJetBTag(DeepJetBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN == 0) return false;


  /* *********** now produce final XCone Jets and write output (especially weight) *********** */
  // store reco jets with and without JEC applied, and also copy uncorrected subjets
  std::vector<TopJet> jets = event.get(h_fatjets);
  if(jets.size() < 2) return false;
  if(!event.isRealData){
    std::vector<GenTopJet> gen33jets = event.get(h_gen33fatjets);
    if(gen33jets.size() < 2) return false;
  }

  if(!event.isRealData) jetprod_gen33->process(event);

  // Here all the Correction is happening
  jetprod_reco_noJEC->process(event);      // first store sum of 'raw' subjets
  copy_jet->process(event);                // copy 'raw' Fatets (with subjets) and name one copy 'noJEC'
  JetCorrections->process(event);          // apply AK4 JEC/JER to subjets of the original Fatjet Collection
  JERSmearing->process(event);             // apply JER smearing on subjets
  jetprod_reco->process(event);            // now store sum of 'jec' subjets
  Correction->process(event);              // apply additional correction (a new 'cor' TopJet Collection is generated)
  jetprod_reco_corrected->process(event);  // finally store sum of 'cor' subjets
  ////
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetAllHadronicSelectionModule)
