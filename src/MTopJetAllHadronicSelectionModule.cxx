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

class MTopJetAllHadronicSelectionModule : public ModuleBASE {

public:
  explicit MTopJetAllHadronicSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners & Correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<JetCleaner> jet_cleaner1;
  std::unique_ptr<JetCleaner> jet_cleaner2;
  std::unique_ptr<JetCorrections_xcone> JetCorrections;
  std::unique_ptr<JER_Smearer_xcone> JERSmearing;
  std::unique_ptr<uhh2::AnalysisModule> Correction;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // Btag efficiency hists
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffHists;

  // selections
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_pupppi;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_noJEC;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco_corrected;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen23;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen33;
  std::unique_ptr<uhh2::AnalysisModule> copy_jet;
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> pv_sel;
  std::unique_ptr<uhh2::Selection> jet_sel;


  Event::Handle<std::vector<TopJet>>h_fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen33fatjets;
  Event::Handle<std::vector<GenTopJet>>h_gen23fatjets;

  // just for testing
  std::unique_ptr<TopPtReweight> ttbar_reweight;

  // store Hist collection as member variables

  bool isMC; //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part

  string corvar ="nominal";

};

MTopJetAllHadronicSelectionModule::MTopJetAllHadronicSelectionModule(uhh2::Context& ctx){

  isMC = (ctx.get("dataset_type") == "MC");
  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700_allHad"  ||
     ctx.get("dataset_version") == "TTbar_Mtt0700to1000_allHad"  ||
     ctx.get("dataset_version") == "TTbar_Mtt1000toInft_allHad" ) isTTbar = true;
  else isTTbar = false;


  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  if(isTTbar) ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  h_fatjets=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  h_gen23fatjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone23TopJets");
  h_gen33fatjets=ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  // combine XCone
  jetprod_reco.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined", "XCone33_lep_Combined", "xconeCHS"));
  jetprod_reco_noJEC.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined_noJEC", "XCone33_lep_Combined_noJEC", "xconeCHS" ));
  jetprod_reco_corrected.reset(new CombineXConeAllHad(ctx, "XCone33_had_Combined_Corrected", "XCone33_lep_Combined_Corrected", "xconeCHS_Corrected"));
  copy_jet.reset(new CopyJets(ctx, "xconeCHS", "xconeCHS_noJEC"));

  if(isMC){
    jetprod_gen23.reset(new CombineXCone23_gen(ctx));
    jetprod_gen33.reset(new CombineXConeAllHad_gen(ctx));
  }
  ////


  // double jecsysfactor = 1.0;
  corvar = ctx.get("JetCorrection_direction","nominal");

  // correct subjets (JEC + additional correction)
  JetCorrections.reset(new JetCorrections_xcone());
  JetCorrections->init(ctx, "xconeCHS");
  Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected", corvar, true));
  JERSmearing.reset(new JER_Smearer_xcone());
  JERSmearing->init(ctx, "xconeCHS", "genXCone33TopJets", "sub");
  //// EVENT SELECTION


  // define IDs
  JetId jetid_cleaner = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  JetId jetid_selection = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4));
  ////

  pv_sel.reset(new NPVSelection(1, -1, PrimaryVertexId(StandardPrimaryVertexId())));


  //// Obj Cleaning

  common.reset(new CommonModules());
  common->set_HTjetid(jetid_cleaner);
  common->switch_jetlepcleaner(true);
  common->switch_metcorrection();
  common->disable_mcpileupreweight(); // do this in PostSel
  common->init(ctx);

  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //
  BTagEffHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",CSVBTag::WP_TIGHT));

  //// set up Hists classes:
  //

}

bool MTopJetAllHadronicSelectionModule::process(uhh2::Event& event){

  /* *********** at least 1 good primary vertex *********** */
  if(!pv_sel->passes(event)) return false;

  /* *********** BTag Effi Hist *********** */
  if(!event.isRealData) BTagEffHists->fill(event);

  /* *********** b-tag counter *********** */
  bool passed_btag = false;
  int jetbtagN(0);
  for(const auto& j : *event.jets) if(CSVBTag(CSVBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN == 0) return false;


  // fill ttbargen class
  if(isTTbar) ttgenprod->process(event);
  ////

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

  if(!event.isRealData){
    jetprod_gen23->process(event);
    jetprod_gen33->process(event);
  }
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
