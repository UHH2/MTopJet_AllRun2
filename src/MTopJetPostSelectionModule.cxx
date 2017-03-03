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
#include "UHH2/common/include/MCWeight.h"

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>

#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/RecoHists_xcone.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>


/*
 *******************************************************************
**************** TO DO ********************************************
*******************************************************************
- scale factors ?
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

  // cleaners
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;


  std::unique_ptr<JetCleaner>                      jet_cleaner1;
  std::unique_ptr<JetCleaner>                      jet_cleaner2;


  // selections
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;


  // get weight (with all SF and weight applied in previous cycle)
  Event::Handle<double>h_weight;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_XCone, h_XCone_noMassCut, h_XCone_reco_gen;


  bool isMC; //define here to use it in "process" part
};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  isMC = (ctx.get("dataset_type") == "MC");

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {

    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";

    throw std::runtime_error(log);
  }

  // get handle for weight
  h_weight=ctx.get_handle<double>("weight");
  ////

 
  //// EVENT SELECTION
  pt_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined", 400));
  mass_sel.reset(new MassCutXCone(ctx));


  // define IDs
  ////

  // scale factors
  // muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID"));
  // muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger"));

  //// Obj Cleaning

  // common.reset(new CommonModules());
  // common->init(ctx);
  //

  //// set up Hists classes:

  h_XCone.reset(new RecoHists_xcone(ctx, "XCone"));
  h_XCone_reco_gen.reset(new RecoGenHists_xcone(ctx, "XCone_reco_gen"));
  h_XCone_noMassCut.reset(new RecoHists_xcone(ctx, "XCone_noMassCut"));

  //

}

bool MTopJetPostSelectionModule::process(uhh2::Event& event){
  //apply weight
  event.weight = event.get(h_weight);
  ////

  /* Events have to pass topjet pt > 400 & Mass_jet1 > Mass_jet2 */
  if(!pt_sel->passes(event)) return false;
  h_XCone_noMassCut->fill(event);
  if(!mass_sel->passes(event)) return false;
  h_XCone->fill(event);
  h_XCone_reco_gen->fill(event);

return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
