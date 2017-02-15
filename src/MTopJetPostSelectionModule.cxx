#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
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

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGenHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/RecoHists_topjet.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>


#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

class MTopJetPostSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  enum lepton { muon, elec };
  lepton channel_;

  // Event Output
  // std::unique_ptr<uhh2::AnalysisModule> output; 
  // Event::Handle<bool> h_reco_sel;
  // Event::Handle<bool> h_2Jet_sel;
  // Event::Handle<bool> h_JetPt_sel;
  // Event::Handle<bool> h_DeltaR_sel;
  // Event::Handle<bool> h_Mass_sel;
  // Event::Handle<int> h_cutflow;

  // cleaners

 
  // selections
 
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;


  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_final_xcone;

};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION

  //// Event Output
  // output.reset(new WriteOutput(ctx));
  // h_reco_sel = ctx.declare_event_output<bool> ("pass_RecoSel");
  // h_2Jet_sel = ctx.declare_event_output<bool> ("pass_2JetSel");
  // h_JetPt_sel = ctx.declare_event_output<bool> ("pass_JetPTSel");
  // h_DeltaR_sel = ctx.declare_event_output<bool> ("pass_DeltaRSel");
  // h_Mass_sel = ctx.declare_event_output<bool> ("pass_MassSel");
  // h_cutflow = ctx.declare_event_output<int> ("cutflow");
  //// COMMON MODULES


  ////

  //// OBJ CLEANING


  //// EVENT SELECTION


  pt_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined", 400));
  mass_sel.reset(new MassCutXCone(ctx));

  ////

  //// set up Hists classes:

  h_final_xcone.reset(new RecoHists_xcone(ctx, "Final_XCone"));

 
}

bool MTopJetPostSelectionModule::process(uhh2::Event& event){


  //--------- Set bools for every selection step --------------------
  //------------------------------------------------------------------

  //==================================================================
  //================= Apply Selection ================================
  //==================================================================

  bool pass_pt = pt_sel->passes(event);
  if(!pass_pt) return false;
  bool pass_mass = mass_sel->passes(event);
  if(!pass_mass) return false;

  h_final_xcone->fill(event);
  //--------------------------------------------------------------------

  // output->process(event);

return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
