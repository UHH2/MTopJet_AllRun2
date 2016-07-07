#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
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
#include <UHH2/common/include/TTbarGenHists.h>

#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/MTopJetSelections.h>
#include <UHH2/MTopJet/include/MTopJetGenSelections.h>
#include <UHH2/MTopJet/include/MTopJetGenHists.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

class MTopJetGenPostSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetGenPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:

  // cleaners
  
  // selections
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;



  // store Hist collection as member variables
  std::unique_ptr<Hists> h_GenHists;


};

MTopJetGenPostSelectionModule::MTopJetGenPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION


  ////

  //// COMMON MODULES

  const std::string ttbar_gen_label("ttbargen");

  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  ////

  //// OBJ CLEANING
 
  //

  // set up Hists classes:
  h_GenHists.reset(new MTopJetGenHists(ctx, "GenHists"));

  // EVENT SELECTION

}

bool MTopJetGenPostSelectionModule::process(uhh2::Event& event){

  //  COMMON MODULES
  ttgenprod->process(event);
  h_GenHists->fill(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetGenPostSelectionModule)
