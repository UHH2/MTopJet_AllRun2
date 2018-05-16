#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

// Hists
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/PseudoXConeHists.h>
//
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

using namespace std;

class PseudoXConeModule : public ModuleBASE {

 public:
  explicit PseudoXConeModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:


  // selections
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;

  uhh2::Event::Handle<std::vector<TopJet>>h_xcone;


  // store Hist collection as member variables
  std::unique_ptr<Hists> h_pseudoxcone, h_pseudoxcone_highpt;
};

PseudoXConeModule::PseudoXConeModule(uhh2::Context& ctx){


  h_xcone=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  ////

  //// EVENT SELECTION REC
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(50, 2.4))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(50, 2.4))));
  met_sel.reset(new METCut  (20, uhh2::infinity));
  muon_sel.reset(new NMuonSelection(1, -1, MuonId(PtEtaCut(40, 2.4 ))));
  elec_sel.reset(new NElectronSelection(1, -1, ElectronId(PtEtaCut(40, 2.4))));
  ////

  h_pseudoxcone.reset(new PseudoXConeHists(ctx, "PseudoXCone"));
  h_pseudoxcone_highpt.reset(new PseudoXConeHists(ctx, "PseudoXCone_highpt"));
}

bool PseudoXConeModule::process(uhh2::Event& event){

  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  const bool pass_jet2 = jet2_sel->passes(event);
  const bool pass_jet1 = jet1_sel->passes(event);
  const bool pass_met = met_sel->passes(event);
  const bool pass_lepsel = (muon_sel->passes(event) || elec_sel->passes(event));
  ///

  if(!(pass_lep1 && pass_jet2 && pass_jet1 && pass_met && pass_lepsel)) return false;

  // fill hists
  h_pseudoxcone->fill(event);
  // go into boosted regime
  std::vector<TopJet> xconejets = event.get(h_xcone);
  if(xconejets.at(0).v4().Pt() < 400) return false;
  //fill hists again
  h_pseudoxcone_highpt->fill(event);
  return true;


}

UHH2_REGISTER_ANALYSIS_MODULE(PseudoXConeModule)
