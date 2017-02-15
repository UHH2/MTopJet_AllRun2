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
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

// Hists
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGenHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
//
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>


class MTopJetPreSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetPreSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:


  // selections
  std::unique_ptr<uhh2::Selection> lumi_sel;

  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::Selection> genflavor_sel;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // store Hist collection as member variables
  // std::unique_ptr<Hists> h_ttbar;
  // std::unique_ptr<Hists> h_PreSel_event, h_PreSel_event2, h_PreSel_elec, h_PreSel_muon, h_PreSel_jets;
};

MTopJetPreSelectionModule::MTopJetPreSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  const bool isMC = (ctx.get("dataset_type") == "MC");


  ////

  //// COMMON MODULES

  if(!isMC) lumi_sel.reset(new LumiSelection(ctx));



  /* GEN M-ttbar selection [TTbar MC "0.<M^{gen}_{ttbar}(GeV)<700.] */
  const std::string ttbar_gen_label("ttbargen");

  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700") genmttbar_sel.reset(new MttbarGenSelection(0., 700.));
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));

  /******************************************************************/

  ////



  //// EVENT SELECTION
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(50, 2.4))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(50, 2.4))));
  met_sel.reset(new METCut  (20, uhh2::infinity));
  muon_sel.reset(new NMuonSelection(1, -1, MuonId(PtEtaCut(40, 2.4 ))));
  elec_sel.reset(new NElectronSelection(1, -1, ElectronId(PtEtaCut(40, 2.4))));
  ////

  //// set up Hists classes:
  // h_PreSel_event.reset(new EventHists(ctx, "01_PreSel_Event"));
  // h_PreSel_event2.reset(new MTopJetHists(ctx, "01_PreSel_Event2"));
  // h_PreSel_elec.reset(new ElectronHists(ctx, "01_PreSel_Elec"));
  // h_PreSel_muon.reset(new MuonHists(ctx, "01_PreSel_Muon"));
  // h_PreSel_jets.reset(new JetHists(ctx, "01_PreSel_Jets"));
  // h_ttbar.reset(new TTbarGenHists(ctx, "TTbarHists"));
}

bool MTopJetPreSelectionModule::process(uhh2::Event& event){

  //// COMMON MODULES

  if(!event.isRealData){

    /* GEN M-ttbar selection */
    ttgenprod->process(event);
    if(!genmttbar_sel->passes(event)) return false;
  }

  /* CMS-certified luminosity sections */
  if(event.isRealData){

    if(!lumi_sel->passes(event)) return false;
  }


  const bool pass_lep1 = ((event.muons->size() >= 1) || (event.electrons->size() >= 1));
  if(!pass_lep1) return false;


  /* JET selection */
  /* 2nd AK4 jet selection */
  const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;

  /* 1st AK4 jet selection */
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;


  /* MET selection */
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;

  /* select min 1 Muon OR min 1 Electron */
  const bool pass_lepsel = (muon_sel->passes(event) || elec_sel->passes(event));
  if(!pass_lepsel) return false;

  // h_PreSel_event->fill(event);
  // h_PreSel_event2->fill(event);
  // h_PreSel_elec->fill(event);
  // h_PreSel_muon->fill(event);
  // h_PreSel_jets->fill(event);
  // h_ttbar->fill(event);
 
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPreSelectionModule)
