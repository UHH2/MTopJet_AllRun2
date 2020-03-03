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

// Hists
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGenHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/common/include/MCWeight.h> // nachträglich eingebaut

//
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

using namespace std;

class MTopJetAllHadronicModule : public ModuleBASE {

public:
  explicit MTopJetAllHadronicModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:

  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;
  std::unique_ptr<AnalysisModule> lumiweight;

  // selections
  std::unique_ptr<uhh2::Selection> lumi_sel;

  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::Selection> genflavor_sel;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> SemiLepDecay;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  //HISTS
  std::unique_ptr<Hists> h_ttbar;

  // handles
  Event::Handle<std::vector<TopJet>>h_fatjets;
  bool isMC;
};

MTopJetAllHadronicModule::MTopJetAllHadronicModule(uhh2::Context& ctx){

  //// CONFIGURATION
  isMC = (ctx.get("dataset_type") == "MC");

  h_fatjets = ctx.get_handle<std::vector<TopJet>>("xconeCHS");

  //HIST-classes
  h_ttbar.reset(new TTbarGenHists(ctx, "TTbar"));

  //// COMMON MODULES
  if(!isMC) lumi_sel.reset(new LumiSelection(ctx));

  /* GEN M-ttbar selection [TTbar MC "0.<M^{gen}_{ttbar}(GeV)<700.] */
  const std::string ttbar_gen_label("ttbargen");

  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700_allHad_2016v3") genmttbar_sel.reset(new MttbarGenSelection(0., 700.));
  else                                                                  genmttbar_sel.reset(new uhh2::AndSelection(ctx));

  /******************************************************************/

  ////

  MuonId muid = AndId<Muon>(MuonID(Muon::Tight), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Summer16_medium_noIso, PtEtaCut(55., 2.4));

  muoSR_cleaner.reset(new     MuonCleaner(muid));
  eleSR_cleaner.reset(new ElectronCleaner(eleid));


  //// EVENT SELECTION REC
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(50, 2.4))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(50, 2.4))));
  met_sel.reset(new METCut  (0, 100));
  muon_sel.reset(new NMuonSelection(0, 0, muid));
  elec_sel.reset(new NElectronSelection(0, 0, eleid));
  ////
  lumiweight.reset(new MCLumiWeight(ctx));

}

bool MTopJetAllHadronicModule::process(uhh2::Event& event){

  bool passed_recsel;

  if(!event.isRealData){
    /* GEN M-ttbar selection */
    ttgenprod->process(event);
    if(!genmttbar_sel->passes(event)) return false;
  }

  /* CMS-certified luminosity sections */
  if(event.isRealData){
    if(!lumi_sel->passes(event)) return false;
  }

  // cut on pT of large jets (R=1.2)
  std::vector<TopJet> jets = event.get(h_fatjets);
  double ptmax = 0;
  double ptcut = 200;
  for(auto jet:jets){
    if(jet.v4().Pt() > ptmax) ptmax = jet.v4().Pt();
  }

  muoSR_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);
  eleSR_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);

  const bool pass_fatjet = (ptmax > ptcut);
  const bool pass_jet2 = jet2_sel->passes(event);
  const bool pass_jet1 = jet1_sel->passes(event);
  const bool pass_met = met_sel->passes(event);
  const bool pass_lepsel = (muon_sel->passes(event) && elec_sel->passes(event));

  if(pass_jet2 && pass_jet1 && pass_met && pass_fatjet && pass_lepsel) passed_recsel = true;
  else passed_recsel = false;

  //FILL HISTS
  h_ttbar->fill(event);

  if(!passed_recsel) return false;
  else               return true;


}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetAllHadronicModule)
