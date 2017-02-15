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
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGenHists.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>


#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>


/*
 *******************************************************************
**************** TO DO ********************************************
*******************************************************************

- MET FILTER
- elec/mu ID
- CLEANER
- COMBINE JETS
- 2D CUT
- HT_lep CUT
- MET CUT
- ==1 mu/elec
- TRIGGER

*******************************************************************
*******************************************************************
*/

class MTopJetSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;

  std::unique_ptr<JetCleaner>                      jet_IDcleaner;
  std::unique_ptr<JetCorrector>                    jet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer>     jetER_smearer;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> jetlepton_cleaner;
  std::unique_ptr<JetCleaner>                      jet_cleaner1;
  std::unique_ptr<JetCleaner>                      jet_cleaner2;

  std::unique_ptr<JetCleaner>                  topjet_IDcleaner;
  std::unique_ptr<TopJetCorrector>             topjet_corrector;
  std::unique_ptr<SubJetCorrector>             topjet_subjet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer> topjetER_smearer;
  std::unique_ptr<TopJetLeptonDeltaRCleaner>   topjetlepton_cleaner;
  std::unique_ptr<TopJetCleaner>               topjet_cleaner;

  // selections
  std::unique_ptr<uhh2::AndSelection> metfilters_sel;
  std::unique_ptr<uhh2::AnalysisModule> jetprod_reco;

  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> triangc_sel;
  // std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> jet_sel;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // store Hist collection as member variables
  std::unique_ptr<Hists> h_PreSel_event,  h_PreSel_elec, h_PreSel_muon, h_PreSel_jets;
  std::unique_ptr<Hists> h_Cleaner_event,  h_Cleaner_elec, h_Cleaner_muon, h_Cleaner_jets;
  std::unique_ptr<Hists> h_Lepton_event,  h_Lepton_elec, h_Lepton_muon, h_Lepton_jets;
  std::unique_ptr<Hists> h_TwoD_event,  h_TwoD_elec, h_TwoD_muon, h_TwoD_jets;
  std::unique_ptr<Hists> h_Jet_event,  h_Jet_elec, h_Jet_muon, h_Jet_jets;
  std::unique_ptr<Hists> h_bTag_event,  h_bTag_elec, h_bTag_muon, h_bTag_jets;
  std::unique_ptr<Hists> h_HTlep_event,  h_HTlep_elec, h_HTlep_muon, h_HTlep_jets;
};

MTopJetSelectionModule::MTopJetSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  // const bool isMC = (ctx.get("dataset_type") == "MC");

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {

    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";

    throw std::runtime_error(log);
  }

  // combine XCone
  jetprod_reco.reset(new CombineXCone33(ctx)); 
  ////

 
 //// EVENT SELECTION

  const std::string& trigger = ctx.get("trigger", "NULL");
  if(trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));
  else                  trigger_sel.reset(new uhh2::AndSelection(ctx));

  /*Only select event with exacly 1 muon or electron */
  if(channel_ == elec){
    muon_sel.reset(new NMuonSelection(0, 0));
    elec_sel.reset(new NElectronSelection(1, 1));
  }
  else if (channel_ == muon){
    muon_sel.reset(new NMuonSelection(1, 1));
    elec_sel.reset(new NElectronSelection(0, 0));
  }

  // jet2_sel.reset(new Jet2Cut(50));
  jet_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(50, 2.4))));

  // met_sel  .reset(new METCut  (MET   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(150, uhh2::infinity));

  twodcut_sel.reset(new TwoDCutALL(0.4, 40));

  if     (channel_ == elec) triangc_sel.reset(new TriangularCuts(M_PI/2., (M_PI/2.)/75.));
  else if(channel_ == muon) triangc_sel.reset(new uhh2::AndSelection(ctx));
  ////


  //// Obj Cleaning
  ElectronId EleId = AndId<Electron>(PtEtaCut(40., 2.4), ElectronID_Spring16_tight_noIso);
  MuonId MuId = AndId<Muon>(PtEtaCut(40., 2.1), MuonIDMedium());
  JetId JetId = JetPFID(JetPFID::WP_LOOSE);

  common.reset(new CommonModules());
  common->set_muon_id(MuId);
  common->set_electron_id(EleId);
  // common->set_jet_id();
  common->init(ctx);

  muoSR_cleaner.reset(new     MuonCleaner(MuId));
  eleSR_cleaner.reset(new ElectronCleaner(EleId));

 
  /*****************************************************/
  /* JEC File for Data Sets!                           */
  /*****************************************************/

  // const JetId jetID(JetPFID(JetPFID::WP_LOOSE));
  // std::vector<std::string> JEC_AK4;
  // if(isMC) JEC_AK4 = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC; //
  // else JEC_AK4 = JERFiles::Spring16_25ns_L123_AK4PFchs_DATA;
  // jet_IDcleaner.reset(new JetCleaner(ctx, jetID));
  // jet_corrector.reset(new JetCorrector(ctx, JEC_AK4));
  // if(isMC) jetER_smearer.reset(new GenericJetResolutionSmearer(ctx));
  // jetlepton_cleaner.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4));
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //

  //// set up Hists classes:
  h_PreSel_event.reset(new MTopJetHists(ctx, "00_PreSel_Event"));
  h_PreSel_elec.reset(new ElectronHists(ctx, "00_PreSel_Elec"));
  h_PreSel_muon.reset(new MuonHists(ctx, "00_PreSel_Muon"));
  h_PreSel_jets.reset(new JetHists(ctx, "00_PreSel_Jets"));

  h_Lepton_event.reset(new MTopJetHists(ctx, "01_Lepton_Event"));
  h_Lepton_elec.reset(new ElectronHists(ctx, "01_Lepton_Elec"));
  h_Lepton_muon.reset(new MuonHists(ctx, "01_Lepton_Muon"));
  h_Lepton_jets.reset(new JetHists(ctx, "01_Lepton_Jets"));

  h_Cleaner_event.reset(new MTopJetHists(ctx, "02_Cleaner_Event"));
  h_Cleaner_elec.reset(new ElectronHists(ctx, "02_Cleaner_Elec"));
  h_Cleaner_muon.reset(new MuonHists(ctx, "02_Cleaner_Muon"));
  h_Cleaner_jets.reset(new JetHists(ctx, "02_Cleaner_Jets"));

  h_TwoD_event.reset(new MTopJetHists(ctx, "03_TwoD_Event"));
  h_TwoD_elec.reset(new ElectronHists(ctx, "03_TwoD_Elec"));
  h_TwoD_muon.reset(new MuonHists(ctx, "03_TwoD_Muon"));
  h_TwoD_jets.reset(new JetHists(ctx, "03_TwoD_Jets"));

  h_Jet_event.reset(new MTopJetHists(ctx, "04_Jet_Event"));
  h_Jet_elec.reset(new ElectronHists(ctx, "04_Jet_Elec"));
  h_Jet_muon.reset(new MuonHists(ctx, "04_Jet_Muon"));
  h_Jet_jets.reset(new JetHists(ctx, "04_Jet_Jets"));

  h_bTag_event.reset(new MTopJetHists(ctx, "05_bTag_Event"));
  h_bTag_elec.reset(new ElectronHists(ctx, "05_bTag_Elec"));
  h_bTag_muon.reset(new MuonHists(ctx, "05_bTag_Muon"));
  h_bTag_jets.reset(new JetHists(ctx, "05_bTag_Jets"));

  h_HTlep_event.reset(new MTopJetHists(ctx, "06_HTlep_Event"));
  h_HTlep_elec.reset(new ElectronHists(ctx, "06_HTlep_Elec"));
  h_HTlep_muon.reset(new MuonHists(ctx, "06_HTlep_Muon"));
  h_HTlep_jets.reset(new JetHists(ctx, "06_HTlep_Jets"));
  //

}

bool MTopJetSelectionModule::process(uhh2::Event& event){

  h_PreSel_event->fill(event);
  h_PreSel_elec->fill(event);
  h_PreSel_muon->fill(event);
  h_PreSel_jets->fill(event);
  ////

  /* Trigger */
  const bool pass_trigger = trigger_sel->passes(event);
  if(!pass_trigger) return false;

  /* Lepton Cleaner and Selection */
  muoSR_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);

  eleSR_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);

  const bool pass_lepsel = (muon_sel->passes(event) && elec_sel->passes(event));
  if(!pass_lepsel) return false;
 
  h_Lepton_event->fill(event);
  h_Lepton_elec->fill(event);
  h_Lepton_muon->fill(event);
  h_Lepton_jets->fill(event);
  ////


  /* Cleaner */
  common->process(event);
  // jet_IDcleaner->process(event);
  // jet_corrector->process(event);
  // if(jetER_smearer.get()) jetER_smearer->process(event);
  // jetlepton_cleaner->process(event);
  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);

  h_Cleaner_event->fill(event);
  h_Cleaner_elec->fill(event);
  h_Cleaner_muon->fill(event);
  h_Cleaner_jets->fill(event);


  /* lepton-2Dcut variables */
  const bool pass_twodcut = twodcut_sel->passes(event); {

    for(auto& muo : *event.muons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);

      muo.set_tag(Muon::twodcut_dRmin, dRmin);
      muo.set_tag(Muon::twodcut_pTrel, pTrel);
    }

    for(auto& ele : *event.electrons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

      ele.set_tag(Electron::twodcut_dRmin, dRmin);
      ele.set_tag(Electron::twodcut_pTrel, pTrel);
    }
  }

  if(!pass_twodcut) return false; // 2D Cut
  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);

  h_TwoD_event->fill(event);
  h_TwoD_elec->fill(event);
  h_TwoD_muon->fill(event);
  h_TwoD_jets->fill(event);
  ////




  /* Jet Selection */
  const bool pass_jetsel = (jet_sel->passes(event));
  if(!pass_jetsel) return false;

  // h_Jet_event->fill(event);
  h_Jet_elec->fill(event);
  h_Jet_muon->fill(event);
  h_Jet_jets->fill(event);

  /* MET selection */
  // const bool pass_met = met_sel->passes(event);
  // if(!pass_met) return false;

  /* Triangular Cut in Electron channel */
  const bool pass_trianc = triangc_sel->passes(event);
  if(!pass_trianc) return false;

  /* b-tag counter */
  int jetbtagN(0);
  for(const auto& j : *event.jets) if(CSVBTag(CSVBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN < 1) return false;

  h_bTag_event->fill(event);
  h_bTag_elec->fill(event);
  h_bTag_muon->fill(event);
  h_bTag_jets->fill(event);

  /* HT lep */
  const bool pass_htlep = htlep_sel->passes(event);
  if(!pass_htlep) return false;

  h_HTlep_event->fill(event);
  h_HTlep_elec->fill(event);
  h_HTlep_muon->fill(event);
  h_HTlep_jets->fill(event);

  /* now produce final XCone Jets */
  jetprod_reco->process(event);


return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetSelectionModule)
