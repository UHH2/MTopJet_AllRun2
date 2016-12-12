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


#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

class MTopJetSelectionModule : public ModuleBASE {

 public:
  explicit MTopJetSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners
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

  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::Selection> genflavor_sel;

  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<uhh2::Selection> muon_sel;
  std::unique_ptr<uhh2::Selection> elec_sel;
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> topjet_sel;
  std::unique_ptr<uhh2::Selection> topjet1_sel;
  std::unique_ptr<uhh2::Selection> topjet2_sel;
  // std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel100;
  std::unique_ptr<uhh2::Selection> htlep_sel150;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  std::unique_ptr<uhh2::Selection> deltarcut_sel;
  std::unique_ptr<uhh2::Selection> deltarak4cut_sel;
  std::unique_ptr<uhh2::Selection> topmass_sel;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // store Hist collection as member variables
  std::unique_ptr<Hists> h_PreSel_event, h_PreSel_elec, h_PreSel_muon, h_PreSel_jets, h_PreSel_topjets, h_PreSel_event2,
    h_twodsel_event, h_twodsel_elec, h_twodsel_muon, h_twodsel_jets, h_twodsel_topjets, h_twodsel_event2,
    h_trianc_event, h_trianc_elec, h_trianc_muon, h_trianc_jets, h_trianc_topjets, h_trianc_event2,
    h_htlep100_event, h_htlep100_elec, h_htlep100_muon, h_htlep100_jets, h_htlep100_topjets, h_htlep100_event2,
    h_htlep150_event, h_htlep150_elec, h_htlep150_muon, h_htlep150_jets, h_htlep150_topjets, h_htlep150_event2,
    // h_topjet_event, h_topjet_elec, h_topjet_muon, h_topjet_jets, h_topjet_topjets, h_topjet_event2,
    // h_toplepcleaner_event, h_toplepcleaner_elec, h_toplepcleaner_muon, h_toplepcleaner_jets, h_toplepcleaner_topjets, h_toplepcleaner_event2,
    // h_toplepdR_event, h_toplepdR_elec, h_toplepdR_muon, h_toplepdR_jets, h_toplepdR_topjets, h_toplepdR_event2,	
  // h_jetlepdR_event, h_jetlepdR_elec, h_jetlepdR_muon, h_jetlepdR_jets, h_jetlepdR_topjets, h_jetlepdR_event2, 
    // h_topmass_event, h_topmass_elec, h_topmass_muon, h_topmass_jets, h_topmass_topjets, h_topmass_event2,
    h_btag_event, h_btag_elec, h_btag_muon, h_btag_jets, h_btag_topjets, h_btag_event2;

  // Event::Handle<float> tt_TMVA_response;// response of TMVA method, dummy value at this step

};

MTopJetSelectionModule::MTopJetSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  //const bool isMC = (ctx.get("dataset_type") == "MC");

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {

    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";

    throw std::runtime_error(log);
  }

  ElectronId eleID;

  //    eleID  = ElectronID_Spring15_25ns_tight_noIso;
  eleID = ElectronID_MVAnotrig_Spring15_25ns_loose; //TEST 

  //// COMMON MODULES


  ////

  //// OBJ CLEANING

  // topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));  
 
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
  // jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));

  // met_sel  .reset(new METCut  (MET   , uhh2::infinity));
  htlep_sel100.reset(new HTlepCut(100, uhh2::infinity));
  htlep_sel150.reset(new HTlepCut(150, uhh2::infinity));

  twodcut_sel.reset(new TwoDCut1(.4, 40.));

  if     (channel_ == elec) triangc_sel.reset(new TriangularCuts(M_PI/2., (M_PI/2.)/75.));
  else if(channel_ == muon) triangc_sel.reset(new uhh2::AndSelection(ctx));
  ////


  //// set up Hists classes:
  h_PreSel_event.reset(new EventHists(ctx, "01_PreSel_Event"));
  h_PreSel_elec.reset(new ElectronHists(ctx, "01_PreSel_Elec"));
  h_PreSel_muon.reset(new MuonHists(ctx, "01_PreSel_Muon"));
  h_PreSel_jets.reset(new JetHists(ctx, "01_PreSel_Jets"));
  h_PreSel_topjets.reset(new TopJetHists(ctx, "01_PreSel_TopJets"));
  h_PreSel_event2.reset(new MTopJetHists(ctx, "01_PreSel_Event2"));

  h_twodsel_event.reset(new EventHists(ctx, "02_twodsel_Event"));
  h_twodsel_elec.reset(new ElectronHists(ctx, "02_twodsel_Elec"));
  h_twodsel_muon.reset(new MuonHists(ctx, "02_twodsel_Muon"));
  h_twodsel_jets.reset(new JetHists(ctx, "02_twodsel_Jets"));
  h_twodsel_topjets.reset(new TopJetHists(ctx, "02_twodsel_TopJets"));
  h_twodsel_event2.reset(new MTopJetHists(ctx, "02_twodsel_Event2"));

  h_trianc_event.reset(new EventHists(ctx, "02_trianc_Event"));
  h_trianc_elec.reset(new ElectronHists(ctx, "02_trianc_Elec"));
  h_trianc_muon.reset(new MuonHists(ctx, "02_trianc_Muon"));
  h_trianc_jets.reset(new JetHists(ctx, "02_trianc_Jets"));
  h_trianc_topjets.reset(new TopJetHists(ctx, "02_trianc_TopJets"));
  h_trianc_event2.reset(new MTopJetHists(ctx, "02_trianc_Event2"));

  h_btag_event.reset(new EventHists(ctx, "03_btag_Event"));
  h_btag_elec.reset(new ElectronHists(ctx, "03_btag_Elec"));
  h_btag_muon.reset(new MuonHists(ctx, "03_btag_Muon"));
  h_btag_jets.reset(new JetHists(ctx, "03_btag_Jets"));
  h_btag_topjets.reset(new TopJetHists(ctx, "03_btag_TopJets"));
  h_btag_event2.reset(new MTopJetHists(ctx, "03_btag_Event2"));

  h_htlep100_event.reset(new EventHists(ctx, "04_htlep100_Event"));
  h_htlep100_elec.reset(new ElectronHists(ctx, "04_htlep100_Elec"));
  h_htlep100_muon.reset(new MuonHists(ctx, "04_htlep100_Muon"));
  h_htlep100_jets.reset(new JetHists(ctx, "04_htlep100_Jets"));
  h_htlep100_topjets.reset(new TopJetHists(ctx, "04_htlep100_TopJets"));
  h_htlep100_event2.reset(new MTopJetHists(ctx, "04_htlep100_Event2"));

  h_htlep150_event.reset(new EventHists(ctx, "05_htlep150_Event"));
  h_htlep150_elec.reset(new ElectronHists(ctx, "05_htlep150_Elec"));
  h_htlep150_muon.reset(new MuonHists(ctx, "05_htlep150_Muon"));
  h_htlep150_jets.reset(new JetHists(ctx, "05_htlep150_Jets"));
  h_htlep150_topjets.reset(new TopJetHists(ctx, "05_htlep150_TopJets"));
  h_htlep150_event2.reset(new MTopJetHists(ctx, "05_htlep150_Event2"));
////

}

bool MTopJetSelectionModule::process(uhh2::Event& event){

  //// COMMON MODULES

  ////

  //// LEPTON selection

  ////

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



  //// JET selection

  /* 2nd AK4 jet selection */
  // const bool pass_jet2 = jet2_sel->passes(event);
  // if(!pass_jet2) return false;

  // /* 1st AK4 jet selection */
  // const bool pass_jet1 = jet1_sel->passes(event);
  // if(!pass_jet1) return false;

  const bool pass_trigger = trigger_sel->passes(event);
  if(!pass_trigger) return false;

  /* Only select event with exacly 1 muon or electron */
  const bool pass_lepsel = (muon_sel->passes(event) && elec_sel->passes(event));
  if(!pass_lepsel) return false;
  h_PreSel_event->fill(event);
  h_PreSel_elec->fill(event);
  h_PreSel_muon->fill(event);
  h_PreSel_jets->fill(event);
  h_PreSel_topjets->fill(event);
  h_PreSel_event2->fill(event);

  ////

  // //// MET selection
  // const bool pass_met = met_sel->passes(event);
  // if(!pass_met) return false;
  // h_metsel_event->fill(event);
  // h_metsel_elec->fill(event);
  // h_metsel_muon->fill(event);
  // h_metsel_jets->fill(event);
  // ////

  //// LEPTON-2Dcut selection
  if(!pass_twodcut) return false;
  h_twodsel_event->fill(event);
  h_twodsel_elec->fill(event);
  h_twodsel_muon->fill(event);
  h_twodsel_jets->fill(event);
  h_twodsel_topjets->fill(event);
  h_twodsel_event2->fill(event);
  ////

  //// Triangular Cut in Electron channel
  const bool pass_trianc = triangc_sel->passes(event);
  if(!pass_trianc) return false;
  h_trianc_event->fill(event);
  h_trianc_elec->fill(event);
  h_trianc_muon->fill(event);
  h_trianc_jets->fill(event);
  h_trianc_topjets->fill(event);
  h_trianc_event2->fill(event);
  ////

  //// b-tag counter
  int jetbtagN(0);
  for(const auto& j : *event.jets) if(CSVBTag(CSVBTag::WP_TIGHT)(j, event)) ++jetbtagN;
  if(jetbtagN < 1) return false;
  h_btag_event->fill(event);
  h_btag_elec->fill(event);
  h_btag_muon->fill(event);
  h_btag_jets->fill(event);
  h_btag_topjets->fill(event);
  h_btag_event2->fill(event);
  ////

  // //// HT_lep selection
  // const bool pass_htlep100 = htlep_sel100->passes(event);
  // if(!pass_htlep100) return false;
  // h_htlep100_event->fill(event);
  // h_htlep100_elec->fill(event);
  // h_htlep100_muon->fill(event);
  // h_htlep100_jets->fill(event);
  // h_htlep100_topjets->fill(event);
  // h_htlep100_event2->fill(event);

  const bool pass_htlep150 = htlep_sel150->passes(event);
  if(!pass_htlep150) return false;
  h_htlep150_event->fill(event);
  h_htlep150_elec->fill(event);
  h_htlep150_muon->fill(event);
  h_htlep150_jets->fill(event);
  h_htlep150_topjets->fill(event);
  h_htlep150_event2->fill(event);
  ////
return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetSelectionModule)
