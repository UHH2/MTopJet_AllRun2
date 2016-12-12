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
  std::unique_ptr<uhh2::AnalysisModule> output; 
  Event::Handle<bool> h_reco_sel;
  Event::Handle<bool> h_2Jet_sel;
  Event::Handle<bool> h_JetPt_sel;
  Event::Handle<bool> h_DeltaR_sel;
  Event::Handle<bool> h_Mass_sel;
  Event::Handle<int> h_cutflow;

  // cleaners
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;

  std::unique_ptr<JetCleaner>                      jet_IDcleaner;
  std::unique_ptr<JetCorrector>                    jet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer>     jetER_smearer;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> topjetlepton_cleaner;
  std::unique_ptr<JetCleaner>                      jet_cleaner1;
  std::unique_ptr<JetCleaner>                      jet_cleaner2;
  std::unique_ptr<JetLeptonCleaner>                jetlep_cleaner;

  std::unique_ptr<JetCleaner>                  topjet_IDcleaner;
  std::unique_ptr<TopJetCorrector>             topjet_corrector;
  std::unique_ptr<SubJetCorrector>             topjet_subjet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer> topjetER_smearer;
  std::unique_ptr<TopJetLeptonDeltaRCleaner>   topjetleptondeltaR_cleaner;
  std::unique_ptr<TopJetCleaner>               topjet_cleaner;
 
  // selections
 
  std::unique_ptr<uhh2::Selection> topjet_sel;
  std::unique_ptr<uhh2::Selection> topjetA_sel;
  std::unique_ptr<uhh2::Selection> topjetB_sel;
  std::unique_ptr<uhh2::Selection> deltarcut_sel;
  std::unique_ptr<uhh2::Selection> topmass_sel;
  std::unique_ptr<uhh2::Selection> met_sel;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // store Hist collection as member variables
  std::unique_ptr<Hists> 
  h_cleaner_event, h_cleaner_elec, h_cleaner_muon, h_cleaner_jets, h_cleaner_topjets, h_cleaner_MTopJetHists, h_cleaner_RecoHists_topjet,
    h_2ndJet_event, h_2ndJet_elec, h_2ndJet_muon, h_2ndJet_jets, h_2ndJet_topjets, h_2ndJet_MTopJetHists, h_2ndJet_RecoHists_topjet,

    h_topjetA_event, h_topjetA_elec, h_topjetA_muon, h_topjetA_jets, h_topjetA_topjets, h_topjetA_MTopJetHists, h_topjetA_RecoHists_topjet,
    h_toplepdRA_event, h_toplepdRA_elec, h_toplepdRA_muon, h_toplepdRA_jets, h_toplepdRA_topjets, h_toplepdRA_MTopJetHists, h_toplepdRA_RecoHists_topjet,
    h_topmassA_event, h_topmassA_elec, h_topmassA_muon, h_topmassA_jets, h_topmassA_topjets, h_topmassA_MTopJetHists, h_topmassA_RecoHists_topjet;

    // h_topjetB_event, h_topjetB_elec, h_topjetB_muon, h_topjetB_jets, h_topjetB_topjets, h_topjetB_MTopJetHists, h_topjetB_RecoHists_topjet,
    // h_toplepdRB_event, h_toplepdRB_elec, h_toplepdRB_muon, h_toplepdRB_jets, h_toplepdRB_topjets, h_toplepdRB_MTopJetHists, h_toplepdRB_RecoHists_topjet,
    // h_topmassB_event, h_topmassB_elec, h_topmassB_muon, h_topmassB_jets, h_topmassB_topjets, h_topmassB_MTopJetHists, h_topmassB_RecoHists_topjet;


};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  //// CONFIGURATION
  const bool isMC = (ctx.get("dataset_type") == "MC");

  //// Event Output
  output.reset(new WriteOutput(ctx));
  h_reco_sel = ctx.declare_event_output<bool> ("pass_RecoSel");
  h_2Jet_sel = ctx.declare_event_output<bool> ("pass_2JetSel");
  h_JetPt_sel = ctx.declare_event_output<bool> ("pass_JetPTSel");
  h_DeltaR_sel = ctx.declare_event_output<bool> ("pass_DeltaRSel");
  h_Mass_sel = ctx.declare_event_output<bool> ("pass_MassSel");
  h_cutflow = ctx.declare_event_output<int> ("cutflow");
  //// COMMON MODULES


  ////

  //// OBJ CLEANING
  std::vector<std::string> JEC_AK4, JEC_AK8;
  if(isMC){
    JEC_AK4 = JERFiles::Fall15_25ns_L123_AK4PFchs_MC;
    JEC_AK8 = JERFiles::Fall15_25ns_L123_AK8PFchs_MC;
  }
  else {
    JEC_AK4 = JERFiles::Fall15_25ns_L123_AK4PFchs_DATA;
    JEC_AK8 = JERFiles::Fall15_25ns_L123_AK8PFchs_DATA;
  }


  topjetlepton_cleaner.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK8, "topjets")); 
  topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(200., 2.4))));

  //// EVENT SELECTION


  topjet_sel.reset(new NTopJetSelection(2, 2, TopJetId(PtEtaCut(200, 2.4))));
  topjetA_sel.reset(new NTopJetSelection(1, 2, TopJetId(PtEtaCut(400, 2.4))));
 
  topmass_sel.reset(new TopJetMassCut());
  
  deltarcut_sel.reset(new deltaRCut(.8));

  ////

  //// set up Hists classes:


  // Hists after Cleaner
  h_cleaner_event.reset(new EventHists(ctx, "00_cleaner_Event"));
  h_cleaner_elec.reset(new ElectronHists(ctx, "00_cleaner_Elec"));
  h_cleaner_muon.reset(new MuonHists(ctx, "00_cleaner_Muon"));
  h_cleaner_jets.reset(new JetHists(ctx, "00_cleaner_Jets"));
  h_cleaner_topjets.reset(new TopJetHists(ctx, "00_cleaner_TopJets"));
  h_cleaner_MTopJetHists.reset(new MTopJetHists(ctx, "00_cleaner_MTopJetHists"));
  h_cleaner_RecoHists_topjet.reset(new RecoHists_topjet(ctx, "00_cleaner_MTopJetHists_RecoJetHists_topjet"));

  // 2 Jets with pt > 200 (SET TO 150??)
  h_2ndJet_event.reset(new EventHists(ctx, "00B_2ndJet_Event"));
  h_2ndJet_elec.reset(new ElectronHists(ctx, "00B_2ndJet_Elec"));
  h_2ndJet_muon.reset(new MuonHists(ctx, "00B_2ndJet_Muon"));
  h_2ndJet_jets.reset(new JetHists(ctx, "00B_2ndJet_Jets"));
  h_2ndJet_topjets.reset(new TopJetHists(ctx, "00B_2ndJet_TopJets"));
  h_2ndJet_MTopJetHists.reset(new MTopJetHists(ctx, "00B_2ndJet_MTopJetHists"));
  h_2ndJet_RecoHists_topjet.reset(new RecoHists_topjet(ctx, "00B_2ndJet_MTopJetHists_RecoJetHists_topjet"));
 
  // Selection A (1st topjet pt > 400)
  h_topjetA_event.reset(new EventHists(ctx, "01A_topjet_Event"));
  h_topjetA_elec.reset(new ElectronHists(ctx, "01A_topjet_Elec"));
  h_topjetA_muon.reset(new MuonHists(ctx, "01A_topjet_Muon"));
  h_topjetA_jets.reset(new JetHists(ctx, "01A_topjet_Jets"));
  h_topjetA_topjets.reset(new TopJetHists(ctx, "01A_topjet_TopJets"));
  h_topjetA_MTopJetHists.reset(new MTopJetHists(ctx, "01A_topjet_MTopJetHists"));
  h_topjetA_RecoHists_topjet.reset(new RecoHists_topjet(ctx, "01A_topjet_MTopJetHists_RecoJetHists_topjet"));

  h_toplepdRA_event.reset(new EventHists(ctx, "02A_toplepdR_Event"));
  h_toplepdRA_elec.reset(new ElectronHists(ctx, "02A_toplepdR_Elec"));
  h_toplepdRA_muon.reset(new MuonHists(ctx, "02A_toplepdR_Muon"));
  h_toplepdRA_jets.reset(new JetHists(ctx, "02A_toplepdR_Jets"));
  h_toplepdRA_topjets.reset(new TopJetHists(ctx, "02A_toplepdR_TopJets"));
  h_toplepdRA_MTopJetHists.reset(new MTopJetHists(ctx, "02A_toplepdR_MTopJetHists"));
  h_toplepdRA_RecoHists_topjet.reset(new RecoHists_topjet(ctx, "02A_toplepdR_RecoJetHists_topjet"));

  h_topmassA_event.reset(new EventHists(ctx, "03A_topmass_Event"));
  h_topmassA_elec.reset(new ElectronHists(ctx, "03A_topmass_Elec"));
  h_topmassA_muon.reset(new MuonHists(ctx, "03A_topmass_Muon"));
  h_topmassA_jets.reset(new JetHists(ctx, "03A_topmass_Jets"));
  h_topmassA_topjets.reset(new TopJetHists(ctx, "03A_topmass_TopJets"));
  h_topmassA_MTopJetHists.reset(new MTopJetHists(ctx, "03A_topmass_MTopJetHists"));
  h_topmassA_RecoHists_topjet.reset(new RecoHists_topjet(ctx, "03A_topmass_RecoJetHists_topjet"));

 
}

bool MTopJetPostSelectionModule::process(uhh2::Event& event){

  //// COMMON MODULES

  ////

  //// LEPTON selection

  ////

  //// JET selection

  ////

  //// Top Jet Selection
  topjetlepton_cleaner->process(event);
  topjet_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets); // Sort TopJets by pT


  h_cleaner_event->fill(event);
  h_cleaner_elec->fill(event);
  h_cleaner_muon->fill(event);
  h_cleaner_jets->fill(event);
  h_cleaner_topjets->fill(event);
  h_cleaner_MTopJetHists->fill(event);
  h_cleaner_RecoHists_topjet->fill(event);

  //--------- Set bools for every selection step --------------------
  bool pass_rec = false;
  const bool pass_topjet = topjet_sel->passes(event);
  const bool pass_topjetA = topjetA_sel->passes(event);
  const bool pass_deltaR = deltarcut_sel->passes(event);
  const bool pass_topmass = topmass_sel->passes(event);
  if(pass_topjet && pass_topjetA && pass_deltaR && pass_topmass) pass_rec = true;
  int cut_index = 0;
  if(pass_topjet){
    cut_index++;
    if(pass_topjetA){
      cut_index++;
      if(pass_deltaR){
	cut_index++;
	if(pass_topmass){
	  cut_index++;
	}
      }
    }
  }
  event.set(h_reco_sel, pass_rec);
  event.set(h_2Jet_sel, pass_topjet);
  event.set(h_JetPt_sel, pass_topjetA);
  event.set(h_DeltaR_sel, pass_deltaR);
  event.set(h_Mass_sel, pass_topmass);
  event.set(h_cutflow, cut_index);

  //------------------------------------------------------------------

  //==================================================================
  //================= Apply Selection ================================
  //==================================================================

  /* 2nd AK8 jet selection */
  if(!pass_topjet) return false;
  h_2ndJet_event->fill(event);
  h_2ndJet_elec->fill(event);
  h_2ndJet_muon->fill(event);
  h_2ndJet_jets->fill(event);
  h_2ndJet_topjets->fill(event);
  h_2ndJet_MTopJetHists->fill(event);
  h_2ndJet_RecoHists_topjet->fill(event);
  
  //--------- Selection with 1 TopJet pT > 400 --------------------
  /* 1st AK8 jet selection */
  if(pass_topjetA){
    h_topjetA_event->fill(event);
    h_topjetA_elec->fill(event);
    h_topjetA_muon->fill(event);
    h_topjetA_jets->fill(event);
    h_topjetA_topjets->fill(event);
    h_topjetA_MTopJetHists->fill(event);
    h_topjetA_RecoHists_topjet->fill(event);
    
    /* delta R (lep, topjet2) < 0.8  */
    if(pass_deltaR){
      h_toplepdRA_event->fill(event);
      h_toplepdRA_elec->fill(event);
      h_toplepdRA_muon->fill(event);
      h_toplepdRA_jets->fill(event);
      h_toplepdRA_topjets->fill(event);
      h_toplepdRA_MTopJetHists->fill(event);
      h_toplepdRA_RecoHists_topjet->fill(event);
      
      /* Top Jet Mass Cut M1 > M2 */
      if(pass_topmass){
	h_topmassA_event->fill(event);
	h_topmassA_elec->fill(event);
	h_topmassA_muon->fill(event);
	h_topmassA_jets->fill(event);
	h_topmassA_topjets->fill(event);
	h_topmassA_MTopJetHists->fill(event);
	h_topmassA_RecoHists_topjet->fill(event);
      }
    }
  }
  

  //--------------------------------------------------------------------


  output->process(event);

return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
