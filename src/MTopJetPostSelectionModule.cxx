#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGen.h>

#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/RecoHists_xcone.h>
#include <UHH2/MTopJet/include/GenHists_xcone.h>
#include <UHH2/MTopJet/include/GenHists_particles.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/SubjetHists_xcone.h>
#include "UHH2/MTopJet/include/RecoGenHists_subjets.h"

#include <vector>

/*
*******************************************************************
**************** TO DO ********************************************
*******************************************************************
- SF apllied to TTbar!
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

  // cleaners & Correctors
  // std::unique_ptr<JetCorrections_xcone> JetCorrections;

  // selections
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> pt200_sel;
  std::unique_ptr<uhh2::Selection> pt300_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;
  std::unique_ptr<uhh2::Selection> pt_gensel;
  std::unique_ptr<uhh2::Selection> mass_gensel;
  std::unique_ptr<uhh2::Selection> pt_gensel23;
  std::unique_ptr<uhh2::Selection> mass_gensel23;
  std::unique_ptr<uhh2::Selection> matched_sub;
  std::unique_ptr<uhh2::Selection> matched_fat;

  // get weight (with all SF and weight applied in previous cycle)
  Event::Handle<double>h_weight;

  // handles for output
  Event::Handle<bool>h_matched;
  Event::Handle<bool>h_recsel;
  Event::Handle<bool>h_gensel23;
  Event::Handle<bool>h_gensel33;
  Event::Handle<bool>h_ttbar;
  Event::Handle<double>h_ttbar_SF;
  Event::Handle<double>h_mass_gen33;
  Event::Handle<double>h_mass_gen23;
  Event::Handle<double>h_mass_rec;
  Event::Handle<std::vector<Jet>>h_recjets_had;
  Event::Handle<std::vector<Particle>>h_genjets23_had;
  Event::Handle<std::vector<Particle>>h_genjets33_had;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // store Hist collection as member variables
  std::unique_ptr<Hists> h_Muon;
  std::unique_ptr<Hists> h_XCone, h_XCone_m, h_XCone_u, h_XCone_m_fat, h_XCone_u_fat, h_XCone_pt200, h_XCone_pt300, h_MTopJet, h_XCone_noMassCut;
  std::unique_ptr<Hists> h_XConeNoJEC, h_XCone_subjets_noJEC, h_XConeNoJEC_noMassCut;
  std::unique_ptr<Hists> h_XCone_subjets;
  std::unique_ptr<Hists> h_XCone_lowPU_subjets, h_XCone_medPU_subjets, h_XCone_highPU_subjets, h_XCone_lowPU, h_XCone_medPU, h_XCone_highPU;
  std::unique_ptr<Hists> h_RecGenHists_lowPU, h_RecGenHists_medPU, h_RecGenHists_highPU, h_RecGenHists_lowPU_noJEC, h_RecGenHists_medPU_noJEC, h_RecGenHists_highPU_noJEC;
  std::unique_ptr<Hists> h_XCone_GEN_RecOnly, h_XCone_GEN_GenOnly, h_XCone_GEN_Both;
  std::unique_ptr<Hists> h_RecGenHists_GenOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly_noJEC;
  std::unique_ptr<Hists> h_RecGenHists_Both;
  std::unique_ptr<Hists> h_RecGenHists_subjets, h_RecGenHists_subjets_noJEC;
  std::unique_ptr<Hists> h_GenParticles_RecOnly, h_GenParticles_GenOnly, h_GenParticles_Both;

  bool isMC; //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part

};

MTopJetPostSelectionModule::MTopJetPostSelectionModule(uhh2::Context& ctx){

  /*************************** CONFIGURATION **********************************************************************************/ 
  isMC = (ctx.get("dataset_type") == "MC");
  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700" || ctx.get("dataset_version") == "TTbar_Mtt0700to1000"  || ctx.get("dataset_version") == "TTbar_Mtt1000toInft") isTTbar = true;
  else  isTTbar = false;

  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  if(isTTbar) ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));

  const std::string& channel = ctx.get("channel", ""); //define Channel
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else {
    std::string log("MTopJetPostSelectionModule::MTopJetPostSelectionModule -- ");
    log += "invalid argument for 'channel' key in xml file (must be 'muon' or 'elec'): \""+channel+"\"";
    throw std::runtime_error(log);
  }

  // get handle for weight
  h_weight=ctx.get_handle<double>("weight");

  // write output
  h_matched = ctx.declare_event_output<bool>("matched");
  h_recsel = ctx.declare_event_output<bool>("passed_recsel");
  h_gensel23 = ctx.declare_event_output<bool>("passed_gensel23");
  h_gensel33 = ctx.declare_event_output<bool>("passed_gensel33");
  h_ttbar = ctx.declare_event_output<bool>("is_TTbar");
  h_ttbar_SF = ctx.declare_event_output<double>("TTbar_SF");
  h_mass_gen23 = ctx.declare_event_output<double>("Mass_Gen23");
  h_mass_gen33 = ctx.declare_event_output<double>("Mass_Gen33");
  h_mass_rec = ctx.declare_event_output<double>("Mass_Rec");
  h_recjets_had = ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
  if(isMC) h_genjets23_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone23_had_Combined");
  if(isMC) h_genjets33_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone33_had_Combined");

  /*************************** Setup Subjet Corrector **********************************************************************************/ 
  // JetCorrections.reset(new JetCorrections_xcone());
  // JetCorrections->init(ctx, "XConeTopJets");
  /*************************** Setup Selections **********************************************************************************/ 
  // RECO Selection
  pt_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined_noJEC", 400));
  pt300_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined_noJEC", 300));
  pt200_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined_noJEC", 200));
  mass_sel.reset(new MassCutXCone(ctx, "XCone33_had_Combined_noJEC", "XCone33_lep_Combined_noJEC"));

  // GEN Selection
  pt_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 400));
  mass_gensel.reset(new MassCut_gen(ctx, "GEN_XCone33_had_Combined", "GEN_XCone33_lep_Combined"));
  pt_gensel23.reset(new LeadingJetPT_gen(ctx, "GEN_XCone23_had_Combined", 400));
  mass_gensel23.reset(new MassCut_gen(ctx, "GEN_XCone23_had_Combined", "GEN_XCone23_lep_Combined"));

  // Selection for matching reco jets to gen particles
  if(isTTbar) matched_sub.reset(new Matching_XCone33(ctx, true));
  if(isTTbar) matched_fat.reset(new Matching_XCone33(ctx, false));

  /*************************** Set up Hists classes **********************************************************************************/ 
  // "true" for XCone Hists means to use jets with JEC applied
  h_XCone.reset(new RecoHists_xcone(ctx, "XCone", true));
  h_XCone_lowPU.reset(new RecoHists_xcone(ctx, "XCone_lowPU", true));
  h_XCone_medPU.reset(new RecoHists_xcone(ctx, "XCone_medPU", true));
  h_XCone_highPU.reset(new RecoHists_xcone(ctx, "XCone_highPU", true));
  h_XConeNoJEC.reset(new RecoHists_xcone(ctx, "XConeNoJEC", false));
  h_XCone_subjets.reset(new SubjetHists_xcone(ctx, "XCone_subjets", true));
  h_XCone_subjets_noJEC.reset(new SubjetHists_xcone(ctx, "XCone_subjets_noJEC", false));
  h_XCone_lowPU_subjets.reset(new SubjetHists_xcone(ctx, "XCone_lowPU_subjets", true));
  h_XCone_medPU_subjets.reset(new SubjetHists_xcone(ctx, "XCone_medPU_subjets", true));
  h_XCone_highPU_subjets.reset(new SubjetHists_xcone(ctx, "XCone_highPU_subjets", true));
  if(isTTbar) h_XCone_m.reset(new RecoHists_xcone(ctx, "XCone_matched", true));
  if(isTTbar) h_XCone_u.reset(new RecoHists_xcone(ctx, "XCone_unmatched", true));
  if(isTTbar) h_XCone_m_fat.reset(new RecoHists_xcone(ctx, "XCone_matched_fat", true));
  if(isTTbar) h_XCone_u_fat.reset(new RecoHists_xcone(ctx, "XCone_unmatched_fat", true));
  h_XCone_pt200.reset(new RecoHists_xcone(ctx, "XCone_pt200", true));
  h_XCone_pt300.reset(new RecoHists_xcone(ctx, "XCone_pt300", true));
  h_MTopJet.reset(new MTopJetHists(ctx, "EventHists"));
  h_XCone_noMassCut.reset(new RecoHists_xcone(ctx, "XCone_noMassCut", true));
  h_XConeNoJEC_noMassCut.reset(new RecoHists_xcone(ctx, "XConeNoJEC_noMassCut", false));
  h_Muon.reset(new MuonHists(ctx, "Muon"));

  if(isMC){
    h_XCone_GEN_GenOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly"));
    h_XCone_GEN_RecOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_RecOnly"));
    h_XCone_GEN_Both.reset(new GenHists_xcone(ctx, "XCone_GEN_Both"));

    h_GenParticles_GenOnly.reset(new GenHists_particles(ctx, "GenParticles_GenOnly"));
    h_GenParticles_RecOnly.reset(new GenHists_particles(ctx, "GenParticles_RecOnly"));
    h_GenParticles_Both.reset(new GenHists_particles(ctx, "GenParticles_Both"));

    h_RecGenHists_subjets.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets", true));
    h_RecGenHists_subjets_noJEC.reset(new RecoGenHists_subjets(ctx, "RecGenHists_subjets_noJEC", false));

    h_RecGenHists_lowPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_lowPU", true));
    h_RecGenHists_medPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_medPU", true));
    h_RecGenHists_highPU.reset(new RecoGenHists_xcone(ctx, "RecGenHists_highPU", true));
    h_RecGenHists_lowPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_lowPU_noJEC", false));
    h_RecGenHists_medPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_medPU_noJEC", false));
    h_RecGenHists_highPU_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_highPU_noJEC", false));

    // "true" for RecGenHists means to use jets with JEC applied
    h_RecGenHists_GenOnly.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly", true));
    h_RecGenHists_RecOnly.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly", true));
    h_RecGenHists_RecOnly_noJEC.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_noJEC", false));
    h_RecGenHists_Both.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both", true));
  }
  /*********************************************************************************************************************************/ 

}

bool MTopJetPostSelectionModule::process(uhh2::Event& event){

  /***************************  some options ***************************************************************************************************************/ 
  bool scale_ttbar = true;
  double SF_tt = 0.75;

  /***************************  some useful variables ******************************************************************************************************/ 
  if(isTTbar) ttgenprod->process(event);

  /***************************  get jets to write mass *****************************************************************************************************/ 
  std::vector<Jet> rec_hadjets = event.get(h_recjets_had);
  double mass_rec = rec_hadjets.at(0).v4().M();
  event.set(h_mass_rec, mass_rec);

  if(isMC){
    std::vector<Particle> gen_hadjets23 = event.get(h_genjets23_had);
    std::vector<Particle> gen_hadjets33 = event.get(h_genjets33_had);
    double mass_gen23 = gen_hadjets23.at(0).v4().M();
    double mass_gen33 = gen_hadjets33.at(0).v4().M();
    event.set(h_mass_gen23, mass_gen23);
    event.set(h_mass_gen33, mass_gen33);
  }
  else{
    event.set(h_mass_gen23, 0.); // set gen mass to 0 for data
    event.set(h_mass_gen33, 0.); // set gen mass to 0 for data
  }

  /***************************  apply weight *****************************************************************************************************/ 
  if(isTTbar && scale_ttbar) event.weight = SF_tt * event.get(h_weight);
  else event.weight = event.get(h_weight);

  /***************************  SubJet Corrector *****************************************************************************************************/ 
  // JetCorrections->process(event);
  /*************************** test with lower pt cut ********************************************************************************************/ 
  if(pt200_sel->passes(event) && mass_sel->passes(event)) h_XCone_pt200->fill(event);
  if(pt300_sel->passes(event) && mass_sel->passes(event)) h_XCone_pt300->fill(event);

  /*************************** Events have to pass topjet pt > 400 & Mass_jet1 > Mass_jet2 */
  bool passed_recsel;
  if(pt_sel->passes(event) && mass_sel->passes(event)) passed_recsel = true;
  else passed_recsel = false;

  /*************************** Selection again on generator level (data events will not pass gen sel but will be stored if they pass rec sel)  ***/ 
  bool passed_gensel33;
  if(isMC && mass_gensel->passes(event) && pt_gensel->passes(event) ) passed_gensel33 = true;
  else passed_gensel33 = false;
  bool passed_gensel23;
  if(isMC && mass_gensel23->passes(event) && pt_gensel23->passes(event) ) passed_gensel23 = true;
  else passed_gensel23 = false;
 
  /*************************** fill hists with reco sel applied ***********************************************************************************/ 

  if(pt_sel->passes(event)) h_XCone_noMassCut->fill(event);
  if(pt_sel->passes(event)) h_XConeNoJEC_noMassCut->fill(event);

  bool is_matched_sub = false;
  bool is_matched_fat = false;
  if(passed_recsel){
    h_XCone_subjets->fill(event);
    h_XCone_subjets_noJEC->fill(event);
    h_XCone->fill(event);
    if(event.pvs->size() <= 10){
      h_XCone_lowPU_subjets->fill(event);
      h_XCone_lowPU->fill(event);
      if(isMC)h_RecGenHists_lowPU->fill(event);
      if(isMC)h_RecGenHists_lowPU_noJEC->fill(event);
    }
    if(event.pvs->size() > 10 && event.pvs->size() <= 20){
      h_XCone_medPU_subjets->fill(event);
      h_XCone_medPU->fill(event);
      if(isMC)h_RecGenHists_medPU->fill(event);
      if(isMC)h_RecGenHists_medPU_noJEC->fill(event);
    }
    if(event.pvs->size() > 20){
      h_XCone_highPU_subjets->fill(event);
      h_XCone_highPU->fill(event);
      if(isMC)h_RecGenHists_highPU->fill(event);
      if(isMC)h_RecGenHists_highPU_noJEC->fill(event);
    }
    h_XConeNoJEC->fill(event);
    h_MTopJet->fill(event);
    h_Muon->fill(event);
    if(isTTbar){
      is_matched_sub = matched_sub->passes(event);
      is_matched_fat = matched_fat->passes(event);
      if(is_matched_sub) h_XCone_m->fill(event);
      else h_XCone_u->fill(event);
      if(is_matched_fat) h_XCone_m_fat->fill(event);
      else h_XCone_u_fat->fill(event);
    }
    if(isMC){
      h_XCone_GEN_RecOnly->fill(event);
      if(isTTbar) h_GenParticles_RecOnly->fill(event);
      h_RecGenHists_RecOnly->fill(event);
      h_RecGenHists_RecOnly_noJEC->fill(event);
      h_RecGenHists_subjets->fill(event);
      h_RecGenHists_subjets_noJEC->fill(event);
    }
  }

  /*************************** fill hists with gen sel applied *************************************************************************************/ 
  if(passed_gensel33){
    h_XCone_GEN_GenOnly->fill(event);
    if(isTTbar) h_GenParticles_GenOnly->fill(event);
    h_RecGenHists_GenOnly->fill(event);
  } 

  /*************************** fill hists with reco+gen selection applied **************************************************************************/ 
  if(passed_recsel && passed_gensel33){
    h_XCone_GEN_Both->fill(event);
    if(isTTbar) h_GenParticles_Both->fill(event);
    h_RecGenHists_Both->fill(event);
  } 

  /*************************** write bools for passing selections **********************************************************************************/ 
  event.set(h_ttbar, isTTbar);
  event.set(h_matched, is_matched_sub);
  event.set(h_ttbar_SF, SF_tt);
  event.set(h_recsel, passed_recsel);
  event.set(h_gensel23, passed_gensel23);
  event.set(h_gensel33, passed_gensel33);

  
  /*************************** only store events that survive one of the selections ****************************************************************/
  if(!passed_recsel && !passed_gensel33 && !passed_gensel23) return false;
  else return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
