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

  // cleaners

  // selections
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> pt200_sel;
  std::unique_ptr<uhh2::Selection> pt300_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;
  std::unique_ptr<uhh2::Selection> pt_gensel;
  std::unique_ptr<uhh2::Selection> mass_gensel;
  std::unique_ptr<uhh2::Selection> pt_gensel23;
  std::unique_ptr<uhh2::Selection> mass_gensel23;
  std::unique_ptr<uhh2::Selection> massbin1;
  std::unique_ptr<uhh2::Selection> massbin2;
  std::unique_ptr<uhh2::Selection> massbin3;
  std::unique_ptr<uhh2::Selection> massbin4;
  std::unique_ptr<uhh2::Selection> massbin5;
  std::unique_ptr<uhh2::Selection> massbin6;
  std::unique_ptr<uhh2::Selection> matched;

  // get weight (with all SF and weight applied in previous cycle)
  Event::Handle<double>h_weight;

  // handles for output
  Event::Handle<bool>h_matched;
  Event::Handle<bool>h_recsel;
  Event::Handle<bool>h_gensel23;
  Event::Handle<bool>h_gensel33;
  Event::Handle<bool>h_ttbar;
  Event::Handle<int>h_massbin;
  Event::Handle<double>h_ttbar_SF;
  Event::Handle<double>h_mass_gen33;
  Event::Handle<double>h_mass_gen23;
  Event::Handle<double>h_mass_rec;
  Event::Handle<std::vector<Jet>>h_recjets_had;
  Event::Handle<std::vector<Particle>>h_genjets23_had;
  Event::Handle<std::vector<Particle>>h_genjets33_had;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


  // store Hist collection as member variables
  std::unique_ptr<Hists> h_XCone, h_XCone_m, h_XCone_u, h_XCone_pt200, h_XCone_pt300, h_MTopJet, h_XCone_noMassCut;
  std::unique_ptr<Hists> h_XCone_GEN_RecOnly, h_XCone_GEN_GenOnly, h_XCone_GEN_Both;
  std::unique_ptr<Hists> h_RecGenHists_GenOnly0, h_RecGenHists_GenOnly1, h_RecGenHists_GenOnly2, h_RecGenHists_GenOnly3, h_RecGenHists_GenOnly4, h_RecGenHists_GenOnly5, h_RecGenHists_GenOnly6;
  std::unique_ptr<Hists> h_RecGenHists_RecOnly0, h_RecGenHists_RecOnly1, h_RecGenHists_RecOnly2, h_RecGenHists_RecOnly3, h_RecGenHists_RecOnly4, h_RecGenHists_RecOnly5, h_RecGenHists_RecOnly6;
  std::unique_ptr<Hists> h_RecGenHists_Both0, h_RecGenHists_Both1, h_RecGenHists_Both2, h_RecGenHists_Both3, h_RecGenHists_Both4, h_RecGenHists_Both5, h_RecGenHists_Both6;
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
  h_massbin = ctx.declare_event_output<int>("massbin");
  h_ttbar_SF = ctx.declare_event_output<double>("TTbar_SF");
  h_mass_gen23 = ctx.declare_event_output<double>("Mass_Gen23");
  h_mass_gen33 = ctx.declare_event_output<double>("Mass_Gen33");
  h_mass_rec = ctx.declare_event_output<double>("Mass_Rec");
  h_recjets_had = ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
  if(isMC) h_genjets23_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone23_had_Combined");
  if(isMC) h_genjets33_had = ctx.get_handle<std::vector<Particle>>("GEN_XCone33_had_Combined");

  /*************************** Setup Selections **********************************************************************************/ 

  // RECO Selection
  pt_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined", 400));
  pt300_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined", 300));
  pt200_sel.reset(new LeadingRecoJetPT(ctx, "XCone33_had_Combined", 200));
  mass_sel.reset(new MassCutXCone(ctx));

  // GEN Selection
  pt_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 400));
  mass_gensel.reset(new MassCut_gen(ctx, "GEN_XCone33_had_Combined", "GEN_XCone33_lep_Combined"));
  pt_gensel23.reset(new LeadingJetPT_gen(ctx, "GEN_XCone23_had_Combined", 400));
  mass_gensel23.reset(new MassCut_gen(ctx, "GEN_XCone23_had_Combined", "GEN_XCone23_lep_Combined"));

  // Selection for Mass binning
  massbin1.reset(new MassCutGen_XCone(ctx,   0, 100));
  massbin2.reset(new MassCutGen_XCone(ctx, 100, 150));
  massbin3.reset(new MassCutGen_XCone(ctx, 150, 200));
  massbin4.reset(new MassCutGen_XCone(ctx, 200, 250));
  massbin5.reset(new MassCutGen_XCone(ctx, 250, 300));
  massbin6.reset(new MassCutGen_XCone(ctx, 300, 500));

  // Selection for matching reco jets to gen particles
  if(isTTbar) matched.reset(new Matching_XCone33(ctx));

  /*************************** Set up Hists classes **********************************************************************************/ 

  h_XCone.reset(new RecoHists_xcone(ctx, "XCone"));
  if(isTTbar) h_XCone_m.reset(new RecoHists_xcone(ctx, "XCone_matched"));
  if(isTTbar) h_XCone_u.reset(new RecoHists_xcone(ctx, "XCone_unmatched"));
  h_XCone_pt200.reset(new RecoHists_xcone(ctx, "XCone_pt200"));
  h_XCone_pt300.reset(new RecoHists_xcone(ctx, "XCone_pt300"));
  h_MTopJet.reset(new MTopJetHists(ctx, "MTopJetHists"));
  h_XCone_noMassCut.reset(new RecoHists_xcone(ctx, "XCone_noMassCut"));

  if(isMC){
    h_XCone_GEN_GenOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_GenOnly"));
    h_XCone_GEN_RecOnly.reset(new GenHists_xcone(ctx, "XCone_GEN_RecOnly"));
    h_XCone_GEN_Both.reset(new GenHists_xcone(ctx, "XCone_GEN_Both"));

    h_GenParticles_GenOnly.reset(new GenHists_particles(ctx, "GenParticles_GenOnly"));
    h_GenParticles_RecOnly.reset(new GenHists_particles(ctx, "GenParticles_RecOnly"));
    h_GenParticles_Both.reset(new GenHists_particles(ctx, "GenParticles_Both"));

    h_RecGenHists_GenOnly0.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_all"));
    h_RecGenHists_GenOnly1.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_000to100"));
    h_RecGenHists_GenOnly2.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_100to150"));
    h_RecGenHists_GenOnly3.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_150to200"));
    h_RecGenHists_GenOnly4.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_200to250"));
    h_RecGenHists_GenOnly5.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_250to300"));
    h_RecGenHists_GenOnly6.reset(new RecoGenHists_xcone(ctx, "RecGenHists_GenOnly_300to500"));

    h_RecGenHists_RecOnly0.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_all"));
    h_RecGenHists_RecOnly1.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_000to100"));
    h_RecGenHists_RecOnly2.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_100to150"));
    h_RecGenHists_RecOnly3.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_150to200"));
    h_RecGenHists_RecOnly4.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_200to250"));
    h_RecGenHists_RecOnly5.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_250to300"));
    h_RecGenHists_RecOnly6.reset(new RecoGenHists_xcone(ctx, "RecGenHists_RecOnly_300to500"));

    h_RecGenHists_Both0.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_all"));
    h_RecGenHists_Both1.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_000to100"));
    h_RecGenHists_Both2.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_100to150"));
    h_RecGenHists_Both3.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_150to200"));
    h_RecGenHists_Both4.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_200to250"));
    h_RecGenHists_Both5.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_250to300"));
    h_RecGenHists_Both6.reset(new RecoGenHists_xcone(ctx, "RecGenHists_Both_300to500"));
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
  /*************************** get massbin (independent from selection) and write it (for purity and stability)  **********************************/ 
  int massbin = 0;
  if(isMC){
    if(massbin1->passes(event)) massbin = 1;
    if(massbin2->passes(event)) massbin = 2;
    if(massbin3->passes(event)) massbin = 3;
    if(massbin4->passes(event)) massbin = 4;
    if(massbin5->passes(event)) massbin = 5;
    if(massbin6->passes(event)) massbin = 6;
  }
  event.set(h_massbin, massbin);


  /*************************** fill hists with reco sel applied ***********************************************************************************/ 
  bool is_matched = false;
  if(passed_recsel){
    h_XCone->fill(event);
    h_MTopJet->fill(event);
    if(isTTbar){
      is_matched = matched->passes(event);
      if(is_matched) h_XCone_m->fill(event);
      else h_XCone_u->fill(event);
    }
    if(!mass_sel->passes(event)){
      h_XCone_noMassCut->fill(event);
    }
    if(isMC){
      h_XCone_GEN_RecOnly->fill(event);
      h_XCone_GEN_RecOnly->fill(event);
      if(isTTbar) h_GenParticles_RecOnly->fill(event);
      h_RecGenHists_RecOnly0->fill(event);
      if(massbin1->passes(event)) h_RecGenHists_RecOnly1->fill(event);
      if(massbin2->passes(event)) h_RecGenHists_RecOnly2->fill(event);
      if(massbin3->passes(event)) h_RecGenHists_RecOnly3->fill(event);
      if(massbin4->passes(event)) h_RecGenHists_RecOnly4->fill(event);
      if(massbin5->passes(event)) h_RecGenHists_RecOnly5->fill(event);
      if(massbin6->passes(event)) h_RecGenHists_RecOnly6->fill(event);
    }
  }

  /*************************** fill hists with gen sel applied *************************************************************************************/ 
  if(passed_gensel33){
    h_XCone_GEN_GenOnly->fill(event);
    h_XCone_GEN_GenOnly->fill(event);
    if(isTTbar) h_GenParticles_GenOnly->fill(event);
    h_RecGenHists_GenOnly0->fill(event);
    if(massbin1->passes(event)) h_RecGenHists_GenOnly1->fill(event);
    if(massbin2->passes(event)) h_RecGenHists_GenOnly2->fill(event);
    if(massbin3->passes(event)) h_RecGenHists_GenOnly3->fill(event);
    if(massbin4->passes(event)) h_RecGenHists_GenOnly4->fill(event);
    if(massbin5->passes(event)) h_RecGenHists_GenOnly5->fill(event);
    if(massbin6->passes(event)) h_RecGenHists_GenOnly6->fill(event);
  } 

  /*************************** fill hists with reco+gen selection applied **************************************************************************/ 
  if(passed_recsel && passed_gensel33){
    h_XCone_GEN_Both->fill(event);
    h_XCone_GEN_Both->fill(event);
    if(isTTbar) h_GenParticles_Both->fill(event);
    h_RecGenHists_Both0->fill(event);
    if(massbin1->passes(event)) h_RecGenHists_Both1->fill(event);
    if(massbin2->passes(event)) h_RecGenHists_Both2->fill(event);
    if(massbin3->passes(event)) h_RecGenHists_Both3->fill(event);
    if(massbin4->passes(event)) h_RecGenHists_Both4->fill(event);
    if(massbin5->passes(event)) h_RecGenHists_Both5->fill(event);
    if(massbin6->passes(event)) h_RecGenHists_Both6->fill(event);
  } 

  /*************************** write bools for passing selections **********************************************************************************/ 
  event.set(h_ttbar, isTTbar);
  event.set(h_matched, is_matched);
  event.set(h_ttbar_SF, SF_tt);
  event.set(h_recsel, passed_recsel);
  event.set(h_gensel23, passed_gensel23);
  event.set(h_gensel33, passed_gensel33);

return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetPostSelectionModule)
