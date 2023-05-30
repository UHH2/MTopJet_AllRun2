#include <iostream>
#include <fstream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>

#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/GenHists_xcone.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/StoreBJet.h>

#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>

#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"


#include <vector>

// ##############################################################################
// ### Definition of class
// ### public: Seen by all functions including class
// ### private: Variables, functions, etc. only seen by class. 
// ###          Outer functions cannot access them! 

class MTopJetGenPostSelectionModule : public ModuleBASE {

public:
  explicit MTopJetGenPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void declare_output(uhh2::Context& ctx);
  void init_handels(uhh2::Context& ctx);
  void init_hists(uhh2::Context& ctx);
  void init_MC_hists(uhh2::Context& ctx);

protected:
  enum lepton { muon, elec };
  lepton channel_;

  // Selection Modules defined in scope of the analysis code
  std::unique_ptr<uhh2::Selection> pt_gensel;

  // Access XCone jets
  Event::Handle<std::vector<GenTopJet>> h_genjets33_had, h_genfatjets;

  // Access events variables stored in input rootfile
  Event::Handle<bool> h_passed_gensel;
  Event::Handle<TTbarGen>h_ttbargen;

  // Handles to store variables in output
  Event::Handle<bool> h_measure_gen;

  // Define Histgrams
  std::unique_ptr<Hists> h_XCone_GEN_noCut, h_XCone_GEN_ptCut;

  // Object construction
  std::unique_ptr<uhh2::AnalysisModule> jetprod_gen;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // global variables for class
  bool year_18;
  bool debug;

};

// ##############################################################################
// ### Constructor; for each xml
// ### Initilize all objects

// Get handles; Variables stored for each event in rootfile e.g. muon pt
void MTopJetGenPostSelectionModule::init_handels(uhh2::Context& ctx){
  // boolean if previous selection on gen level to select ttbar events passed 
  // See MTopJetSelection or MTopJetPreSelection
  h_passed_gensel = ctx.get_handle<bool>("passed_gensel_2");
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");
}

// Store handles in output rootfile
void MTopJetGenPostSelectionModule::declare_output(uhh2::Context& ctx){
  h_measure_gen = ctx.declare_event_output<bool>("passed_measurement_gen");
}


// Histograms only stored for MC e.g. gen info is not accesible in data.
void MTopJetGenPostSelectionModule::init_MC_hists(uhh2::Context& ctx){
  // GenHists_xcone definiton also in MTopJet/src folder
  // Give name e.g. XCone_GEN_noCut as it should be stored in output file
  // Hist classes will be stored as `directories` and contains the defined histograms
  h_XCone_GEN_noCut.reset(new GenHists_xcone(ctx, "XCone_GEN_noCut"));
  h_XCone_GEN_ptCut.reset(new GenHists_xcone(ctx, "XCone_GEN_ptCut"));
}

// Histograms for MC and data, so only rec level
void MTopJetGenPostSelectionModule::init_hists(uhh2::Context& ctx){

}

// Main constructor. add previous functions in here.
// Not necessary, but when constructor gets to large (too many lines),
// the compile complains and crashes.
MTopJetGenPostSelectionModule::MTopJetGenPostSelectionModule(uhh2::Context& ctx){

  // /////////////
  // Set variables

  // Only set true if you want to test the code and/or the code crashes. 
  // Only true if run locally and not via HTCondor.
  // Set in xml file, so one does not has to compile everytime
  debug = string2bool(ctx.get("Debug","false")); // look for Debug, expect false if not found

  if(debug) cout << "--- Start Module - CTX ---" << endl;

  // Example, not important now, but when you analyse multiple years.
  year_18 = false;

  // //////////////////
  // construct gen jets 
  if(debug) cout << "\t--- Construct XCone jets" << endl;

  // gat ttbar gen info; Important for clustering, ignore for now
  const std::string ttbar_gen_label("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  h_ttbargen=ctx.get_handle<TTbarGen>("ttbargen");

  jetprod_gen.reset(new CombineXCone33_gen(ctx, true)); // keep true for now, since you look only at ttbar

  // //////////////
  // Access genjets
  if(debug) cout << "\t--- Access genejts" << endl;

  // fatjets: XCone jets with R=1.2 (both hadronic and leptonic)
  // `combined jets`: Combined subjets (R=0.4) within fatjet; Used for measurement
  // 33 means both fatjets have 3 subjets. 
  h_genfatjets = ctx.get_handle<std::vector<GenTopJet>>("genXCone33TopJets");

  // GEN_XCone33_had_Combined produces in CombineXCone33_gen, so define afterwards
  // Only hadronic now
  h_genjets33_had = ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");

  // /////////////////
  // Define Selections
  if(debug) cout << "\t--- Define selections" << endl;

  // Defined in GenSelections.cxx
  pt_gensel.reset(new LeadingJetPT_gen(ctx, "GEN_XCone33_had_Combined", 400));

  // /////////////////////////
  // Initiate input and output
  if(debug) cout << "\t--- Initiate input and output" << endl;
  init_handels(ctx);
  init_MC_hists(ctx);

  // /////////////
  // Handle output
  if(debug) cout << "\t--- Declare output" << endl;

  ctx.undeclare_all_event_output(); // Delete everything to save space
  declare_output(ctx); // Declare new output

}

// ##############################################################################
// ### Process; Event loop
// ### Selection, Corrections, fill of hists for each event               

bool MTopJetGenPostSelectionModule::process(uhh2::Event& event){

  if(debug) cout << "--- New Event ---" << endl;

  // /////////////////
  // Construct objects
  if(debug) cout << "\t--- Construct Obects" << endl;

  ttgenprod->process(event);
  jetprod_gen->process(event);

  // /////////////////////////
  // Get Objects from rootfile
  if(debug) cout << "\t--- Get handles" << endl;

  std::vector<GenTopJet> gen_hadjets33 = event.get(h_genjets33_had);

  // ///////////////////////////////
  // Run selection for single events
  if(debug) cout << "\t--- Run Selection" << endl;

  // fill all events
  h_XCone_GEN_noCut->fill(event);

  if(debug) cout << "\t\t--- before passed gensel" << endl;
  bool passed_gensel33 = event.get(h_passed_gensel);
  if(debug) cout << "\t\t--- after passed gensel" << endl;
  bool pass_selection_gen;

  // Check of events passes selection
  // Use function with pointer: pointer->func()
  // Pointer points to memory location of object;
  //    --- Faster and more efficient, because no object has to be copied
  //    --- Be careful, since function changes object you give as an arguments
  //        --- p1->func(p2, o1); p2 can also be changed. o1 will be copied and be a new object (p=pointer, o=object) 
  pass_selection_gen = passed_gensel33 && pt_gensel->passes(event); // Pass both: logical and && (in c++)

  // fill events after selection
  if(pass_selection_gen) h_XCone_GEN_ptCut->fill(event);

  // Fill handles for output
  if(debug) cout << "\t--- Set output" << endl;
  event.set(h_measure_gen, pass_selection_gen); // either 0 or 1, but filled for all events

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetGenPostSelectionModule)
