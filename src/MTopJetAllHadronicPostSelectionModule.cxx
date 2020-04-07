#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>

#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/MCWeight.h>

#include <UHH2/MTopJet/include/MTopJetHists.h>
#include <UHH2/MTopJet/include/CombineXCone.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>
#include <UHH2/MTopJet/include/RecoSelections.h>
#include <UHH2/MTopJet/include/RecoSelections_topjet.h>
#include <UHH2/MTopJet/include/GenSelections.h>
#include <UHH2/MTopJet/include/RecoHists_xcone.h>
#include <UHH2/MTopJet/include/RecoHists_puppi.h>
#include <UHH2/MTopJet/include/RecoHists_topjet.h>
#include <UHH2/MTopJet/include/PDFHists.h>

#include <UHH2/MTopJet/include/GenHists_xcone.h>
#include <UHH2/MTopJet/include/RecoGenHists_xcone.h>
#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/AnalysisOutput.h>
#include <UHH2/MTopJet/include/SubjetHists_xcone.h>
#include <UHH2/MTopJet/include/RecoGenHists_allHad.h>
#include <UHH2/MTopJet/include/RecoGenHists_ak4.h>
#include <UHH2/MTopJet/include/CorrectionHists_allHad.h>
#include <UHH2/MTopJet/include/CorrectionFactor.h>
#include <UHH2/MTopJet/include/tt_width_reweight.h>
#include <UHH2/MTopJet/include/JetCorrections_xcone.h>
#include <UHH2/MTopJet/include/ControlHists.h>

#include <vector>

class MTopJetAllHadronicPostSelectionModule : public ModuleBASE {

public:
  explicit MTopJetAllHadronicPostSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:

  std::unique_ptr<CommonModules> common;
  std::unique_ptr<uhh2::Selection> njet_had;
  std::unique_ptr<uhh2::Selection> njet_lep;
  std::unique_ptr<uhh2::Selection> pt_sel1;
  std::unique_ptr<uhh2::Selection> pt_sel2;
  // use top pt reweight module instead of derived sf?
  Event::Handle<std::vector<TopJet>>h_recjets_had;
  Event::Handle<std::vector<GenTopJet>>h_genjets33_had;

  std::unique_ptr<AnalysisModule> PUreweight;
  std::unique_ptr<AnalysisModule> lumiweight;

  // store Hist collection as member variables
  std::unique_ptr<Hists> h_CorrectionHists;


  std::unique_ptr<Hists> h_XCone_cor, h_XCone_jec, h_XCone_raw;
  std::unique_ptr<Hists> h_XCone_cor_SF, h_XCone_jec_SF, h_XCone_raw_SF;
  std::unique_ptr<Hists> h_XCone_cor_subjets, h_XCone_jec_subjets, h_XCone_raw_subjets;
  std::unique_ptr<Hists> h_XCone_cor_subjets_SF, h_XCone_jec_subjets_SF, h_XCone_raw_subjets_SF;

  std::unique_ptr<Hists> h_RecGenHists_rec0_gen400;
  std::unique_ptr<Hists> h_RecGenHists_rec200_gen0;
  std::unique_ptr<Hists> h_RecGenHists_rec300_gen0;
  std::unique_ptr<Hists> h_RecGenHists_rec400_gen0;

  std::unique_ptr<Hists> h_RecGenHists_allHad_jet1, h_RecGenHists_allHad_jet2, h_RecGenHists_allHad_jet3;
  std::unique_ptr<Hists> h_RecGenHists_allHad, h_RecGenHists_allHad_noJEC, h_RecGenHists_allHad_corrected;
  std::unique_ptr<Hists> h_RecGenHists_allHad_lowPU, h_RecGenHists_allHad_medPU, h_RecGenHists_allHad_highPU;
  std::unique_ptr<Hists> h_RecGenHists_allHad_noJEC_lowPU, h_RecGenHists_allHad_noJEC_medPU, h_RecGenHists_allHad_noJEC_highPU;

  bool isMC;    //define here to use it in "process" part
  bool isTTbar; //define here to use it in "process" part

  std::unique_ptr<uhh2::AnalysisModule> BTagScaleFactors;


};

MTopJetAllHadronicPostSelectionModule::MTopJetAllHadronicPostSelectionModule(uhh2::Context& ctx){

  /*
  .██████ ████████ ██   ██
  ██         ██     ██ ██
  ██         ██      ███
  ██         ██     ██ ██
  .██████    ██    ██   ██
  */
  /*************************** CONFIGURATION **********************************************************************************/
  TString dataset_version = (TString) ctx.get("dataset_version");
  if(dataset_version.Contains("TTbar")) isTTbar = true;
  else  isTTbar = false;
  // PU reweighting
  PUreweight.reset(new MCPileupReweight(ctx, "central"));
  lumiweight.reset(new MCLumiWeight(ctx));
  // ttbar gen
  const std::string ttbar_gen_label("ttbargen");
  h_recjets_had = ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");
  h_genjets33_had = ctx.get_handle<std::vector<GenTopJet>>("GEN_XCone33_had_Combined");
  /*************************** Setup Subjet Corrector **********************************************************************************/
  // Correction.reset(new CorrectionFactor(ctx, "xconeCHS_Corrected"));
  /*************************** Setup Selections **********************************************************************************/
  // RECO Selection
  // chose: XCone33_had_Combined_Corrected,  XCone33_had_Combined_noJEC,  XCone33_had_Combined
  const std::string& jet_label_had = "XCone33_had_Combined";
  const std::string& jet_label_lep = "XCone33_lep_Combined";
  njet_had.reset(new NJetXCone(ctx, jet_label_had, 1));
  njet_lep.reset(new NJetXCone(ctx, jet_label_lep, 1));

  pt_sel1.reset(new LeadingRecoJetPT(ctx, jet_label_had, 400));
  pt_sel2.reset(new LeadingRecoJetPT(ctx, jet_label_lep, 400));

  // Scale factors
  BTagScaleFactors.reset(new MCBTagScaleFactor(ctx,BTag::DEEPJET,BTag::WP_TIGHT,"jets","central"));
  /*************************** Set up Hists classes **********************************************************************************/

  // XCone Combined Jet
  h_XCone_raw.reset(new RecoHists_xcone(ctx, "XCone_raw", "raw"));
  h_XCone_cor.reset(new RecoHists_xcone(ctx, "XCone_cor", "cor"));
  h_XCone_jec.reset(new RecoHists_xcone(ctx, "XCone_jec", "jec"));
  h_XCone_raw_SF.reset(new RecoHists_xcone(ctx, "XCone_raw_SF", "raw"));
  h_XCone_cor_SF.reset(new RecoHists_xcone(ctx, "XCone_cor_SF", "cor"));
  h_XCone_jec_SF.reset(new RecoHists_xcone(ctx, "XCone_jec_SF", "jec"));

  // XCone Subjets
  h_XCone_jec_subjets.reset(new SubjetHists_xcone(ctx, "XCone_jec_subjets", "jec"));
  h_XCone_raw_subjets.reset(new SubjetHists_xcone(ctx, "XCone_raw_subjets", "raw"));
  h_XCone_cor_subjets.reset(new SubjetHists_xcone(ctx, "XCone_cor_subjets", "cor"));
  h_XCone_jec_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_jec_subjets_SF", "jec"));
  h_XCone_raw_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_raw_subjets_SF", "raw"));
  h_XCone_cor_subjets_SF.reset(new SubjetHists_xcone(ctx, "XCone_cor_subjets_SF", "cor"));

  h_CorrectionHists.reset(new CorrectionHists_allHad(ctx, "CorrectionHists"));
  h_RecGenHists_allHad_jet1.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_jet1", "jec", "first", 0, 400));
  h_RecGenHists_allHad_jet2.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_jet2", "jec", "second", 0, 400));
  h_RecGenHists_allHad_jet3.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_jet3", "jec", "third", 0, 400));

  h_RecGenHists_allHad.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad", "jec", "all", 0, 400));
  h_RecGenHists_allHad_corrected.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_corrected", "cor", "all", 0, 400));
  h_RecGenHists_allHad_lowPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_lowPU", "jec", "all", 0, 400));
  h_RecGenHists_allHad_medPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_medPU", "jec", "all", 0, 400));
  h_RecGenHists_allHad_highPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_highPU", "jec", "all", 0, 400));
  h_RecGenHists_allHad_noJEC.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_noJEC", "raw", "all", 0, 400));
  h_RecGenHists_allHad_noJEC_lowPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_noJEC_lowPU", "raw", "all", 0, 400));
  h_RecGenHists_allHad_noJEC_medPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_noJEC_medPU", "raw", "all", 0, 400));
  h_RecGenHists_allHad_noJEC_highPU.reset(new RecoGenHists_allHad(ctx, "RecGenHists_allHad_noJEC_highPU", "raw", "all", 0, 400));

  h_RecGenHists_rec0_gen400.reset(new RecoGenHists_allHad(ctx, "RecGenHists_rec0_gen400", "jec", "all", 0, 400));
  h_RecGenHists_rec200_gen0.reset(new RecoGenHists_allHad(ctx, "RecGenHists_rec200_gen0", "jec", "all", 200, 0));
  h_RecGenHists_rec300_gen0.reset(new RecoGenHists_allHad(ctx, "RecGenHists_rec300_gen0", "jec", "all", 300, 0));
  h_RecGenHists_rec400_gen0.reset(new RecoGenHists_allHad(ctx, "RecGenHists_rec400_gen0", "jec", "all", 400, 0));

  common.reset(new CommonModules());
  common->disable_jersmear();
  common->init(ctx);
  // undeclare event output (jet collections etc) to get small root files
  // ctx.undeclare_all_event_output();
}

bool MTopJetAllHadronicPostSelectionModule::process(uhh2::Event& event){

  /*
  ██████  ██████   ██████   ██████ ███████ ███████ ███████
  ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
  ██████  ██████  ██    ██ ██      █████   ███████ ███████
  ██      ██   ██ ██    ██ ██      ██           ██      ██
  ██      ██   ██  ██████   ██████ ███████ ███████ ███████
  */
  lumiweight->process(event);

  if( !(njet_had->passes(event)) ) return false;
  if( !(njet_lep->passes(event)) ) return false;
  ////
  /***************************  apply weight *****************************************************************************************************/
  // bool scale_ttbar = true;           // match MC and data cross-section (for plots only)?
  // double SF_tt = 0.75;

  /** PU Reweighting *********************/
  PUreweight->process(event);

  /** b-tagging *********************/
  int jetbtagN(0);
  for(const auto& j : *event.jets){
    if(BTag(BTag::DEEPJET, BTag::WP_TIGHT)(j, event)) ++jetbtagN;
  }
  if(jetbtagN < 1) return false;
  BTagScaleFactors->process(event);

  /** pT Selection *********************/
  // if(!pt_sel1->passes(event) && !pt_sel2->passes(event)) return false;

  /*************************** Pile Up bools  ***************************************************************************************************/
  bool lowPU = (event.pvs->size() <= 10);
  bool midPU = (event.pvs->size() > 10 && event.pvs->size() <= 20);
  bool highPU = (event.pvs->size() > 20);

  /*************************** fill hists with reco sel applied ***********************************************************************************/
  if(pt_sel1->passes(event)){
    h_XCone_raw->fill(event);
    h_XCone_jec->fill(event);
    h_XCone_cor->fill(event);

    h_XCone_raw_subjets->fill(event);
    h_XCone_jec_subjets->fill(event);
    h_XCone_cor_subjets->fill(event);

    // if(scale_ttbar) event.weight *= SF_tt;
    // Since ttbar_SF is excluded, the next Hists with _SF are unnecessary
    // h_XCone_raw_SF->fill(event);
    // h_XCone_jec_SF->fill(event);
    // h_XCone_cor_SF->fill(event);
    // h_XCone_raw_subjets_SF->fill(event);
    // h_XCone_jec_subjets_SF->fill(event);
    // h_XCone_cor_subjets_SF->fill(event);

  }
  if(lowPU){
    h_RecGenHists_allHad_lowPU->fill(event);
    h_RecGenHists_allHad_noJEC_lowPU->fill(event);
  }
  if(midPU){
    h_RecGenHists_allHad_medPU->fill(event);
    h_RecGenHists_allHad_noJEC_medPU->fill(event);
  }
  if(highPU){
    h_RecGenHists_allHad_highPU->fill(event);
    h_RecGenHists_allHad_noJEC_highPU->fill(event);
  }

  h_RecGenHists_allHad_jet1->fill(event);
  h_RecGenHists_allHad_jet2->fill(event);
  h_RecGenHists_allHad_jet3->fill(event);

  h_RecGenHists_allHad->fill(event);
  h_RecGenHists_allHad_noJEC->fill(event);
  h_RecGenHists_allHad_corrected->fill(event);
  h_CorrectionHists->fill(event);
  h_RecGenHists_rec0_gen400->fill(event);
  h_RecGenHists_rec200_gen0->fill(event);
  h_RecGenHists_rec300_gen0->fill(event);
  h_RecGenHists_rec400_gen0->fill(event);

  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(MTopJetAllHadronicPostSelectionModule)
