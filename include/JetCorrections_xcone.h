#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Jet.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/core/include/LorentzVector.h>

using namespace uhh2;
using namespace std;



std::vector<TopJet> set_JEC_factor(std::vector<TopJet> jets);

class JetCorrections_xcone : public uhh2::AnalysisModule {

 public:
  JetCorrections_xcone();
  void init(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<GenericSubJetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  std::unique_ptr<GenericJetResolutionSmearer> JER_Smearer;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;
  bool isMC;
  Event::Handle<std::vector<TopJet>>h_topjets;
  Event::Handle<std::vector<GenTopJet>>h_GENtopjets;
};
