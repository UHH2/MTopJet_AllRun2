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
  void init(uhh2::Context & ctx, const std::string& jet_collection);
  virtual bool process(uhh2::Event & event) override;

 private:
   std::unique_ptr<GenericSubJetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;
  bool isMC;
  Event::Handle<std::vector<TopJet>>h_topjets;

};



// define JER Smearer for XCone Jets
class JER_Smearer : public uhh2::AnalysisModule {

 public:
  explicit JER_Smearer(uhh2::Context&, const std::string& recj="jets", const std::string& genj="genjets", const bool allow_met_smear=true,
                                       const JERSmearing::SFtype1& JER_sf=JERSmearing::SF_13TeV_2016);
  virtual ~JER_Smearer() {}

  virtual bool process(uhh2::Event&) override;

  void apply_smearing(std::vector<Jet>&, const std::vector<Particle>&, LorentzVector&);

 private:
  uhh2::Event::Handle<std::vector<TopJet> >    h_rectopjets_;
  uhh2::Event::Handle<std::vector<GenTopJet> > h_gentopjets_;

  int direction = 0; // -1 = down, +1 = up, 0 = nominal
  JERSmearing::SFtype1 JER_SFs_;
};


