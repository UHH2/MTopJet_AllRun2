#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/core/include/Utils.h>

#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

using namespace std;


namespace uhh2 {

  class SubjetQuality_gen : public Selection {

  public:
    explicit SubjetQuality_gen(Context&, const std::string &, float, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_jets;
    float ptmin, etamax;
  };

  class GenMuonCount : public Selection {

  public:
    explicit GenMuonCount(Context&);
    virtual bool passes(const Event&) override;

  protected:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;

  };

  class GenElecCount : public Selection {

  public:
    explicit GenElecCount(Context&);
    virtual bool passes(const Event&) override;

  protected:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class GenMuonSel : public Selection {

  public:
    explicit GenMuonSel(Context&, double);
    virtual bool passes(const Event&) override;

  protected:
    double ptmin;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class GenElecSel : public Selection {

  public:
    explicit GenElecSel(Context&, double);
    virtual bool passes(const Event&) override;

  protected:
    double ptmin;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class TTbarSemilep : public Selection {

  public:
    explicit TTbarSemilep(Context&);
    virtual bool passes(const Event&) override;

  protected:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;

  };

  class TTbarSemilep_herwig : public Selection {

  public:
    explicit TTbarSemilep_herwig(Context&);
    virtual bool passes(const Event&) override;
  };

  class TopHadPT : public Selection {

  public:
    explicit TopHadPT(Context&, float ptmin=300);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float ptmin_;

  };

  class Matching : public Selection {

  public:
    explicit Matching(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

  class Matching_HOTVR : public Selection {

  public:
    explicit Matching_HOTVR(Context&, const std::string &, double);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float rho_;
  };

  class Matching_XCone : public Selection {

  public:
    explicit Matching_XCone(Context&, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class Matching_XCone23 : public Selection {

  public:
    explicit Matching_XCone23(Context&, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class Matching_XCone33 : public Selection {

  public:
    explicit Matching_XCone33(Context&, bool);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_fatjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    bool subjet_matching;
  };

  class Matching_XCone33GEN : public Selection {

  public:
    explicit Matching_XCone33GEN(Context&, const std::string &, bool);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_fatjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    bool subjet_matching;
  };

  class Matching_XCone_botlep_lep : public Selection {

  public:
    explicit Matching_XCone_botlep_lep(Context&, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class Matching_top : public Selection {

  public:
    explicit Matching_top(Context&, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };


  class LeadingJetPT : public Selection {

  public:
    explicit LeadingJetPT(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    float ptcut_;
  };

  class LeadingJetPT_gen : public Selection {

  public:
    explicit LeadingJetPT_gen(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_jets;
    float ptcut_;
  };


  class LeadingTopJetPT : public Selection {

  public:
    explicit LeadingTopJetPT(Context&, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    float ptcut_;
  };

  class MassCut : public Selection {

  public:
    explicit MassCut(Context&, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class MassCut_gen : public Selection {

  public:
    explicit MassCut_gen(Context&, const std::string &, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_hadjets;
    uhh2::Event::Handle<std::vector<GenTopJet>> h_lepjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class MassCut_top : public Selection {

  public:
    explicit MassCut_top(Context&);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
  };

  class DeltaPhiCut : public Selection {

  public:
    explicit DeltaPhiCut(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

  class DeltaRCut : public Selection {

  public:
    explicit DeltaRCut(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

  class DeltaRCut_HOTVR : public Selection {

  public:
    explicit DeltaRCut_HOTVR(Context&, const std::string &, double);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float rho_;
  };


  class DeltaRCut_top : public Selection {

  public:
    explicit DeltaRCut_top(Context&, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

  class NGenJets : public Selection {

  public:
    explicit NGenJets(Context&, const std::string &, float min_pt = 0, float min = 0, float max = 9999);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    float min_pt_;
    float min_;
    float max_;
  };

  class NGenTopJets : public Selection {

  public:
    explicit NGenTopJets(Context&, float min_pt = 0, float min = 0, float max = 9999);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    float min_pt_;
    float min_;
    float max_;
  };

  class GenJetLeptonCleaner : public uhh2::AnalysisModule {

  public:
    explicit GenJetLeptonCleaner(Context&, const std::string &, float);
    virtual bool process(Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

  class GenTopJetLeptonCleaner : public uhh2::AnalysisModule {

  public:
    explicit GenTopJetLeptonCleaner(Context&, float);
    virtual bool process(Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };


  class MassCutGen1 : public Selection{

  public:
    explicit MassCutGen1(Context&, const std::string &, float, float);
    virtual bool passes(const Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    float M_min_;
    float M_max_;
  };

  class MassCutGen_XCone : public Selection{

  public:
    explicit MassCutGen_XCone(Context&, float, float);
    virtual bool passes(const Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_jets;
    float M_min_;
    float M_max_;
  };

  class MassCutGen1_top : public Selection{

  public:
    explicit MassCutGen1_top(Context&, float, float);
    virtual bool passes(const Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;
    float M_min_;
    float M_max_;
  };

}
