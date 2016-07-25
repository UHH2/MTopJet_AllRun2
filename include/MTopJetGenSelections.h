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


namespace uhh2 {

  class TTbarSemilep : public Selection { 

  public: 
    explicit TTbarSemilep(Context&);
    virtual bool passes(const Event&) override;

  protected:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;

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

 class LeadingJetPT : public Selection { 

  public: 
    explicit LeadingJetPT(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
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

 class DeltaRCut : public Selection { 

  public: 
    explicit DeltaRCut(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    float jetradius_;
  };

 class NGenJets : public Selection { 

  public: 
    explicit NGenJets(Context&, const std::string &, float min = 0, float max = 9999);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
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

}
