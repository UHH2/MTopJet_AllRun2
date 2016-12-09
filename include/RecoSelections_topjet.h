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

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>




namespace uhh2 {


 class LeadingRecoJetPT_topjet : public Selection { 

  public: 
    explicit LeadingRecoJetPT_topjet(Context&, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
    float ptcut_;
  };


 class MassCutReco_topjet : public Selection { 

  public: 
    explicit MassCutReco_topjet(Context&);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
  };

 class DeltaRCutReco_topjet : public Selection { 

  public: 
    explicit DeltaRCutReco_topjet(Context&, float);
    virtual bool passes(const Event&) override;

  private:
     uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
   float jetradius_;
  };

 class NRecoJets_topjet : public Selection { 

  public: 
    explicit NRecoJets_topjet(Context&, float min_pt = 0, float min = 0, float max = 9999);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
    float min_pt_;
    float min_;
    float max_;
  };

 class RecoJetLeptonCleaner_topjet : public uhh2::AnalysisModule { 

  public: 
    explicit RecoJetLeptonCleaner_topjet(Context&);
    virtual bool process(Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
  };
 

  class MassCut_topjet : public Selection {

   public:
    virtual bool passes(const Event&) override;

  };

  ////

  class deltaRCut_topjet : public Selection {

   public:
    explicit deltaRCut_topjet(float max_deltaR=100);
    virtual bool passes(const Event&) override;

   private:
    float max_deltaR_;
  };

  ////


  class Jet2Cut_topjet : public Selection {
  public:
    explicit Jet2Cut_topjet(float min_pt=0);
    virtual bool passes(const Event&) override;   

   private:
    float min_pt_; 
  };
  ////

}
