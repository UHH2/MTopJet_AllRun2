#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>
#include <vector>
#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

class CorrectionHists_subjets: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  CorrectionHists_subjets(uhh2::Context & ctx, const std::string & dirname,  const std::string & type);

    virtual void fill(const uhh2::Event & ev) override;

protected:

    /* std::vector<TH1F*> pt_reso; */
    std::vector<std::vector<TH1F*>> pt_reso;
    std::vector<std::vector<TH1F*>> pt_rec;
    std::vector<std::vector<TH1F*>> pt_gen;
    std::vector<std::vector<TH1F*>> event_count;
    std::vector<std::vector<TH2F*>> pt_eta;

    uhh2::Event::Handle<std::vector<TopJet>>h_recjets;
    uhh2::Event::Handle<std::vector<TopJet>>h_recjets_noJEC;
    uhh2::Event::Handle<std::vector<TopJet>>h_hadjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_genjets;
    uhh2::Event::Handle<std::vector<GenTopJet>>h_hadgenjets;

    std::vector<int> ptgen_binning;
    std::vector<int> ptrec_binning;
    std::vector<double> etarec_binning;

};
