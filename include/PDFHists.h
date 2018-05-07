#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include <math.h>
#include <sstream>
#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

class PDFHists: public uhh2::Hists {
  public:
    PDFHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

  protected:
    std::vector< std::string > hist_names1;
    std::vector< std::string > hist_names2;
    std::vector< std::string > hist_names3;
    bool isTTbar;
    uhh2::Event::Handle<std::vector<Jet>>h_hadjets;

};
