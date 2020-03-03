#pragma once
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include <math.h>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <random>

using namespace std;
using namespace uhh2;


class CountingEventHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
  CountingEventHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:
  TH1F *Events;

};
