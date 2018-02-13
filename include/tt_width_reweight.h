#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>
#include <iostream>

using namespace uhh2;
using namespace std;

class tt_width_reweight{
 public:

  explicit tt_width_reweight(uhh2::Context &, double);
  bool process(uhh2::Event & );
  double get_factor(uhh2::Event & ); 

 protected:
  uhh2::Event::Handle<TTbarGen>h_ttbargen;


 private:
  double Breit_Wigner(double width, double mass, double energy);
  double sm_mass = 172.5;
  double sm_width = 1.324;
  double width_factor;
};
