#include "UHH2/MTopJet/include/MTopJetGenSelections.h"


uhh2::TTbarSemilep::TTbarSemilep(uhh2::Context& ctx):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")) {}

bool uhh2::TTbarSemilep::passes(const uhh2::Event& event){

    const auto & ttbargen = event.get(h_ttbargen);
    bool semilep = false;

    if(ttbargen.DecayChannel() == TTbarGen::e_muhad || ttbargen.DecayChannel() == TTbarGen::e_ehad){
      semilep = true;
    }
    return semilep;
}
