#include <UHH2/MTopJet/include/tt_width_reweight.h>


tt_width_reweight::tt_width_reweight(uhh2::Context & ctx, double width_factor_):
  h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")),
  width_factor(width_factor_)
{}

bool tt_width_reweight::process(uhh2::Event & event){
  const auto & ttbargen = event.get(h_ttbargen);
  GenParticle tophad, toplep;
  if(ttbargen.IsTopHadronicDecay()) tophad = ttbargen.Top();
  else tophad = ttbargen.Antitop();

  // nur had. top betrachten?

  double mass = tophad.v4().M();

  double factor = Breit_Wigner(width_factor* sm_width, sm_mass, mass) / Breit_Wigner(sm_width, sm_mass, mass);
  double new_weight = event.weight * factor;
  event.weight = new_weight;
  return true;
}

double tt_width_reweight::get_factor(uhh2::Event & event){
  const auto & ttbargen = event.get(h_ttbargen);
  GenParticle tophad, toplep;
  if(ttbargen.IsTopHadronicDecay()) tophad = ttbargen.Top();
  else tophad = ttbargen.Antitop();
  double mass = tophad.v4().M();
  double factor = Breit_Wigner(width_factor* sm_width, sm_mass, mass) / Breit_Wigner(sm_width, sm_mass, mass);
  return factor;
}


double tt_width_reweight::Breit_Wigner(double width, double peak, double mass){
  double gamma = sqrt(peak*peak*(peak*peak+width*width));
  double k = (2*sqrt(2)*peak*width*gamma) / (M_PI*sqrt(peak*peak+gamma));
  double value = k/((mass*mass-peak*peak)*(mass*mass-peak*peak)+peak*peak*width*width);
  return value;
}
