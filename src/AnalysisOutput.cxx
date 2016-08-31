#include "UHH2/MTopJet/include/AnalysisOutput.h"

using namespace std;

WriteOutput::WriteOutput(uhh2::Context & ctx):
  h_Mass(ctx.declare_event_output<double>("1st_TopJet_Mass")),
  h_Number(ctx.declare_event_output<int>("Top_Jet_Number")) {}

bool WriteOutput::process(uhh2::Event & event){

  int number = event.topjets->size();
  double mass = -1;
  if((event.topjets->size())>0){
    const Particle* Top1 = &event.topjets->at(0);
    mass = Top1->v4().M();
  }

  event.set(h_Mass, mass);
  event.set(h_Number, number);

  return true;
}
