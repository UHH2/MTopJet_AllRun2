#include <UHH2/MTopJet/include/AnalysisOutput.h>

using namespace std;

WriteOutput::WriteOutput(uhh2::Context & ctx):
  h_weight(ctx.declare_event_output<double>("weight")) {}

bool WriteOutput::process(uhh2::Event & event){

  double weight = event.weight;

  event.set(h_weight, weight);

  return true;
}
