#include "UHH2/MTopJet/include/PartonShowerWeight.h"

PartonShowerWeight::PartonShowerWeight(uhh2::Context & ctx, std::string sys){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";

  weightIndex = -1;
  if(sys == "ISRup_sqrt2") weightIndex = 2;
  else if(sys == "ISRup_2") weightIndex = 6;
  else if(sys == "ISRup_4") weightIndex = 10;

  else if(sys == "ISRdown_sqrt2") weightIndex = 4;
  else if(sys == "ISRdown_2") weightIndex = 8;
  else if(sys == "ISRdown_4") weightIndex = 12;

  else if(sys == "FSRup_sqrt2") weightIndex = 3;
  else if(sys == "FSRup_2") weightIndex = 7;
  else if(sys == "FSRup_4") weightIndex = 11;

  else if(sys == "FSRdown_sqrt2") weightIndex = 5;
  else if(sys == "FSRdown_2") weightIndex = 9;
  else if(sys == "FSRdown_4") weightIndex = 13;
  else
    std::cout << "No parton shower variation sepecified -> module will not have an effect" << std::endl;

  if(!is_mc){
    weightIndex = -1;
    std::cout << "No MC sample. Parton shower weights will have no effec on data" << std::endl;
  }
}

bool PartonShowerWeight::process(uhh2::Event & event){

  if( weightIndex == -1) return true;
  if( event.genInfo->weights().size() == 1) {
    std::cout << "no parton shower weights stored. Nothing to be done." << std::endl;
    return true;
  }

  double centralWeight = event.genInfo->weights().at(0);
  double PSweight = event.genInfo->weights().at(weightIndex);

  double factor = PSweight/centralWeight;
  event.weight *= factor;

  return true;
}
