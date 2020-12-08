#include <UHH2/MTopJet/include/StoreBJet.h>

StoreBJet::StoreBJet(uhh2::Context & ctx, const std::string & collectionname):
h_bjets(ctx.declare_event_output<std::vector<Jet>>(collectionname)){}

bool StoreBJet::process(uhh2::Event & event){

  vector<Jet> bjets;

  double maxdeepjet = -1;
  int i_maxdeepjet = -1;
  bool bjet_found = false;
  for(unsigned int i=0; i<event.jets->size(); i++){
    if(event.jets->at(i).btag_DeepJet() > maxdeepjet){
      maxdeepjet = event.jets->at(i).btag_DeepJet();
      i_maxdeepjet = i;
      bjet_found = true;
    }
  }

  if(bjet_found) bjets.push_back(event.jets->at(i_maxdeepjet));
  event.set(h_bjets, bjets);
  
  return true;
}
