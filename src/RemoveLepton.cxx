#include <UHH2/MTopJet/include/RemoveLepton.h>

RemoveLepton::RemoveLepton(uhh2::Context & ctx, const std::string & name_jet):
h_topjets(ctx.get_handle<std::vector<TopJet>>(name_jet)){}

bool RemoveLepton::process(uhh2::Event & event){
  std::vector<TopJet>   topjets = event.get(h_topjets);
  std::vector<Particle> leptons;
  for(unsigned int i=0; i<event.muons->size(); i++){
    leptons.push_back(event.muons->at(i));
  }
  for(unsigned int i=0; i<event.electrons->size(); i++){
    leptons.push_back(event.electrons->at(i));
  }
  if(topjets.size() == 0 || leptons.size() == 0) return true;

  for(unsigned int i=0; i<topjets.size(); i++){
    std::vector<Jet> new_subjets = topjets[i].subjets();
    for(unsigned int j=0; j<leptons.size(); j++){
      // loop over leptons and subtract v4 from topjet if dR < 1.2
      if(deltaR(topjets[i], leptons[j]) < 1.2){
        LorentzVector new_v4_top = topjets[i].v4() - leptons[j].v4();
        if(new_v4_top.pt() > 5 && deltaR(new_v4_top, topjets[i].v4()) > M_PI/2){
          cout << "Warning: subtracting lepton flipped jet direction" << endl;
          topjets[i].set_v4(LorentzVector());
        }
        else topjets[i].set_v4(new_v4_top);

        // if there is also a subjet within 0.4, find subjet that is closed to
        // the lepton and subtract the lepton from this one subjet only
        double dRmin = 0.4;
        int index = 100;
        for(unsigned int k=0; k<new_subjets.size(); k++){
          if(deltaR(new_subjets[k], leptons[j]) < dRmin){
            dRmin = deltaR(new_subjets[k], leptons[j]);
            index = k;
          }
        }
        if(index != 100){
          LorentzVector new_v4_sub = new_subjets[index].v4() - leptons[j].v4();
          if(new_v4_sub.pt() > 5 && deltaR(new_v4_sub, new_subjets[index].v4()) > M_PI/2){
            cout << "Warning: subtracting lepton flipped jet direction" << endl;
            new_subjets[index].set_v4(LorentzVector());
          }
          else new_subjets[index].set_v4(new_v4_sub);
        }
      }
    }
    topjets[i].set_subjets(new_subjets);
  }
  event.set(h_topjets, topjets);

  return true;
}

/*
.██████  ███████ ███    ██
██       ██      ████   ██
██   ███ █████   ██ ██  ██
██    ██ ██      ██  ██ ██
.██████  ███████ ██   ████
*/


RemoveLeptonGen::RemoveLeptonGen(uhh2::Context & ctx, const std::string & name_jet):
h_topjets(ctx.get_handle<std::vector<GenTopJet>>(name_jet)),
h_ttbargen(ctx.get_handle<TTbarGen>("ttbargen")){}

bool RemoveLeptonGen::process(uhh2::Event & event){
  std::vector<GenTopJet> topjets = event.get(h_topjets);
  GenParticle lepton;
  std::vector<GenParticle> leptons;
  const auto & ttbargen = event.get(h_ttbargen);
  if(ttbargen.IsSemiLeptonicDecay()) {
    lepton = ttbargen.ChargedLepton();
    leptons.push_back(lepton);
  }

  if(topjets.size() == 0 || leptons.size() == 0) return true;

  // for GenTopJet there is no method 'set_subjets'
  // So, you need to create new GenTopJets without subjets and add the
  // corrected subjets in the end
  std::vector<GenTopJet> new_topjets;
  for(unsigned int i=0; i<topjets.size(); i++){
    std::vector<GenJet> new_subjets = topjets[i].subjets();
    GenTopJet new_topjet;
    new_topjet.set_v4(topjets[i].v4());
    // loop over leptons and subtract v4 from topjet if dR < 1.2
    if(deltaR(new_topjet, lepton) < 1.2){
      LorentzVector new_v4_top = new_topjet.v4() - lepton.v4();
      if(new_v4_top.pt() > 5 && deltaR(new_v4_top, new_topjet.v4()) > M_PI/2){
        cout << "Warning: subtracting lepton flipped jet direction" << endl;
        new_topjet.set_v4(LorentzVector());
      }
      else new_topjet.set_v4(new_v4_top);

      // if there is also a subjet within 0.4, find subjet that is closed to
      // the lepton and subtract the lepton from this one subjet only
      double dRmin = 0.4;
      int index = 100;
      for(unsigned int k=0; k<new_subjets.size(); k++){
        if(deltaR(new_subjets[k], lepton) < dRmin){
          dRmin = deltaR(new_subjets[k], lepton);
          index = k;
        }
      }
      if(index != 100){
        LorentzVector new_v4_sub = new_subjets[index].v4() - lepton.v4();
        if(new_v4_sub.pt() > 5 && deltaR(new_v4_sub, new_subjets[index].v4()) > M_PI/2){
          cout << "Warning: subtracting lepton flipped jet direction" << endl;
          new_subjets[index].set_v4(LorentzVector());
        }
        else new_subjets[index].set_v4(new_v4_sub);
      }
    }
    for(auto subjet: new_subjets) new_topjet.add_subjet(subjet);
    new_topjets.push_back(new_topjet);
  }

  // cout << "----------------------------------------------" << endl;
  // for(auto top:new_topjets){
  //   cout << top.pt() << endl;
  //   for(auto sub: top.subjets()) cout << sub.pt() << endl;
  //   cout << "----------------" << endl;
  // }

  event.set(h_topjets, new_topjets);
  return true;
}
