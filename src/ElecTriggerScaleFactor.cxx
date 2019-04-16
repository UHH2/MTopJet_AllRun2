#include <UHH2/MTopJet/include/ElecTriggerScaleFactor.h>

#include <limits>
#include <cmath>
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <stdexcept>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

#include "UHH2/core/include/Hists.h"
using namespace uhh2;
using namespace std;

////////////////////////////////////////////////////////

ElectronTriggerSF::ElectronTriggerSF(Context & ctx){
  sysdirection=ctx.get("ElTrigger_variation");
  auto dataset_type = ctx.get("dataset_type");
  is_mc = dataset_type == "MC";
  if (!is_mc) {
    cout << "Warning: ElectronTriggerSF will not have an effect on "
         <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
}

bool ElectronTriggerSF::process(uhh2::Event & event){
  if(!is_mc) return true;

  if(event.electrons->size() == 0) return true;
  Particle elec = event.electrons->at(0);
  double elec_pt=elec.pt();
  double add_unc=0.01;

  if(elec_pt<=140){
    sf=0.9802;
    if(sysdirection == "up")sf=sqrt(pow(0.9834,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.977,2)-pow(add_unc,2));
  }
  else if(elec_pt<=160){
    sf=0.9873;
    if(sysdirection == "up")sf=sqrt(pow(0.9909,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9835,2)-pow(add_unc,2));
  }
  else if(elec_pt<=180){
    sf=0.9867;
    if(sysdirection == "up")sf=sqrt(pow(0.9911,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9818,2)-pow(add_unc,2));
  }
  else if(elec_pt<=200){
    sf =0.9905;
    if(sysdirection == "up")sf=sqrt(pow(0.9959,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9845,2)-pow(add_unc,2));
  }
  else if(elec_pt<=250){
    sf =0.979421;
    if(sysdirection == "up")sf=sqrt(pow(0.9852,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9731,2)-pow(add_unc,2));
  }
  else if(elec_pt<=400){
    sf=0.996993;
    if(sysdirection == "up")sf=sqrt(pow(0.9989,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9939,2)-pow(add_unc,2));

  }
  else{
    sf=1;
    if(sysdirection == "up")sf=sqrt(pow(1.001,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9837,2)-pow(add_unc,2));
  }
  double weight= event.weight *sf;
  event.weight = weight;

  return true;
}
