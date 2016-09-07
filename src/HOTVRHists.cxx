#include "UHH2/MTopJet/include/HOTVRHists.h"
#include "UHH2/MTopJet/include/JetCluster.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/PFParticle.h"

#include "UHH2/common/include/TTbarGen.h"

#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace uhh2;

HOTVRHists::HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const std::string & jetname, double rho): Hists(ctx, dirname){
  // book all histograms here
  Reff_HOTVR = book<TH1F>("Reff_HOTVR", "R_{eff} of leading Jet", 40, 0, 2.0);
  Mass_Reff_HOTVR = book<TH2F>("Mass_Reff_HOTVR", "x=M_Jet1 y=R_eff_2ndJet", 50, 0, 500., 40, 0, 2.0);
  Reff2_HOTVR = book<TH1F>("Reff2_HOTVR", "R_{eff} of 2nd Jet", 40, 0, 2.0);
  Mass2_Reff2_HOTVR = book<TH2F>("Mass2_Reff2_HOTVR", "x=M_Jet2 y=R_eff_2ndJet", 50, 0, 500., 40, 0, 2.0);
  // handle for clustered jets
  h_jets=ctx.get_handle<std::vector<Jet>>(jetname);
  rho_=rho;
}



void HOTVRHists::fill(const Event & event){


  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<Jet> jets = event.get(h_jets);
  TLorentzVector jet1_v4, jet2_v4;
  Jet jet1,jet2;
  if(jets.size()>0) jet1 = jets.at(0);
  if(jets.size()>1) jet2 = jets.at(1);

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------


 
  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of 2 jets and lepton -------------------=-----------------
  //---------------------------------------------------------------------------------------
   if(jets.size() > 1){
     jet1_v4.SetPxPyPzE(jets.at(0).v4().Px(), jets.at(0).v4().Py(), jets.at(0).v4().Pz(), jets.at(0).v4().E());
     jet2_v4.SetPxPyPzE(jets.at(1).v4().Px(), jets.at(1).v4().Py(), jets.at(1).v4().Pz(), jets.at(1).v4().E()); //v4 of first jet
   }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;
  ////


  if(jets.size() > 0){
  double pt1 = jet1_v4.Pt();
  double reff = rho_/pt1;
  if(reff < 0.1) reff = 0.1;
  if(reff > 1.5) reff = 1.5;
  Reff_HOTVR->Fill(reff, weight);
  Mass_Reff_HOTVR->Fill(jet1_v4.M(), reff, weight);
  }

  if(jets.size() > 1){
  double pt2 = jet2_v4.Pt();
  double reff2 = rho_/pt2;
  if(reff2 < 0.1) reff2 = 0.1;
  if(reff2 > 1.5) reff2 = 1.5;
  Reff2_HOTVR->Fill(reff2, weight);
  Mass2_Reff2_HOTVR->Fill(jet2_v4.M(), reff2, weight);
  }
  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------
  jet1_v4.Delete();
  jet2_v4.Delete();
  //---------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------- 

}


