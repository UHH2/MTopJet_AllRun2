#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/GenJet.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/PFParticle.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/common/include/Utils.h"

#include <math.h>
#include <vector>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"

/*
Calling a function needs time. Ist copies all necessary objetcs and moves them
place of the function (not an accurate description!). Small functions need less
time than bigger ones. But the copy process does not speed up. Inline tells the
compiler to place the function inside the code to reduce the time needed.
Only useful for small functions.

Inline functions are defined in the header inside a class or as single functions
*/

//------------------------------------------------------------------------------
inline TLorentzVector lorentz_to_tlorentz(const LorentzVector v4){
  TLorentzVector lorentz_v4;
  double px, py, pz, Energy;
  px = v4.Px();
  py = v4.Py();
  pz = v4.Pz();
  Energy = v4.E();
  lorentz_v4.SetPxPyPzE(px, py, pz, Energy);
  return lorentz_v4;
}

//------------------------------------------------------------------------------
inline LorentzVector tlorentz_to_lorentz(const TLorentzVector v4){
  LorentzVector lorentz_v4;
  double px, py, pz, Energy;
  px = v4.Px();
  py = v4.Py();
  pz = v4.Pz();
  Energy = v4.E();
  lorentz_v4.SetPxPyPzE(px, py, pz, Energy);
  return lorentz_v4;
}

//------------------------------------------------------------------------------
inline TLorentzVector jet_to_tlorentz(Jet jet){
  TLorentzVector lorentz_v4;
  double px, py, pz, Energy;
  px = jet.v4().Px();
  py = jet.v4().Py();
  pz = jet.v4().Pz();
  Energy  = jet.v4().E();
  lorentz_v4.SetPxPyPzE(px, py, pz, Energy);
  return lorentz_v4;
}

//------------------------------------------------------------------------------
inline TopJet jets_to_topjet(std::vector<Jet> jets){

  // Get momentum and energy of jets
  double px=0, py=0, pz=0, Energy=0;
  for(unsigned int i=0; i < jets.size(); ++i){
    px += jets[i].v4().Px();
    py += jets[i].v4().Py();
    pz += jets[i].v4().Pz();
    Energy  += jets[i].v4().E();
  }

  // Set new TLorentzVector for topjet from jets
  TLorentzVector topjet_v4;
  topjet_v4.SetPxPyPzE(px, py, pz, Energy);

  // Set TopJet from TLV
  TopJet topjet;
  topjet.set_pt(topjet_v4.Pt());
  topjet.set_eta(topjet_v4.Eta());
  topjet.set_phi(topjet_v4.Phi());
  topjet.set_energy(topjet_v4.E());

  return topjet;
}

//------------------------------------------------------------------------------
inline TopJet tlorentz_to_topjet(std::vector<TLorentzVector> v4){

  // Get momentum and energy of jets
  double px=0, py=0, pz=0, Energy=0;
  for(unsigned int i=0; i < v4.size(); ++i){
    px += v4[i].Px();
    py += v4[i].Py();
    pz += v4[i].Pz();
    Energy  += v4[i].E();
  }

  // Set new TLorentzVector for topjet from jets
  TLorentzVector topjet_v4;
  topjet_v4.SetPxPyPzE(px, py, pz, Energy);

  // Set TopJet from TLV
  TopJet topjet;
  topjet.set_pt(topjet_v4.Pt());
  topjet.set_eta(topjet_v4.Eta());
  topjet.set_phi(topjet_v4.Phi());
  topjet.set_energy(topjet_v4.E());

  return topjet;
}

//------------------------------------------------------------------------------
inline TopJet tlorentz_to_topjet(TLorentzVector v4){

  // Get momentum and energy of jets
  double px=v4.Px();
  double py=v4.Py();
  double pz=v4.Pz();
  double E=v4.E();

  // Set new TLorentzVector for topjet from jets
  TLorentzVector topjet_v4;
  topjet_v4.SetPxPyPzE(px, py, pz, E);

  // Set TopJet from TLV
  TopJet topjet;
  topjet.set_pt(topjet_v4.Pt());
  topjet.set_eta(topjet_v4.Eta());
  topjet.set_phi(topjet_v4.Phi());
  topjet.set_energy(topjet_v4.E());

  return topjet;
}

// -----------------------------------------------------------------------------
inline TVector3 get_parallel_component(TVector3 v1, TVector3 v2){
  //Get the parallel component of v1 to v2
  // Formula for the vector projection:
  // vec{v1_parallel} = [](vec(v1)*vec(v2))/norm(v2)^2] * vec(v2)
  v2.SetMag(1); // Redefine to unit vector - v2/norm(v2)
  double product = v1 * v2; // Scalar Product
  double mag_v2 = v2.Mag();  // Get magnitued
  // get parralel component of v1 to v2 - the vector v2 is already normilized.
  TVector3 v1_par = (product/mag_v2) * v2;
  return v1_par;
}

// -----------------------------------------------------------------------------
// This function could be easily done in the code itself.
// To save some comments and make the code more obvious, it is done in an extra function
inline TVector3 get_perpendicular_component(TVector3 v, TVector3 v_par){return(v - v_par);}
