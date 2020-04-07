#include <UHH2/MTopJet/include/Vector_utils.h>

//------------------------------------------------------------------------------
TLorentzVector lorentz_to_tlorentz(const LorentzVector v4){
  TLorentzVector lorentz_v4;
  double px, py, pz, Energy;
  px = v4.Px();
  py = v4.Py();
  pz = v4.Pz();
  Energy = v4.E();
  lorentz_v4.SetPxPyPzE(px, py, pz, Energy);
  return lorentz_v4;
}

// -----------------------------------------------------------------------------
TVector3 get_parallel_component(TVector3 v1, TVector3 v2){
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
TVector3 get_perpendicular_component(TVector3 v, TVector3 v_par){return(v - v_par);}
