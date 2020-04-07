#pragma once
#include "UHH2/core/include/LorentzVector.h"

#include <math.h>
#include <vector>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"

TLorentzVector lorentz_to_tlorentz(const LorentzVector v4); // Lorentz to TLorentzVector
TVector3 get_parallel_component(TVector3 v1, TVector3 v2); // Parallel component of v1 w.r.t v2
TVector3 get_perpendicular_component(TVector3 v, TVector3 v_par); // Perpendicular component
