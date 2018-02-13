#include "unfold.h"

using namespace std;

unfolding::unfolding(TH1D* input, TH1D* background, TH1D* signal, TH2* migration_matrix, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan){
  cout << "Do unfolding" << endl;
  cout << "   -->  You selected: ";
  if(do_lcurve) cout << "L-Curve Scan" << endl;
  if(!do_lcurve) cout << "Rho Scan" << endl;


 // preserve the area
  TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;

  // choice of regularisation scheme:
  TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  // TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
  // TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;

  // density flags
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;

  // detailed steering for regularisation
  const char *REGULARISATION_DISTRIBUTION=0;
  const char *REGULARISATION_AXISSTEERING="*[B]";

  const char *SCAN_DISTRIBUTION="measurement_gen";
  const char *SCAN_AXISSTEERING=0;

  // define scaling to account for different cross section in data and MC
  double scale = (input->Integral() - background->Integral()) / signal->Integral();
  cout << "   -->  fB = " << scale << endl;

  // turn lcurve on/off
  cout << "   -->  Regularisation Mode: ";
  if(regMode == TUnfold::kRegModeSize) cout<< "Size";
  else if(regMode == TUnfold::kRegModeDerivative) cout<< "Derivative";
  else if(regMode == TUnfold::kRegModeCurvature) cout<< "Curvature";
  cout << endl;


  // set up matrix of migrations
  TUnfoldDensity unfold(migration_matrix,TUnfold::kHistMapOutputHoriz,
                        regMode,constraintMode,densityFlags,
                        binning_gen,binning_rec,
  			REGULARISATION_DISTRIBUTION,
  			REGULARISATION_AXISSTEERING);


  unfold.SetInput(input, scale);
  unfold.SubtractBackground(background, "BG", 1.0, 0.0);

  TSpline *rhoLogTau=0;
  TSpline *logTauX=0,*logTauY=0;

  // L-curve scan
  if(do_lcurve) unfold.ScanLcurve(nscan,0.,0.,&lcurve,&logTauX,&logTauY);

  // rho scan
  if(!do_lcurve) unfold.ScanTau(nscan,0.0001,100.,&rhoLogTau, TUnfoldDensity::kEScanTauRhoAvgSys, SCAN_DISTRIBUTION,SCAN_AXISSTEERING, &lcurve,&logTauX,&logTauY);

  tau=unfold.GetTau();
  double logTau=TMath::Log10(tau);
  lcurveX=logTauX->Eval(logTau);
  lcurveY=logTauY->Eval(logTau);

  // write hists
  output = unfold.GetOutput("output",0,"measurement_gen","pt[C]",kTRUE);
  output_all = unfold.GetOutput("output_all",0,0,0,kFALSE);
  CovM = unfold.GetRhoIJtotal("CovM", 0, "measurement_gen","pt[C]",kTRUE);
  ProbM = unfold.GetProbabilityMatrix("ProbM", 0, kTRUE);

}

TH1* unfolding::get_output(bool merged){
  if(merged) return output;
  else return output_all;
}

TH2* unfolding::get_prob_matrix(){
  return ProbM;
}

TH2* unfolding::get_cov_matrix(){
  return CovM;
}

TGraph* unfolding::get_lcurve(){
  return lcurve;
}

double unfolding::get_best_point(TString xy){
  if(xy == "x" || xy == "X") return lcurveX;
  else return lcurveY;
}

double unfolding::get_tau(){
  return tau;
}
