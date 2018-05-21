#include "unfold.h"

using namespace std;


unfolding::unfolding(TH1D* input, vector<TH1D*> backgrounds,  vector<TString>bgr_name, TH1D* signal, TH2* migration_matrix, vector<TH2*>sys_matrix, vector<TString>sys_name, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan){
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
  double bkg_events = 0;
  for(auto bkg: backgrounds){
    bkg_events +=  bkg->Integral();
  }
  double scale = (input->Integral() - bkg_events) / signal->Integral();
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

  for(unsigned int i=0; i<bgr_name.size(); i++){
    unfold.SubtractBackground(backgrounds[i], bgr_name[i], 1.0, 0.0);
  }


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

  // write hists before including sys uncertainties
  output = unfold.GetOutput("output",0,"measurement_gen","pt[C]",kTRUE);
  output_all = unfold.GetOutput("output_all",0,0,0,kFALSE);
  CorM = unfold.GetRhoIJtotal("CorM", 0, "measurement_gen","pt[C]",kTRUE);
  ProbM = unfold.GetProbabilityMatrix("ProbM", 0, kTRUE);
  CovM = unfold.GetEmatrixTotal("CovM", 0, "measurement_gen","pt[C]", kTRUE); // this has to stay, beacuse

  // treat sys uncertainties
  for(unsigned int i=0; i<sys_name.size(); i++){
    unfold.AddSysError(sys_matrix[i], sys_name[i], TUnfold::kHistMapOutputHoriz, TUnfoldDensity::kSysErrModeMatrix);
    SysDelta.push_back(unfold.GetDeltaSysSource(sys_name[i], "delta_"+sys_name[i],0,"measurement_gen","pt[C]",kTRUE));
  }
  for(unsigned int i=0; i<SysDelta.size(); i++){
    cout << " ====================== " << sys_name[i] << " ====================== " << endl;
    SysCov.push_back(CreateCovMatrixFromDelta(SysDelta[i]));
  }

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

TH2* unfolding::get_cor_matrix(){
  return CorM;
}

vector<TH1*> unfolding::get_delta(){
  return SysDelta;
}

vector<TH2*> unfolding::get_sys_matrix(){
  return SysCov;
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

TH2* unfolding::CreateCovMatrixFromDelta(TH1* delta){
  TH2* cov = (TH2*) CovM->Clone();
  cov->Reset();
  int nbins = delta->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    for(int j=1; j<=nbins; j++){
      double entry = delta->GetBinContent(i) * delta->GetBinContent(j);
      cout << i << ", " << j << ", " << entry << endl;
      cov->SetBinContent(i,j,entry);
    }
  }
  return cov;
}
