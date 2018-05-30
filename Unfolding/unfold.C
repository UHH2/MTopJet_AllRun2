#include "unfold.h"

using namespace std;


unfolding::unfolding(TH1D* input, vector<TH1D*> backgrounds,  vector<TString>bgr_name, TH1D* signal, TH2* migration_matrix, vector< vector<TH2*> > sys_matrix, vector< vector<TString> > sys_name, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, bool do_lcurve, int nscan){
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

  // handle backgrounds
  for(unsigned int i=0; i<bgr_name.size(); i++){
    double scale_error = 0;
    if(bgr_name[i] == "WJets") scale_error = 0.19;
    else if(bgr_name[i] == "SingleTop") scale_error = 0.23;
    else if(bgr_name[i] == "other") scale_error = 1.0;
    else cout << "!!!!!ATTENTION!!!!! Background not known (from unfold.C)" << endl;
    unfold.SubtractBackground(backgrounds[i], bgr_name[i], 1.0, scale_error);
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
  output = unfold.GetOutput("",0,"measurement_gen","pt[C]",kTRUE);
  output_all = unfold.GetOutput("",0,0,0,kFALSE);
  CorM = unfold.GetRhoIJtotal("", 0, "measurement_gen","pt[C]",kTRUE);
  ProbM = unfold.GetProbabilityMatrix("", 0, kTRUE);

  // Statistical uncertainties of input distribution
  CovInputStat = unfold.GetEmatrixInput("", 0, "measurement_gen","pt[C]", kTRUE);
  // Statistical ucertainties of matrix
  CovMatrixStat = unfold.GetEmatrixSysUncorr("", 0, "measurement_gen","pt[C]", kTRUE);

  for(unsigned int i=0; i<bgr_name.size(); i++){
    // Statistical uncertainties from Background
    CovBgrStat.push_back(unfold.GetEmatrixSysBackgroundUncorr(bgr_name[i], "", 0, "measurement_gen","pt[C]", kTRUE));
    // Sys uncertainties from background scale (delta)
    BgrDelta.push_back(unfold.GetDeltaSysBackgroundScale(bgr_name[i], "",0,"measurement_gen","pt[C]",kTRUE));
  }
  for(unsigned int i=0; i<BgrDelta.size(); i++){
    // Sys uncertainties from background scale (cov matrix)
    CovBgrScale.push_back(CreateCovMatrixFromDelta(BgrDelta[i]));
  }


  // treat sys uncertainties
  for(unsigned int i=0; i<sys_name.size(); i++){
    vector<TH1*> dummy1;
    SysDelta.push_back(dummy1);
    for(unsigned int j=0; j<sys_name[i].size(); j++){
      unfold.AddSysError(sys_matrix[i][j], sys_name[i][j], TUnfold::kHistMapOutputHoriz, TUnfoldDensity::kSysErrModeMatrix);
      SysDelta[i].push_back(unfold.GetDeltaSysSource(sys_name[i][j], "",0,"measurement_gen","pt[C]",kTRUE));
    }
  }
  for(unsigned int i=0; i<SysDelta.size(); i++){
    vector<TH2*> dummy2;
    SysCov.push_back(dummy2);
    for(unsigned int j=0; j<SysDelta[i].size(); j++){
      SysCov[i].push_back(CreateCovMatrixFromDelta(SysDelta[i][j]));
    }
  }

  // CovTotal = unfold.GetEmatrixTotal("", 0, "measurement_gen","pt[C]", kTRUE);



}


TH1* unfolding::get_output(bool merged){
  if(merged) return output;
  else return output_all;
}

TH2* unfolding::get_prob_matrix(){
  return ProbM;
}

TH2* unfolding::GetInputStatCov(){
  return CovInputStat;
}

TH2* unfolding::GetMatrixStatCov(){
  return CovMatrixStat;
}

vector<TH2*> unfolding::GetBgrStatCov(){
  return CovBgrStat;
}

vector<TH2*> unfolding::GetBgrScaleCov(){
  return CovBgrScale;
}

// TH2* unfolding::GetTotalCov(){
//   return CovTotal;
// }

vector< vector<TH2*> > unfolding::GetSysCov(){
  return SysCov;
}

vector< vector<TH1*> > unfolding::get_sys_delta(){
  return SysDelta;
}

TH2* unfolding::get_cor_matrix(){
  return CorM;
}

vector<TH1*> unfolding::get_bgr_delta(){
  return BgrDelta;
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
  TH2* cov = (TH2*) CovInputStat->Clone();
  cov->Reset();
  int nbins = delta->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    for(int j=1; j<=nbins; j++){
      double entry = delta->GetBinContent(i) * delta->GetBinContent(j);
      cov->SetBinContent(i,j,entry);
    }
  }
  return cov;
}
