#include "unfold.h"

using namespace std;


unfolding::unfolding(TH1D* input, vector<TH1D*> backgrounds,  vector<TString>bgr_name, TH1D* signal, TH2* migration_matrix, vector< vector<TH2*> > sys_matrix, vector< vector<TString> > sys_name, TUnfoldBinning *binning_rec, TUnfoldBinning *binning_gen, TUnfoldBinning *binning_gen20, bool do_lcurve, int nscan, double tau_){
  cout << "Do unfolding" << endl;
  cout << "   -->  You selected: ";
  if(do_lcurve) cout << "L-Curve Scan" << endl;
  if(!do_lcurve) cout << "Rho Scan" << endl;


  // get larger binning
  measurement_gen20 = binning_gen20->FindNode("measurement_gen20");

  /*
   █████  ██████  ███████  █████       ██████  ██████  ███    ██ ███████ ████████ ██████   █████  ██ ███    ██ ████████
  ██   ██ ██   ██ ██      ██   ██     ██      ██    ██ ████   ██ ██         ██    ██   ██ ██   ██ ██ ████   ██    ██
  ███████ ██████  █████   ███████     ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ███████ ██ ██ ██  ██    ██
  ██   ██ ██   ██ ██      ██   ██     ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██   ██ ██ ██  ██ ██    ██
  ██   ██ ██   ██ ███████ ██   ██      ██████  ██████  ██   ████ ███████    ██    ██   ██ ██   ██ ██ ██   ████    ██
  */

  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
  // TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;

  /*
  ██████  ███████  ██████  ██    ██ ██       █████  ██████  ██ ███████  █████  ████████ ██  ██████  ███    ██
  ██   ██ ██      ██       ██    ██ ██      ██   ██ ██   ██ ██ ██      ██   ██    ██    ██ ██    ██ ████   ██
  ██████  █████   ██   ███ ██    ██ ██      ███████ ██████  ██ ███████ ███████    ██    ██ ██    ██ ██ ██  ██
  ██   ██ ██      ██    ██ ██    ██ ██      ██   ██ ██   ██ ██      ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
  ██   ██ ███████  ██████   ██████  ███████ ██   ██ ██   ██ ██ ███████ ██   ██    ██    ██  ██████  ██   ████
  */

  // TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  // TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;

  /*
██████  ███████ ███    ██ ███████ ██ ████████ ██    ██     ███████ ██       █████   ██████  ███████
██   ██ ██      ████   ██ ██      ██    ██     ██  ██      ██      ██      ██   ██ ██       ██
██   ██ █████   ██ ██  ██ ███████ ██    ██      ████       █████   ██      ███████ ██   ███ ███████
██   ██ ██      ██  ██ ██      ██ ██    ██       ██        ██      ██      ██   ██ ██    ██      ██
██████  ███████ ██   ████ ███████ ██    ██       ██        ██      ███████ ██   ██  ██████  ███████
*/
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
  // TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  // TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidthAndUser;
  // TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeUser;

  /*
  ███████  ██████  █████  ███    ██     ███    ███  ██████  ██████  ███████
  ██      ██      ██   ██ ████   ██     ████  ████ ██    ██ ██   ██ ██
  ███████ ██      ███████ ██ ██  ██     ██ ████ ██ ██    ██ ██   ██ █████
       ██ ██      ██   ██ ██  ██ ██     ██  ██  ██ ██    ██ ██   ██ ██
  ███████  ██████ ██   ██ ██   ████     ██      ██  ██████  ██████  ███████
  */

  // TUnfoldDensity::EScanTauMode scanMode = TUnfoldDensity::kEScanTauRhoAvg;
  // TUnfoldDensity::EScanTauMode scanMode = TUnfoldDensity::kEScanTauRhoMax;
  TUnfoldDensity::EScanTauMode scanMode = TUnfoldDensity::kEScanTauRhoAvgSys;
  // TUnfoldDensity::EScanTauMode scanMode = TUnfoldDensity::kEScanTauRhoSquareAvg;
  // TUnfoldDensity::EScanTauMode scanMode = TUnfoldDensity::kEScanTauRhoSquareAvgSys;



  // Torbens Setup
  // const char *REGULARISATION_DISTRIBUTION=0;
  // const char *REGULARISATION_AXISSTEERING="*[B]";
  //
  // const char *SCAN_DISTRIBUTION="measurement_gen";
  // const char *SCAN_AXISSTEERING=0;

  // latest standard
  // const char *REGULARISATION_DISTRIBUTION=0;
  // const char *REGULARISATION_AXISSTEERING=0;
  //
  // const char *SCAN_DISTRIBUTION=0;
  // const char *SCAN_AXISSTEERING=0;

  const char *REGULARISATION_DISTRIBUTION=0;
  const char *REGULARISATION_AXISSTEERING=0;
  const char *SCAN_DISTRIBUTION=0;
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


  // set up TUnfold Class
  TUnfoldDensity unfold(migration_matrix,TUnfold::kHistMapOutputHoriz,
    regMode,constraintMode,densityFlags,
    binning_gen,binning_rec,
    REGULARISATION_DISTRIBUTION,
    REGULARISATION_AXISSTEERING);

    /*
  ██ ███    ██ ██████  ██    ██ ████████
  ██ ████   ██ ██   ██ ██    ██    ██
  ██ ██ ██  ██ ██████  ██    ██    ██
  ██ ██  ██ ██ ██      ██    ██    ██
  ██ ██   ████ ██       ██████     ██
  */

  unfold.SetInput(input, scale);

  /*
  ██████   █████   ██████ ██   ██  ██████  ██████   ██████  ██    ██ ███    ██ ██████  ███████
  ██   ██ ██   ██ ██      ██  ██  ██       ██   ██ ██    ██ ██    ██ ████   ██ ██   ██ ██
  ██████  ███████ ██      █████   ██   ███ ██████  ██    ██ ██    ██ ██ ██  ██ ██   ██ ███████
  ██   ██ ██   ██ ██      ██  ██  ██    ██ ██   ██ ██    ██ ██    ██ ██  ██ ██ ██   ██      ██
  ██████  ██   ██  ██████ ██   ██  ██████  ██   ██  ██████   ██████  ██   ████ ██████  ███████
  */

  for(unsigned int i=0; i<bgr_name.size(); i++){
    double scale_error = 0;
    if(bgr_name[i] == "WJets") scale_error = 0.19;
    else if(bgr_name[i] == "SingleTop") scale_error = 0.23;
    else if(bgr_name[i] == "other") scale_error = 1.0;
    else cout << "!!!!!ATTENTION!!!!! Background not known (from unfold.C)" << endl;
    unfold.SubtractBackground(backgrounds[i], bgr_name[i], 1.0, scale_error);
  }


  /*
██    ██ ███    ██ ███████  ██████  ██      ██████
██    ██ ████   ██ ██      ██    ██ ██      ██   ██
██    ██ ██ ██  ██ █████   ██    ██ ██      ██   ██
██    ██ ██  ██ ██ ██      ██    ██ ██      ██   ██
 ██████  ██   ████ ██       ██████  ███████ ██████
*/


  TSpline *logTauX=0,*logTauY=0;

  double tau_min = 0.00001;
  double tau_max = 0.9;

  // L-curve scan
  if(do_lcurve && nscan != 0){
    unfold.ScanLcurve(nscan,tau_min,tau_max,&lcurve,&logTauX,&logTauY);
    // get tau value and position on l-curve
    tau=unfold.GetTau();
    double logTau=TMath::Log10(tau);
    lcurveX=logTauX->Eval(logTau);
    lcurveY=logTauY->Eval(logTau);
  }

  // directly put tau value here
  if(nscan == 0){
    tau = tau_;
    unfold.DoUnfold(tau);
  }

  // rho scan
  if(!do_lcurve && nscan != 0){
    unfold.ScanTau(nscan,tau_min,tau_max,&rhoLogTau, scanMode, SCAN_DISTRIBUTION,SCAN_AXISSTEERING, &lcurve,&logTauX,&logTauY);
    // get tau value and position on l-curve
    tau=unfold.GetTau();
    double logTau=TMath::Log10(tau);
    lcurveX=logTauX->Eval(logTau);
    lcurveY=logTauY->Eval(logTau);
  }


  // axis steering for outputs
  const char *axissteer = "pt[C]";

  // create bin map
  // vector<Int_t> binmap = {1,1,2,2,3,3,4,4,5,5};
  vector<Int_t> binmap = {1,2,2,3,3,4,5,-1,-1};



  // write hists before including sys uncertainties
  output_fine = unfold.GetOutput("",0,"measurement_gen",axissteer,kTRUE);
  output = ConvertToNewBinning(output_fine, binmap);
  output_all = unfold.GetOutput("",0,0,0,kFALSE);
  CorM = unfold.GetRhoIJtotal("", 0, "measurement_gen",axissteer,kTRUE);
  ProbM = unfold.GetProbabilityMatrix("", 0, kTRUE);

  // Get Bias Distribution
  TH1* bias_fine = unfold.GetBias( "",0,"measurement_gen",axissteer,kTRUE);
  bias = ConvertToNewBinning(bias_fine, binmap);

  // Statistical uncertainties of input distribution
  TH2* CovInputStat_fine = unfold.GetEmatrixInput("", 0, "measurement_gen",axissteer, kTRUE);
  CovInputStat = ConvertToNewBinning(CovInputStat_fine, binmap);

  // Statistical ucertainties of matrix
  TH2* CovMatrixStat_fine = unfold.GetEmatrixSysUncorr("", 0, "measurement_gen",axissteer, kTRUE);
  CovMatrixStat = ConvertToNewBinning(CovMatrixStat_fine, binmap);

  for(unsigned int i=0; i<bgr_name.size(); i++){
    // Statistical uncertainties from Background
    TH2* CovBgrStat_fine = unfold.GetEmatrixSysBackgroundUncorr(bgr_name[i], "", 0, "measurement_gen",axissteer, kTRUE);
    CovBgrStat.push_back(ConvertToNewBinning(CovBgrStat_fine, binmap));
    // Sys uncertainties from background scale (delta)
    TH1* BgrDelta_fine = unfold.GetDeltaSysBackgroundScale(bgr_name[i], "",0,"measurement_gen",axissteer,kTRUE);
    BgrDelta.push_back(ConvertToNewBinning(BgrDelta_fine, binmap));
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
      TH1* SysDelta_fine = unfold.GetDeltaSysSource(sys_name[i][j], "",0,"measurement_gen",axissteer,kTRUE);
      SysDelta[i].push_back(ConvertToNewBinning(SysDelta_fine, binmap));
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

TH1* unfolding::get_output_fine(){
  return output_fine;
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

TH1* unfolding::GetBiasDistribution(){
  return bias;
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

TSpline* unfolding::GetRhoTau(){
  return rhoLogTau;
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

TH1* unfolding::ConvertToNewBinning(TH1* oldhist, vector<int> binmap){
  TH1* newhist = measurement_gen20->CreateHistogram(oldhist->GetTitle(),kTRUE,0,0,"pt[C]");
  int nbins = newhist->GetXaxis()->GetNbins();
  int nbinsold = oldhist->GetXaxis()->GetNbins();
  if(binmap.size() != nbinsold) cout << "not correct amount of entries in binmap" << endl;
  for(int newbin=1; newbin<=nbins; newbin++){
    double sum=0, err=0, err2=0;
    for(int i=0; i<binmap.size(); i++){
      int oldbin = i+1;
      if(newbin == binmap[i]){
        sum += oldhist->GetBinContent(oldbin);
        err2 += pow(oldhist->GetBinError(oldbin),2);
      }
    }
    err = sqrt(err2);
    newhist->SetBinContent(newbin, sum);
    newhist->SetBinError(newbin, err);
  }
  // for(int bin=1; bin<=nbins; bin++){
  //   int oldbin = bin*2 - 1;
  //   double sum = 0, err =0;
  //   sum = oldhist->GetBinContent(oldbin) + oldhist->GetBinContent(oldbin+1);
  //   err = sqrt(oldhist->GetBinError(oldbin)*oldhist->GetBinError(oldbin) + oldhist->GetBinError(oldbin+1)*oldhist->GetBinError(oldbin+1));
  //   newhist->SetBinContent(bin, sum);
  //   newhist->SetBinError(bin, err);
  // }
  return newhist;
}

TH2* unfolding::ConvertToNewBinning(TH2* oldcov, vector<int> binmap){
  TH2* newcov = measurement_gen20->CreateErrorMatrixHistogram(oldcov->GetTitle(), kTRUE, 0, 0, "pt[C]");
  int nbins = newcov->GetXaxis()->GetNbins();
  int nbinsold = oldcov->GetXaxis()->GetNbins();
  if(binmap.size() != nbinsold) cout << "not correct amount of entries in binmap" << endl;
  for(int xbin=1; xbin<=nbins; xbin++){
    for(int ybin=1; ybin<=nbins; ybin++){
      // cout << "------------------------------" << endl;
      // cout << "x = " << xbin << endl;
      // cout << "y = " << ybin << endl;
      double new_entry = 0;
      for(int i=0; i<binmap.size(); i++){
        if(xbin == binmap[i]){
          int oldxbin = i+1;
          for(int j=0; j<binmap.size(); j++){
            if(ybin == binmap[j]){
              int oldybin = j+1;
              // cout << "x' = " << oldxbin << ", y' = " << oldybin << endl;
              new_entry += oldcov->GetBinContent(oldxbin,   oldybin);
            }
          }
        }
      }
      newcov->SetBinContent(xbin, ybin, new_entry);
    }
  }

  // for(int xbin=1; xbin<=nbins; xbin++){
  //   int oldxbin = xbin*2 - 1;
  //   for(int ybin=1; ybin<=nbins; ybin++){
  //     int oldybin = ybin*2 - 1;
  //     double new_entry = 0;
  //     new_entry += oldcov->GetBinContent(oldxbin,   oldybin);
  //     new_entry += oldcov->GetBinContent(oldxbin,   oldybin+1);
  //     new_entry += oldcov->GetBinContent(oldxbin+1, oldybin);
  //     new_entry += oldcov->GetBinContent(oldxbin+1, oldybin+1);
  //     newcov->SetBinContent(xbin, ybin, new_entry);
  //   }
  // }
  return newcov;
}
