#include "chi2fit.h"

chi2fit::chi2fit(TH1D* data, TH2D* cov_matrix, std::vector<TH1D*> masspoints, std::vector<double> masses, double lower, double upper, bool width){
  width_ = width;
  SetData(data);
  SetCovMatrix(cov_matrix);
  SetMasspoints(masspoints);
  SetMasses(masses);
  SetRange(lower, upper);
  SetBins(masspoints);
}


void chi2fit::SetMasspoints(std::vector<TH1D*> masspoints){
  for(unsigned int i = 0; i<masspoints.size(); i++) {
      TH1D *mass = (TH1D*)masspoints[i]->Clone();
      masspoints_.push_back(mass);
  }
  return;
}

void chi2fit::SetData(TH1D* data){
  data_ = (TH1D*)data->Clone("data_");
  return;
}

void chi2fit::SetCovMatrix(TH2D* cov_matrix){
  cov_ = (TH2D*)cov_matrix->Clone("cov_");
  return;
}


void chi2fit::SetRange(double lower, double upper){
  lower_ = data_->GetXaxis()->FindBin(lower);
  upper_ = data_->GetXaxis()->FindBin(upper);
  return;
}

void chi2fit::SetBins(std::vector<TH1D*> masspoints){
  // when calculating the chi2, one degree of freedom is lost by normalisation
  // therefore covariance matrix is only invertable if one excludes one bin from chi2
  // this bin is selected randomly (chi2 should not change depending on the bin)
  // for every masspoint, the bins that should be used are written into a vector
  bool skip;
  int rand_bin = rand()%(upper_-lower_+1)+lower_;
  std::cout << "Exclude bin " << rand_bin << " from chi2" << std::endl;
  for(unsigned int j=0; j<masspoints.size(); j++){
    std::vector<int> bins_this;
    for(int i=lower_; i<=upper_; i++){
      if(i == rand_bin) skip = true;
      else skip = false;
      if(!skip) bins_this.push_back(i);
    }
    bins_.push_back(bins_this);
  }
  return;
}


void chi2fit::SetMasses(std::vector<double> masses){
  masses_.ResizeTo(masses.size());
  for(unsigned int i=0; i<masses.size(); i++) masses_[i] = masses[i];
  return;
}

void chi2fit::CalculateChi2(){
  int j = 0;
  const int NMasses = masspoints_.size();
  chi2_.ResizeTo(NMasses);
  for(auto &masspoint: masspoints_){
    double chi2 = ComputeChi2(masspoint, bins_[j]);
    chi2_[j] = chi2;
    j++;
  }

  TGraph* chi2Hist = new TGraph(masses_, chi2_);
  TF1*f1 = new TF1("f1","pol2",0,500);
  chi2Hist->Fit("f1","R");
  fit_ = chi2Hist->GetFunction("f1");

  return;
}

double chi2fit::GetMass(){
  double minY = fit_->GetMinimum();
  double mass = fit_->GetX(minY, 140, 200);
  return mass;
}

double chi2fit::GetMin(){
  double minY = fit_->GetMinimum();
  return minY;
}

double chi2fit::GetUncertainty(){
  double minY = fit_->GetMinimum();
  double minX = fit_->GetX(minY, 140, 200);
  double minX_sigmaup = fit_->GetX(minY+1, minX, 200);
  // double minX_sigmadown = fit_->GetX(minY+1, 140, minX);
  double uncert = minX_sigmaup-minX;
  return uncert;
}

std::vector<double> chi2fit::GetChi2Values(){
  std::vector<double> chi2;
  int N = chi2_.GetNoElements();
  for(int i=0; i<N; i++) chi2.push_back(chi2_[i]);
  return chi2;
}


TF1* chi2fit::GetChi2Fit(){
  return fit_;
}


//
// compute chi^2
//
double chi2fit::ComputeChi2(TH1D* MC, std::vector<int> bins){
  //double chi2 = 0;
  TH1D* data_sc = (TH1D*)data_->Clone("data_sc");
  TH2D* cov_sc = (TH2D*)cov_->Clone("cov_sc");
  TH1D* MC_sc = (TH1D*)MC->Clone("MC_sc");

  // bins
  int N =  bins.size();

  //vectors and matrices
  TVectorD vdata(N);
  TVectorD vMC(N);
  TMatrixDSym mat(N);

  for (unsigned int i=0; i < bins.size(); ++i){
    vdata[i] = data_sc->GetBinContent(bins.at(i));
    vMC[i] = MC_sc->GetBinContent(bins.at(i));
    for (unsigned int j=0; j < bins.size(); ++j){
      mat[i][j] = cov_sc->GetBinContent(bins.at(i),bins.at(j));
      // cout << "cov[" << i << "," << j << "] = " << cov->GetBinContent(bins.at(i),bins.at(j)) << endl;
    }
  }

  //invert the matrix
  // TDecompSVD lu(mat);
  TDecompLU lu(mat);

  TMatrixD imat = TMatrixD(mat);
  lu.Invert(imat);

  TMatrixD test(mat, TMatrixD::kMult, imat);
  TMatrixD unit(TMatrixD::kUnit, mat);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 10e-10);
  if (!ok){
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<N; ++j)  printf("% 5.2f ", mat[i][j] );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("(%3d,%3d)   %9.3f",i,i,mat[i][i]);
      if ( mat[i][i] < 1.e-10 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // calculate chi2
  double chi2 = 0;

  TVectorD diff = vdata - vMC;
  TVectorD right = imat * diff;
  chi2 = diff * right;
  std::cout << "Chi2 = " << chi2 << std::endl;
  return chi2;
}
