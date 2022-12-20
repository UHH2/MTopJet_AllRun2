#include "../include/CentralInclude.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Fit Data to one hypothesis

class chi2value{
public:
  chi2value(TH1F* data, TH2F* cov_matrix, TH1F* hypothesis, double lower, double upper, bool width);

  void CalculateChi2();

  TF1* GetChi2Fit();
  double GetChi2Value();
  double GetMin();
  double GetUncertainty();
  TGraph* GetChi2Graph();

private:
  void SetMasses(double masses);
  double ComputeChi2(TH1F* MC, std::vector<int> bins);
  void SetBins(TH1F* hypothesis);
  void SetRange(double lower, double upper);
  void SetHypothesis(TH1F* hypothesis);
  void SetData(TH1F* data);
  void SetCovMatrix(TH2F* cov_matrix);


  TH1F* hypothesis_;
  TH1F* data_;
  int lower_, upper_;
  bool width_;
  TH2F* cov_;
  std::vector<int>  bins_;
  double chi2_;
  TF1 * fit_;

};

chi2value::chi2value(TH1F* data, TH2F* cov_matrix, TH1F* hypothesis, double lower, double upper, bool width){
  width_ = width;
  SetData(data);
  SetCovMatrix(cov_matrix);
  SetHypothesis(hypothesis);
  SetRange(lower, upper);
  SetBins(hypothesis);
}


void chi2value::SetHypothesis(TH1F* hypothesis){
  hypothesis_ = (TH1F*)hypothesis->Clone();
  return;
}

void chi2value::SetData(TH1F* data){
  data_ = (TH1F*)data->Clone("data_");
  return;
}

void chi2value::SetCovMatrix(TH2F* cov_matrix){
  cov_ = (TH2F*)cov_matrix->Clone("cov_");
  return;
}


void chi2value::SetRange(double lower, double upper){
  lower_ = data_->GetXaxis()->FindBin(lower);
  upper_ = data_->GetXaxis()->FindBin(upper-1);
  cout << "lower bin = " << lower_ << ", upper bin = " << upper_ << endl;
  return;
}

void chi2value::CalculateChi2(){
  chi2_ = ComputeChi2(hypothesis_, bins_);
  return;
}

double chi2value::GetChi2Value(){
  return chi2_;
}

void chi2value::SetBins(TH1F* hypothesis){
  // when calculating the chi2, one degree of freedom is lost by normalisation
  // therefore covariance matrix is only invertable if one excludes one bin from chi2
  // this bin is selected randomly (chi2 should not change depending on the bin)
  // for every masspoint, the bins that should be used are written into a vector
  bool skip;
  int rand_bin = rand()%(upper_-lower_+1)+lower_;
  // int rand_bin = 20;
  std::cout << "Exclude bin " << rand_bin << " from chi2" << std::endl;

  std::vector<int> bins_this;
  for(int i=lower_; i<=upper_; i++){
    if(i == rand_bin) skip = true;
    else skip = false;
    if(!skip) bins_this.push_back(i);
  }
  bins_ = bins_this;

  return;
}

//
// compute chi^2
//
double chi2value::ComputeChi2(TH1F* MC, std::vector<int> bins){
  //double chi2 = 0;
  TH1F* data_sc = (TH1F*)data_->Clone("data_sc");
  TH2F* cov_sc = (TH2F*)cov_->Clone("cov_sc");
  TH1F* MC_sc = (TH1F*)MC->Clone("MC_sc");

  // bins
  int N =  bins.size();

  //vectors and matrices
  TVectorD vdata(N);
  TVectorD vMC(N);
  TMatrixDSym mat(N);

  // vector<bool> skip;
  // for(int i=0; i<180; i++) skip.push_pack(true);

  for (unsigned int i=0; i < bins.size(); ++i){
    vdata[i] = data_sc->GetBinContent(bins.at(i));
    vMC[i] = MC_sc->GetBinContent(bins.at(i));
    for (unsigned int j=0; j < bins.size(); ++j){
      mat[i][j] = cov_sc->GetBinContent(bins.at(i),bins.at(j));
      // cout << "cov[" << i << "," << j << "] = " << cov_sc->GetBinContent(bins.at(i),bins.at(j)) << endl;
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
      if ( mat[i][i] < 1.e-20 ) printf("   <<<<<<<< \n");
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Fit Data to multiple hypothesis

class chi2fit{
public:
  chi2fit(TH1F* data, TH2F* cov_matrix, std::vector<TH1F*> masspoints, std::vector<double> masses, double lower, double upper, bool width);

  void CalculateChi2();

  TF1* GetChi2Fit();
  std::vector<double> GetChi2Values();
  double GetMass();
  double GetMin();
  double GetUncertainty();
  TGraph* GetChi2Graph();

private:
  void SetMasses(std::vector<double> masses);
  double ComputeChi2(TH1F* MC, std::vector<int> bins);
  void SetBins(std::vector<TH1F*> masspoints);
  void SetRange(double lower, double upper);
  void SetMasspoints(std::vector<TH1F*> masspoints);
  void SetData(TH1F* data);
  void SetCovMatrix(TH2F* cov_matrix);


  std::vector<TH1F*> masspoints_;
  TH1F* data_;
  int lower_, upper_;
  bool width_;
  TH2F* cov_;
  std::vector< std::vector<int> > bins_;
  TVectorD chi2_;
  TVectorD masses_;
  TF1 * fit_;

};

chi2fit::chi2fit(TH1F* data, TH2F* cov_matrix, std::vector<TH1F*> masspoints, std::vector<double> masses, double lower, double upper, bool width){
  width_ = width;
  SetData(data);
  SetCovMatrix(cov_matrix);
  SetMasspoints(masspoints);
  SetMasses(masses);
  SetRange(lower, upper);
  SetBins(masspoints);
}


void chi2fit::SetMasspoints(std::vector<TH1F*> masspoints){
  for(unsigned int i = 0; i<masspoints.size(); i++) {
    TH1F *mass = (TH1F*)masspoints[i]->Clone();
    masspoints_.push_back(mass);
  }
  return;
}

void chi2fit::SetData(TH1F* data){
  data_ = (TH1F*)data->Clone("data_");
  return;
}

void chi2fit::SetCovMatrix(TH2F* cov_matrix){
  cov_ = (TH2F*)cov_matrix->Clone("cov_");
  return;
}


void chi2fit::SetRange(double lower, double upper){
  lower_ = data_->GetXaxis()->FindBin(lower);
  upper_ = data_->GetXaxis()->FindBin(upper-1);
  cout << "lower bin = " << lower_ << ", upper bin = " << upper_ << endl;
  return;
}

void chi2fit::SetBins(std::vector<TH1F*> masspoints){
  // when calculating the chi2, one degree of freedom is lost by normalisation
  // therefore covariance matrix is only invertable if one excludes one bin from chi2
  // this bin is selected randomly (chi2 should not change depending on the bin)
  // for every masspoint, the bins that should be used are written into a vector
  bool skip;
  int rand_bin = rand()%(upper_-lower_+1)+lower_;
  // int rand_bin = 20;
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


TGraph * chi2fit::GetChi2Graph(){
  TGraph* chi_hist = new TGraph(masses_,chi2_);
  return chi_hist;
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
double chi2fit::ComputeChi2(TH1F* MC, std::vector<int> bins){
  //double chi2 = 0;
  TH1F* data_sc = (TH1F*)data_->Clone("data_sc");
  TH2F* cov_sc = (TH2F*)cov_->Clone("cov_sc");
  TH1F* MC_sc = (TH1F*)MC->Clone("MC_sc");

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
      // cout << "cov[" << i << "," << j << "] = " << cov_sc->GetBinContent(bins.at(i),bins.at(j)) << endl;
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
      if ( mat[i][i] < 1.e-20 ) printf("   <<<<<<<< \n");
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
