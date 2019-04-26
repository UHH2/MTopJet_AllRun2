#include "chi2sys.h"

chi2sys::chi2sys(TH1* data_, TH2* cov_, TH1D* truth_){
  // set all member variables
  data  = (TH1D*) data_->Clone();
  DELTA = (TH1D*) data_->Clone();
  truth = (TH1D*)truth_->Clone();
  oldcov = (TH2D*)cov_->Clone();
  Nbins =  truth->GetSize() - 2;

  // also create a delta matrix with zeros everywhere
  DELTAcov = (TH2D*)oldcov->Clone();
  CreateDeltaCov(0.0);

  // calculate chi2/ndf
  double chi2ndf = GetChi2NDFValue(oldcov);

  // if this vaule is smaller one, you are done
  // else create a new shift, update cov and calculate chi2/ndf again
  bool ChiSmallerOne = false;
  if(chi2ndf < 1) ChiSmallerOne = true;
  factor = 0.0;           // set start value
  double stepsize = 0.02; // set stepsize
  while(!ChiSmallerOne){
    CreateDeltaCov(factor);
    TH2D* newcov = (TH2D*)oldcov->Clone();
    newcov->Add(DELTAcov);
    double newchi2ndf = GetChi2NDFValue(newcov);
    if(newchi2ndf < 1) ChiSmallerOne = true;
    else               factor += stepsize;
  }
}

double chi2sys::GetChi2NDFValue(TH2D* inputcov){
  //vectors and matrices
  TVectorD vdata(Nbins);
  TVectorD vtruth(Nbins);
  TMatrixDSym mcov(Nbins);

  for (unsigned int i=0; i < Nbins; ++i){
    vdata[i] = data->GetBinContent(i+1);
    vtruth[i] = truth->GetBinContent(i+1);
    for (unsigned int j=0; j < Nbins; ++j){
      mcov[i][j] = inputcov->GetBinContent(i+1,j+1);
      // cout << "inputcov[" << i+1 << "," << j+1 << "] = " << inputcov->GetBinContent(i+1,j+1) << endl;
    }
  }

  //invert the mcov
  // TDecompSVD lu(mcov);
  TDecompLU lu(mcov);

  TMatrixD imcov = TMatrixD(mcov);
  lu.Invert(imcov);

  TMatrixD test(mcov, TMatrixD::kMult, imcov);
  TMatrixD unit(TMatrixD::kUnit, mcov);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 10e-10);
  if (!ok){
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<Nbins; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<Nbins; ++j)  printf("% 5.2f ", mcov[i][j] );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<Nbins; ++i){
      printf("(%3d,%3d)   %9.3f",i,i,mcov[i][i]);
      if ( mcov[i][i] < 1.e-10 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // calculate chi2
  double chi2 = 0;
  double chi2_ndf = 0;

  TVectorD diff = vdata - vtruth;
  TVectorD right = imcov * diff;
  chi2 = diff * right;
  chi2_ndf = chi2/Nbins;
  // std::cout << "Chi2 = " << chi2 << std::endl;
  std::cout << "Chi2/ndf = " << chi2_ndf << std::endl;

  return chi2_ndf;
}

void chi2sys::CreateDeltaCov(double factor){
  for (unsigned int i=1; i <= Nbins; ++i){
    double delta_i = data->GetBinContent(i) - truth->GetBinContent(i);
    DELTA->SetBinContent(i, delta_i*factor);
    for (unsigned int j=1; j <= Nbins; ++j){
      double delta_j = data->GetBinContent(j) - truth->GetBinContent(j);
      double CovEntry = delta_i * factor * delta_j * factor;
      DELTAcov->SetBinContent(i,j,CovEntry);
    }
  }
  return;
}
