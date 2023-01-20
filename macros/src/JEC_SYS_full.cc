#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"
#include "../include/Chi2Fit.h"
#include "../include/tdrstyle_all.h"

#include "TRandom2.h"
#include "TF3.h"
#include "TError.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"

// #include "../CovMatrices/JMS/CollectCovHeaders.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>
#include "TSystem.h"

using namespace std;

typedef map<int, TGraph2DErrors*> MapI2DGe;
typedef map<TString, MapI2DGe> MapSI2DGe;

typedef map<TString, TMatrixD> MapM;
typedef map<TString, MapM> MapMM;
typedef map<TString, MapMM> MapMMM;

typedef map<TString, TF2*> MapF2;
typedef map<TString, MapF2> MapFF2;
typedef map<int, TF2*> MapIF2;
typedef map<TString, MapIF2> MapSIF2;

// ---------------------------------------------------------------------------
// Declare functions

// MapHHH SetUpMap(VecTS hists);
// TH1F* SetError(TH1F* data, TH2F* cov);
// void PrintValues(VecD& nominal, VecD& uu, VecD& ud, VecD& du, VecD& dd, TString name, double &sigma_cor_center, double& rho, double& sigma_cor, double& sigma_jec, VecD& PoU_COR, VecD& PoU_JEC, double& sigma_jms);
// void AddTwoMaps(MapHHH& map1, MapHHH& map2, int option);
// bool sortcolx( const VecD& v1, const VecD& v2 );
// MapHH GetHistograms(TString process, TString h_name);
// TH1F* RebinMap(MapHH& process, int& width);
// TH2F* CreateMatrix(TH1F* hist);
// MapH NormalizeMap(MapH& process);
// TMatrixD GetCovMatrix(TH1F* map, const TString& save, const TString& process);
// TMatrixD GetCovMatrixNorm(TMatrixD mtemp, TH1F* hist, const TString& save, const TString& process);
// TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel=0, TString add = "");
// TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
// vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist);
// void SubtractBackgroundsMap(MapHH& m_data, MapHH& m_st, MapHH& m_wjets, MapHH& m_other);
//
// vector<vector<double>> FindXY_alt(TF2 *function, double zfix, double xmin, double xmax, double ymin, double ymax, int steps = 1000, double accuracy = 0, bool count=false);
//
// TCanvas *c_sys = new TCanvas("sys", "sys", 600, 600); // to Plot all systematics in on pdf
// TString cut = "";

bool debug = false;
bool fast = false; // Reads Cov Matrix for data and SYS from root file
TString save, year, channel;
TString save_path, save_general, save_nfs;
int space;
int precision = 10; // For double to string function - dtos()
int cout_precision = 6;

// inline TString wbin(TString h){return h.ReplaceAll(cut, "");}

TH1F* BuildSYSHist(TH1F* sys_, TH1F* nom_, TString name);
TH2F* BuildCovMatrix(TH1F* hist_, TString name);
TH2F* BuildSYSCovMatrix(TH1F* sys_, TString name);
TH2F* NormCovMatrix(TH2F* cov_, TH1F* hist_, TString name);
TH2F* BuildTotalCovMatrix(vector<TH2F*> add_, TString name);
TH1F* TrimHistogram(TH1F* hist_, double r_low, double r_high, TString name);
double GetChi2(TH1F* data_, TH1F* pred_, TH2F* cov_, TString name, bool norm);
void PrintCovMatrix(TH2F* cov_, double factor, TString name);
void PrintHist(TH1F* hist_, double factor, TString name);

TH1F* BuildSYSHist(TH1F* sys_, TH1F* nom_, TString name){
  TH1F* sys = (TH1F*) sys_->Clone();
  TH1F* nom = (TH1F*) nom_->Clone();
  int n_sys = sys->GetNbinsX();
  int n_nom = nom->GetNbinsX();
  if(n_sys != n_nom){
    throw runtime_error("<E> Systematic does not have same dimension as nominal hist for "+name);
  }
  TH1F* unc = (TH1F*) sys_->Clone();
  for(int i=1; i<=n_sys; i++){
    double c_sys = sys->GetBinContent(i);
    double c_nom = nom->GetBinContent(i);
    double content = c_nom-c_sys;
    unc->SetBinContent(i, content);
  }
  if (debug){
    cout << endl << " === " << name << endl;
    cout << "sys | "; for(int i=1; i<=n_sys; i++) cout << setw(13) << sys->GetBinContent(i); cout << endl;
    cout << "nom | "; for(int i=1; i<=n_sys; i++) cout << setw(13) << nom->GetBinContent(i); cout << endl;
    cout << "unc | "; for(int i=1; i<=n_sys; i++) cout << setw(13) << unc->GetBinContent(i); cout << endl;
  }
  return unc;
}

TH2F* BuildCovMatrix(TH1F* hist_, TString name){
  TH1F* hist = (TH1F*) hist_->Clone();
  int n_hist = hist->GetNbinsX();
  double low = hist->GetXaxis()->GetBinLowEdge(1);
  double high = hist->GetXaxis()->GetBinUpEdge(n_hist);
  TH2F *cov = new TH2F(name, name, n_hist,low,high,n_hist,low,high);
  for(int i=1; i<=n_hist; i++){
    double error = hist->GetBinError(i);
    cov->SetBinContent(i,i, error*error);
  }
  if(debug) PrintCovMatrix(cov, 1, name);
  return cov;
}


TH2F* BuildSYSCovMatrix(TH1F* sys_, TString name){
  TH1F* sys = (TH1F*) sys_->Clone();
  int n_sys = sys->GetNbinsX();
  double low = sys->GetXaxis()->GetBinLowEdge(1);
  double high = sys->GetXaxis()->GetBinUpEdge(n_sys);
  TH2F *cov = new TH2F(name, name, n_sys,low,high,n_sys,low,high);
  if(debug) cout << setw(21) << name << " - " << setw(7) << n_sys << setw(7) << low << setw(7) << high << setw(7) << n_sys << setw(7) << low << setw(7) << high << endl;
  for(int i=1; i<=n_sys; i++){
    double i_sys = sys->GetBinContent(i);
    for(int j=1; j<=n_sys; j++){
      double j_sys = sys->GetBinContent(j);
      cov->SetBinContent(i,j, i_sys*j_sys);
    }
  }
  if(debug) PrintCovMatrix(cov, 1, name);
  return cov;
}

TH2F* NormCovMatrix(TH2F* cov_, TH1F* hist_, TString name){
  TH2F* cov = (TH2F*) cov_->Clone();
  TH2F* cov_norm = (TH2F*) cov_->Clone();
  TH2F* hist = (TH2F*) hist_->Clone();
  int n_cols = cov->GetNbinsX();
  int n_rows = cov->GetNbinsY();
  int n_bins = hist->GetNbinsX();
  double integral = hist->Integral();
  double integral2 = pow(integral,2);
  if(n_cols != n_rows) throw runtime_error("<E> Cov. matrix is not sysmmetric - "+name);
  if(n_cols != n_bins) throw runtime_error("<E> Cov. matrix and hist have not the same dimension - "+name);
  for(int i=1; i <= n_cols; i++){
    double content_i = hist->GetBinContent(i);
    double subtracted_i = integral - content_i;
    for(int j=1; j <= n_cols; j++){
      double content_j = hist->GetBinContent(j);
      double subtracted_j = integral - content_j;
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= n_cols; k++){
        for(int l=1; l <= n_cols; l++){
          old_entry = cov->GetBinContent(k,l);
          if(i==k) derivation_i =   (subtracted_i) / integral2;
          else     derivation_i = - (content_i)    / integral2;
          if(j==l) derivation_j =   (subtracted_j) / integral2;
          else     derivation_j = - (content_j)    / integral2;
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      cov_norm->SetBinContent(i, j, sum);
    }
  }
  return cov_norm;
}

TH2F* BuildTotalCovMatrix(vector<TH2F*> add_, TString name){
  TH2F* cov = (TH2F*) add_[0]->Clone();
  int n_cols = cov->GetNbinsX();
  int n_rows = cov->GetNbinsY();
  int n_mat = add_.size();
  if(n_cols != n_rows) {
    printf("First cov matrix: cols %3i vs. rows %3i\n",n_cols,n_rows);
    throw runtime_error("<E> First cov. matrix is not sysmmetric - "+name);
  }
  if(debug) printf("Add %3i matrices\n",n_mat);
  for(int i=1; i < n_mat; i++){
    TH2F* add = (TH2F*) add_[i]->Clone();
    int n_cols_i = add->GetNbinsX();
    int n_rows_i = add->GetNbinsY();
    if(n_cols_i != n_rows_i){
      printf("From added cov matrix from bin %3i: cols %3i vs. rows %3i\n",i,n_cols_i,n_rows_i);
      throw runtime_error("<E> Added cov. matrix is not sysmmetric - "+name);
    }
    if(n_cols_i != n_rows_i){
      printf("Cov matrices in bin %3i: added cols %3i vs. cov cols %3i\n",i,n_cols_i,n_cols);
      throw runtime_error("<E> Covariance matrices have not the same dimensions - "+name);
    }
    cov->Add(add, 1);
  }
  return cov;
}

TH1F* TrimHistogram(TH1F* hist_, double r_low, double r_high, TString name){
  TH1F* hist = (TH1F*) hist_->Clone();
  int n_hist = hist->GetNbinsX();
  vector<int> bins;
  for(int i=1; i<=n_hist; i++){
    double center = hist->GetBinCenter(i);
    if(center > r_low){
      if(center < r_high){
        bins.push_back(i);
      }
    }
  }
  int N = bins.size();
  int e_low = hist->GetXaxis()->GetBinLowEdge(bins[0]);
  int e_high = hist->GetXaxis()->GetBinUpEdge(bins[N-1]);
  if(e_low != r_low){
    printf("range: %3d vs. edge: %3d\n",r_low,e_low);
    throw runtime_error("<E> Lower boundaries do not fit - "+name);
  }
  if(e_high != r_high){
    printf("range: %3d vs. edge: %3d\n",r_high,e_high);
    throw runtime_error("<E> Higher boundaries do not fit - "+name);
  }
  if(debug) cout << "Consider" << setw(3) << N << "from" << setw(5) << e_low << "to" << setw(5) << e_high << endl;
  TH1F* hist_trim = new TH1F(name, name, N, e_low, e_high);
  for(int i=0; i<N;i++){
    hist_trim->SetBinContent(i+1, hist->GetBinContent(bins[i]));
    hist_trim->SetBinError(i+1, hist->GetBinError(bins[i]));
  }
  return hist_trim;
}

double GetChi2Cut(TH1F* data_, TH1F* pred_, TH2F* cov_, double xmin, double xmax, TString name, bool norm){
  TH1F* data = (TH1F*) data_->Clone();
  TH1F* pred = (TH1F*) pred_->Clone();
  TH2F* cov = (TH2F*) cov_->Clone();
  int dim_d = data->GetNbinsX();
  int dim_p = pred->GetNbinsX();
  int c_row = cov->GetNbinsX();
  int c_col = cov->GetNbinsY();
  if(dim_d != dim_p){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Data and prediciton have not the same dimension - "+name);
  }
  if(c_row != c_col){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Covariance matrix is not symmetric - "+name);
  }
  if(c_row != dim_d){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Covariance matrix has not same dimension as data - "+name);
  }
  if(debug) printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);

  int b_xmin = data->FindBin(xmin);
  int b_xmax = data->FindBin(xmax);
  int bin_skip = rand()%(b_xmin-1)+1; // not remove bins in peak
  vector<int> bins = {};
  vector<int> bins_trim = {};
  for(unsigned int i=1; i<=dim_d; i++){
    bool skip = data->GetBinContent(i)==0 || pred->GetBinContent(i)==0 || (i == bin_skip && norm);
    if(!skip){
      bins.push_back(i);
      if(i>b_xmin && i<b_xmax) bins_trim.push_back(i);
    }
  }
  int N = bins.size();
  if(debug) printf("Use %3i bins\n",N);

  TMatrixDSym mat(N);
  for (unsigned int i=0; i < N; ++i){
    for (unsigned int j=0; j < N; ++j){
      mat[i][j] = cov->GetBinContent(bins.at(i),bins.at(j));
    }
  }

  //invert the matrix
  TDecompLU lu(mat);
  TMatrixD imat = TMatrixD(mat);
  lu.Invert(imat);
  TMatrixD test(mat, TMatrixD::kMult, imat);
  TMatrixD unit(TMatrixD::kUnit, mat);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 1.e-8);
  if (!ok){
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<N; ++j)  printf("%7.2f ", mat[i][j] );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("(%3d,%3d)   %9.8f",i,i,mat[i][i]*1e8);
      if ( mat[i][i] < 1.e-8 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // trim
  int N_trim = bins_trim.size();
  TVectorD vdata(N_trim);
  TVectorD vMC(N_trim);
  for(int i=0;i<N_trim;i++){
    vdata[i] = data->GetBinContent(bins_trim.at(i));
    vMC[i] = pred->GetBinContent(bins_trim.at(i));
  }
  TMatrixD imat_trim = TMatrixD(N_trim, N_trim);
  for(int i=0;i<N_trim;i++){
    for(int j=0;j<N_trim;j++){
      imat_trim[i][j] = imat[bins_trim.at(i)-1][bins_trim.at(j)-1];
    }
  }

  int n_d = vdata.GetNrows();
  int n_m = vMC.GetNrows();
  int n_r = imat_trim.GetNrows();
  int n_c = imat_trim.GetNcols();
  int n_b = N_trim;
  if(n_d!=n_m) cout << "FUCK d and m" << endl;
  if(n_r!=n_c) cout << "FUCK r and c" << endl;
  if(n_d!=n_r) cout << "FUCK d and r" << endl;
  if(n_d!=n_b) cout << "FUCK d and b" << endl;

  // calculate chi2
  double chi2 = 0;
  TVectorD diff = vdata - vMC;
  TVectorD right = imat_trim * diff;
  chi2 = diff * right;
  double chi2ndf = chi2/N;
  cout << setw(20) << name;
  if(norm) printf(": chi2 = %5.2f with ndf = %2i lead to chi2/ndf = %5.2f (skip bin %2i)\n",chi2,N,chi2ndf,bin_skip);
  else printf(": chi2 = %5.2f with ndf = %2i lead to chi2/ndf = %5.2f\n",chi2,N,chi2ndf);
  return chi2;
}


double GetChi2(TH1F* data_, TH1F* pred_, TH2F* cov_, TString name, bool norm){
  TH1F* data = (TH1F*) data_->Clone();
  TH1F* pred = (TH1F*) pred_->Clone();
  TH2F* cov = (TH2F*) cov_->Clone();
  int dim_d = data->GetNbinsX();
  int dim_p = pred->GetNbinsX();
  int c_row = cov->GetNbinsX();
  int c_col = cov->GetNbinsY();
  if(dim_d != dim_p){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Data and prediciton have not the same dimension - "+name);
  }
  if(c_row != c_col){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Covariance matrix is not symmetric - "+name);
  }
  if(c_row != dim_d){
    printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);
    throw runtime_error("<E> Covariance matrix has not same dimension as data - "+name);
  }
  if(debug) printf("data: %3d vs. prediction: %3d vs. columns: %3d vs. rows: %3d\n",dim_d,dim_p,c_row,c_col);

  int bin_skip = rand()%dim_d+1;
  vector<int> bins = {};
  for(unsigned int i=1; i<=dim_d; i++){
    bool skip = data->GetBinContent(i)==0 || pred->GetBinContent(i)==0 || (i == bin_skip && norm);
    if(!skip) bins.push_back(i);
  }
  int N = bins.size();
  if(debug) printf("Use %3i bins\n",N);

  TVectorD vdata(N);
  TVectorD vMC(N);
  TMatrixDSym mat(N);
  for (unsigned int i=0; i < N; ++i){
    vdata[i] = data->GetBinContent(bins.at(i));
    vMC[i] = pred->GetBinContent(bins.at(i));
    for (unsigned int j=0; j < N; ++j){
      mat[i][j] = cov->GetBinContent(bins.at(i),bins.at(j));
    }
  }

  //invert the matrix
  TDecompLU lu(mat);
  TMatrixD imat = TMatrixD(mat);
  lu.Invert(imat);
  TMatrixD test(mat, TMatrixD::kMult, imat);
  TMatrixD unit(TMatrixD::kUnit, mat);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 10e-8);
  if (!ok){
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<N; ++j)  printf("% 7.2f ", mat[i][j]*1e6 );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("(%3d,%3d)   %9.3f",i,i,mat[i][i]*1e6);
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
  double chi2ndf = chi2/N;
  cout << setw(20) << name;
  if(norm) printf(": chi2 = %5.2f with ndf = %2i lead to chi2/ndf = %5.2f (skip bin %2i)\n",chi2,N,chi2ndf,bin_skip);
  else printf(": chi2 = %5.2f with ndf = %2i lead to chi2/ndf = %5.2f\n",chi2,N,chi2ndf);
  return chi2;
}

void PrintCovMatrix(TH2F* cov_, double factor, TString name){
  TH2F* cov = (TH2F*) cov_->Clone();
  int n_cols = cov->GetNbinsX();
  int n_rows = cov->GetNbinsY();
  if(n_cols != n_rows) {
    printf("First cov matrix: cols %3i vs. rows %3i\n",n_cols,n_rows);
    throw runtime_error("<E> First cov. matrix is not sysmmetric - "+name);
  }
  printf("\n\n------------------------------- "+name+" ----------------------------- \n");
  // printf("    | ");
  // for (int i=1; i<=n_cols; ++i) printf("%9i",i); cout << endl;
  for (int i=1; i<=n_cols; ++i){
    printf("%3d | ",i);
    for (int j=1; j<=n_cols; ++j)  printf("% 7.2f ", cov->GetBinContent(i,j)*factor );
    printf("\n");
  }
}

void PrintHist(TH1F* hist_, double factor, TString name){
  TH2F* hist = (TH2F*) hist_->Clone();
  int N = hist->GetNbinsX();
  printf("\n\n------------------------------- "+name+" ----------------------------- \n");
  for (int i=1; i<=N; ++i) cout << setw(15) << hist->GetBinContent(i);
  cout << endl;
}

int main(int argc, char* argv[]){

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  // Input ---------------------------------------------------------------------
  year = "combine";
  channel = "combine";

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptTitle(0);
  SetupGlobalStyle();

  // year ---------------------------------------------------------------------
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16  = true;
  else if(strcmp(year, "2017")==0)     is17  = true;
  else if(strcmp(year, "2018")==0)     is18  = true;
  else if(strcmp(year, "combine")==0)  isAll = true;
  // else cerr("Give me the correct year please (2016, 2017, 2018 or combine)");

  // Bins ----------------------------------------------------------------------
  int bin_width = 1;

  // Create Directories --------------------------------------------------------
  if(debug) cout << "Create Directiories ..." << endl;

  save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_peak/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/WMass/"+year+"/"+channel+"/";
  save = save_nfs;
  CreateSavePath((string) save_nfs+"SYS");

  // ================================================================================
  // === I load all histograms; This code was created after studies ended. Only one
  // === option is available. If necessary include other options as well (no bins etc.)

  cout << "Start collecting Histograms ..." << endl;

  TString h = "mW__";
  TFile *file = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/macros/RecoLevelPlots.root");
  TH1F* h_data = (TH1F*) file->Get(h+"DATA");
  TH1F* h_ttbar = (TH1F*) file->Get(h+"TTbar");
  TH1F* h_st = (TH1F*) file->Get(h+"SingleTop");
  TH1F* h_wjets = (TH1F*) file->Get(h+"WJets");
  TH1F* h_other = (TH1F*) file->Get(h+"other");

  cout << "\t ... Subtract bkg from data" << endl;
  TH1F* h_bkg = (TH1F*) h_st->Clone();
  h_bkg->Add(h_wjets, 1);
  h_bkg->Add(h_other, 1);
  TH1F* h_data_sub = (TH1F*) h_data->Clone();
  h_data_sub->Add(h_bkg, -1);

  for(unsigned int i=1;i<h_data_sub->GetNbinsX();i++){
    cout << setw(2) << i << setw(15) << h_data_sub->GetBinContent(i) << setw(15) << h_data_sub->GetBinError(i) << endl;
  }

  cout << setw(25) << "h_data" << setw(20) << h_data->GetBinError(5) << setw(20) << h_data->GetBinContent(5) << endl;
  cout << setw(25) << "h_ttbar" << setw(20) << h_ttbar->GetBinError(5) << setw(20) << h_ttbar->GetBinContent(5) << endl;
  cout << setw(25) << "h_st" << setw(20) << h_st->GetBinError(5) << setw(20) << h_st->GetBinContent(5) << endl;
  cout << setw(25) << "h_wjets" << setw(20) << h_wjets->GetBinError(5) << setw(20) << h_wjets->GetBinContent(5) << endl;
  cout << setw(25) << "h_other" << setw(20) << h_other->GetBinError(5) << setw(20) << h_other->GetBinContent(5) << endl;
  cout << setw(25) << "h_bkg" << setw(20) << h_bkg->GetBinError(5) << setw(20) << h_bkg->GetBinContent(5) << endl;
  cout << setw(25) << "h_data_sub" << setw(20) << h_data_sub->GetBinError(5) << setw(20) << h_data_sub->GetBinContent(5) << endl;


  cout << "\t ... Get systematics" << endl;
  TH1F* h_jec_up = (TH1F*) file->Get(h+"TTbar__JEC__plus");
  TH1F* h_jms_up = (TH1F*) file->Get(h+"TTbar__JMS__plus");
  TH1F* h_jmsf_up = (TH1F*) file->Get(h+"TTbar__JMS_flavor__plus");
  TH1F* h_cor_up = (TH1F*) file->Get(h+"TTbar__COR__plus");
  TH1F* h_jer_up = (TH1F*) file->Get(h+"TTbar__JER__plus");
  TH1F* h_pre_up = (TH1F*) file->Get(h+"TTbar__prefire__plus");
  TH1F* h_pu_up = (TH1F*) file->Get(h+"TTbar__pu__plus");
  TH1F* h_btag_up = (TH1F*) file->Get(h+"TTbar__btag__plus");

  TH1F* h_jec_down = (TH1F*) file->Get(h+"TTbar__JEC__minus");
  TH1F* h_jms_down = (TH1F*) file->Get(h+"TTbar__JMS__minus");
  TH1F* h_jmsf_down = (TH1F*) file->Get(h+"TTbar__JMS_flavor__minus");
  TH1F* h_cor_down = (TH1F*) file->Get(h+"TTbar__COR__minus");
  TH1F* h_jer_down = (TH1F*) file->Get(h+"TTbar__JER__minus");
  TH1F* h_pre_down = (TH1F*) file->Get(h+"TTbar__prefire__minus");
  TH1F* h_pu_down = (TH1F*) file->Get(h+"TTbar__pu__minus");
  TH1F* h_btag_down = (TH1F*) file->Get(h+"TTbar__btag__minus");

  int mDim = h_data->GetNbinsX();
  int ndf_peak = 9;
  int ndf = 35;

  // ================================================================================
  // === Build all necessary TObjects

  cout << "Start building covariance matrix ..." << endl;
  TH2F* m_cov_tt = BuildCovMatrix(h_ttbar, "tt");
  TH2F* m_cov_data = BuildCovMatrix(h_data_sub, "data");

  cout << "\t ... build uncertainty hists" << endl;
  TH1F* h_sys_jec = BuildSYSHist(h_jec_up, h_ttbar, "jec");
  TH1F* h_sys_cor = BuildSYSHist(h_cor_up, h_ttbar, "cor");
  TH1F* h_sys_jer = BuildSYSHist(h_jer_up, h_ttbar, "jer");
  TH1F* h_sys_jms = BuildSYSHist(h_jms_up, h_ttbar, "jms");
  TH1F* h_sys_jmsf = BuildSYSHist(h_jmsf_up, h_ttbar, "jmsf");
  TH1F* h_sys_pre = BuildSYSHist(h_pre_up, h_ttbar, "pre");
  TH1F* h_sys_pu = BuildSYSHist(h_pu_up, h_ttbar, "pu");
  TH1F* h_sys_btag = BuildSYSHist(h_btag_up, h_ttbar, "btag");

  cout << "\t ... build sys cov matrix" << endl;
  TH2F* m_cov_sys_jec = BuildSYSCovMatrix(h_sys_jec, "jec");
  TH2F* m_cov_sys_cor = BuildSYSCovMatrix(h_sys_cor, "cor");
  TH2F* m_cov_sys_jer = BuildSYSCovMatrix(h_sys_jer, "jer");
  TH2F* m_cov_sys_jms = BuildSYSCovMatrix(h_sys_jms, "jms");
  TH2F* m_cov_sys_jmsf = BuildSYSCovMatrix(h_sys_jmsf, "jmsf");
  TH2F* m_cov_sys_pre = BuildSYSCovMatrix(h_sys_pre, "pre");
  TH2F* m_cov_sys_pu = BuildSYSCovMatrix(h_sys_pu, "pu");
  TH2F* m_cov_sys_btag = BuildSYSCovMatrix(h_sys_btag, "btag");

  cout << "Start building normalized covariance matrix ..." << endl;
  TH2F* m_cov_tt_norm = NormCovMatrix(m_cov_tt, h_ttbar, "tt_norm");
  TH2F* m_cov_data_norm = NormCovMatrix(m_cov_data, h_data_sub, "data_norm");

  cout << "\t ... normalize hists" << endl;
  TH1F *h_data_norm = Normalize(h_data_sub);
  TH1F *h_ttbar_norm = Normalize(h_ttbar);
  TH1F *h_jec_up_norm = Normalize(h_jec_up);
  TH1F *h_jms_up_norm = Normalize(h_jms_up);
  TH1F *h_jmsf_up_norm = Normalize(h_jmsf_up);
  TH1F *h_cor_up_norm = Normalize(h_cor_up);
  TH1F *h_jer_up_norm = Normalize(h_jer_up);
  TH1F *h_pre_up_norm = Normalize(h_pre_up);
  TH1F *h_pu_up_norm = Normalize(h_pu_up);
  TH1F *h_btag_up_norm = Normalize(h_btag_up);

  cout << "\t ... build norm uncertainty hists" << endl;
  TH1F* h_sys_jec_norm = BuildSYSHist(h_jec_up_norm, h_ttbar_norm, "jec_norm");
  TH1F* h_sys_cor_norm = BuildSYSHist(h_cor_up_norm, h_ttbar_norm, "cor_norm");
  TH1F* h_sys_jer_norm = BuildSYSHist(h_jer_up_norm, h_ttbar_norm, "jer_norm");
  TH1F* h_sys_jms_norm = BuildSYSHist(h_jms_up_norm, h_ttbar_norm, "jms_norm");
  TH1F* h_sys_jmsf_norm = BuildSYSHist(h_jmsf_up_norm, h_ttbar_norm, "jmsf_norm");
  TH1F* h_sys_pre_norm = BuildSYSHist(h_pre_up_norm, h_ttbar_norm, "pre_norm");
  TH1F* h_sys_pu_norm = BuildSYSHist(h_pu_up_norm, h_ttbar_norm, "pu_norm");
  TH1F* h_sys_btag_norm = BuildSYSHist(h_btag_up_norm, h_ttbar_norm, "btag_norm");

  vector<TH2F*> m_cov_sys_all = {m_cov_sys_jec, m_cov_sys_cor, m_cov_sys_jer, m_cov_sys_jms, m_cov_sys_jmsf, m_cov_sys_pre, m_cov_sys_pu, m_cov_sys_btag};
  TH2F* m_cov_sys_total= BuildTotalCovMatrix(m_cov_sys_all, "Full sys");
  TH2F* m_cov_total = BuildTotalCovMatrix({m_cov_sys_total,m_cov_data,m_cov_tt}, "Full total");


  cout << "\t ... build norm sys cov matrix" << endl;
  TH2F* m_cov_sys_jec_norm = BuildSYSCovMatrix(h_sys_jec_norm, "jec_norm");
  TH2F* m_cov_sys_cor_norm = BuildSYSCovMatrix(h_sys_cor_norm, "cor_norm");
  TH2F* m_cov_sys_jer_norm = BuildSYSCovMatrix(h_sys_jer_norm, "jer_norm");
  TH2F* m_cov_sys_jms_norm = BuildSYSCovMatrix(h_sys_jms_norm, "jms_norm");
  TH2F* m_cov_sys_jmsf_norm = BuildSYSCovMatrix(h_sys_jmsf_norm, "jmsf_norm");
  TH2F* m_cov_sys_pre_norm = BuildSYSCovMatrix(h_sys_pre_norm, "pre_norm");
  TH2F* m_cov_sys_pu_norm = BuildSYSCovMatrix(h_sys_pu_norm, "pu_norm");
  TH2F* m_cov_sys_btag_norm = BuildSYSCovMatrix(h_sys_btag_norm, "btag_norm");

  vector<TH2F*> m_cov_sys_norm_all = {m_cov_sys_jec_norm, m_cov_sys_cor_norm, m_cov_sys_jer_norm, m_cov_sys_jms_norm, m_cov_sys_jmsf_norm, m_cov_sys_pre_norm, m_cov_sys_pu_norm, m_cov_sys_btag_norm};
  TH2F* m_cov_sys_norm_total= BuildTotalCovMatrix(m_cov_sys_norm_all, "Full sys");
  TH2F* m_cov_total_norm = BuildTotalCovMatrix({m_cov_sys_norm_total,m_cov_data_norm,m_cov_tt_norm}, "Full total");

  // ================================================================================
  // === Sub study to cross check results

  cout << "Start building trimmed hists ..." << endl;
  double r_low = 65;
  double r_high = 110;

  cout << "\t ... trim non-normalized" << endl;
  TH1F* h_tt_trim = TrimHistogram(h_ttbar, r_low, r_high, "h_tt_trim");
  TH1F* h_data_trim = TrimHistogram(h_data_sub, r_low, r_high, "h_data_trim");
  TH2F* m_cov_tt_trim = BuildCovMatrix(h_tt_trim, "m_tt_trim");
  TH2F* m_cov_data_trim = BuildCovMatrix(h_data_trim, "m_data_trim");

  cout << "\t ... trim non-normalized sys" << endl;
  TH1F* h_jec_trim = TrimHistogram(h_jec_up, r_low, r_high, "h_jec_trim");
  TH1F* h_jms_trim = TrimHistogram(h_jms_up, r_low, r_high, "h_jms_trim");
  TH1F* h_jmsf_trim = TrimHistogram(h_jmsf_up, r_low, r_high, "h_jmsf_trim");
  TH1F* h_cor_trim = TrimHistogram(h_cor_up, r_low, r_high, "h_cor_trim");
  TH1F* h_jer_trim = TrimHistogram(h_jer_up, r_low, r_high, "h_jer_trim");
  TH1F* h_pre_trim = TrimHistogram(h_pre_up, r_low, r_high, "h_pre_trim");
  TH1F* h_pu_trim = TrimHistogram(h_pu_up, r_low, r_high, "h_pu_trim");
  TH1F* h_btag_trim = TrimHistogram(h_btag_up, r_low, r_high, "h_btag_trim");

  cout << "\t ... build uncertainty hists" << endl;
  TH1F* h_sys_jec_trim = BuildSYSHist(h_jec_trim, h_tt_trim, "h_sys_jec_trim");
  TH1F* h_sys_cor_trim = BuildSYSHist(h_cor_trim, h_tt_trim, "h_sys_cor_trim");
  TH1F* h_sys_jer_trim = BuildSYSHist(h_jer_trim, h_tt_trim, "h_sys_jer_trim");
  TH1F* h_sys_jms_trim = BuildSYSHist(h_jms_trim, h_tt_trim, "h_sys_jms_trim");
  TH1F* h_sys_jmsf_trim = BuildSYSHist(h_jmsf_trim , h_tt_trim, "h_sys_jmsf_trim");
  TH1F* h_sys_pre_trim = BuildSYSHist(h_pre_trim, h_tt_trim, "h_sys_pre_trim");
  TH1F* h_sys_pu_trim = BuildSYSHist(h_pu_trim, h_tt_trim, "h_sys_pu_trim");
  TH1F* h_sys_btag_trim = BuildSYSHist(h_btag_trim, h_tt_trim, "h_sys_btag_trim");

  cout << "\t ... build sys cov matrix" << endl;
  TH2F* m_cov_sys_jec_trim = BuildSYSCovMatrix(h_sys_jec_trim, "m_sys_jec_trim");
  TH2F* m_cov_sys_cor_trim = BuildSYSCovMatrix(h_sys_cor_trim, "m_sys_cor_trim");
  TH2F* m_cov_sys_jer_trim = BuildSYSCovMatrix(h_sys_jer_trim, "m_sys_jer_trim");
  TH2F* m_cov_sys_jms_trim = BuildSYSCovMatrix(h_sys_jms_trim, "m_sys_jms_trim");
  TH2F* m_cov_sys_jmsf_trim = BuildSYSCovMatrix(h_sys_jmsf_trim, "m_sys_jmsf_trim");
  TH2F* m_cov_sys_pre_trim = BuildSYSCovMatrix(h_sys_pre_trim, "m_sys_pre_trim");
  TH2F* m_cov_sys_pu_trim = BuildSYSCovMatrix(h_sys_pu_trim, "m_sys_pu_trim");
  TH2F* m_cov_sys_btag_trim = BuildSYSCovMatrix(h_sys_btag_trim, "m_sys_btag_trim");
  vector<TH2F*> m_cov_sys_trim_all = {m_cov_sys_jec_trim, m_cov_sys_cor_trim, m_cov_sys_jer_trim, m_cov_sys_jms_trim, m_cov_sys_jmsf_trim, m_cov_sys_pre_trim, m_cov_sys_pu_trim, m_cov_sys_btag_trim};
  TH2F* m_cov_sys_total_trim = BuildTotalCovMatrix(m_cov_sys_trim_all, "trim");
  TH2F* m_cov_total_trim = BuildTotalCovMatrix({m_cov_sys_total_trim,m_cov_data_trim,m_cov_tt_trim}, "trim");
  TH2F* m_cov_total_trim_norm_added = NormCovMatrix(m_cov_total_trim, h_data_trim, "trim_norm_added");

  cout << "\t ... normalize trim" << endl;
  TH1F* h_tt_trim_norm = Normalize(h_tt_trim);
  TH1F* h_data_trim_norm = Normalize(h_data_trim);
  TH2F* m_cov_tt_trim_norm = NormCovMatrix(m_cov_tt_trim, h_tt_trim, "m_t_trim_norm");
  TH2F* m_cov_data_trim_norm = NormCovMatrix(m_cov_data_trim, h_data_trim, "m_data_trim_norm");
  TH2F* m_cov_tt_trim_norm_v2 = BuildCovMatrix(h_tt_trim_norm, "m_tt_trim_v2");

  cout << "\t ... normalize trim sys" << endl;
  TH1F *h_jec_trim_norm = Normalize(h_jec_trim);
  TH1F *h_jms_trim_norm = Normalize(h_jms_trim);
  TH1F *h_jmsf_trim_norm = Normalize(h_jmsf_trim);
  TH1F *h_cor_trim_norm = Normalize(h_cor_trim);
  TH1F *h_jer_trim_norm = Normalize(h_jer_trim);
  TH1F *h_pre_trim_norm = Normalize(h_pre_trim);
  TH1F *h_pu_trim_norm = Normalize(h_pu_trim);
  TH1F *h_btag_trim_norm = Normalize(h_btag_trim);

  cout << "\t ... build sys hist from trim sys" << endl;
  TH1F* h_sys_jec_trim_norm = BuildSYSHist(h_jec_trim_norm, h_tt_trim_norm, "h_sys_jec_trim");
  TH1F* h_sys_cor_trim_norm = BuildSYSHist(h_cor_trim_norm, h_tt_trim_norm, "h_sys_cor_trim");
  TH1F* h_sys_jer_trim_norm = BuildSYSHist(h_jer_trim_norm, h_tt_trim_norm, "h_sys_jer_trim");
  TH1F* h_sys_jms_trim_norm = BuildSYSHist(h_jms_trim_norm, h_tt_trim_norm, "h_sys_jms_trim");
  TH1F* h_sys_jmsf_trim_norm = BuildSYSHist(h_jmsf_trim_norm, h_tt_trim_norm, "h_sys_jmsf_trim");
  TH1F* h_sys_pre_trim_norm = BuildSYSHist(h_pre_trim_norm, h_tt_trim_norm, "h_sys_pre_trim");
  TH1F* h_sys_pu_trim_norm = BuildSYSHist(h_pu_trim_norm, h_tt_trim_norm, "h_sys_pu_trim");
  TH1F* h_sys_btag_trim_norm = BuildSYSHist(h_btag_trim_norm, h_tt_trim_norm, "h_sys_btag_trim");

  cout << "\t ... build normalized trim cov matrix from sys" << endl;
  TH2F* m_cov_sys_jec_trim_norm = BuildSYSCovMatrix(h_sys_jec_trim_norm, "m_sys_jec_trim_norm");
  TH2F* m_cov_sys_cor_trim_norm = BuildSYSCovMatrix(h_sys_cor_trim_norm, "m_sys_cor_trim_norm");
  TH2F* m_cov_sys_jer_trim_norm = BuildSYSCovMatrix(h_sys_jer_trim_norm, "m_sys_jer_trim_norm");
  TH2F* m_cov_sys_jms_trim_norm = BuildSYSCovMatrix(h_sys_jms_trim_norm, "m_sys_jms_trim_norm");
  TH2F* m_cov_sys_jmsf_trim_norm = BuildSYSCovMatrix(h_sys_jmsf_trim_norm, "m_sys_jmsf_trim_norm");
  TH2F* m_cov_sys_pre_trim_norm = BuildSYSCovMatrix(h_sys_pre_trim_norm, "m_sys_pre_trim_norm");
  TH2F* m_cov_sys_pu_trim_norm = BuildSYSCovMatrix(h_sys_pu_trim_norm, "m_sys_pu_trim_norm");
  TH2F* m_cov_sys_btag_trim_norm = BuildSYSCovMatrix(h_sys_btag_trim_norm, "m_sys_btag_trim_norm");

  vector<TH2F*> v_m_cov_sys_trim_norm = {m_cov_sys_jec_trim_norm, m_cov_sys_cor_trim_norm, m_cov_sys_jer_trim_norm, m_cov_sys_jms_trim_norm, m_cov_sys_jmsf_trim_norm, m_cov_sys_pre_trim_norm, m_cov_sys_pu_trim_norm, m_cov_sys_btag_trim_norm};
  TH2F* m_cov_sys_trim_norm = BuildTotalCovMatrix(v_m_cov_sys_trim_norm, "sys_trim_norm");
  TH2F* m_cov_total_trim_norm = BuildTotalCovMatrix({m_cov_sys_trim_norm, m_cov_data_trim_norm, m_cov_tt_trim_norm}, "total_trim_norm");

  // ================================================================================
  // === Get Chi2

  cout << endl;
  double chi2_cut = GetChi2Cut(h_data_sub, h_ttbar, m_cov_total, 65, 110, "Full cut", false);
  double chi2_cut_norm = GetChi2Cut(h_data_norm, h_ttbar_norm, m_cov_total_norm, 65, 110, "Full cut norm", true);
  double chi2 = GetChi2(h_data_sub, h_ttbar, m_cov_total, "Full", false);
  double chi2_norm = GetChi2(h_data_norm, h_ttbar_norm, m_cov_total_norm, "Full norm", true);
  double chi2_trim = GetChi2(h_data_trim, h_tt_trim, m_cov_total_trim, "Peak", false);
  double chi2_trim_norm_added = GetChi2(h_data_trim_norm, h_tt_trim_norm, m_cov_total_trim_norm_added, "Peak norm cov total", true);
  double chi2_trim_norm = GetChi2(h_data_trim_norm, h_tt_trim_norm, m_cov_total_trim_norm, "Peak normed", true);

  // ================================================================================
  // === Hists ans covs from dennis

  // TFile *in = new TFile("/afs/desy.de/user/p/paaschal/WorkingArea/output.root");
  // TH2F* c_dennis = (TH2F*) in->Get("cov");
  // TH1F* h_data_dennis =  (TH1F*) in->Get("data");
  // TH1F* h_ttbar_dennis =  (TH1F*) in->Get("ttbar");
  //
  // TH2F* c_dennis_norm = (TH2F*) in->Get("cov_norm");
  // TH1F* h_data_dennis_norm =  (TH1F*) in->Get("data_norm");
  // TH1F* h_ttbar_dennis_norm =  (TH1F*) in->Get("ttbar_norm");
  //
  // double chi2_dennis = GetChi2(h_data_dennis, h_ttbar_dennis, c_dennis, "dennis", false);
  // double chi2_dennis_norm = GetChi2(h_data_dennis_norm, h_ttbar_dennis_norm, c_dennis_norm, "dennis_norm", true);


  cout << "\nEnding programm ..." << endl;
  return 0;
}

//   // // ================================================================================
//   // // === Draw with ratio
//   //
//   // TH1F* ratio_data    = GetRatio(h_data_norm,    h_ttbar_norm, false);
//   // TH1F* ratio_ttbar   = GetRatio(h_ttbar_norm,   h_ttbar_norm, true);
//   // TH1F* ratio_JECup   = GetRatio(h_JECup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_JECdown = GetRatio(h_JECdown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_CORup   = GetRatio(h_CORup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_CORdown = GetRatio(h_CORdown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_JMSup   = GetRatio(h_JMSup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_JMSdown = GetRatio(h_JMSdown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_JMSflavorup   = GetRatio(h_JMSflavorup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_JMSflavordown = GetRatio(h_JMSflavordown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_JERup   = GetRatio(h_JERup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_JERdown = GetRatio(h_JERdown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_prefireup   = GetRatio(h_prefireup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_prefiredown = GetRatio(h_prefiredown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_puup   = GetRatio(h_puup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_pudown = GetRatio(h_pudown_norm, h_ttbar_norm, false);
//   // TH1F* ratio_btagup   = GetRatio(h_btagup_norm,   h_ttbar_norm, false);
//   // TH1F* ratio_btagdown = GetRatio(h_btagdown_norm, h_ttbar_norm, false);
//   //
//   // // TH1F* h_ratio_sys_up = (TH1F*) ratio_JECup->Clone();
//   // // h_ratio_sys_up->Add(ratio_JECup, 1);
//   // // h_ratio_sys_up->Add(ratio_JECup, 1);
//   // // h_ratio_sys_up->Add(ratio_JECup, 1);
//   // // h_ratio_sys_up->Add(ratio_JECup, 1);
//   //
//   // SetupCanvas(true);
//   // // void CosmeticsLocal(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
//   // Cosmetics(h_data_norm,    "#it{m}_{W} [GeV]", kBlack, kSolid,  8, 2, true);
//   // Cosmetics(h_ttbar_norm,   "#it{m}_{W} [GeV]", kBlack, kSolid,  0, 2, true);
//   // Cosmetics(h_JECup_norm,   "#it{m}_{W} [GeV]", kAzure-3, kDashed,  0, 2, true);
//   // Cosmetics(h_JECdown_norm, "#it{m}_{W} [GeV]", kAzure-3, kSolid, 0, 2, true);
//   // Cosmetics(h_JERup_norm,   "#it{m}_{W} [GeV]", kGreen-7, kDashed,  0, 2, true);
//   // Cosmetics(h_JERdown_norm, "#it{m}_{W} [GeV]", kGreen-7, kSolid, 0, 2, true);
//   // Cosmetics(h_CORup_norm,   "#it{m}_{W} [GeV]", kYellow+1, kDashed, 0, 2, true);
//   // Cosmetics(h_CORdown_norm, "#it{m}_{W} [GeV]", kYellow+1,   kSolid,  0, 2, true);
//   // Cosmetics(h_JMSup_norm,   "#it{m}_{W} [GeV]", kRed+1, kDashed, 0, 2, true);
//   // Cosmetics(h_JMSdown_norm, "#it{m}_{W} [GeV]", kRed+1,   kSolid,  0, 2, true);
//   // Cosmetics(h_JMSflavorup_norm,   "#it{m}_{W} [GeV]", kOrange-7, kDashed, 0, 2, true);
//   // Cosmetics(h_JMSflavordown_norm, "#it{m}_{W} [GeV]", kOrange-7,   kSolid,  0, 2, true);
//   // Cosmetics(h_prefireup_norm,   "#it{m}_{W} [GeV]", kGreen+1, kDashed, 0, 2, true);
//   // Cosmetics(h_prefiredown_norm, "#it{m}_{W} [GeV]", kGreen+1,   kSolid,  0, 2, true);
//   // Cosmetics(h_puup_norm,   "#it{m}_{W} [GeV]", kGreen+1, kDashed, 0, 2, true);
//   // Cosmetics(h_pudown_norm, "#it{m}_{W} [GeV]", kGreen+1,   kSolid,  0, 2, true);
//   // Cosmetics(h_btagup_norm,   "#it{m}_{W} [GeV]", kGreen+1, kDashed, 0, 2, true);
//   // Cosmetics(h_btagdown_norm, "#it{m}_{W} [GeV]", kGreen+1,   kSolid,  0, 2, true);
//   //
//   // // void RatioCosmeticsLocal(TH1* hist,TString xtitle, int color, int style, int width, int marker)
//   // RatioCosmetics(ratio_ttbar,   "#it{m}_{W} [GeV]", kBlack,  kSolid,  2, 0);
//   // RatioCosmetics(ratio_JECup,   "#it{m}_{W} [GeV]", kAzure-3, kSolid,  2, 0);
//   // RatioCosmetics(ratio_JECdown, "#it{m}_{W} [GeV]", kAzure-3, kDashed,  2, 0);
//   // RatioCosmetics(ratio_JERup,   "#it{m}_{W} [GeV]", kGreen-7, kSolid, 2, 0);
//   // RatioCosmetics(ratio_JERdown, "#it{m}_{W} [GeV]", kGreen-7, kDashed, 2, 0);
//   // RatioCosmetics(ratio_CORup,   "#it{m}_{W} [GeV]", kYellow+1, kSolid, 2, 0);
//   // RatioCosmetics(ratio_CORdown, "#it{m}_{W} [GeV]", kYellow+1, kDashed, 2, 0);
//   // RatioCosmetics(ratio_JMSup,   "#it{m}_{W} [GeV]", kRed+1, kSolid, 2, 0);
//   // RatioCosmetics(ratio_JMSdown, "#it{m}_{W} [GeV]", kRed+1, kDashed, 2, 0);
//   // RatioCosmetics(ratio_JMSflavorup,   "#it{m}_{W} [GeV]", kOrange-7, kSolid, 2, 0);
//   // RatioCosmetics(ratio_JMSflavordown, "#it{m}_{W} [GeV]", kOrange-7, kDashed, 2, 0);
//   // RatioCosmetics(ratio_prefireup,   "#it{m}_{W} [GeV]", kGreen+1, kSolid, 2, 0);
//   // RatioCosmetics(ratio_prefiredown, "#it{m}_{W} [GeV]", kGreen+1, kDashed, 2, 0);
//   // RatioCosmetics(ratio_puup,   "#it{m}_{W} [GeV]", kGreen+1, kSolid, 2, 0);
//   // RatioCosmetics(ratio_pudown, "#it{m}_{W} [GeV]", kGreen+1, kDashed, 2, 0);
//   // RatioCosmetics(ratio_btagup,   "#it{m}_{W} [GeV]", kGreen+1, kSolid, 2, 0);
//   // RatioCosmetics(ratio_btagdown, "#it{m}_{W} [GeV]", kGreen+1, kDashed, 2, 0);
//   // RatioCosmetics(ratio_data,    "#it{m}_{W} [GeV]", kBlack,   kSolid,  2, 8);
//   //
//   // h_ttbar_norm->GetYaxis()->SetTitleSize(0.08);
//   // h_ttbar_norm->GetYaxis()->SetTitleOffset(1.05);
//   // h_ttbar_norm->GetYaxis()->SetLabelOffset(0.01);
//   // ratio_ttbar->GetYaxis()->SetTitleSize(0.15);
//   // ratio_ttbar->GetYaxis()->SetTitleOffset(0.58);
//   // ratio_ttbar->GetYaxis()->SetLabelOffset(0.011);
//   //
//   // m_rp1_top->cd();
//   // h_ttbar_norm->GetYaxis()->SetRangeUser(0.001, h_ttbar_norm->GetMaximum()*1.2);
//   // h_ttbar_norm->Draw("hist ][");
//   // h_JECup_norm->Draw("hist same ][");
//   // h_JECdown_norm->Draw("hist same ][");
//   // h_CORup_norm->Draw("hist same ][");
//   // h_CORdown_norm->Draw("hist same ][");
//   // h_JMSup_norm->Draw("hist same ][");
//   // h_JMSdown_norm->Draw("hist same ][");
//   // h_JMSflavorup_norm->Draw("hist same ][");
//   // h_JMSflavordown_norm->Draw("hist same ][");
//   // h_JERup_norm->Draw("hist same ][");
//   // h_JERdown_norm->Draw("hist same ][");
//   // // h_prefireup_norm->Draw("hist same ][");
//   // // h_prefiredown_norm->Draw("hist same ][");
//   // h_puup_norm->Draw("hist same ][");
//   // h_pudown_norm->Draw("hist same ][");
//   // // h_btagup_norm->Draw("hist same ][");
//   // // h_btagdown_norm->Draw("hist same ][");
//   // h_ttbar_norm->Draw("axis same ][");
//   // h_data_norm->Draw("pe same ][");
//   // gPad->RedrawAxis();
//   // DrawLegend({h_data_norm, h_ttbar_norm, h_JECup_norm, h_JECdown_norm, h_CORup_norm, h_CORdown_norm}, "wmass");
//   // // leg->SetTextSize(0.065);
//   //
//   // // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
//   // DrawLumi(138., true, true, false, "wmass");
//   // // if(bin.Contains("hl")||bin.Contains("hh")) WMassLabel(0.68, 0.87, bin);
//   // // else                                       WMassLabel(0.22, 0.78, bin);
//   //
//   // m_rp1->cd();
//   // ratio_ttbar->Draw("hist ][");
//   // ratio_JECup->Draw("hist same ][");
//   // ratio_CORup->Draw("hist same ][");
//   // ratio_JERup->Draw("hist same ][");
//   // ratio_JMSup->Draw("hist same ][");
//   // ratio_JMSflavorup->Draw("hist same ][");
//   // // ratio_prefireup->Draw("hist same ][");
//   // ratio_puup->Draw("hist same ][");
//   // // ratio_btagup->Draw("hist same ][");
//   // ratio_JECdown->Draw("hist same ][");
//   // ratio_CORdown->Draw("hist same ][");
//   // ratio_JERdown->Draw("hist same ][");
//   // ratio_JMSdown->Draw("hist same ][");
//   // ratio_JMSflavordown->Draw("hist same ][");
//   // // ratio_prefiredown->Draw("hist same ][");
//   // ratio_pudown->Draw("hist same ][");
//   // // ratio_btagdown->Draw("hist same ][");
//   // ratio_data->Draw("pe same ][");
//   // gPad->RedrawAxis();
//   //
//   // // m_can->Print(save_afs+"wmass_"+bin+".pdf");
//   // m_can->Print(save_nfs+"wmass.pdf");
//
//   return 0;
//
// }
//
//
// // --------------------------------------------------------------------------------------------------------------------------------------------
// // --------------------------------------------------------------------------------------------------------------------------------------------
// // --------------------------------------------------------------------------------------------------------------------------------------------
// // --------------------------------------------------------------------------------------------------------------------------------------------
// // --------------------------------------------------------------------------------------------------------------------------------------------
// // Define functions
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
// MapHH GetHistograms(TString process, TString h_name){
//   VecTS collection_muon, collection_elec;
//
//   if(debug) cout << "\t\t ... " << process << " - " << h_name << endl;
//
//   if(process.EqualTo("data")){ collection_muon = data_muon; collection_elec = data_elec;}
//   else if(process.EqualTo("ttbar")){ collection_muon = ttbar_muon; collection_elec = ttbar_elec;}
//   else if(process.EqualTo("wjets")){ collection_muon = wjets_muon; collection_elec = wjets_elec;}
//   else if(process.EqualTo("st")){ collection_muon = st_muon; collection_elec = st_elec;}
//   else if(process.EqualTo("other")){ collection_muon = other_muon; collection_elec = other_elec;}
//   else if(process.EqualTo("JECup")){ collection_muon = jec_up_muon; collection_elec = jec_up_elec;}
//   else if(process.EqualTo("JECdown")){ collection_muon = jec_down_muon; collection_elec = jec_down_elec;}
//   else if(process.EqualTo("JMSup")){ collection_muon = jms_uu_muon; collection_elec = jms_uu_elec;}
//   else if(process.EqualTo("JMSdown")){ collection_muon = jms_dd_muon; collection_elec = jms_dd_elec;}
//   else if(process.EqualTo("CORup")){ collection_muon = cor_up_muon; collection_elec = cor_up_elec;}
//   else if(process.EqualTo("CORdown")){ collection_muon = cor_down_muon; collection_elec = cor_down_elec;}
//   else if(process.EqualTo("JERup")){ collection_muon = jer_up_muon; collection_elec = jer_up_elec;}
//   else if(process.EqualTo("JERdown")){ collection_muon = jer_down_muon; collection_elec = jer_down_elec;}
//   else if(process.EqualTo("hdampup")){ collection_muon = hdamp_up_muon; collection_elec = hdamp_up_elec;}
//   else if(process.EqualTo("hdampdown")){ collection_muon = hdamp_down_muon; collection_elec = hdamp_down_elec;}
//   else if(process.EqualTo("tuneup")){ collection_muon = tune_up_muon; collection_elec = tune_up_elec;}
//   else if(process.EqualTo("tunedown")){ collection_muon = tune_down_muon; collection_elec = tune_down_elec;}
//   else if(process.EqualTo("CRgm")){ collection_muon = cr_gm_muon; collection_elec = cr_gm_elec;}
//   else if(process.EqualTo("CRqb")){ collection_muon = cr_qb_muon; collection_elec = cr_qb_elec;}
//   else if(process.EqualTo("fsrup")){ collection_muon = fsr_up_2_muon; collection_elec = fsr_up_2_elec;}
//   else if(process.EqualTo("fsrdown")){ collection_muon = fsr_down_2_muon; collection_elec = fsr_down_2_elec;}
//   else if(process.EqualTo("isrup")){ collection_muon = isr_up_2_muon; collection_elec = isr_up_2_elec;}
//   else if(process.EqualTo("isrdown")){ collection_muon = isr_down_2_muon; collection_elec = isr_down_2_elec;}
//   else throw runtime_error("Check the process ´"+process+"´ to obtain the histograms");
//
//   MapHH map;
//   vector<TH1F*> h_muon = get_all_hists(collection_muon, h_name);
//   vector<TH1F*> h_elec = get_all_hists(collection_elec, h_name);
//   vector<TH1F*> h_combine = combine_channels(h_muon, h_elec);
//
//   for(auto nhist:h_muon) nhist->SetTitle(h_name);
//   for(auto nhist:h_elec) nhist->SetTitle(h_name);
//   for(auto nhist:h_combine) nhist->SetTitle(h_name);
//
//   map["muon"]["2016"] = h_muon[0];
//   map["muon"]["2017"] = h_muon[1];
//   map["muon"]["2018"] = h_muon[2];
//   map["muon"]["combine"] = h_muon[3];
//
//   map["elec"]["2016"] = h_elec[0];
//   map["elec"]["2017"] = h_elec[1];
//   map["elec"]["2018"] = h_elec[2];
//   map["elec"]["combine"] = h_elec[3];
//
//   map["combine"]["2016"] = h_combine[0];
//   map["combine"]["2017"] = h_combine[1];
//   map["combine"]["2018"] = h_combine[2];
//   map["combine"]["combine"] = h_combine[3];
//
//   return map;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
// TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize){
//   int nbins = data->GetXaxis()->GetNbins();
//   TH1F* result = (TH1F*) data->Clone();
//   for(unsigned int i=0; i<bgr.size(); i++){
//     result->Add(bgr[i], -1);
//   }
//   for(int bin=1; bin<=nbins; bin++){
//     double syserror2 = 0;
//     for(unsigned int i=0; i<bgr.size(); i++){
//       syserror2 += pow(bgr[i]->GetBinContent(bin) * syssize[i], 2);
//     }
//     double olderror2 = pow(result->GetBinError(bin), 2);
//     // cout << bin << "\t" << sqrt(syserror2) << "\t" << sqrt(olderror2) << "\t" << sqrt(syserror2+olderror2) << endl;
//     result->SetBinError(bin, sqrt(syserror2+olderror2));
//   }
//   for(int bin=1; bin<=nbins; bin++){
//     if(result->GetBinContent(bin)<0){
//       result->SetBinContent(bin, 0);
//       result->SetBinError(bin, 0);
//     }
//   }
//
//   return result;
// }
//
// void SubtractBackgroundsMap(MapHH& m_data, MapHH& m_st, MapHH& m_wjets, MapHH& m_other){
//   for(auto channel: m_data){
//     TString c = channel.first;
//     for(auto year: channel.second){
//       TString y = year.first;
//       vector<TH1F*> bkgs = {m_st[c][y], m_wjets[c][y], m_other[c][y]};
//       vector<double> bkgrate = {0.23, 0.19, 1.00};
//       TH1F* hist = SubtractBackgrounds(m_data[c][y], bkgs, bkgrate);
//       m_data[c][y] = hist;
//     }
//   }
// }
//
// void AddTwoMaps(MapHHH& map1, MapHHH& map2, int option){
//   // MapHHH map_out = map1;
//   for(auto hist: map1){
//     TString h = hist.first;
//     for(auto channel: map1[h]){
//       TString c = channel.first;
//       for(auto year: map1[h][c]){
//         TString y = year.first;
//         // if(y.Contains("combine")&&c.Contains("combine")) cout << h << " - " << map1[h][c][y]->GetBinContent(90) << " & " << map2[h][c][y]->GetBinContent(90);
//         map1[h][c][y]->Add(map2[h][c][y], option);
//         // if(y.Contains("combine")&&c.Contains("combine")) cout << " -> " << map1[h][c][y]->GetBinContent(90) << endl;;
//       }
//     }
//   }
//   // return map_out;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
//
// TMatrixD GetCovMatrix(TH1F* hist, const TString& save, const TString& process){
//
//   int dim = hist->GetNbinsX();
//
//   TString c = channel; TString y = year;
//   TMatrixD mtemp = TMatrixD(dim, dim);
//   mtemp = GetCovMatrix(hist, debug);
//
//   TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, c+y+process);
//   DrawCov(htemp, save_nfs+"SYS/"+process+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
//
//   return mtemp;
// }
//
// TMatrixD GetCovMatrixNorm(TMatrixD mtemp, TH1F* hist, const TString& save, const TString& process){
//   auto start = high_resolution_clock::now(); // Calculation time - start
//
//   TString c = channel; TString y = year;
//   bool onlyDiag = process.Contains("data")?false:true;
//   TMatrixD norm = TMatrixD(NormCovMatrixAlt(hist, mtemp, false, onlyDiag));
//   TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, c+y+process+"norm");
//   DrawCov(htemp_norm, save_nfs+"SYS/"+process+"_norm_"+c+"_"+y, "#it{m}_{W}", 0.3);
//
//   auto stop = high_resolution_clock::now();  // Calculation time - stop
//   auto duration = duration_cast<seconds>(stop - start);
//   if(!onlyDiag) cout << process << " took " << RED << duration.count() << "s" << RESET << endl;
//
//   return norm;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
//
// TH1F* RebinMap(MapHH& process, int& width){
//   TString c = channel; TString y = year;
//   TH1F* hist = (TH1F*) process[c][y]->Rebin(width);
//   return hist;
// }
//
// MapH NormalizeMap(MapH& process){
//   MapH map;
//   for(auto hist: process){
//     TString h = hist.first;
//     map[h] = Normalize(((TH1F*) process[h])); // already sets normalized error
//   }
//   return map;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
//
// vector<bool> GetSkipBins(MapHHH &map, MapVI& peak, TString& hist){
//   vector<bool> skipBin;
//   VecI v = peak[hist];
//   for(unsigned int i=1; i<=map[hist][channel][year]->GetNbinsX(); i++){
//     if(find(v.begin(), v.end(), i) != v.end()) skipBin.push_back(false);
//     else skipBin.push_back(true);
//   }
//   return skipBin;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
// bool sortcolx( const VecD& v1, const VecD& v2 ) {
//   return v1[0] < v2[0];
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
//
// TMatrixD GetSYSCov(TH1F* hist, TH1F* var, TString sysname, int panel, TString add){
//   int nbins = hist->GetXaxis()->GetNbins();
//   TH1F* sys = (TH1F*) hist->Clone();
//   sys->Reset();
//   vector<double> binUnc;
//   for(int i=1; i<=nbins; i++){
//     double diff = hist->GetBinContent(i) - var->GetBinContent(i);
//     sys->SetBinContent(i, diff);
//     binUnc.push_back(diff);
//   }
//   TMatrix matrix = TMatrixD(nbins, nbins);
//   for(int i=0; i<nbins; i++){
//     for(int j=0; j<nbins; j++){
//       matrix[i][j] = binUnc[i]*binUnc[j];
//     }
//   }
//   TH2D* h_norm = TMatrixDtoTH2D(matrix, nbins, 0, 180, "cov_"+sysname+"_"+year+to_string(panel)+add);
//   DrawCov(h_norm, save_nfs+"/SYS/cov_"+sysname+"_"+year+"_"+add, to_string(panel)+add);
//
//   TCanvas *c = new TCanvas(sysname+add, "", 600, 600);
//   gPad->SetLeftMargin(0.2);
//   gPad->SetBottomMargin(0.2);
//   sys->SetTitle(sysname);
//   sys->GetXaxis()->SetTitle("m_{W}");
//   sys->GetYaxis()->SetTitle("Uncertainty");
//   sys->GetYaxis()->SetRangeUser(-100, 100);
//   sys->SetLineColor(kRed);
//   sys->SetFillColor(kRed);
//   sys->Draw("HIST ][");
//   // c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
//   // if(saveNFS) c->SaveAs(save_nfs+"/"+year+"/SYS/SYS_"+sysname+".pdf");
//   c->SaveAs(save_nfs+"/SYS/SYS_"+sysname+"_"+add+".pdf");
//
//   // c_sys->cd(panel);
//   sys->Draw("HIST ][");
//
//   delete c;
//   return matrix;
// }
//
// // ======================================================================================================
// // ===                                                                                                ===
// // ======================================================================================================
//
// TH1F* SetError(TH1F* data, TH2F* cov){
//   TH1F* hist = (TH1F*) data->Clone();
//   int nbins = hist->GetXaxis()->GetNbins();
//   for(int i=1; i<=nbins; i++){
//     double error = sqrt(cov->GetBinContent(i,i));
//     hist->SetBinError(i, error);
//   }
//   return hist;
// }
