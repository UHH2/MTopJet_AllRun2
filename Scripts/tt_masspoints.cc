#include "TMath.h"
#include <iostream>
#include <iomanip>  
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include <TTree.h>
#include <TH1.h>
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include <TStyle.h>

using namespace std;

TH2D* calculate_cov(TH1F* hist, int lower, int upper, bool binwidth);
double ComputeChi2(TH1F* data, TH1F* MC, TH2D* cov, std::vector<int> bins);
void show_matrix(TH2D* matrix, int lower, int upper);

void tt_masspoints()
{
  bool use_data = false;
  bool binwidth = false;

  TString data_file_string;
  if(use_data) data_file_string = "uhh2.AnalysisModuleRunner.DATA.DATA.root";
  else         data_file_string = "uhh2.AnalysisModuleRunner.MC.TTbar.root";

  TFile *data = new TFile(data_file_string);
  TFile *bgr = new TFile("uhh2.AnalysisModuleRunner.MC.Background_only.root");

  TFile *TT = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *TT_1665 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1665.root");
  TFile *TT_1695 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695_ext2.root");
  TFile *TT_1715 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1715.root");
  TFile *TT_1735 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1735.root");
  TFile *TT_1755 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
  TFile *TT_1785 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_mtop1785.root");


  TH1F * JetMass_data = (TH1F*)data->Get("XCone_cor/M_jet1");

  if(use_data){
    TH1F * JetMass_bgr = (TH1F*)bgr->Get("XCone_cor/M_jet1");
    JetMass_data->Add(JetMass_bgr, -1); // subtract background
  }

  TH1F * JetMass_1665 = (TH1F*)TT_1665->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1695 = (TH1F*)TT_1695->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1715 = (TH1F*)TT_1715->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1725 = (TH1F*)TT->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1735 = (TH1F*)TT_1735->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1755 = (TH1F*)TT_1755->Get("XCone_cor/M_jet1");
  TH1F * JetMass_1785 = (TH1F*)TT_1785->Get("XCone_cor/M_jet1");

  std::vector<TH1F*> mass_samples;
  if(use_data) mass_samples = {JetMass_1665, JetMass_1695, JetMass_1715, JetMass_1725, JetMass_1735, JetMass_1755, JetMass_1785};
  else         mass_samples = {JetMass_1665, JetMass_1695, JetMass_1715, JetMass_1735, JetMass_1755, JetMass_1785};

  double from = 140;
  double to = 260;
  int lower = JetMass_data->GetXaxis()->FindBin(from);
  int upper = JetMass_data->GetXaxis()->FindBin(to);

  // calculate covariance matrix with normalisation applied
  // careful: do not use scaled distribution here!
  TH2D*cov = calculate_cov(JetMass_data, lower, upper, binwidth);



  // select bins that enter chi2 calculation
  // leave one bin out (chi2 should not change)
  int j = 0;
  std::vector< std::vector<int> > bins;

  // normalise here
  for(auto &hist: mass_samples){
    if(binwidth) hist->Scale(1/hist->Integral(lower, upper), "width");
    else         hist->Scale(1/hist->Integral(lower, upper));

  }
  if(binwidth) JetMass_data->Scale(1/JetMass_data->Integral(lower, upper),"width");
  else         JetMass_data->Scale(1/JetMass_data->Integral(lower, upper));
  //

  bool skip;
  int rand_bin = rand()%(upper-lower+1)+lower; // select random bin that is not used
  for(auto &masspoint: mass_samples){
    std::vector<int> bins_this;
    for(int i=lower; i<=upper; i++){
      if(i == rand_bin) skip = true;
      else skip = false;
      if(!skip) bins_this.push_back(i);
    }
    bins.push_back(bins_this);
    j++;
  }

  std::vector<double> chi_vec;
  j = 0;
  for(auto &masspoint: mass_samples){
    double chi2 = ComputeChi2(JetMass_data, masspoint, cov, bins[j]);
    cout << chi2 << endl;
    chi_vec.push_back(chi2);
    j++;
  }

  const int N_samples = mass_samples.size();
  Double_t x[7] = {166.5, 169.5, 171.5, 173.5, 175.5, 178.5};
  Double_t y[7];
  for(int i=0; i<N_samples; i++){
    y[i] = chi_vec[i];
  }
 
  TGraph* chi_hist = new TGraph(N_samples,x,y);
  TF1 * fit;
  TF1*f1 = new TF1("f1","pol2",0,500);
  chi_hist->Fit("f1","R");
  fit = chi_hist->GetFunction("f1");
  double minY = fit->GetMinimum();
  double minX = fit->GetX(minY, 160, 180);
  double minX_sigmaup = fit->GetX(minY+1, minX, 180);
  double minX_sigmadown = fit->GetX(minY+1, 160, minX);
  cout << endl << "=========================== " << endl;
  cout << "| value = " << std::setprecision(5) << minX << " +- " << std::setprecision(4) << (minX_sigmaup-minX) << "|" << endl; 
  cout << "=========================== " << endl << endl;

  double xlow, xhigh;
  double ylow, yhigh;
  xlow = from;
  xhigh = to;
  // xlow = 0;
  // xhigh = 1000;
  ylow = 0;
  yhigh = JetMass_data->GetMaximum() * 1.2;
  
  for(auto &i: mass_samples){
    i->SetLineWidth(3);
    i->SetTitle(" ");
    i->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
    i->GetYaxis()->SetTitle("a.u.");
    i->GetYaxis()->SetRangeUser(ylow, yhigh);
    i->GetXaxis()->SetRangeUser(xlow, xhigh);
    i->GetXaxis()->SetNdivisions(505);
    i->GetYaxis()->SetNdivisions(505);
    i->GetYaxis()->SetTitleOffset(1.5);
    i->SetLineWidth(4);
  }

  JetMass_1715->SetLineColor(kRed);
  JetMass_1735->SetLineColor(kAzure+7);

  JetMass_1665->SetLineColor(kRed);
  JetMass_1725->SetLineColor(13);
  JetMass_1785->SetLineColor(kAzure+7);

  JetMass_data->SetLineColor(kBlack);
  JetMass_data->SetMarkerColor(kBlack);
  JetMass_data->SetMarkerStyle(8);
  JetMass_data->SetMarkerSize(1); 
  JetMass_data->GetYaxis()->SetTitleOffset(1.5);
  JetMass_data->GetXaxis()->SetNdivisions(505);
  JetMass_data->GetYaxis()->SetNdivisions(505);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3); 
  if(use_data){ 
    for(auto i: {JetMass_1665, JetMass_1725, JetMass_1785}) i->Draw("HIST SAME");
  }
  else{
    for(auto i: {JetMass_1715, JetMass_1735})i->Draw("HIST SAME");
  }

  JetMass_data->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.55,0.55,0.86,0.8);
  if(use_data) leg->AddEntry(JetMass_data,"data","pl");
  else         leg->AddEntry(JetMass_data,"pseudo data","pl");
  if(use_data){
    leg->AddEntry(JetMass_1665,"t#bar{t} (m_{top} = 166.5 GeV)","l");
    leg->AddEntry(JetMass_1725,"t#bar{t} (m_{top} = 172.5 GeV)","l");
    leg->AddEntry(JetMass_1785,"t#bar{t} (m_{top} = 178.5 GeV)","l");
  }
  else{
    leg->AddEntry(JetMass_1715,"t#bar{t} (m_{top} = 171.5 GeV)","l");
    leg->AddEntry(JetMass_1735,"t#bar{t} (m_{top} = 173.5 GeV)","l");
  }
  leg->Draw("");
  if(use_data) A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_norm.pdf");
  else         A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_norm_pseudo.pdf");

  TCanvas *B = new TCanvas("B", "B", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3); 
  chi_hist->SetTitle(" ");
  chi_hist->GetXaxis()->SetTitle("m_{top} [GeV]");
  chi_hist->GetYaxis()->SetTitle("#chi^{2}");
  chi_hist->GetYaxis()->SetTitleOffset(1.5);
  chi_hist->GetXaxis()->SetNdivisions(505);
  chi_hist->GetYaxis()->SetNdivisions(505);
  chi_hist->SetMarkerStyle(20);
  chi_hist->SetMarkerSize(1.5);
  chi_hist->SetLineColor(1);
  chi_hist->Draw("AP");
  if(use_data) B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_chi2.pdf");
  else         B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_chi2_pseudo.pdf");


}





double ComputeChi2(TH1F* data, TH1F* MC, TH2D* cov, std::vector<int> bins){
  //double chi2 = 0;
  TH1F* data_sc = (TH1F*)data->Clone("data_sc");
  TH2D* cov_sc = (TH2D*)cov->Clone("cov_sc");
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
    cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << endl;
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

  return chi2;
}


TH2D* calculate_cov(TH1F* hist, int lower, int upper, bool binwidth){
  int nbins = upper;
  TH2D* old_cov = new TH2D("old_cov", "old_cov", nbins, 1, upper+1, nbins, 1, upper+1);
  TH2D* new_cov = new TH2D("new_cov", "new_cov", nbins, 1, upper+1, nbins, 1, upper+1); 
  for(int i=1; i <= upper; i++){
    double error = hist->GetBinError(i);
    if(i < lower ) old_cov->Fill(i,i,0); 
    else  old_cov->Fill(i,i,pow(error,2));
  }
  double integral = hist->Integral(lower, upper);
  for(int i=1; i <= upper; i++){
    for(int j=1; j <= upper; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= upper; k++){
	for(int l=1; l <= upper; l++){
	  old_entry = old_cov->GetBinContent(k,l);
	  double binwidth_i = hist->GetBinWidth(i);
	  double binwidth_j = hist->GetBinWidth(j);
	  if(!binwidth){
	    binwidth_i = 1.0;
	    binwidth_j = 1.0;
	  }
	  if(i==k) derivation_i = (integral - hist->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
	  else     derivation_i = - (hist->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
	  if(j==l) derivation_j = (integral - hist->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
	  else     derivation_j = - (hist->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
	  sum += derivation_i * derivation_j * old_entry;
	}
      }
      new_cov->Fill(i, j, sum);
    }
  }

  // show_matrix(old_cov, 1, upper);
  // show_matrix(new_cov, 1, upper);

  return new_cov;
}


void show_matrix(TH2D* matrix, int lower, int upper){

  for(int i=lower; i <= upper; i++){
    printf("%3d:  ",i);
    for(int j=lower; j <= upper; j++){
      printf("% 5.2f ",  matrix->GetBinContent(j,i));
    }
    printf("\n");
  }
  cout << "================================================" << endl;

  return;
}
