#pragma once
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <vector>
#include <iostream>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <TColor.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>


using namespace std;

TH2D* GetCovFromStat(TH1F* hist, int firstbin, int lastbin, bool binwidth);
TH2D* NormaliseCov(TH2D* old_cov, TH1F* central, int firstbin, int lastbin, bool binwidth);
TH2D* create_cov(TH1F* errorsUP, TH1F* errorsDOWN, int firstbin, int lastbin);
double ComputeChi2(TH1F* data, TH1F* MC, TH2D* cov_data, TH2D* cov_mc, std::vector<int> bins);
void show_matrix(TH2D* matrix, int firstbin, int lastbin);
TH1F* DataWithCorrectErrors(TH1F* JetMass_data, TH2D* cov, int firstbin, int lastbin);
TH1F* GetNormErrorFromSYS(TH1F* Central, TH1F* Variation, int firstbin, int lastbin);
TH2D* GetCovFromErrors(TH1F* error, int firstbin, int lastbin);
TH1F* GetSymmetricError(TH1F* up, TH1F* down);
