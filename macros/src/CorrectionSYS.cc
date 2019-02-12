#include "TMath.h"
#include <iostream>
#include <iomanip>  
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include <TTree.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include <TStyle.h>
#include <TVectorD.h>
#include <TMatrix.h>
#include <TDecompLU.h>
#include <vector>
#include <string> 
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>

// compile with:
// g++ -o CorrectionSYS CorrectionSYS.cc `root-config --cflags --glibs`

using namespace std;

int main(){

  Float_t xbins[] = {0, 80, 130, 180, 250, 350, 500};
  Float_t ybins[] = {-4, -1.5, -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 4};
  int no_ptbins = 6; // rows
  int no_etabins = 12; // columns

  double Parameter[no_ptbins][no_etabins], err[no_ptbins][no_etabins];

  vector<TF1*> FitCentral;
  for(int i=0; i<no_etabins; i++){
    FitCentral.push_back( new TF1(" ","pol2",0,500) );
  }


  // read out parameter values (and fill Central Fit Function)
  ifstream correction_file;
  correction_file.open ("CorrectionFactors_new.txt");
  int eta;
  double p0, p1, p2;
  vector<double> p0_vec;
  vector<double> p1_vec;
  vector<double> p2_vec;
  while(correction_file >> eta >> p0 >> p1 >> p2){
    FitCentral[eta]->SetParameter(0, p0);
    FitCentral[eta]->SetParameter(1, p1);
    FitCentral[eta]->SetParameter(2, p2);
    p0_vec.push_back(p0);
    p1_vec.push_back(p1);
    p2_vec.push_back(p2);
  }

  // read out errors on parameters
  ifstream error_file;
  error_file.open ("CorrectionErrors_new.txt");
  double e0, e1, e2;
  vector<double> e0_vec;
  vector<double> e1_vec;
  vector<double> e2_vec;
  while(error_file >> eta >> e0 >> e1 >> e2){
    e0_vec.push_back(e0);
    e1_vec.push_back(e1);
    e2_vec.push_back(e2);
  }

  vector<double> p0up, p1up, p2up, p0down, p1down, p2down, p0none, p1none, p2none;
  for(int i=0; i<no_etabins; i++){
    p0up.push_back(p0_vec[i] + e0_vec[i]);
    p1up.push_back(p1_vec[i] + e1_vec[i]);
    p2up.push_back(p2_vec[i] + e2_vec[i]);
    p0down.push_back(p0_vec[i]);
    p1down.push_back(p1_vec[i]);
    p2down.push_back(p2_vec[i]);
    p0none.push_back(p0_vec[i]);
    p1none.push_back(p1_vec[i]);
    p2none.push_back(p2_vec[i]);  
  }

  // construct functions for all permutations of up, down, none
  vector< vector<TF1*> > VariedFits;
  vector<TF1*> f;
  for(int j=0; j<6; j++){
    f.push_back( new TF1("f1","pol2",0,500) );
  }

  for(int i=0; i<no_etabins; i++){
    VariedFits.push_back( f );
  }
  
  for(int i_eta=0; i_eta<no_etabins; i_eta++){
    VariedFits[i_eta][0]->SetParameter(0, p0up[i_eta]);
    VariedFits[i_eta][0]->SetParameter(1, p1none[i_eta]);
    VariedFits[i_eta][0]->SetParameter(2, p2none[i_eta]);

    VariedFits[i_eta][1]->SetParameter(0, p0none[i_eta]);
    VariedFits[i_eta][1]->SetParameter(1, p1up[i_eta]);
    VariedFits[i_eta][1]->SetParameter(2, p2none[i_eta]);

    VariedFits[i_eta][2]->SetParameter(0, p0none[i_eta]);
    VariedFits[i_eta][2]->SetParameter(1, p1none[i_eta]);
    VariedFits[i_eta][2]->SetParameter(2, p2up[i_eta]);

    VariedFits[i_eta][3]->SetParameter(0, p0down[i_eta]);
    VariedFits[i_eta][3]->SetParameter(1, p1none[i_eta]);
    VariedFits[i_eta][3]->SetParameter(2, p2none[i_eta]);

    VariedFits[i_eta][4]->SetParameter(0, p0none[i_eta]);
    VariedFits[i_eta][4]->SetParameter(1, p1down[i_eta]);
    VariedFits[i_eta][4]->SetParameter(2, p2none[i_eta]);

    VariedFits[i_eta][5]->SetParameter(0, p0none[i_eta]);
    VariedFits[i_eta][5]->SetParameter(1, p1none[i_eta]);
    VariedFits[i_eta][5]->SetParameter(2, p2down[i_eta]);

  }

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3); 
  FitCentral[0]->SetLineWidth(4);
  FitCentral[0]->SetLineColor(kRed-4);
  FitCentral[0]->Draw();
  FitCentral[0]->SetTitle(" ");

  for(auto &i: VariedFits[0]){
    i->SetLineWidth(1);
    i->SetLineColor(13);
    i->Draw("SAME");
  }
  FitCentral[0]->Draw("SAME");

  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_SYS/Variations.pdf");

  return 0;
}
