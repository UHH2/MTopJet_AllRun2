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

double xmin = 0;
double xmax = 0;

vector<double> CalculateParameterUncertainty(TGraphErrors* data_, TString function, TFitResultPtr fitresult, vector<double> parameters, TString updown);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/


int main(int argc, char* argv[]){

  // true masses
  vector<double> truth = {169.5, 171.5, 172.5, 173.5, 175.5};


  // read txt files and get measured masses with uncert.
  vector<double> measured_stat;
  vector<double> error_stat;
  string directory = "/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/";
  vector<string> massdir = {"pseudo1695", "pseudo1715", "pseudo1", "pseudo1735", "pseudo1755"};
  for(auto mdir: massdir){
    std::ifstream infile(directory+mdir+"/Mass.txt");
    double m, e;
    while (infile >> m >> e){
      cout << "stat only:        " << m << " +- " << e << endl;
      measured_stat.push_back(m);
      error_stat.push_back(e);

    }
  }

  TGraphErrors* masses_stat = new TGraphErrors(5, &truth[0], &measured_stat[0], 0, &error_stat[0]);

  // diagonal line where measured = generated
  xmin = 167;
  xmax = 178;
  double ymin = xmin;
  double ymax = xmax;

  TLine *diag_line = new TLine(xmin, ymin, xmax, ymax);
  diag_line->SetLineColor(kRed);
  diag_line->SetLineWidth(3);
  diag_line->SetLineStyle(7);

  /*
  ██████  ██       ██████  ████████
  ██   ██ ██      ██    ██    ██
  ██████  ██      ██    ██    ██
  ██      ██      ██    ██    ██
  ██      ███████  ██████     ██
  */


  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  masses_stat->SetTitle(" ");
  masses_stat->GetXaxis()->SetLimits(xmin, xmax);
  masses_stat->GetHistogram()->SetMinimum(ymin);
  masses_stat->GetHistogram()->SetMaximum(ymax);
  masses_stat->GetXaxis()->SetTitle("m_{top} truth [GeV]");
  masses_stat->GetYaxis()->SetTitle("m_{top} measured [GeV]");
  masses_stat->GetYaxis()->SetTitleOffset(1.5);
  masses_stat->GetXaxis()->SetNdivisions(505);
  masses_stat->GetYaxis()->SetNdivisions(505);
  masses_stat->SetMarkerStyle(20);
  masses_stat->SetMarkerSize(1.3);
  masses_stat->SetLineColor(1);
  masses_stat->Draw("AP");
  // Area->SetFillColor(16);
  // Area->Draw("f SAME");
  // fit3->SetLineColor(kAzure+7);
  // fit3->SetLineWidth(3);
  // fit3->SetLineStyle(1);
  diag_line->Draw("SAME");
  // fit3->Draw("SAME");
  masses_stat->Draw("P SAME");
  // up->Draw("l SAME");
  // down->Draw("l SAME");
  TLegend *leg = new TLegend(0.6,0.15,0.88,0.38);
  leg->AddEntry(masses_stat, "extracted m_{top}^{MC} (stat only)", "l");
  leg->AddEntry(diag_line, "perfect measurement", "l");
  // leg->AddEntry(fit3, "calibration fit", "l");
  // leg->AddEntry(Area, "fit uncertainty", "f");
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/MassCalibratrion_stat.pdf");

}
