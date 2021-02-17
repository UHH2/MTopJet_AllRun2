#include "../include/CentralInclude.h"

using namespace std;


/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/


int main(int argc, char* argv[]){

  TString save_path = get_save_path();
  TString channel;
  TString year;

  if(argc != 3){
    cout << "[ERROR] you have to specify the channel (muon/elec) and year!" << endl;
    return 0;
  }

  if(strcmp(argv[1], "muon") == 0)         channel = "muon";
  else if(strcmp(argv[1], "elec") == 0)    channel = "elec";
  else if(strcmp(argv[1], "combine") == 0) channel = "combine";
  else{
    cout << "[ERROR] Channel not known, select muon, elec or combine!" << endl;
    return 0;
  }

  year = argv[2];

  cout << "CHANNEL: " << argv[1] << endl;
  cout << "YEAR   : " << argv[2] << endl;


  // true masses
  // vector<double> truth = {169.5, 171.5, 172.5, 173.5, 175.5};
  vector<double> truth = {169.5, 171.5, 173.5, 175.5};


  // read txt files and get measured masses with uncert.
  vector<double> measured_stat;
  vector<double> error_stat;
  vector<double> measured;
  vector<double> error;
  TString directory = save_path+"/Plots/Unfolding_Run2/"+year+"/";
  // vector<string> massdir = {"Pseudo1695/", "Pseudo1715/", "Pseudo1/", "Pseudo1735/", "Pseudo1755/"};
  vector<TString> massdir = {"Pseudo1695/", "Pseudo1715/", "Pseudo1735/", "Pseudo1755/"};
  int i=0;
  for(auto mdir: massdir){
    std::ifstream infile(directory+mdir+channel+"/Mass.txt");
    double m, e, m2, e2;
    while (infile >> m >> e >> m2 >> e2){
      cout << "------------------------------------------" << endl;
      cout << "truth:            " << truth[i] << endl;
      cout << "stat only:        " << m << " +- " << e;
      if(fabs(truth[i]-m) > e) cout << "  <-- NOT COMPATIBLE";
      cout << endl;
      cout << "with mass uncert: " << m2 << " +- " << e2;
      if(fabs(truth[i]-m2) > e2) cout << "  <-- NOT COMPATIBLE";
      cout << endl;
      measured_stat.push_back(m);
      error_stat.push_back(e);
      measured.push_back(m2);
      error.push_back(e2);
    }
    i++;
  }
  cout << "------------------------------------------" << endl;

  TGraphErrors* masses_stat = new TGraphErrors(5, &truth[0], &measured_stat[0], 0, &error_stat[0]);
  TGraphErrors* masses = new TGraphErrors(5, &truth[0], &measured[0], 0, &error[0]);

  // diagonal line where measured = generated
  double xmin = 167;
  double xmax = 178;
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
  masses_stat->GetXaxis()->SetTitle("True #it{m}_{t} [GeV]");
  masses_stat->GetYaxis()->SetTitle("Measured #it{m}_{t} [GeV]");
  masses_stat->GetYaxis()->SetTitleOffset(1.1);
  masses_stat->GetXaxis()->SetTitleOffset(0.9);
  masses_stat->GetYaxis()->SetTitleSize(0.05);
  masses_stat->GetXaxis()->SetTitleSize(0.05);
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
  CMSSimLabel(true, 0.2, 0.85);

  TLegend *leg = new TLegend(0.5,0.15,0.88,0.38);
  leg->AddEntry(masses_stat, "Extracted #it{m}_{t} (stat only)", "pl");
  leg->AddEntry(diag_line, "Perfect measurement", "l");
  // leg->AddEntry(fit3, "calibration fit", "l");
  // leg->AddEntry(Area, "fit uncertainty", "f");
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs(save_path+"/Plots/Unfolding_Run2/MassCalibratrion_"+year+"_"+channel+"_stat.pdf");

  TCanvas *B = new TCanvas("B", "B", 600, 600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  masses->SetTitle(" ");
  masses->GetXaxis()->SetLimits(xmin, xmax);
  masses->GetHistogram()->SetMinimum(ymin);
  masses->GetHistogram()->SetMaximum(ymax);
  masses->GetXaxis()->SetTitle("True #it{m}_{t} [GeV]");
  masses->GetYaxis()->SetTitle("Measured #it{m}_{t} [GeV]");
  masses->GetYaxis()->SetTitleOffset(1.1);
  masses->GetXaxis()->SetTitleOffset(0.9);
  masses->GetYaxis()->SetTitleSize(0.05);
  masses->GetXaxis()->SetTitleSize(0.05);  masses->GetXaxis()->SetNdivisions(505);
  masses->GetYaxis()->SetNdivisions(505);
  masses->SetMarkerStyle(20);
  masses->SetMarkerSize(1.3);
  masses->SetLineColor(1);
  masses->Draw("AP");
  // Area->SetFillColor(16);
  // Area->Draw("f SAME");
  // fit3->SetLineColor(kAzure+7);
  // fit3->SetLineWidth(3);
  // fit3->SetLineStyle(1);
  diag_line->Draw("SAME");
  // fit3->Draw("SAME");
  masses->Draw("P SAME");
  // up->Draw("l SAME");
  // down->Draw("l SAME");
  CMSSimLabel(true, 0.2, 0.85);

  TLegend *leg2 = new TLegend(0.5,0.15,0.88,0.38);
  leg2->AddEntry(masses, "Extracted #it{m}_{t}^{MC}", "pl");
  leg2->AddEntry(diag_line, "Perfect measurement", "l");
  // leg->AddEntry(fit3, "calibration fit", "l");
  // leg->AddEntry(Area, "fit uncertainty", "f");
  leg2->Draw();
  gPad->RedrawAxis();
  B->SaveAs(save_path+"/Plots/Unfolding_Run2/MassCalibratrion_"+year+"_"+channel+".pdf");

}
