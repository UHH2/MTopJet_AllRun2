#include "../include/CentralInclude.h"

// compile with:
// g++ -o SubjetFlavor SubjetFlavor.cc `root-config --cflags --glibs`

using namespace std;


int main(int argc, char* argv[]){

  TString mode;
  if(argc >1 && strcmp(argv[1], "A") == 0) mode = "A";
  else if(argc >1 && strcmp(argv[1], "B") == 0) mode = "B";
  else if(argc >1 && strcmp(argv[1], "C") == 0) mode = "C";
  else mode = "A";
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get lumi plots -------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TString dir = "GenParticles_GenOnly/";
  vector<TString> names = {"fraction_ud_", "fraction_s_", "fraction_c_", "fraction_b_", "fraction_g_"};
  for(unsigned int i=0; i<names.size(); i++){
    names[i] += mode;
  }
  names.push_back("fraction_nomatch_A"); // there is only version A but B and C would be the same
  vector<TString> flavors = {"ud", "s", "c", "b", "g", "no match"};
  vector<TH1F*> hists;

  bool skip_unmatched = true;

  TH1F *nMatches = (TH1F*)TT_file->Get(dir+"N_matches");

  for(auto name: names){
    if(skip_unmatched && name == "fraction_nomatch_A") continue;
    TString histname = dir + name;
    TH1F *h = (TH1F*)TT_file->Get(histname);
    h->Rebin(2);
    hists.push_back(h);
  }

  const int nbins = hists[0]->GetSize()-2;   // also do not consider highest bins


  // first calculate sum in every bin
  vector<double> sum;
  for(unsigned int i=1; i<=nbins; i++){
    double s = 0;
    for(unsigned int k=0; k<hists.size(); k++){
      s += hists[k]->GetBinContent(i);
    }
    sum.push_back(s);
  }

  // now read x and y values
  vector<vector<double>> xvals;
  vector<vector<double>> yvals;
  for(unsigned int j=0; j<hists.size(); j++){
    vector<double> x,y;
    for(unsigned int i=1; i<=nbins; i++){
      x.push_back(hists[j]->GetXaxis()->GetBinCenter(i));
      y.push_back(hists[j]->GetBinContent(i)/sum[i-1]); // careful: sum[0] is for bin 1
    }
    xvals.push_back(x);
    yvals.push_back(y);
  }

  // Now write into TGraph for each flavor
  vector<TGraph*> fractions;
  for(unsigned int i=0; i<hists.size(); i++){
    TGraph* g = new TGraph(xvals[i].size(), &xvals[i][0], &yvals[i][0]);
    fractions.push_back(g);
  }

  // write output to file
  TString rootname = "files/FlavorFractions_"+mode+".root";
  TFile * f_fractions = new TFile(rootname,"RECREATE");;
  for(unsigned int i=0; i<fractions.size(); i++){
    TString name = flavors[i] + "_fraction";
    if(i == 5) name = "no_match";
    fractions[i]->Write(name);
  }
  f_fractions->Close();


  // Set Stile and Plot
  vector<Color_t> colors = {kRed+1, kBlue+1, kGreen+1, 798, kAzure+7, 13};
  vector<Style_t> markers = {21, 22, 23, 33, 29, 20};

  for(unsigned int i=0; i<fractions.size(); i++){
    fractions[i]->SetTitle(" ");
    fractions[i]->GetXaxis()->SetTitle("subjet p_{T}");
    fractions[i]->GetYaxis()->SetTitle("flavor fraction");
    fractions[i]->GetYaxis()->SetTitleSize(0.06);
    fractions[i]->GetXaxis()->SetTitleSize(0.05);
    fractions[i]->GetXaxis()->SetTitleOffset(0.9);
    fractions[i]->GetYaxis()->SetTitleOffset(1.1);
    fractions[i]->GetXaxis()->SetNdivisions(505);
    fractions[i]->GetYaxis()->SetNdivisions(505);
    fractions[i]->SetLineWidth(4);
    fractions[i]->SetLineColor(colors[i]);
    fractions[i]->SetMarkerColor(colors[i]);
    fractions[i]->SetMarkerStyle(markers[i]);
  }

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  TCanvas *c = new TCanvas("c", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  // gPad->SetLogx();
  // gPad->SetLogy();
  fractions[0]->GetXaxis()->SetRangeUser(0, 600);
  fractions[0]->GetYaxis()->SetRangeUser(0.0, 0.7);
  fractions[0]->Draw("AP");
  TLegend* leg = new TLegend(0.65,0.65,0.85,0.85);
  leg->SetNColumns(2);
  for(unsigned int i=0; i<fractions.size(); i++){
    fractions[i]->Draw("P SAME");
    leg->AddEntry(fractions[i], flavors[i], "p");
  }
  leg->Draw();
  c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/FlavorFractions/Flavor_fraction_"+mode+".pdf");

  TCanvas *d = new TCanvas("d", " ", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  nMatches->SetTitle(" ");
  nMatches->GetXaxis()->SetTitle("number of particles matched to subjet");
  nMatches->GetYaxis()->SetTitle("a.u.");
  nMatches->GetYaxis()->SetTitleSize(0.06);
  nMatches->GetXaxis()->SetTitleSize(0.05);
  nMatches->GetZaxis()->SetTitleSize(0.05);
  nMatches->GetXaxis()->SetTitleOffset(0.9);
  nMatches->GetYaxis()->SetTitleOffset(1.1);
  nMatches->GetZaxis()->SetTitleOffset(0.9);
  nMatches->GetXaxis()->SetNdivisions(505);
  nMatches->GetYaxis()->SetNdivisions(505);
  nMatches->GetXaxis()->SetRangeUser(0, 4);
  nMatches->Scale(1/nMatches->Integral());
  nMatches->SetFillColor(13);
  nMatches->SetLineColor(1);
  nMatches->Draw("HIST");
  d->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/FlavorFractions/Flavor_nMatches.pdf");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  return 0;
}
