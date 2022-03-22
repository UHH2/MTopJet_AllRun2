#include "../include/CentralInclude.h"
#include "../include/tdrstyle_all.h"
#include "../include/HistogramUtils.h"

// #include <ranges>

using namespace std;

Double_t lumi_plot;
TString hist;
TString year, year_v;
TString fdir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/";
TString uhh2MC = "uhh2.AnalysisModuleRunner.MC.";
bool debug = false;

VecTS process = {"other", "WJets", "SingleTop", "TTbar"};

MapVecTS hists = {
  {"XCone_cor",                        {"M_jet1_"}        },
  {"XCone_cor_SF",                     {"M_jet1_"}        },
  {"XCone_cor_subjets_SF",             {"pt_had_subjets"} },
  {"comparison_topjet_xcone_pass_rec", {"wmass_match"}    },
  {"JetMassScaleHists",                {"hadjet_jms_mass"}},
  {"PreSel04_Event",                   {"MET"}            },
  {"PreSel04_Muon",                    {"pt"}             },
  {"MuonHists",                        {"pt"}             },
  {"PreSel04_Elec",                    {"pt"}             },
  {"ElecHists",                        {"pt"}             },
};

MapTS xaxis = {
  {"M_jet1_", "#it{m}_{jet}"},
  {"wmass_match", "#it{m}_{W}"},
  {"hadjet_jms_mass", "#it{m}_{jet}"},
  {"pt", "#it{p}_{T}"},
  {"MET", "missing #it{p}_{T}"},
};

MapI colors = {
  {"TTbar", kRed+1},
  {"SingleTop", kBlue+2},
  {"WJets", kGreen+1},
  {"other", kOrange}
};

MapFF files;
MapHHH hists_data, hists_mc;
MapHHH ratios_data, ratios_mc;
Map3Stack stacks;

void DrawText(TString name);
void DrawControl(TH1F* data, THStack* stack, TH1F* rData, TH1F* rMC, TString xaxis, TString h_class, TString h);
TH1F* GetDataHistograms(TFile* data, TString h_class, TString h);
TH1F* GetMCHistograms(TFile* file, TString h_class, TString h);
TH1F* GetRatioStack(TH1F* h1, TH1F* h2);
THStack* StackHistograms(MapF files, TString h_class, TString h);

int main(int argc, char* argv[]){
  cout << "===================================================" << endl << endl;

  gErrorIgnoreLevel = kWarning;

  year = argv[1];
  year_v = (year.EqualTo("2016"))?"_2016v3":(year.EqualTo("2017"))?"_2017v2":(year.EqualTo("2018"))?"_2018":"WRONG_YEAR";
  lumi_plot = (year.EqualTo("2016"))?35.9:(year.EqualTo("2017"))?41.5:(year.EqualTo("2018"))?59.74:00.00;

  // ===========================================================================
  //                        Define Hists
  // ===========================================================================

  // ===========================================================================
  // === Get Files
  cout << "Collect Files ... " << endl;

  for(auto c: {"muon", "elec"}){
    files[c];
    for(auto p: process){
      files[c][p] = new TFile(fdir+c+"/"+uhh2MC+p+year_v+".root");
    }
  }

  files["muon"]["Data"] = new TFile(fdir+"muon/uhh2.AnalysisModuleRunner.DATA.DATA"+year_v+".root");
  files["elec"]["Data"] = new TFile(fdir+"elec/uhh2.AnalysisModuleRunner.DATA.DATA"+year_v+".root");

  // ===========================================================================
  // === Get hists
  cout << "Get Histograms ... " << endl;

  for(TString c: {"muon", "elec"}){
    hists_data[c]; hists_mc[c];
    for (auto const& [h_class, vec_h] : hists)
    {
      hists_data[c][h_class]; hists_mc[c][h_class];
      for(auto const& h: vec_h){
        hists_data[c][h_class][h] = GetDataHistograms(files[c]["Data"], h_class, h);
        stacks[c][h_class][h] = StackHistograms(files[c], h_class, h);

        // === For ratio plot - Find different solution! ============
        TH1F* h_full = (TH1F*) files[c]["TTbar"]->Get(h_class+"/"+h);
        for(TString p: process){
          if(p.EqualTo("TTbar")) continue;
          h_full->Add((TH1F*) files[c][p]->Get(h_class+"/"+h), 1);
        }
        hists_mc[c][h_class][h] = h_full;
        // ==========================================================

      }
    }
  }

  // ===========================================================================
  // === Get Ratios
  cout << "Get Ratios ... " << endl;

  for(TString c: {"muon", "elec"}){
    ratios_data[c]; ratios_data[c];
    for (auto const& [h_class, vec_h] : hists)
    {
      ratios_data[c][h_class]; ratios_data[c][h_class];
      for(auto const& h: vec_h){
        ratios_mc[c][h_class][h] = GetRatioStack(hists_mc[c][h_class][h], hists_mc[c][h_class][h]);
        ratios_data[c][h_class][h] = GetRatioStack(hists_data[c][h_class][h], hists_mc[c][h_class][h]);
      }
    }
  }

  // ===========================================================================
  //                        Draw Plots
  // ===========================================================================

  SetupGlobalStyle();

  for(TString c: {"muon", "elec"}){
    for (auto const& [h_class, vec_h] : hists)
    {
      for(auto const& h: vec_h){
        DrawControl(hists_data[c][h_class][h], stacks[c][h_class][h], ratios_data[c][h_class][h], ratios_mc[c][h_class][h], xaxis[h], h_class, h);
      }
    }
  }

  return 0;
}

// ==================================================================================================
TH1F* GetDataHistograms(TFile* data, TString h_class, TString h){
  TH1F* hist = (TH1F*) data->Get(h_class+"/"+h);
  hist->SetTitle("");
  hist->SetMarkerStyle(8);
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);

  hist->GetXaxis()->SetTitleSize(0); // x axis in ratio
  hist->GetXaxis()->SetLabelSize(0);

  hist->GetYaxis()->SetTitle("Events");
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.15);

  return hist;
}

// ==================================================================================================
TH1F* GetMCHistograms(TFile* file, TString h_class, TString h){
  TH1F* hist = (TH1F*) file->Get(h_class+"/"+h);
  return hist;
}

// ==================================================================================================
THStack* StackHistograms(MapF files, TString h_class, TString h){
  THStack* stack = new THStack(h_class+"/"+h, "");
  for(TString p: process){
    // files contains data, but data not in process
    TH1F* hist = (TH1F*) files[p]->Get(h_class+"/"+h);
    hist->SetFillColor(colors[p]);
    hist->SetLineColor(colors[p]);
    hist->SetLineWidth(0);
    stack->Add(hist);
  }
  return stack;
}

// ==================================================================================================
void DrawControl(TH1F* data, THStack* stack, TH1F* rData, TH1F* rMC, TString xaxis, TString h_class, TString h)
{
  if(debug) cout << "Draw Control Hists ("+h+") ..." << endl;
  gStyle->SetOptStat(0); // turn off stats box

  TString unique_path = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/PostSel/"+year+"/"+h_class+"_"+h+".pdf";
  double bins = data->GetNbinsX();

  TLegend *leg = new TLegend(0.7,0.6,0.85,0.85);
  leg->SetBorderSize(0);
  TCanvas *A = new TCanvas(unique_path, unique_path, 600, 600);

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  TH1F* errorUp = (TH1F*) rMC->Clone(); TH1F* errorDown = (TH1F*) rMC->Clone();
  for(unsigned int i=1; i<=rMC->GetNbinsX(); i++){
    errorUp->SetBinContent(i, rMC->GetBinContent(i)+rMC->GetBinError(i));
    errorDown->SetBinContent(i, rMC->GetBinContent(i)-rMC->GetBinError(i));
  }
  errorUp->SetFillColorAlpha(kGray+2, 0.5);
  errorUp->SetLineWidth(0);
  errorDown->SetFillColorAlpha(kWhite, 1);
  errorDown->SetLineWidth(0);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();

  data->Draw("P");
  data->GetYaxis()->SetRangeUser(0.1, data->GetMaximum()*1.5);
  stack->Draw("same hist");
  data->Draw("same PE"); // redraw
  gPad->RedrawAxis();

  leg->AddEntry(data, "Data", "pl");

  // dummy hists for legend
  MapH dummies;
  vector<TString> flip = process;
  reverse(flip.begin(), flip.end());
  for (auto p: flip) {
    TH1F* dummy = new TH1F();
    dummy->SetFillColor(colors[p]);
    dummy->SetLineColor(colors[p]);
    leg->AddEntry(dummy, p, "f");
  }
  leg->SetTextSize(0.04);
  leg->Draw();

  if(debug) cout << "\t ... Plot info (CMS+Lumi)" << endl;
  CMSLabel(0.22, 0.85, "Work in Progress");
  LumiInfo(lumi_plot, false, 0.9, 0.96);

  // Ratio Plot kommt in unteres pad
  A->cd();
  if(debug) cout << "\t ... Pad2" << endl;
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  A->SetTickx(0);
  A->SetTicky(0);

  TLine* lMC = new TLine(-2.5, 1, 2.5, 1);
  lMC->SetLineWidth(1);
  lMC->SetLineColor(kBlack);

  rData->SetMarkerStyle(8);
  rData->SetLineColor(kBlack);
  rMC->SetFillStyle(0);
  rMC->SetLineStyle(kSolid);
  rMC->SetLineWidth(0);
  rMC->GetXaxis()->SetLabelSize(0.15);
  rMC->GetYaxis()->SetLabelSize(0.13);
  rMC->GetXaxis()->SetTitleSize(0.1);
  rMC->GetYaxis()->SetRangeUser(0.4, 1.6);
  rMC->GetXaxis()->SetTitle(xaxis);
  rMC->GetYaxis()->SetTitle("#frac{DATA}{MC}");
  rMC->GetXaxis()->SetTitleOffset(0.8);
  rMC->GetXaxis()->SetTitleSize(0.2);
  rMC->GetYaxis()->SetTitleSize(0.15);
  rMC->GetYaxis()->SetNdivisions(505);
  rMC->GetYaxis()->SetTitleOffset(0.4);
  rMC->GetYaxis()->CenterTitle();
  rMC->Draw("Axis");
  errorUp->Draw("Hist same");
  errorDown->Draw("Hist same");
  lMC->Draw("same");
  rData->Draw("SAME PE");

  gPad->RedrawAxis();
  A->SaveAs(unique_path);

  delete A;
  leg->Clear();
}

// ==================================================================================================
void SetupCanvasForEPS()
{
  // set up a canvas for single EPS files
  // optimised plots for including in theses or publications and documents
  // different possibilities
  // ratio/no ratio plots

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 400;
  CanHeight = 400;

  // set up the canvas
  TCanvas* m_can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //
  Float_t y1, y2, y3;                           //  y3 +-------------+
  y3 = 0.99;                                    //     |             |
  y2 = y3-yplot;                                //     |     pad1    |
  y1 = y2-yratio;                               //  y2 |-------------|
  Float_t x1, x2;                               //     |     rp1     |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
  //
  // No Pad 2!


  TPad* m_rp1_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  TPad* m_rp1 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  TPad* m_pad1 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);

  TPad* m_rp2_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  TPad* m_rp2 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  TPad* m_pad2 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);


  m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
  m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);

  m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.0);  m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
  m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.0);  m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);
  m_rp1->SetTopMargin(0.0);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19);  m_rp1->SetRightMargin(0.05);
  m_rp2->SetTopMargin(0.0);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19);  m_rp2->SetRightMargin(0.05);

  if (debug){
    m_rp1_top->SetFillColor(kYellow);
    m_rp2_top->SetFillColor(kOrange);
    // if (bPlotRatio){
    //   m_rp1->SetFillColor(kGray);
    //   m_rp2->SetFillColor(kGray);
    // }
  }

  m_pad1->Draw();
  m_pad2->Draw();

  m_rp1_top->Draw();
  m_rp2_top->Draw();

  // m_rp1->Draw();
  // m_rp2->Draw();

  return;

}

// ==================================================================================================
TH1F* GetRatioStack(TH1F* h1, TH1F* h2){
  CompareHistStructure(h1, h2);
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetNbinsX();
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 0);
      ratio->SetBinError(i, 10);
    }
    else{
      double r = N1/N2;
      double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
