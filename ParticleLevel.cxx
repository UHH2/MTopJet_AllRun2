#include <iostream>
#include <vector>

#include "PlottingStyle.h"

using namespace std;

// void SetHist(TH1F *hist);
// void SetupGlobalStyle();

int ParticleLevel(){



  //==================================
  //important settings:
  //linestyles, markerstyles, colors, historgramnames, legendentries, axeslegends
  //==================================




  vector<TString> FileNames;
  FileNames.push_back("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  FileNames.push_back("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  FileNames.push_back("/nfs/dust/cms/user/paaschal/MTopJet/PostSel/2016/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
  FileNames.push_back("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2016v3.root");
  FileNames.push_back("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2017v2.root");
  FileNames.push_back("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_2018.root");


  TString histname = "Mass_HadJet33";
  TString dirall = "XCone_GEN_GenOnly/";
  TString dirmatch = "XCone_GEN_GenOnly_matched/";
  TString dirunmatch = "XCone_GEN_GenOnly_unmatched/";
  TString LegendEntries[] = {"t#bar{t} total", "t#bar{t} fully merged", "t#bar{t} not merged"};

  int colors[] = {TColor::GetColor("#000000"),
                  TColor::GetColor("#008ae6"),
                  TColor::GetColor("#ff8000") };
  int linestyles[] = {kSolid,kSolid,kDashed};
  //int linestyles[] = {kSolid, 7, 3};
  int markerstyles[] = {0,0,0};
  TString y_label = "Events / 20 GeV";
  TString x_label = "#it{m}_{jet} [GeV]";

  //==================================
  //read input histograms
  //==================================
  vector<TH1F*> hists;
  bool isFirstFile = true;
  TH1F *h_all, *h_match, *h_nomatch;
  for(auto FileName: FileNames){
    TFile* file = new TFile(FileName);
    TH1F *h1 = (TH1F*)file->Get(dirall+histname);
    TH1F *h2 = (TH1F*)file->Get(dirmatch+histname);
    TH1F *h3 = (TH1F*)file->Get(dirunmatch+histname);
    if(isFirstFile){
      h_all = h1;
      h_match = h2;
      h_nomatch = h3;
    }
    else{
      h_all->Add(h1);
      h_match->Add(h2);
      h_nomatch->Add(h3);
    }
    isFirstFile = false;
  }
  hists.push_back(h_all);
  hists.push_back(h_match);
  hists.push_back(h_nomatch);


  SetupGlobalStyle();


  //==================================
  //legend settings
  //==================================

  float top = 0.92;
  float xleft = 0.58;
  float xright = 0.92;
  float bottom = 0.68;

  TLegend *leg = new TLegend(xleft,bottom,xright,top, NULL,"brNDC");

  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);

  TLegend *leg2 = new TLegend(xleft,0.55,xright,bottom, NULL,"brNDC");

  leg2->SetFillColor(0);
  leg2->SetLineColor(1);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.03);
  //==================================
  //canvas settings
  //==================================
  TCanvas *can = new TCanvas();
  //  SetupCanvas(can);

  Int_t CanWidth;
  Int_t CanHeight;
  // CanWidth = 400;//default
  CanWidth = 440;
  CanHeight = 400;

  can->SetCanvasSize( CanWidth, CanHeight);
  Float_t yplot = 0.65;
  Float_t yratio = 0.34;


  Float_t y1, y2, y3;                           //  y3 +-------------+
  y3 = 0.99;                                    //     |             |
  y2 = y3-yplot;                                //     |     pad1    |
  y1 = y2-yratio;                               //  y2 |-------------|
  Float_t x1, x2;                               //     |     rp1     |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
                                                //
                                                // No Pad 2!

  TPad *pad = new TPad("pad", "Control Plots 2", x1, y1, x2, y3);
  pad->SetTopMargin(0.055);
  pad->SetBottomMargin(0.16);
  // pad->SetLeftMargin(0.19);//default
  pad->SetLeftMargin(0.19);
  //  pad->SetRightMargin(0.05);//default
  pad->SetRightMargin(0.045);

  pad->Draw();


  //==================================
  //histogram settings and plotting
  //==================================

  for(unsigned int i = 0 ;  i< hists.size(); i++){
    pad->cd();
    SetHist(hists.at(i));

    // normalisation to new NNLO+NNLL cross sections
    // TString hn = hists.at(i)->GetName();
    // if (!hn.Contains("data")){
    //   double fac = 252.9/245.8;
    //   hists.at(i)->Scale(fac);
    //   cout << "Scaling by factor: " << fac << endl;
    // }

    //individual settings
    hists.at(i)->GetYaxis()->SetTitleOffset(1.8);
    hists.at(i)->Rebin(2);

    hists.at(i)->SetLineStyle(linestyles[i]);
    hists.at(i)->SetLineWidth(3);
    hists.at(i)->SetLineColor(colors[i]);
    hists.at(i)->SetMarkerStyle(markerstyles[i]);
    hists.at(i)->SetMarkerColor(colors[i]);

    hists.at(i)->SetTitle(";"+x_label+";" +y_label);

    hists.at(i)->SetAxisRange(0.,400 , "X");
    hists.at(i)->SetAxisRange(0,52000 ,"Y");

    TString legin = "";
    TString drawopt = "";
    if(markerstyles[i]>0){
      legin = "lpe";
      drawopt = "E1";
    }
    else{
      legin = "l";
      drawopt = "Hist";
    }

    if(LegendEntries[i]!= "_")leg->AddEntry(hists.at(i), LegendEntries[i], legin);


    if(i==0) hists.at(i)->Draw(drawopt);
    else  hists.at(i)->Draw(drawopt + "same");
  }

  // leg2->AddEntry((TObject*)0, "\ell + jets", "");
  //leg2->AddEntry((TObject*)0, "R_{jet} = 1.2", "");
  //leg2->AddEntry((TObject*)0, "N_{sub} = 3, R_{sub} = 0.4", "");
  //leg2->AddEntry((TObject*)0, "p_{T} > 400 GeV", "");

  leg->Draw();
  //leg2->Draw();

  // CMSSimLabel(false);
  //CMSPrelimSimLabel(false);

  LumiInfo(137, false);

  XConeLabel(0.59, 0.62);


  cout << "Numbers: total = " << hists.at(0)->Integral() << endl;
  cout << "         fully-merged = " << hists.at(1)->Integral() << endl;
  cout << "         unmerged = " << hists.at(2)->Integral() << endl;

  gPad->RedrawAxis();

  can->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2/ParticleLevel.pdf");

  return 0;
}
