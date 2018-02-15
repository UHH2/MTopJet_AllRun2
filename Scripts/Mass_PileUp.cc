void Mass_PileUp()
{
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TT_file = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get lumi plots -------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TH2F * TT_ = (TH2F*)TT_file->Get("XCone_cor_SF/Mass_Vertices");
  TH1F * mean = new TH1F("mean","number of primary vertices", 10, 0, 50);
  TH2F * TT = new TH2F("TT"," ", 50, 0, 50, 50, 0 ,500);

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  double TT_count = TT->Integral();

  // TT->Scale(1/TT_count);

  // TT->GetZaxis()->SetRangeUser(0.0000001, 0.1);


  TT->SetTitle("t #bar{t}");
  TT->GetXaxis()->SetTitle("number of primary vertices");
  TT->GetYaxis()->SetTitle("m_{jet}");
  TT->GetZaxis()->SetTitle("#Delta N / N");
  TT->GetYaxis()->SetTitleSize(0.06);
  TT->GetXaxis()->SetTitleSize(0.05);
  TT->GetZaxis()->SetTitleSize(0.05);
  TT->GetXaxis()->SetTitleOffset(0.9);
  TT->GetYaxis()->SetTitleOffset(1.1);
  TT->GetZaxis()->SetTitleOffset(0.9);

  TT->GetXaxis()->SetNdivisions(505);
  TT->GetYaxis()->SetNdivisions(505);

  mean->SetTitle("t #bar{t}");
  mean->GetXaxis()->SetTitle("number of primary vertices");
  mean->GetYaxis()->SetTitle("mean m_{jet}");
  mean->GetYaxis()->SetRangeUser(140, 220);

  mean->GetXaxis()->SetTitleOffset(1.1);
  mean->GetYaxis()->SetTitleOffset(1.5);
  mean->GetXaxis()->SetNdivisions(505);
  mean->GetYaxis()->SetNdivisions(505);

  for(int i=1; i<=50; i++){
    double a = 0, b = 0, c = 0;
    for(int j=1; j<=50; j++){
      a += TT_->GetBinContent(i,j) * TT_->GetYaxis()->GetBinCenter(j);
      c += TT_->GetBinContent(i,j);
    }
    if(c != 0) mean->SetBinContent(i, a/c);
    for(int j=1; j<=50; j++){
      b = TT_->GetBinContent(i,j);
      if(c !=0 ) TT->Fill(i, j*10, b/c);
    }
  }

  mean->SetMarkerStyle(20);
  mean->SetMarkerSize(0.8);
  mean->SetLineColor(1);

  line = new TLine(0, 173, 50, 173);
  line->SetLineColor(2);
  line->SetLineWidth(3);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
 

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  TCanvas *b = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogz();
  TT->Draw("COLZ");
  TT->Draw("BOX SAME");
  line->Draw("SAME");
  b->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Mass_PileUp_2D.pdf");

  TCanvas *a = new TCanvas();
  mean->Draw("E1");
  a->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Mass_PileUp_1D.pdf");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
