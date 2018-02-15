void resolution_XCone_GenSelOnly()
{

  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.Background.root");

 
  // get plots
  TH1F * reso0 = (TH1F*)All->Get("RecGenHists_GenOnly_all/MassResolution");
  TH1F * reso1 = (TH1F*)All->Get("RecGenHists_GenOnly_000to100/MassResolution");
  TH1F * reso2 = (TH1F*)All->Get("RecGenHists_GenOnly_100to150/MassResolution");
  TH1F * reso3 = (TH1F*)All->Get("RecGenHists_GenOnly_150to200/MassResolution");
  TH1F * reso4 = (TH1F*)All->Get("RecGenHists_GenOnly_200to250/MassResolution");
  TH1F * reso5 = (TH1F*)All->Get("RecGenHists_GenOnly_250to300/MassResolution");
  TH1F * reso6 = (TH1F*)All->Get("RecGenHists_GenOnly_300to500/MassResolution");

  TH1F * reso_pt0 = (TH1F*)All->Get("RecGenHists_GenOnly_all/PtResolution");
  TH1F * reso_pt1 = (TH1F*)All->Get("RecGenHists_GenOnly_000to100/PtResolution");
  TH1F * reso_pt2 = (TH1F*)All->Get("RecGenHists_GenOnly_100to150/PtResolution");
  TH1F * reso_pt3 = (TH1F*)All->Get("RecGenHists_GenOnly_150to200/PtResolution");
  TH1F * reso_pt4 = (TH1F*)All->Get("RecGenHists_GenOnly_200to250/PtResolution");
  TH1F * reso_pt5 = (TH1F*)All->Get("RecGenHists_GenOnly_250to300/PtResolution");
  TH1F * reso_pt6 = (TH1F*)All->Get("RecGenHists_GenOnly_300to500/PtResolution");

   // Histogramme fuellen
  Float_t bins[] = {0, 100, 150, 200, 250, 300, 500};
  TH1F * resolution1 = new TH1F("resolution1","Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution2 = new TH1F("resolution2","RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt1 = new TH1F("resolution_pt1","Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt2 = new TH1F("resolution_pt2","RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);

  //-----------------------------------------------------------------------
  //------------------- mass reso -----------------------------------------
  //-----------------------------------------------------------------------
  reso0->Fit("gaus");
  TF1 *fit0 = reso0->GetFunction("gaus");

  reso1->Fit("gaus");
  TF1 *fit1 = reso1->GetFunction("gaus");
  Double_t p1a = fit1->GetParameter(1);
  Double_t p1b = fit1->GetParameter(2);
  Double_t e1a = fit1->GetParError(1);
  Double_t e1b = fit1->GetParError(2);
 
  reso2->Fit("gaus");
  TF1 *fit2 = reso2->GetFunction("gaus");
  Double_t p2a = fit2->GetParameter(1);
  Double_t p2b = fit2->GetParameter(2);
  Double_t e2a = fit2->GetParError(1);
  Double_t e2b = fit2->GetParError(2);

  reso3->Fit("gaus");
  TF1 *fit3 = reso3->GetFunction("gaus");
  Double_t p3a = fit3->GetParameter(1);
  Double_t p3b = fit3->GetParameter(2);
  Double_t e3a = fit3->GetParError(1);
  Double_t e3b = fit3->GetParError(2);

  reso4->Fit("gaus");
  TF1 *fit4 = reso4->GetFunction("gaus");
  Double_t p4a = fit4->GetParameter(1);
  Double_t p4b = fit4->GetParameter(2);
  Double_t e4a = fit4->GetParError(1);
  Double_t e4b = fit4->GetParError(2);

  reso5->Fit("gaus");
  TF1 *fit5 = reso5->GetFunction("gaus");
  Double_t p5a = fit5->GetParameter(1);
  Double_t p5b = fit5->GetParameter(2);
  Double_t e5a = fit5->GetParError(1);
  Double_t e5b = fit5->GetParError(2);

  reso6->Fit("gaus");
  TF1 *fit6 = reso6->GetFunction("gaus");
  Double_t p6a = fit6->GetParameter(1);
  Double_t p6b = fit6->GetParameter(2);
  Double_t e6a = fit6->GetParError(1);
  Double_t e6b = fit6->GetParError(2);

  TAxis *xaxis1 = resolution1->GetXaxis();
  
  resolution1->Fill(xaxis1->GetBinCenter(1), p1a);
  resolution1->Fill(xaxis1->GetBinCenter(2), p2a);
  resolution1->Fill(xaxis1->GetBinCenter(3), p3a);
  resolution1->Fill(xaxis1->GetBinCenter(4), p4a);
  resolution1->Fill(xaxis1->GetBinCenter(5), p5a);
  resolution1->Fill(xaxis1->GetBinCenter(6), p6a);
  resolution1->SetBinError(1, e1a);
  resolution1->SetBinError(2, e2a);
  resolution1->SetBinError(3, e3a);
  resolution1->SetBinError(4, e4a);
  resolution1->SetBinError(5, e5a);
  resolution1->SetBinError(6, e6a);

  TAxis *xaxis2 = resolution2->GetXaxis();
  
  resolution2->Fill(xaxis2->GetBinCenter(1), p1b);
  resolution2->Fill(xaxis2->GetBinCenter(2), p2b);
  resolution2->Fill(xaxis2->GetBinCenter(3), p3b);
  resolution2->Fill(xaxis2->GetBinCenter(4), p4b);
  resolution2->Fill(xaxis2->GetBinCenter(5), p5b);
  resolution2->Fill(xaxis2->GetBinCenter(6), p6b);
  resolution2->SetBinError(1, e1b);
  resolution2->SetBinError(2, e2b);
  resolution2->SetBinError(3, e3b);
  resolution2->SetBinError(4, e4b);
  resolution2->SetBinError(5, e5b);
  resolution2->SetBinError(6, e6b);
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------


  //-----------------------------------------------------------------------
  //------------------- pt reso -------------------------------------------
  //-----------------------------------------------------------------------

  reso_pt0->Fit("gaus");
  TF1 *fit_pt0 = reso_pt0->GetFunction("gaus");
 
  reso_pt1->Fit("gaus");
  TF1 *fit_pt1 = reso_pt1->GetFunction("gaus");
  Double_t p1a_pt = fit_pt1->GetParameter(1);
  Double_t p1b_pt = fit_pt1->GetParameter(2);
  Double_t e1a_pt = fit_pt1->GetParError(1);
  Double_t e1b_pt = fit_pt1->GetParError(2);

  reso_pt2->Fit("gaus");
  TF1 *fit_pt2 = reso_pt2->GetFunction("gaus");
  Double_t p2a_pt = fit_pt2->GetParameter(1);
  Double_t p2b_pt = fit_pt2->GetParameter(2);
  Double_t e2a_pt = fit_pt2->GetParError(1);
  Double_t e2b_pt = fit_pt2->GetParError(2);

  reso_pt3->Fit("gaus");
  TF1 *fit_pt3 = reso_pt3->GetFunction("gaus");
  Double_t p3a_pt = fit_pt3->GetParameter(1);
  Double_t p3b_pt = fit_pt3->GetParameter(2);
  Double_t e3a_pt = fit_pt3->GetParError(1);
  Double_t e3b_pt = fit_pt3->GetParError(2);

  reso_pt4->Fit("gaus");
  TF1 *fit_pt4 = reso_pt4->GetFunction("gaus");
  Double_t p4a_pt = fit_pt4->GetParameter(1);
  Double_t p4b_pt = fit_pt4->GetParameter(2);
  Double_t e4a_pt = fit_pt4->GetParError(1);
  Double_t e4b_pt = fit_pt4->GetParError(2);

  reso_pt5->Fit("gaus");
  TF1 *fit_pt5 = reso_pt5->GetFunction("gaus");
  Double_t p5a_pt = fit_pt5->GetParameter(1);
  Double_t p5b_pt = fit_pt5->GetParameter(2);
  Double_t e5a_pt = fit_pt5->GetParError(1);
  Double_t e5b_pt = fit_pt5->GetParError(2);

  reso_pt6->Fit("gaus");
  TF1 *fit_pt6 = reso_pt6->GetFunction("gaus");
  Double_t p6a_pt = fit_pt6->GetParameter(1);
  Double_t p6b_pt = fit_pt6->GetParameter(2);
  Double_t e6a_pt = fit_pt6->GetParError(1);
  Double_t e6b_pt = fit_pt6->GetParError(2);

  TAxis *xaxis_pt1 = resolution_pt1->GetXaxis();
  
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(1), p1a_pt);
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(2), p2a_pt);
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(3), p3a_pt);
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(4), p4a_pt);
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(5), p5a_pt);
  resolution_pt1->Fill(xaxis_pt1->GetBinCenter(6), p6a_pt);
  resolution_pt1->SetBinError(1, e1a_pt);
  resolution_pt1->SetBinError(2, e2a_pt);
  resolution_pt1->SetBinError(3, e3a_pt);
  resolution_pt1->SetBinError(4, e4a_pt);
  resolution_pt1->SetBinError(5, e5a_pt);
  resolution_pt1->SetBinError(6, e6a_pt);

  TAxis *xaxis_pt2 = resolution_pt2->GetXaxis();
  
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(1), p1b_pt);
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(2), p2b_pt);
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(3), p3b_pt);
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(4), p4b_pt);
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(5), p5b_pt);
  resolution_pt2->Fill(xaxis_pt2->GetBinCenter(6), p6b_pt);
  resolution_pt2->SetBinError(1, e1b_pt);
  resolution_pt2->SetBinError(2, e2b_pt);
  resolution_pt2->SetBinError(3, e3b_pt);
  resolution_pt2->SetBinError(4, e4b_pt);
  resolution_pt2->SetBinError(5, e5b_pt);
  resolution_pt2->SetBinError(6, e6b_pt);
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kTRUE);
  resolution1->GetYaxis()->SetRangeUser(-0.1, 0.6);
  resolution2->GetYaxis()->SetRangeUser(-0.1, 0.5);
  reso0->GetXaxis()->SetTitle("(M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen} ");
  reso0->GetYaxis()->SetTitle("# Events");
  resolution1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  resolution1->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  resolution2->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  reso0->SetLineWidth(2);
  resolution1->SetLineWidth(2);
  resolution2->SetLineWidth(2);

  resolution_pt1->GetYaxis()->SetRangeUser(-0.3, 0.1);
  resolution_pt2->GetYaxis()->SetRangeUser(0, 0.1);
  reso_pt0->GetXaxis()->SetTitle("(p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}");
  reso_pt0->GetYaxis()->SetTitle("# Events");
  resolution_pt1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  resolution_pt1->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  resolution_pt2->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  reso_pt0->SetLineWidth(2);
  resolution_pt1->SetLineWidth(2);
  resolution_pt2->SetLineWidth(2);

  //turn on stats for gaussian fit
  reso0->SetStats(kTRUE);
  reso_pt0->SetStats(kTRUE);
  resolution1->SetStats(kFALSE);
  resolution2->SetStats(kFALSE);
  resolution_pt1->SetStats(kFALSE);
  resolution_pt2->SetStats(kFALSE);


  // Draw Histogram
  TCanvas *A = new TCanvas();
  reso_pt1->Draw("E");

  TCanvas *C = new TCanvas();
  reso0->Draw("E");

  TCanvas *C_pt = new TCanvas();
  reso_pt0->Draw("E");
  

  TCanvas *all = new TCanvas();
  all->Divide(2,2);
  all->cd(1);
  resolution1->Draw("E");
  all->cd(2);
  resolution_pt1->Draw("E");
  all->cd(3);
  resolution2->Draw("E");
  all->cd(4);
  resolution_pt2->Draw("E");
}

