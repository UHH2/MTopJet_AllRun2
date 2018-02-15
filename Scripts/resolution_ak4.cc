void resolution_ak4()
{

  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");


  // get plots with JEC applied
  TH1F * reso_pt[7];
  reso_pt[0] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution");
  reso_pt[1] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_1");
  reso_pt[2] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_2");
  reso_pt[3] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_3");
  reso_pt[4] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_4");
  reso_pt[5] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_5");
  reso_pt[6] = (TH1F*)All->Get("RecGenHists_ak4/PtResolution_6");
  // get plots without JEC applied
  TH1F * reso_pt_noJEC[7];
  reso_pt_noJEC[0] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution");
  reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_1");
  reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_2");
  reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_3");
  reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_4");
  reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_5");
  reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_ak4_noJEC/PtResolution_6");

   // Histogramme fuellen
  Float_t bins[] = {0, 50, 100, 200, 300, 400, 500};
  TH1F * resolution_pt1 = new TH1F("resolution_pt1","AK4, Mean((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen}) ",6, bins);
  TH1F * resolution_pt2 = new TH1F("resolution_pt2","AK4, RMS((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen}) ",6, bins);
  TH1F * resolution_pt1_noJEC = new TH1F("resolution_pt1_noJEC","AK4, Mean((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen}) ",6, bins);
  TH1F * resolution_pt2_noJEC = new TH1F("resolution_pt2_noJEC","AK4, RMS((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen}) ",6, bins);
  

  //-----------------------------------------------------------------------
  //------------------- pt reso -------------------------------------------
  //-----------------------------------------------------------------------
  double mean_pt[7], rms_pt[7], mean_error_pt[7], rms_error_pt[7];
  TF1 *fit_pt[7];

  for(unsigned int i=1; i<=6; i++){
    reso_pt[i]->Fit("gaus");
    fit_pt[i] = reso_pt[i]->GetFunction("gaus");
    mean_pt[i] = fit_pt[i]->GetParameter(1);
    rms_pt[i] = fit_pt[i]->GetParameter(2);
    mean_error_pt[i] = fit_pt[i]->GetParError(1);
    rms_error_pt[i] = fit_pt[i]->GetParError(2);
  }
  TAxis *xaxis_pt1 = resolution_pt1->GetXaxis();
  TAxis *xaxis_pt2 = resolution_pt2->GetXaxis();
  for(unsigned int i=1; i<=6; i++){
    resolution_pt1->Fill(xaxis_pt1->GetBinCenter(i), mean_pt[i]);
    resolution_pt1->SetBinError(i, mean_error_pt[i]);
    resolution_pt2->Fill(xaxis_pt2->GetBinCenter(i), rms_pt[i]);
    resolution_pt2->SetBinError(i, rms_error_pt[i]);
  }
  //-----------------------------------------------------------------------
  //------------------- pt reso (NO JEC) ----------------------------------
  //-----------------------------------------------------------------------
  double mean_pt_noJEC[7], rms_pt_noJEC[7], mean_error_pt_noJEC[7], rms_error_pt_noJEC[7];
  TF1 *fit_pt_noJEC[7];

  for(unsigned int i=1; i<=6; i++){
    reso_pt_noJEC[i]->Fit("gaus");
    fit_pt_noJEC[i] = reso_pt_noJEC[i]->GetFunction("gaus");
    mean_pt_noJEC[i] = fit_pt_noJEC[i]->GetParameter(1);
    rms_pt_noJEC[i] = fit_pt_noJEC[i]->GetParameter(2);
    mean_error_pt_noJEC[i] = fit_pt_noJEC[i]->GetParError(1);
    rms_error_pt_noJEC[i] = fit_pt_noJEC[i]->GetParError(2);
  }
  TAxis *xaxis_pt1_noJEC = resolution_pt1_noJEC->GetXaxis();
  TAxis *xaxis_pt2_noJEC = resolution_pt2_noJEC->GetXaxis();
  for(unsigned int i=1; i<=6; i++){
    resolution_pt1_noJEC->Fill(xaxis_pt1_noJEC->GetBinCenter(i), mean_pt_noJEC[i]);
    resolution_pt1_noJEC->SetBinError(i, mean_error_pt_noJEC[i]);
    resolution_pt2_noJEC->Fill(xaxis_pt2_noJEC->GetBinCenter(i), rms_pt_noJEC[i]);
    resolution_pt2_noJEC->SetBinError(i, rms_error_pt_noJEC[i]);
  }
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);

  resolution_pt1->GetYaxis()->SetRangeUser(-0.3, 0.3);
  resolution_pt2->GetYaxis()->SetRangeUser(0, 0.3);
  resolution_pt1->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_pt1->GetYaxis()->SetTitle("Mean((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  resolution_pt2->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_pt2->GetYaxis()->SetTitle("RMS((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  resolution_pt1->SetLineWidth(2);
  resolution_pt2->SetLineWidth(2);

  resolution_pt1_noJEC->GetYaxis()->SetRangeUser(-0.3, 0.3);
  resolution_pt2_noJEC->GetYaxis()->SetRangeUser(0, 0.3);
  resolution_pt1_noJEC->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_pt1_noJEC->GetYaxis()->SetTitle("Mean((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  resolution_pt2_noJEC->GetXaxis()->SetTitle("p_{T}^{gen}");
  resolution_pt2_noJEC->GetYaxis()->SetTitle("RMS((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  resolution_pt1_noJEC->SetLineWidth(2);
  resolution_pt2_noJEC->SetLineWidth(2);
  resolution_pt1_noJEC->SetLineColor(2);
  resolution_pt2_noJEC->SetLineColor(2);

  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------

  TCanvas *c2 = new TCanvas();
  resolution_pt1->Draw("E");
  resolution_pt1_noJEC->Draw("SAME E");
  leg2 = new TLegend(0.65,0.85,0.85,0.65);
  leg2->AddEntry(resolution_pt1,"JEC applied","l");
  leg2->AddEntry(resolution_pt1_noJEC,"no JEC applied","l");
  leg2->Draw();
  gPad->RedrawAxis();
  c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_ak4_pt_mean.pdf"); 

  TCanvas *c4 = new TCanvas();
  resolution_pt2->Draw("E");
  resolution_pt2_noJEC->Draw("SAME E");
  leg4 = new TLegend(0.65,0.85,0.85,0.65);
  leg4->AddEntry(resolution_pt2,"JEC applied","l");
  leg4->AddEntry(resolution_pt2_noJEC,"no JEC applied","l");
  leg4->Draw();
  gPad->RedrawAxis();
  c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_ak4_pt_rms.pdf"); 
}

