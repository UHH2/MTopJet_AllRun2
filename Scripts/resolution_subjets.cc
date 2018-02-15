void resolution_subjets()
{

  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");

  bool use_rec = false;
  n_ptbin = 11;
  TH1F * reso_pt[n_ptbin];
  TH1F * reso_pt_noJEC[n_ptbin];
  TH1F * reso_pt_corrected[n_ptbin];

  std::string dir_jec = "RecGenHists_subjets/";
  std::string dir_raw = "RecGenHists_subjets_noJEC/";
  std::string dir_cor = "RecGenHists_subjets_corrected/";

  std::string name_jec;
  std::string name_raw;
  std::string name_cor;

  for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){
    if(use_rec){
      name_jec = dir_jec + "PtResolution_" + "rec" + std::to_string(ptbin);
      name_raw = dir_raw + "PtResolution_" + "rec" + std::to_string(ptbin);
      name_cor = dir_cor + "PtResolution_" + "rec" + std::to_string(ptbin);
    }
    else{
      name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
      name_raw = dir_raw + "PtResolution_" + std::to_string(ptbin);
      name_cor = dir_cor + "PtResolution_" + std::to_string(ptbin);
    }
    reso_pt[ptbin] = (TH1F*)All->Get(name_jec.c_str());
    reso_pt_noJEC[ptbin] = (TH1F*)All->Get(name_raw.c_str());
    reso_pt_corrected[ptbin] = (TH1F*)All->Get(name_cor.c_str());
  }

  // Histogramme fuellen
  Float_t bins[] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 600};
  int n_bins = sizeof(bins)/sizeof(bins[0]) - 1;
  TH1F * resolution_pt1 = new TH1F("resolution_pt1"," ",n_bins, bins);
  TH1F * resolution_pt2 = new TH1F("resolution_pt2"," ",n_bins, bins);
  TH1F * resolution_pt1_noJEC = new TH1F("resolution_pt1_noJEC"," ",n_bins, bins);
  TH1F * resolution_pt2_noJEC = new TH1F("resolution_pt2_noJEC"," ",n_bins, bins);
  TH1F * resolution_pt1_corrected = new TH1F("resolution_pt1_corrected"," ",n_bins, bins);
  TH1F * resolution_pt2_corrected = new TH1F("resolution_pt2_corrected"," ",n_bins, bins); 

  //-----------------------------------------------------------------------
  //------------------- pt reso -------------------------------------------
  //-----------------------------------------------------------------------
  double mean_pt[11], rms_pt[11], mean_error_pt[11], rms_error_pt[11];
  TF1 *fit_pt[11];

  for(unsigned int i=1; i<=n_bins; i++){
    reso_pt[i]->Fit("gaus");
    fit_pt[i] = reso_pt[i]->GetFunction("gaus");
    mean_pt[i] = fit_pt[i]->GetParameter(1);
    rms_pt[i] = fit_pt[i]->GetParameter(2);
    mean_error_pt[i] = fit_pt[i]->GetParError(1);
    rms_error_pt[i] = fit_pt[i]->GetParError(2);
  }
  TAxis *xaxis_pt1 = resolution_pt1->GetXaxis();
  TAxis *xaxis_pt2 = resolution_pt2->GetXaxis();
  for(unsigned int i=1; i<=n_bins; i++){
    resolution_pt1->Fill(xaxis_pt1->GetBinCenter(i), mean_pt[i]);
    resolution_pt1->SetBinError(i, mean_error_pt[i]);
    resolution_pt2->Fill(xaxis_pt2->GetBinCenter(i), rms_pt[i]);
    resolution_pt2->SetBinError(i, rms_error_pt[i]);
  }
  //-----------------------------------------------------------------------
  //------------------- pt reso (NO JEC) ----------------------------------
  //-----------------------------------------------------------------------
  double mean_pt_noJEC[11], rms_pt_noJEC[11], mean_error_pt_noJEC[11], rms_error_pt_noJEC[11];
  TF1 *fit_pt_noJEC[11];

  for(unsigned int i=1; i<=n_bins; i++){
    reso_pt_noJEC[i]->Fit("gaus");
    fit_pt_noJEC[i] = reso_pt_noJEC[i]->GetFunction("gaus");
    mean_pt_noJEC[i] = fit_pt_noJEC[i]->GetParameter(1);
    rms_pt_noJEC[i] = fit_pt_noJEC[i]->GetParameter(2);
    mean_error_pt_noJEC[i] = fit_pt_noJEC[i]->GetParError(1);
    rms_error_pt_noJEC[i] = fit_pt_noJEC[i]->GetParError(2);
  }
  TAxis *xaxis_pt1_noJEC = resolution_pt1_noJEC->GetXaxis();
  TAxis *xaxis_pt2_noJEC = resolution_pt2_noJEC->GetXaxis();
  for(unsigned int i=1; i<=n_bins; i++){
    resolution_pt1_noJEC->Fill(xaxis_pt1_noJEC->GetBinCenter(i), mean_pt_noJEC[i]);
    resolution_pt1_noJEC->SetBinError(i, mean_error_pt_noJEC[i]);
    resolution_pt2_noJEC->Fill(xaxis_pt2_noJEC->GetBinCenter(i), rms_pt_noJEC[i]);
    resolution_pt2_noJEC->SetBinError(i, rms_error_pt_noJEC[i]);
  }
  //-----------------------------------------------------------------------
  //------------------- pt reso corrected ---------------------------------
  //-----------------------------------------------------------------------
  double mean_pt_corrected[11], rms_pt_corrected[11], mean_error_pt_corrected[11], rms_error_pt_corrected[11];
  TF1 *fit_pt_corrected[11];

  for(unsigned int i=1; i<=n_bins; i++){
    reso_pt_corrected[i]->Fit("gaus");
    fit_pt_corrected[i] = reso_pt_corrected[i]->GetFunction("gaus");
    mean_pt_corrected[i] = fit_pt_corrected[i]->GetParameter(1);
    rms_pt_corrected[i] = fit_pt_corrected[i]->GetParameter(2);
    mean_error_pt_corrected[i] = fit_pt_corrected[i]->GetParError(1);
    rms_error_pt_corrected[i] = fit_pt_corrected[i]->GetParError(2);
  }
  TAxis *xaxis_pt1_corrected = resolution_pt1_corrected->GetXaxis();
  TAxis *xaxis_pt2_corrected = resolution_pt2_corrected->GetXaxis();
  for(unsigned int i=1; i<=n_bins; i++){
    resolution_pt1_corrected->Fill(xaxis_pt1_corrected->GetBinCenter(i), mean_pt_corrected[i]);
    resolution_pt1_corrected->SetBinError(i, mean_error_pt_corrected[i]);
    resolution_pt2_corrected->Fill(xaxis_pt2_corrected->GetBinCenter(i), rms_pt_corrected[i]);
    resolution_pt2_corrected->SetBinError(i, rms_error_pt_corrected[i]);
  }
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------

  zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);

  resolution_pt1->GetYaxis()->SetRangeUser(-0.3, 0.3);
  resolution_pt2->GetYaxis()->SetRangeUser(0, 0.3);

  resolution_pt1->GetXaxis()->SetNdivisions(505);
  resolution_pt2->GetXaxis()->SetNdivisions(505);
  resolution_pt1->GetYaxis()->SetNdivisions(505);
  resolution_pt2->GetYaxis()->SetNdivisions(505);

  resolution_pt1->GetXaxis()->SetTitleSize(0.05);
  resolution_pt1->GetYaxis()->SetTitleSize(0.06);
  resolution_pt1->GetXaxis()->SetTitleOffset(0.9);
  resolution_pt1->GetYaxis()->SetTitleOffset(1.1);

  resolution_pt2->GetXaxis()->SetTitleSize(0.05);
  resolution_pt2->GetYaxis()->SetTitleSize(0.06);
  resolution_pt2->GetXaxis()->SetTitleOffset(0.9);
  resolution_pt2->GetYaxis()->SetTitleOffset(1.1);

  if(!use_rec) resolution_pt1->GetXaxis()->SetTitle("p_{T}^{gen}");
  else resolution_pt1->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_pt1->GetYaxis()->SetTitle("Mean((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  if(!use_rec) resolution_pt2->GetXaxis()->SetTitle("p_{T}^{gen}");
  else resolution_pt2->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_pt2->GetYaxis()->SetTitle("RMS((p_{T}^{rec} - p_{T}^{gen}) /p_{T}^{gen})");
  resolution_pt1->SetLineWidth(4);
  resolution_pt2->SetLineWidth(4);
  resolution_pt1->SetLineColor(kRed+1);
  resolution_pt2->SetLineColor(kRed+1);

  resolution_pt1_noJEC->SetLineWidth(4);
  resolution_pt2_noJEC->SetLineWidth(4);
  resolution_pt1_noJEC->SetLineColor(kOrange+1);
  resolution_pt2_noJEC->SetLineColor(kOrange+1);


  resolution_pt1_corrected->SetLineWidth(4);
  resolution_pt2_corrected->SetLineWidth(4);
  resolution_pt1_corrected->SetLineColor(kAzure+7);
  resolution_pt2_corrected->SetLineColor(kAzure+7);
  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------

  bool after = true;

  TCanvas *c2 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  resolution_pt1->Draw("SAME E");
  zero_line->Draw("SAME");
  resolution_pt1->Draw("SAME E");
  resolution_pt1_noJEC->Draw("SAME E");
  if(after) resolution_pt1_corrected->Draw("SAME E");
  leg2 = new TLegend(0.5,0.85,0.85,0.65);
  leg2->AddEntry(resolution_pt1,"JEC applied","l");
  leg2->AddEntry(resolution_pt1_noJEC,"no JEC applied","l");
  if(after) leg2->AddEntry(resolution_pt1_corrected,"corrected","l");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  gPad->RedrawAxis();
  if(after){
    if(!use_rec) c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_mean_after.pdf"); 
    else c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_mean_rec_after.pdf");
  }
  else{
    if(!use_rec) c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_mean_before.pdf"); 
    else c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_mean_rec_before.pdf"); 
  }

  TCanvas *c4 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  resolution_pt2->Draw("SAME E");
  resolution_pt2_noJEC->Draw("SAME E");
  if(after) resolution_pt2_corrected->Draw("SAME E");
  leg4 = new TLegend(0.5,0.85,0.85,0.65);
  leg4->AddEntry(resolution_pt2,"JEC applied","l");
  leg4->AddEntry(resolution_pt2_noJEC,"no JEC applied","l");
  if(after) leg4->AddEntry(resolution_pt2_corrected,"corrected","l");
  leg4->SetTextSize(0.05);
  leg4->Draw();
  gPad->RedrawAxis();
  if(after){
    if(!use_rec) c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_rms_after.pdf"); 
    else c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_rms_rec_after.pdf");
  } 
  else{
    if(!use_rec) c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_rms_before.pdf"); 
    else c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_Subjets/pt_rms_rec_before.pdf");
  } 

}

