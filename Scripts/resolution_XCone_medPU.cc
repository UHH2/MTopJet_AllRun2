void resolution_XCone_medPU(bool mass)
{

  // use mass or pt bins?
  // bool mass = false;

  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
 

  // get plots with JEC applied
  TH1F * reso[7];
  TH1F * reso_pt[7];
  TH1F * reso_noJEC[7];
  TH1F * reso_pt_noJEC[7];

  if(mass){
    reso[0] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution");
    reso[1] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass1");
    reso[2] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass2");
    reso[3] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass3");
    reso[4] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass4");
    reso[5] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass5");
    reso[6] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_mass6");

    reso_pt[0] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution");
    reso_pt[1] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass1");
    reso_pt[2] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass2");
    reso_pt[3] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass3");
    reso_pt[4] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass4");
    reso_pt[5] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass5");
    reso_pt[6] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_mass6");

    reso_noJEC[0] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution");
    reso_noJEC[1] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass1");
    reso_noJEC[2] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass2");
    reso_noJEC[3] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass3");
    reso_noJEC[4] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass4");
    reso_noJEC[5] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass5");
    reso_noJEC[6] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_mass6");

    reso_pt_noJEC[0] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution");
    reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass1");
    reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass2");
    reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass3");
    reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass4");
    reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass5");
    reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_mass6");
  }
  else{
    reso[0] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution");
    reso[1] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt1");
    reso[2] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt2");
    reso[3] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt3");
    reso[4] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt4");
    reso[5] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt5");
    reso[6] = (TH1F*)All->Get("RecGenHists_medPU/MassResolution_pt6");

    reso_pt[0] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution");
    reso_pt[1] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt1");
    reso_pt[2] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt2");
    reso_pt[3] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt3");
    reso_pt[4] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt4");
    reso_pt[5] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt5");
    reso_pt[6] = (TH1F*)All->Get("RecGenHists_medPU/PtResolution_pt6");

    reso_noJEC[0] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution");
    reso_noJEC[1] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt1");
    reso_noJEC[2] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt2");
    reso_noJEC[3] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt3");
    reso_noJEC[4] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt4");
    reso_noJEC[5] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt5");
    reso_noJEC[6] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/MassResolution_pt6");

    reso_pt_noJEC[0] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution");
    reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt1");
    reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt2");
    reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt3");
    reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt4");
    reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt5");
    reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_medPU_noJEC/PtResolution_pt6");
  }

   // Histogramme fuellen
  Float_t bins[7];
  if(mass){
    bins[0] = 0;
    bins[1] = 100;
    bins[2] = 150;
    bins[3] = 200;
    bins[4] = 250;
    bins[5] = 300;
    bins[6] = 500;
  }
  else{
    bins[0] = 400;
    bins[1] = 450;
    bins[2] = 500;
    bins[3] = 600;
    bins[4] = 700;
    bins[5] = 800;
    bins[6] = 2000;
  }
  TH1F * resolution1 = new TH1F("resolution1","Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution2 = new TH1F("resolution2","RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt1 = new TH1F("resolution_pt1","Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt2 = new TH1F("resolution_pt2","RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);
  TH1F * resolution1_noJEC = new TH1F("resolution1_noJEC","Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution2_noJEC = new TH1F("resolution2_noJEC","RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt1_noJEC = new TH1F("resolution_pt1_noJEC","Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);
  TH1F * resolution_pt2_noJEC = new TH1F("resolution_pt2_noJEC","RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}) ",6, bins);
  //-----------------------------------------------------------------------
  //------------------- mass reso -----------------------------------------
  //-----------------------------------------------------------------------
  double mean[7], rms[7], mean_error[7], rms_error[7];
  TF1 * fit[7];

  for(unsigned int i=1; i<=6; i++){
    reso[i]->Fit("gaus");
    fit[i] = reso[i]->GetFunction("gaus");
    mean[i] = fit[i]->GetParameter(1);
    rms[i] = fit[i]->GetParameter(2);
    mean_error[i] = fit[i]->GetParError(1);
    rms_error[i] = fit[i]->GetParError(2);
  }

  TAxis *xaxis1 = resolution1->GetXaxis();
  TAxis *xaxis2 = resolution2->GetXaxis();

  for(unsigned int i=1; i<=6; i++){
    resolution1->Fill(xaxis1->GetBinCenter(i), mean[i]);
    resolution1->SetBinError(i, mean_error[i]);
    resolution2->Fill(xaxis2->GetBinCenter(i), rms[i]);
    resolution2->SetBinError(i, rms_error[i]);
  }
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
  //------------------- mass reso (NO JEC)  -------------------------------
  //-----------------------------------------------------------------------
  double mean_noJEC[7], rms_noJEC[7], mean_error_noJEC[7], rms_error_noJEC[7];
  TF1 *fit_noJEC[7];

  for(unsigned int i=1; i<=6; i++){
    reso_noJEC[i]->Fit("gaus");
    fit_noJEC[i] = reso_noJEC[i]->GetFunction("gaus");
    mean_noJEC[i] = fit_noJEC[i]->GetParameter(1);
    rms_noJEC[i] = fit_noJEC[i]->GetParameter(2);
    mean_error_noJEC[i] = fit_noJEC[i]->GetParError(1);
    rms_error_noJEC[i] = fit_noJEC[i]->GetParError(2);
  }

  TAxis *xaxis1_noJEC = resolution1_noJEC->GetXaxis();
  TAxis *xaxis2_noJEC = resolution2_noJEC->GetXaxis();

  for(unsigned int i=1; i<=6; i++){
    resolution1_noJEC->Fill(xaxis1_noJEC->GetBinCenter(i), mean_noJEC[i]);
    resolution1_noJEC->SetBinError(i, mean_error_noJEC[i]);
    resolution2_noJEC->Fill(xaxis2_noJEC->GetBinCenter(i), rms_noJEC[i]);
    resolution2_noJEC->SetBinError(i, rms_error_noJEC[i]);
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
  resolution1->GetYaxis()->SetRangeUser(-0.1, 0.5);
  if(mass)resolution1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution1->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution1->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution1->SetLineWidth(2);

  resolution2->GetYaxis()->SetRangeUser(0.0, 0.3);
  if(mass)resolution2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution2->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution2->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2->GetYaxis()->SetTitleOffset(1.1);
  resolution2->SetLineWidth(2);

  resolution_pt1->GetYaxis()->SetRangeUser(-0.1, 0.1);
  if(mass)resolution_pt1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution_pt1->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution_pt1->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt1->GetYaxis()->SetTitleOffset(1.2);
  resolution_pt1->SetLineWidth(2);

  resolution_pt2->GetYaxis()->SetRangeUser(0.0, 0.1);
  if(mass)resolution_pt2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution_pt2->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution_pt2->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2->GetYaxis()->SetTitleOffset(1.2);
  resolution_pt2->SetLineWidth(2);

  resolution1_noJEC->GetYaxis()->SetRangeUser(-0.1, 0.5);
  if(mass)resolution1_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution1_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution1_noJEC->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution1_noJEC->SetLineWidth(2);
  resolution1_noJEC->SetLineColor(2);

  resolution2_noJEC->GetYaxis()->SetRangeUser(0.0, 0.3);
  if(mass)resolution2_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution2_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution2_noJEC->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2_noJEC->SetLineWidth(2);
  resolution2_noJEC->SetLineColor(2);

  resolution_pt1_noJEC->GetYaxis()->SetRangeUser(-0.1, 0.1);
  if(mass)resolution_pt1_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution_pt1_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution_pt1_noJEC->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt1_noJEC->SetLineWidth(2);
  resolution_pt1_noJEC->SetLineColor(2);

  resolution_pt2_noJEC->GetYaxis()->SetRangeUser(0.0, 0.1);
  if(mass)resolution_pt2_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else resolution_pt2_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  resolution_pt2_noJEC->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2_noJEC->SetLineWidth(2);
  resolution_pt2_noJEC->SetLineColor(2);

  // reso_pt[0]->SetLineColor(810);
  // reso_pt[0]->SetFillColor(810);
  // reso_pt[0]->SetTitle("JEC applied");
  // reso_pt[0]->GetXaxis()->SetRangeUser(-1, 1);
  // reso_pt[0]->GetXaxis()->SetTitle("(p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}");
  // reso_pt[0]->GetYaxis()->SetTitle("events");

  // reso_pt_noJEC[0]->SetLineColor(810);
  // reso_pt_noJEC[0]->SetFillColor(810);
  // reso_pt_noJEC[0]->SetTitle("no JEC applied");
  // reso_pt_noJEC[0]->GetXaxis()->SetRangeUser(-1, 1);
  // reso_pt_noJEC[0]->GetXaxis()->SetTitle("(p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen}");
  // reso_pt_noJEC[0]->GetYaxis()->SetTitle("events");

  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas();
  resolution1->Draw("E");
  resolution1_noJEC->Draw("SAME E");
  leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(resolution1,"JEC applied","l");
  leg->AddEntry(resolution1_noJEC,"no JEC applied","l");
  leg->Draw();
  gPad->RedrawAxis();
  if(mass)c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_mass_mean_massbin.pdf"); 
  else c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_mass_mean_ptbin.pdf"); 

  TCanvas *c2 = new TCanvas();
  resolution_pt1->Draw("E");
  resolution_pt1_noJEC->Draw("SAME E");
  leg2 = new TLegend(0.65,0.85,0.85,0.65);
  leg2->AddEntry(resolution_pt1,"JEC applied","l");
  leg2->AddEntry(resolution_pt1_noJEC,"no JEC applied","l");
  leg2->Draw();
  gPad->RedrawAxis();
  if(mass)c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_pt_mean_massbin.pdf"); 
  else c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_pt_mean_ptbin.pdf"); 

  TCanvas *c3 = new TCanvas();
  resolution2->Draw("E");
  resolution2_noJEC->Draw("SAME E");
  leg3 = new TLegend(0.65,0.85,0.85,0.65);
  leg3->AddEntry(resolution2,"JEC applied","l");
  leg3->AddEntry(resolution2_noJEC,"no JEC applied","l");
  leg3->Draw();
  gPad->RedrawAxis();
  if(mass)c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_mass_rms_massbin.pdf"); 
  else c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_mass_rms_ptbin.pdf"); 

  TCanvas *c4 = new TCanvas();
  resolution_pt2->Draw("E");
  resolution_pt2_noJEC->Draw("SAME E");
  leg4 = new TLegend(0.65,0.85,0.85,0.65);
  leg4->AddEntry(resolution_pt2,"JEC applied","l");
  leg4->AddEntry(resolution_pt2_noJEC,"no JEC applied","l");
  leg4->Draw();
  gPad->RedrawAxis();
  if(mass)c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_pt_rms_massbin.pdf"); 
  else c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/resolution_medPU_pt_rms_ptbin.pdf"); 


}

