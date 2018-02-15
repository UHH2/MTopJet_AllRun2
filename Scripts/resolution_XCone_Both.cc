void resolution_XCone_Both()
{

  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");

 
  // get plots
  TH1F * reso0 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution");
  TH1F * reso1 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass1");
  TH1F * reso2 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass2");
  TH1F * reso3 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass3");
  TH1F * reso4 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass4");
  TH1F * reso5 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass5");
  TH1F * reso6 = (TH1F*)All->Get("RecGenHists_Both_corrected/MassResolution_mass6");


   // Histogramme fuellen
  Float_t bins[] = {0, 150, 175, 200, 250, 300, 500};
  TH1F * resolution1 = new TH1F("resolution1","Mass shift",6, bins);
  TH1F * resolution2 = new TH1F("resolution2","Mass resolution",6, bins);

  Float_t bins2[] = {0, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480};
  int n_bins = 12;
  TH1F * bin_width = new TH1F("bin_width","bin width",n_bins, bins2);
  TH1F * bin_width2 = new TH1F("bin_width2","bin width",n_bins, bins2);

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


  for(int bin = 1; bin <= n_bins; bin++){
    double bin_center = bin_width->GetXaxis()->GetBinCenter(bin);
    int rms_bin = resolution2->GetXaxis()->FindBin(bin_center);
    double rms_value = resolution2->GetBinContent(rms_bin);
    double rms_err = resolution2->GetBinError(rms_bin);

    bin_width->SetBinContent(bin, rms_value * bin_center);
    bin_width->SetBinError(bin, rms_err * bin_center);

    bin_width2->SetBinContent(bin, rms_value * bin_center * 2);
    bin_width2->SetBinError(bin, rms_err * bin_center * 2);

    cout << "bin        = " << bin << endl;
    cout << "bin center = " << bin_center << endl;
    cout << "rms bin    = " << rms_bin << endl;
    cout << "rms value  = " << rms_value << endl;
    cout << "bin width  = " << rms_value * bin_center << endl;
    cout << "===============================-=========" << endl;



  }
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------


  

  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);


  TCanvas *c_mean = new TCanvas("mean","",600,600);
  gPad->SetLeftMargin(0.15); 
  resolution1->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution1->GetXaxis()->SetTitle("m_{gen}");
  resolution1->GetYaxis()->SetTitleSize(0.06);
  resolution1->GetXaxis()->SetTitleSize(0.05);
  resolution1->GetXaxis()->SetTitleOffset(0.9);
  resolution1->GetYaxis()->SetTitleOffset(1.1);
  resolution1->GetXaxis()->SetNdivisions(505);
  resolution1->GetYaxis()->SetNdivisions(505);
  resolution1->GetYaxis()->SetTitle("Mean( (m_{rec} - m_{gen}) / m_{gen} )");
  resolution1->SetLineColor(kBlack);
  resolution1->SetMarkerColor(kBlack);
  resolution1->SetMarkerStyle(8);
  resolution1->SetMarkerSize(1);
  resolution1->Draw("E1");
  c_mean->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Binning/MassResolution_mean.pdf");

  TCanvas *c_rms = new TCanvas("rms","",600,600);
  gPad->SetLeftMargin(0.15); 
  resolution2->GetYaxis()->SetRangeUser(0, 0.1);
  resolution2->GetXaxis()->SetTitle("m_{gen}");
  resolution2->GetYaxis()->SetTitleSize(0.06);
  resolution2->GetXaxis()->SetTitleSize(0.05);
  resolution2->GetXaxis()->SetTitleOffset(0.9);
  resolution2->GetYaxis()->SetTitleOffset(1.1);
  resolution2->GetYaxis()->SetTitle("RMS( (m_{rec} - m_{gen}) / m_{gen} )");
  resolution2->GetXaxis()->SetNdivisions(505);
  resolution2->GetYaxis()->SetNdivisions(505);
  resolution2->SetLineColor(kBlack);
  resolution2->SetMarkerColor(kBlack);
  resolution2->SetMarkerStyle(8);
  resolution2->SetMarkerSize(1);
  resolution2->Draw("E1");
  c_rms->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Binning/MassResolution_rms.pdf");

  TCanvas *c_bins = new TCanvas("bins","",600,600);
  gPad->SetLeftMargin(0.15); 
  bin_width->GetXaxis()->SetTitle("m_{gen}");
  bin_width->GetYaxis()->SetTitle("bin width [GeV]");
  bin_width->GetYaxis()->SetTitleSize(0.06);
  bin_width->GetXaxis()->SetTitleSize(0.05);
  bin_width->GetXaxis()->SetTitleOffset(0.9);
  bin_width->GetYaxis()->SetTitleOffset(1.1);
  bin_width->GetXaxis()->SetNdivisions(505);
  bin_width->GetYaxis()->SetNdivisions(505);
  bin_width->SetLineColor(kBlack);
  bin_width->SetMarkerColor(kBlack);
  bin_width->SetMarkerStyle(8);
  bin_width->SetMarkerSize(1);
  bin_width->Draw("E1");
  c_bins->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Binning/Bin_width.pdf");

  TCanvas *c_bins2 = new TCanvas("bins2","",600,600);
  gPad->SetLeftMargin(0.15); 
  bin_width2->GetXaxis()->SetTitle("m_{gen}");
  bin_width2->GetYaxis()->SetTitle("bin width [GeV]");
  bin_width2->GetYaxis()->SetTitleSize(0.06);
  bin_width2->GetXaxis()->SetTitleSize(0.05);
  bin_width2->GetXaxis()->SetTitleOffset(0.9);
  bin_width2->GetYaxis()->SetTitleOffset(1.1);
  bin_width2->GetXaxis()->SetNdivisions(505);
  bin_width2->GetYaxis()->SetNdivisions(505);
  bin_width2->SetLineColor(kBlack);
  bin_width2->SetMarkerColor(kBlack);
  bin_width2->SetMarkerStyle(8);
  bin_width2->SetMarkerSize(1);
  bin_width2->Draw("E1");
  c_bins2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Binning/Bin_width2.pdf");


}

