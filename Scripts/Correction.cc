void Correction()
{
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *file = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get plots ------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  Float_t xbins[] = {0, 80, 130, 180, 250, 350, 500};
  Float_t ybins[] = {-4, -1.5, -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 4};

  TH2F * resolution_mean = new TH2F("resolution_mean","Mean",6, xbins, 12, ybins);
  TH2F * resolution_rms = new TH2F("resolution_mean","RMS",6, xbins, 12, ybins);
  TH2F * resolution_mean_err = new TH2F("resolution_mean_err","Mean Error",6, xbins, 12, ybins);
  TH2F * resolution_rms_err = new TH2F("resolution_mean_err","RMS Error",6, xbins, 12, ybins);

  int no_ptbins = 6; // rows
  int no_etabins = 12; // columns
  std::string dir = "CorrectionHists/PtReso_";
  std::string filename;
  TH1F * reso[no_ptbins][no_etabins];
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      filename = dir + std::to_string(pt_bin) + std::to_string(eta_bin) ;
      reso[pt_bin][eta_bin] = (TH1F*)file->Get(filename.c_str());
    }
  }

  
  double m, s;
  double mean[no_ptbins][no_etabins], rms[no_ptbins][no_etabins], factor[no_ptbins][no_etabins];
  double mean_err[no_ptbins][no_etabins], rms_err[no_ptbins][no_etabins];
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      m = reso[pt_bin][eta_bin]->GetMean();
      s = reso[pt_bin][eta_bin]->GetRMS();
      reso[pt_bin][eta_bin]->GetXaxis()->SetRangeUser(m-2*s, m+2*s); // only consider events inside 2sigma from mean
      mean[pt_bin][eta_bin] = reso[pt_bin][eta_bin]->GetMean();
      mean_err[pt_bin][eta_bin] = reso[pt_bin][eta_bin]->GetMeanError();
      rms[pt_bin][eta_bin] = reso[pt_bin][eta_bin]->GetRMS();
      rms_err[pt_bin][eta_bin] = reso[pt_bin][eta_bin]->GetRMSError();
      factor[pt_bin][eta_bin] = 1/mean[pt_bin][eta_bin];
    }
  }

  // now bin in eta and fit in pt direction
  TH1F * factor_binned[no_etabins];
  std::string title = "factor_binned_";
  std::string name;
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    name = title + std::to_string(eta_bin);
    factor_binned[eta_bin] = new TH1F(name.c_str(), "", 6, xbins);
    factor_binned[eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}");
    factor_binned[eta_bin]->GetYaxis()->SetTitle("correction factor");
 }
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
      factor_binned[eta_bin]->SetBinContent(pt_bin+1, factor[pt_bin][eta_bin]);
      factor_binned[eta_bin]->SetBinError(pt_bin+1, (log(factor[pt_bin][eta_bin]) * rms[pt_bin][eta_bin]));
    }
  }
  // fit and write function parameters in file
  TF1 *fit[no_etabins];
  ofstream correction;
  correction.open ("CorrectionFactors.txt");
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    factor_binned[eta_bin]->Fit("pol2");
    fit[eta_bin] = factor_binned[eta_bin]->GetFunction("pol2");
    correction << fixed << setprecision(10) << eta_bin << " " << fit[eta_bin]->GetParameter(0) << " " << fit[eta_bin]->GetParameter(1) << " " << fit[eta_bin]->GetParameter(2) << endl;
    cout << fixed << setprecision(10) << eta_bin << " " << (double)fit[eta_bin]->GetParameter(0) << " " << (double)fit[eta_bin]->GetParameter(1) << " " << (double)fit[eta_bin]->GetParameter(2) << " " << endl;
  }

  // write values from 2D Plot in table 
  ofstream table;
  table.open ("CorrectionFactors_table.txt");
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      table << pt_bin << " " << eta_bin << " " << factor[pt_bin][eta_bin] << endl;
      // cout << pt_bin << " " << eta_bin << " " << factor << endl;
    }
  }


  TAxis *xaxis1 = resolution_mean->GetXaxis();
  TAxis *yaxis1 = resolution_mean->GetYaxis();
  TAxis *xaxis1_err = resolution_mean_err->GetXaxis();
  TAxis *yaxis1_err = resolution_mean_err->GetYaxis();
  TAxis *xaxis2 = resolution_rms->GetXaxis();
  TAxis *yaxis2 = resolution_rms->GetYaxis();
  TAxis *xaxis2_err = resolution_rms_err->GetXaxis();
  TAxis *yaxis2_err = resolution_rms_err->GetYaxis();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      resolution_mean->Fill(xaxis1->GetBinCenter(pt_bin+1), yaxis1->GetBinCenter(eta_bin+1), mean[pt_bin][eta_bin]);
      resolution_mean_err->Fill(xaxis1_err->GetBinCenter(pt_bin+1), yaxis1_err->GetBinCenter(eta_bin+1), mean_err[pt_bin][eta_bin]);
      resolution_rms->Fill(xaxis2->GetBinCenter(pt_bin+1), yaxis2->GetBinCenter(eta_bin+1), rms[pt_bin][eta_bin]);
      resolution_rms_err->Fill(xaxis2_err->GetBinCenter(pt_bin+1), yaxis2_err->GetBinCenter(eta_bin+1), rms_err[pt_bin][eta_bin]);
    }
  }
  // resolution_mean->Fill(xaxis1->GetBinCenter(0), yaxis1->GetBinCenter(3), mean[0][3]);

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  resolution_mean->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_mean->GetYaxis()->SetTitle("#eta^{rec}");
  resolution_mean->GetZaxis()->SetTitle("MEAN(p_{T}^{rec} / p_{T}^{gen})");
  resolution_mean->GetXaxis()->SetTitleSize(0.05);
  resolution_mean->GetYaxis()->SetTitleSize(0.05);
  resolution_mean->GetZaxis()->SetTitleSize(0.05);
  resolution_mean->GetXaxis()->SetTitleOffset(0.9);
  resolution_mean->GetYaxis()->SetTitleOffset(0.8);
  resolution_mean->GetZaxis()->SetTitleOffset(0.9);
  resolution_rms->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_rms->GetYaxis()->SetTitle("#eta^{rec}");
  resolution_rms->GetZaxis()->SetTitle("RMS(p_{T}^{rec} / p_{T}^{gen})");
  resolution_rms->GetXaxis()->SetTitleSize(0.05);
  resolution_rms->GetYaxis()->SetTitleSize(0.05);
  resolution_rms->GetZaxis()->SetTitleSize(0.05);
  resolution_rms->GetXaxis()->SetTitleOffset(0.9);
  resolution_rms->GetYaxis()->SetTitleOffset(0.8);
  resolution_rms->GetZaxis()->SetTitleOffset(0.9);
  resolution_mean_err->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_mean_err->GetYaxis()->SetTitle("#eta^{rec}");
  resolution_mean_err->GetZaxis()->SetTitle("MEAN ERR(p_{T}^{rec} / p_{T}^{gen})");
  resolution_mean_err->GetXaxis()->SetTitleSize(0.05);
  resolution_mean_err->GetYaxis()->SetTitleSize(0.05);
  resolution_mean_err->GetZaxis()->SetTitleSize(0.05);
  resolution_mean_err->GetXaxis()->SetTitleOffset(0.9);
  resolution_mean_err->GetYaxis()->SetTitleOffset(0.8);
  resolution_mean_err->GetZaxis()->SetTitleOffset(0.9);
  resolution_rms_err->GetXaxis()->SetTitle("p_{T}^{rec}");
  resolution_rms_err->GetYaxis()->SetTitle("#eta^{rec}");
  resolution_rms_err->GetZaxis()->SetTitle("RMS ERR(p_{T}^{rec} / p_{T}^{gen})");
  resolution_rms_err->GetXaxis()->SetTitleSize(0.05);
  resolution_rms_err->GetYaxis()->SetTitleSize(0.05);
  resolution_rms_err->GetZaxis()->SetTitleSize(0.05);
  resolution_rms_err->GetXaxis()->SetTitleOffset(0.9);
  resolution_rms_err->GetYaxis()->SetTitleOffset(0.8);
  resolution_rms_err->GetZaxis()->SetTitleOffset(0.9);



  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
 
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TCanvas *A = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  resolution_mean->Draw("COLZ");
  resolution_mean->Draw("text:same");
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Mean_numbers.pdf"); 

  TCanvas *B = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  resolution_rms->Draw("COLZ");
  resolution_rms->Draw("text:same");
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/RMS_numbers.pdf"); 

  TCanvas *C = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  resolution_mean_err->Draw("COLZ");
  resolution_mean_err->Draw("text:same");
  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Mean_Error_numbers.pdf"); 

  TCanvas *D = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  resolution_rms_err->Draw("COLZ");
  resolution_rms_err->Draw("text:same");
  D->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/RMS_Error_numbers.pdf"); 

  TCanvas *E = new TCanvas();
  E->Divide(3,4);
  E->SetCanvasSize(800, 1200);
  E->SetWindowSize(800, 1200);
  TLegend *leg[12];
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    E->cd(eta_bin+1);
    gPad->SetLeftMargin(0.15);  
    leg[eta_bin] = new TLegend(0.4,0.25,0.85,0.4);
    leg[eta_bin]->AddEntry(factor_binned[eta_bin],"correction factor","pl");
    leg[eta_bin]->AddEntry(fit[eta_bin],"fit","l");
    std::stringstream title;
    title << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1];
    factor_binned[eta_bin]->SetTitle(title.str().c_str());
    factor_binned[eta_bin]->GetYaxis()->SetTitleSize(0.06);
    factor_binned[eta_bin]->GetXaxis()->SetTitleSize(0.05);
    factor_binned[eta_bin]->GetXaxis()->SetTitleOffset(0.9);
    factor_binned[eta_bin]->GetYaxis()->SetTitleOffset(1.1);
    factor_binned[eta_bin]->GetXaxis()->SetNdivisions(505);
    factor_binned[eta_bin]->GetYaxis()->SetNdivisions(505);
    factor_binned[eta_bin]->SetMarkerStyle(20);
    factor_binned[eta_bin]->SetMarkerSize(0.8);
    factor_binned[eta_bin]->SetLineColor(1);
    factor_binned[eta_bin]->Draw("E1");
    leg[eta_bin]->Draw();
    gPad->RedrawAxis();
  }
  E->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Fits.pdf"); 

  TCanvas *F = new TCanvas();
  F->SetCanvasSize(600, 600);
  gPad->SetLeftMargin(0.15);  
  factor_binned[7]->Draw("E1");
  leg[7]->Draw();
  gPad->RedrawAxis();
  F->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Fits_example.pdf"); 
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
