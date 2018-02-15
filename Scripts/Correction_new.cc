void Correction_new()
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

  TH2F * h_ratio_mean = new TH2F("h_ratio_mean","Mean",6, xbins, 12, ybins);
  TH2F * h_ratio_mean_err = new TH2F("h_ratio_mean_err","Mean Error",6, xbins, 12, ybins);

  TH2F * h_ptrec_mean = new TH2F("h_ptrec_mean","Mean",6, xbins, 12, ybins);
  TH2F * h_ptrec_mean_err = new TH2F("h_ptrec_mean_err","Mean Error",6, xbins, 12, ybins);
  
  TH2F * h_count = new TH2F("h_count","Event Count",6, xbins, 12, ybins);


  int no_ptbins = 6; // rows
  int no_etabins = 12; // columns

  // get histograms from root file
  std::string dir = "CorrectionHists/";
  std::string name_ratio = "PtReso_";
  std::string name_ptrec = "PtRec_";
  std::string name_count = "Count_";


  std::string filename_ratio, filename_ptrec, filename_count;
  TH1F * ratio[no_ptbins][no_etabins];
  TH1F * ptrec[no_ptbins][no_etabins];
  TH1F * count[no_ptbins][no_etabins];

  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      filename_ratio = dir + name_ratio + std::to_string(pt_bin) + std::to_string(eta_bin) ;
      filename_ptrec = dir + name_ptrec + std::to_string(pt_bin) + std::to_string(eta_bin) ;
      filename_count = dir + name_count + std::to_string(pt_bin) + std::to_string(eta_bin) ;

      ratio[pt_bin][eta_bin] = (TH1F*)file->Get(filename_ratio.c_str());
      ptrec[pt_bin][eta_bin] = (TH1F*)file->Get(filename_ptrec.c_str());
      count[pt_bin][eta_bin] = (TH1F*)file->Get(filename_count.c_str());
    }
  }
  ////

  // get mean values and writ to array  
  double ratio_mean[no_ptbins][no_etabins], ratio_mean_err[no_ptbins][no_etabins];
  double ratio_rms[no_ptbins][no_etabins], ratio_rms_err[no_ptbins][no_etabins];
  double ptrec_mean[no_ptbins][no_etabins], ptrec_mean_err[no_ptbins][no_etabins];
  double ptrec_rms[no_ptbins][no_etabins], ptrec_rms_err[no_ptbins][no_etabins];
  double factor[no_ptbins][no_etabins], factor_err[no_ptbins][no_etabins];
  double n_events[no_ptbins][no_etabins];
  TF1 *ratio_fit[no_ptbins][no_etabins];
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){

      // get mean from ratio with gaussian fit
      double upper = ratio[pt_bin][eta_bin]->GetMean() + ratio[pt_bin][eta_bin]->GetRMS();
      double lower = ratio[pt_bin][eta_bin]->GetMean() - ratio[pt_bin][eta_bin]->GetRMS();
      TF1*ratio_func = new TF1("ratio_func","gaus",lower,upper);
      ratio[pt_bin][eta_bin]->Fit("ratio_func","R");
      ratio_fit[pt_bin][eta_bin] = ratio[pt_bin][eta_bin]->GetFunction("ratio_func");
      ratio_mean[pt_bin][eta_bin] = ratio_fit[pt_bin][eta_bin]->GetParameter(1);
      ratio_mean_err[pt_bin][eta_bin] = ratio_fit[pt_bin][eta_bin]->GetParError(1);
      
      // get correction factor
      factor[pt_bin][eta_bin] = 1/ratio_mean[pt_bin][eta_bin];
      factor_err[pt_bin][eta_bin] = abs(-ratio_mean_err[pt_bin][eta_bin]/(ratio_mean[pt_bin][eta_bin] * ratio_mean[pt_bin][eta_bin]));
      
      // get mean from ptrec with RootMean (because ptrec is only needed to get bin center in fit right)
      double upper_ = ptrec[pt_bin][eta_bin]->GetMean() + 2*ptrec[pt_bin][eta_bin]->GetRMS();
      double lower_ = ptrec[pt_bin][eta_bin]->GetMean() - 2*ptrec[pt_bin][eta_bin]->GetRMS();
      cout << "-------------------------------" << endl;
      cout << "Mean before = " << ptrec[pt_bin][eta_bin]->GetMean() << endl;
      cout << "lower = " << lower_ << endl;
      cout << "upper = " << upper_ << endl;
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetRangeUser(lower_, upper_);
      ptrec_mean[pt_bin][eta_bin] = ptrec[pt_bin][eta_bin]->GetMean();
      ptrec_mean_err[pt_bin][eta_bin] = ptrec[pt_bin][eta_bin]->GetMeanError();
      cout << "Mean after = " << ptrec_mean[pt_bin][eta_bin] << endl;

      // get event number in each bin
      n_events[pt_bin][eta_bin] = count[pt_bin][eta_bin]->Integral(); 
    }
  }
  ////

  // fill Hists with values from arrays
  TAxis *xaxis_ratio = h_ratio_mean->GetXaxis();
  TAxis *yaxis_ratio = h_ratio_mean->GetYaxis();
  TAxis *xaxis_ratio_err = h_ratio_mean_err->GetXaxis();
  TAxis *yaxis_ratio_err = h_ratio_mean_err->GetYaxis();
  TAxis *xaxis_ptrec = h_ptrec_mean->GetXaxis();
  TAxis *yaxis_ptrec = h_ptrec_mean->GetYaxis();
  TAxis *xaxis_ptrec_err = h_ptrec_mean_err->GetXaxis();
  TAxis *yaxis_ptrec_err = h_ptrec_mean_err->GetYaxis();
  TAxis *xaxis_count = h_count->GetXaxis();
  TAxis *yaxis_count = h_count->GetYaxis();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      h_ratio_mean->Fill(xaxis_ratio->GetBinCenter(pt_bin+1), yaxis_ratio->GetBinCenter(eta_bin+1), ratio_mean[pt_bin][eta_bin]);
      h_ratio_mean_err->Fill(xaxis_ratio->GetBinCenter(pt_bin+1), yaxis_ratio->GetBinCenter(eta_bin+1), ratio_mean_err[pt_bin][eta_bin]);
      h_ptrec_mean->Fill(xaxis_ptrec->GetBinCenter(pt_bin+1), yaxis_ptrec->GetBinCenter(eta_bin+1), ptrec_mean[pt_bin][eta_bin]);
      h_ptrec_mean_err->Fill(xaxis_ptrec->GetBinCenter(pt_bin+1), yaxis_ptrec->GetBinCenter(eta_bin+1), ptrec_mean_err[pt_bin][eta_bin]);
      h_count->Fill(xaxis_count->GetBinCenter(pt_bin+1), yaxis_count->GetBinCenter(eta_bin+1), n_events[pt_bin][eta_bin]);
    }
  }
  ////

  // now create graph Correction_factor vs <ptrec> in every eta bin
  Double_t g_factor[no_ptbins], g_factor_err[no_ptbins];
  Double_t g_ptrec[no_ptbins], g_ptrec_err[no_ptbins];
  TGraphErrors* factor_pt[no_etabins];
  TF1 * fit[no_etabins];
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
      g_factor[pt_bin] = factor[pt_bin][eta_bin];
      g_factor_err[pt_bin] = factor_err[pt_bin][eta_bin];
      g_ptrec[pt_bin] = ptrec_mean[pt_bin][eta_bin];
      g_ptrec_err[pt_bin] = ptrec_mean_err[pt_bin][eta_bin];
    }
    factor_pt[eta_bin] = new TGraphErrors(no_ptbins,g_ptrec,g_factor, g_ptrec_err, g_factor_err );
    TF1*f1 = new TF1("f1","pol2",0,500);
    f1->SetParLimits(0,0.1,1);
    f1->SetParLimits(1,0.0001,0.001);
    f1->SetParLimits(2,-0.000001,-0.0000001);

    factor_pt[eta_bin]->Fit("f1","R");
    fit[eta_bin] = factor_pt[eta_bin]->GetFunction("f1");
  }


  ofstream correction_file;
  correction_file.open ("CorrectionFactors_new.txt");
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    correction_file << fixed << setprecision(10) << eta_bin << " " 
		    << fit[eta_bin]->GetParameter(0) << " " 
		    << fit[eta_bin]->GetParameter(1) << " " 
		    << fit[eta_bin]->GetParameter(2) << endl;
  }
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
  h_ratio_mean->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ratio_mean->GetYaxis()->SetTitle("#eta^{rec}");
  h_ratio_mean->GetZaxis()->SetTitle("MEAN(p_{T}^{rec} / p_{T}^{gen})");
  h_ratio_mean->GetXaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetYaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetZaxis()->SetTitleSize(0.05);
  h_ratio_mean->GetXaxis()->SetTitleOffset(0.9);
  h_ratio_mean->GetYaxis()->SetTitleOffset(0.8);
  h_ratio_mean->GetZaxis()->SetTitleOffset(0.9);
  h_ratio_mean->Draw("COLZ");
  h_ratio_mean->Draw("text:same");
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Ratio_Mean.pdf"); 

  TCanvas *B = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  h_ptrec_mean->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ptrec_mean->GetYaxis()->SetTitle("#eta^{rec}");
  h_ptrec_mean->GetZaxis()->SetTitle("MEAN(p_{T}^{rec})");
  h_ptrec_mean->GetXaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetYaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetZaxis()->SetTitleSize(0.05);
  h_ptrec_mean->GetXaxis()->SetTitleOffset(0.9);
  h_ptrec_mean->GetYaxis()->SetTitleOffset(0.8);
  h_ptrec_mean->GetZaxis()->SetTitleOffset(0.9);
  h_ptrec_mean->Draw("COLZ");
  h_ptrec_mean->Draw("text:same");
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Ptrec_Mean.pdf"); 

  TCanvas *C = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  h_ratio_mean_err->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ratio_mean_err->GetYaxis()->SetTitle("#eta^{rec}");
  h_ratio_mean_err->GetZaxis()->SetTitle("MEAN ERR(p_{T}^{rec} / p_{T}^{gen})");
  h_ratio_mean_err->GetXaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetYaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetZaxis()->SetTitleSize(0.05);
  h_ratio_mean_err->GetXaxis()->SetTitleOffset(0.9);
  h_ratio_mean_err->GetYaxis()->SetTitleOffset(0.8);
  h_ratio_mean_err->GetZaxis()->SetTitleOffset(0.9);
  h_ratio_mean_err->Draw("COLZ");
  h_ratio_mean_err->Draw("text:same");
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Ratio_MEAN_ERR.pdf"); 

  TCanvas *D = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  h_ptrec_mean_err->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_ptrec_mean_err->GetYaxis()->SetTitle("#eta^{rec}");
  h_ptrec_mean_err->GetZaxis()->SetTitle("MEAN ERR(p_{T}^{rec})");
  h_ptrec_mean_err->GetXaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetYaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetZaxis()->SetTitleSize(0.05);
  h_ptrec_mean_err->GetXaxis()->SetTitleOffset(0.9);
  h_ptrec_mean_err->GetYaxis()->SetTitleOffset(0.8);
  h_ptrec_mean_err->GetZaxis()->SetTitleOffset(0.9);
  h_ptrec_mean_err->Draw("COLZ");
  h_ptrec_mean_err->Draw("text:same");
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Ptrec_MEAN_ERR.pdf"); 

  TCanvas *E = new TCanvas();
  E->Divide(3,4);
  E->SetCanvasSize(800, 1200);
  E->SetWindowSize(800, 1200);
  TLegend *leg[12];
  for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
    E->cd(eta_bin+1);
    gPad->SetLeftMargin(0.15);  
    leg[eta_bin] = new TLegend(0.4,0.25,0.85,0.4);
    leg[eta_bin]->AddEntry(factor_pt[eta_bin],"correction factor","pl");
    leg[eta_bin]->AddEntry(fit[eta_bin],"fit","l");
    std::stringstream title;
    title << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1];
    factor_pt[eta_bin]->SetTitle(title.str().c_str());
    factor_pt[eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}");
    factor_pt[eta_bin]->GetYaxis()->SetTitle("correction factor");
    factor_pt[eta_bin]->GetXaxis()->SetTitleSize(0.05);
    factor_pt[eta_bin]->GetXaxis()->SetTitleOffset(0.9);
    factor_pt[eta_bin]->GetXaxis()->SetNdivisions(505);
    factor_pt[eta_bin]->GetYaxis()->SetTitleSize(0.06);
    factor_pt[eta_bin]->GetYaxis()->SetTitleOffset(1.1);
    factor_pt[eta_bin]->GetYaxis()->SetNdivisions(505);
    factor_pt[eta_bin]->SetMarkerStyle(20);
    factor_pt[eta_bin]->SetMarkerSize(0.8);
    factor_pt[eta_bin]->SetLineColor(1);
    factor_pt[eta_bin]->Draw("AP");
    fit[eta_bin]->Draw("SAME");
    leg[eta_bin]->Draw();
    gPad->RedrawAxis();
  }
  E->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Fits.pdf"); 

  
  TCanvas *F = new TCanvas();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    F->Clear(); 
    F->Divide(3,4);
    F->SetCanvasSize(800, 1200);
    F->SetWindowSize(800, 1200);
    TLegend *leg_ratio[12];
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      F->cd(eta_bin+1);
      gPad->SetLeftMargin(0.15);  
      leg_ratio[eta_bin] = new TLegend(0.2,0.6,0.45,0.75);
      leg_ratio[eta_bin]->AddEntry(ratio[pt_bin][eta_bin],"p_{T}^{rec}/p_{T}^{gen}","f");
      leg_ratio[eta_bin]->AddEntry(ratio_fit[pt_bin][eta_bin],"fit","l");
      std::stringstream title2;
      title2 << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1] << " | " << xbins[pt_bin] << " < p_{T}^{gen} < " << xbins[pt_bin+1];
      ratio[pt_bin][eta_bin]->SetTitle(title2.str().c_str());
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}/p_{T}^{gen}");
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitle("events");
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitleSize(0.05);
      ratio[pt_bin][eta_bin]->GetXaxis()->SetTitleOffset(0.9);
      ratio[pt_bin][eta_bin]->GetXaxis()->SetNdivisions(505);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitleSize(0.06);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetTitleOffset(1.1);
      ratio[pt_bin][eta_bin]->GetYaxis()->SetNdivisions(505);
      ratio[pt_bin][eta_bin]->SetFillColor(kGray);
      ratio[pt_bin][eta_bin]->SetLineColor(1);
      ratio[pt_bin][eta_bin]->Draw("HIST");
      ratio_fit[pt_bin][eta_bin]->Draw("SAME");
      leg_ratio[eta_bin]->Draw();
      gPad->RedrawAxis();
    }
    std::string name = "/afs/desy.de/user/s/schwarzd/Plots/Correction/Fits_ratio_ptbin";
    std::string ending = std::to_string(pt_bin) + ".pdf";
    F->SaveAs((name+ending).c_str());
  }


  TCanvas *F2 = new TCanvas();
  for(int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    F2->Clear();
    F2->Divide(3,4);
    F2->SetCanvasSize(800, 1200);
    F2->SetWindowSize(800, 1200);
    TLegend *leg_ptrec[12];
    for(int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      F2->cd(eta_bin+1);
      gPad->SetLeftMargin(0.15);
      TLine *ptrec_line = new TLine(ptrec_mean[pt_bin][eta_bin], 0, ptrec_mean[pt_bin][eta_bin], ptrec[pt_bin][eta_bin]->GetMaximum());
      ptrec_line->SetLineColor(kRed);
      ptrec_line->SetLineWidth(3);
      ptrec_line->SetLineStyle(7);
      std::stringstream title3;
      title3 << ybins[eta_bin] << " < #eta < " << ybins[eta_bin+1] << " | " << xbins[pt_bin] << " < p_{T}^{gen} < " << xbins[pt_bin+1];
      ptrec[pt_bin][eta_bin]->SetTitle(title3.str().c_str());
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitle("p_{T}^{rec}");
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitle("events");
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitleSize(0.05);
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetTitleOffset(0.9);
      ptrec[pt_bin][eta_bin]->GetXaxis()->SetNdivisions(505);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitleSize(0.06);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetTitleOffset(1.1);
      ptrec[pt_bin][eta_bin]->GetYaxis()->SetNdivisions(505);
      ptrec[pt_bin][eta_bin]->SetFillColor(kGray);
      ptrec[pt_bin][eta_bin]->SetLineColor(1);
      ptrec[pt_bin][eta_bin]->Draw("HIST");
      ptrec_line->Draw("SAME");
      gPad->RedrawAxis();
    }
    std::string name = "/afs/desy.de/user/s/schwarzd/Plots/Correction/ptrec_mean_ptbin";
    std::string ending = std::to_string(pt_bin) + ".pdf";
    F2->SaveAs((name+ending).c_str());
  }

  TCanvas *G = new TCanvas("G", "G", 600, 600);
  gPad->SetLeftMargin(0.15);  
  factor_pt[7]->Draw("AP");
  fit[7]->Draw("SAME");
  leg[7]->Draw();
  gPad->RedrawAxis();
  G->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Fits_example.pdf"); 

  TCanvas *G2 = new TCanvas("G2", "G2", 600, 600);
  gPad->SetLeftMargin(0.15);  
  TLine *example_line = new TLine(ptrec_mean[2][3], 0, ptrec_mean[2][3], ptrec[2][3]->GetMaximum());
  example_line->SetLineColor(kRed);
  example_line->SetLineWidth(3);
  example_line->SetLineStyle(7);
  ptrec[2][3]->Draw("HIST");
  example_line->Draw("SAME");
  gPad->RedrawAxis();
  G2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/ptrec_mean_example.pdf"); 

  TCanvas *H = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);   
  h_count->GetXaxis()->SetTitle("p_{T}^{gen}");
  h_count->GetYaxis()->SetTitle("#eta^{rec}");
  h_count->GetZaxis()->SetTitle("Events");
  h_count->GetXaxis()->SetTitleSize(0.05);
  h_count->GetYaxis()->SetTitleSize(0.05);
  h_count->GetZaxis()->SetTitleSize(0.05);
  h_count->GetXaxis()->SetTitleOffset(0.9);
  h_count->GetYaxis()->SetTitleOffset(0.8);
  h_count->GetZaxis()->SetTitleOffset(0.9);
  h_count->Draw("COLZ");
  h_count->Draw("text:same");
  H->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Correction/Event_Count.pdf"); 
  
}
