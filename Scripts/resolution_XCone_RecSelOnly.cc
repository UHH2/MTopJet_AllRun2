void resolution_XCone_RecSelOnly(bool mass, bool genbins)
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
  TH1F * reso_corrected[7];
  TH1F * reso_pt_corrected[7];

  if(mass){
    reso[0] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution");
    reso[1] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass1");
    reso[2] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass2");
    reso[3] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass3");
    reso[4] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass4");
    reso[5] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass5");
    reso[6] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_mass6");

    reso_pt[0] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution");
    reso_pt[1] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass1");
    reso_pt[2] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass2");
    reso_pt[3] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass3");
    reso_pt[4] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass4");
    reso_pt[5] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass5");
    reso_pt[6] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_mass6");

    reso_noJEC[0] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution");
    reso_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass1");
    reso_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass2");
    reso_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass3");
    reso_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass4");
    reso_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass5");
    reso_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_mass6");

    reso_pt_noJEC[0] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution");
    reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass1");
    reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass2");
    reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass3");
    reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass4");
    reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass5");
    reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_mass6");

    reso_corrected[0] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution");
    reso_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass1");
    reso_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass2");
    reso_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass3");
    reso_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass4");
    reso_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass5");
    reso_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_mass6");

    reso_pt_corrected[0] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution");
    reso_pt_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass1");
    reso_pt_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass2");
    reso_pt_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass3");
    reso_pt_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass4");
    reso_pt_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass5");
    reso_pt_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_mass6");
  }
  else{
    reso[0] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution");
    reso_pt[0] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution");
    reso_noJEC[0] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution");
    reso_pt_noJEC[0] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution");
    reso_corrected[0] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution");
    reso_pt_corrected[0] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution");

    if(genbins){
      reso[1] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt1");
      reso[2] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt2");
      reso[3] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt3");
      reso[4] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt4");
      reso[5] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt5");
      reso[6] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt6");

      reso_pt[1] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt1");
      reso_pt[2] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt2");
      reso_pt[3] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt3");
      reso_pt[4] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt4");
      reso_pt[5] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt5");
      reso_pt[6] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt6");

      reso_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt1");
      reso_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt2");
      reso_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt3");
      reso_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt4");
      reso_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt5");
      reso_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt6");

      reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt1");
      reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt2");
      reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt3");
      reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt4");
      reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt5");
      reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt6");

      reso_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt1");
      reso_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt2");
      reso_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt3");
      reso_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt4");
      reso_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt5");
      reso_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt6");

      reso_pt_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt1");
      reso_pt_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt2");
      reso_pt_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt3");
      reso_pt_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt4");
      reso_pt_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt5");
      reso_pt_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt6");
    }
    else{
      reso[1] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt1_rec");
      reso[2] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt2_rec");
      reso[3] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt3_rec");
      reso[4] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt4_rec");
      reso[5] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt5_rec");
      reso[6] = (TH1F*)All->Get("RecGenHists_RecOnly/MassResolution_pt6_rec");

      reso_pt[1] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt1_rec");
      reso_pt[2] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt2_rec");
      reso_pt[3] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt3_rec");
      reso_pt[4] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt4_rec");
      reso_pt[5] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt5_rec");
      reso_pt[6] = (TH1F*)All->Get("RecGenHists_RecOnly/PtResolution_pt6_rec");

      reso_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt1_rec");
      reso_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt2_rec");
      reso_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt3_rec");
      reso_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt4_rec");
      reso_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt5_rec");
      reso_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/MassResolution_pt6_rec");

      reso_pt_noJEC[1] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt1_rec");
      reso_pt_noJEC[2] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt2_rec");
      reso_pt_noJEC[3] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt3_rec");
      reso_pt_noJEC[4] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt4_rec");
      reso_pt_noJEC[5] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt5_rec");
      reso_pt_noJEC[6] = (TH1F*)All->Get("RecGenHists_RecOnly_noJEC/PtResolution_pt6_rec");

      reso_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt1_rec");
      reso_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt2_rec");
      reso_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt3_rec");
      reso_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt4_rec");
      reso_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt5_rec");
      reso_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/MassResolution_pt6_rec");

      reso_pt_corrected[1] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt1_rec");
      reso_pt_corrected[2] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt2_rec");
      reso_pt_corrected[3] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt3_rec");
      reso_pt_corrected[4] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt4_rec");
      reso_pt_corrected[5] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt5_rec");
      reso_pt_corrected[6] = (TH1F*)All->Get("RecGenHists_RecOnly_corrected/PtResolution_pt6_rec");
    }

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
  TH1F * resolution1 = new TH1F("resolution1"," ",6, bins);
  TH1F * resolution2 = new TH1F("resolution2"," ",6, bins);
  TH1F * resolution_pt1 = new TH1F("resolution_pt1"," ",6, bins);
  TH1F * resolution_pt2 = new TH1F("resolution_pt2"," ",6, bins);
  TH1F * resolution1_noJEC = new TH1F("resolution1_noJEC"," ",6, bins);
  TH1F * resolution2_noJEC = new TH1F("resolution2_noJEC"," ",6, bins);
  TH1F * resolution_pt1_noJEC = new TH1F("resolution_pt1_noJEC"," ",6, bins);
  TH1F * resolution_pt2_noJEC = new TH1F("resolution_pt2_noJEC"," ",6, bins);
  TH1F * resolution1_corrected = new TH1F("resolution1_corrected"," ",6, bins);
  TH1F * resolution2_corrected = new TH1F("resolution2_corrected"," ",6, bins);
  TH1F * resolution_pt1_corrected = new TH1F("resolution_pt1_corrected"," ",6, bins);
  TH1F * resolution_pt2_corrected = new TH1F("resolution_pt2_corrected"," ",6, bins);
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
  //------------------- mass reso (corrected) -----------------------------
  //-----------------------------------------------------------------------
  double mean_corrected[7], rms_corrected[7], mean_error_corrected[7], rms_error_corrected[7];
  TF1 *fit_corrected[7];

  for(unsigned int i=1; i<=6; i++){
    reso_corrected[i]->Fit("gaus");
    fit_corrected[i] = reso_corrected[i]->GetFunction("gaus");
    mean_corrected[i] = fit_corrected[i]->GetParameter(1);
    rms_corrected[i] = fit_corrected[i]->GetParameter(2);
    mean_error_corrected[i] = fit_corrected[i]->GetParError(1);
    rms_error_corrected[i] = fit_corrected[i]->GetParError(2);
  }

  TAxis *xaxis1_corrected = resolution1_corrected->GetXaxis();
  TAxis *xaxis2_corrected = resolution2_corrected->GetXaxis();

  for(unsigned int i=1; i<=6; i++){
    resolution1_corrected->Fill(xaxis1_corrected->GetBinCenter(i), mean_corrected[i]);
    resolution1_corrected->SetBinError(i, mean_error_corrected[i]);
    resolution2_corrected->Fill(xaxis2_corrected->GetBinCenter(i), rms_corrected[i]);
    resolution2_corrected->SetBinError(i, rms_error_corrected[i]);
  }
  //-----------------------------------------------------------------------
  //------------------- pt reso (corrected) -------------------------------
  //-----------------------------------------------------------------------
  double mean_pt_corrected[7], rms_pt_corrected[7], mean_error_pt_corrected[7], rms_error_pt_corrected[7];
  TF1 *fit_pt_corrected[7];

  for(unsigned int i=1; i<=6; i++){
    reso_pt_corrected[i]->Fit("gaus");
    fit_pt_corrected[i] = reso_pt_corrected[i]->GetFunction("gaus");
    mean_pt_corrected[i] = fit_pt_corrected[i]->GetParameter(1);
    rms_pt_corrected[i] = fit_pt_corrected[i]->GetParameter(2);
    mean_error_pt_corrected[i] = fit_pt_corrected[i]->GetParError(1);
    rms_error_pt_corrected[i] = fit_pt_corrected[i]->GetParError(2);
  }

  TAxis *xaxis_pt1_corrected = resolution_pt1_corrected->GetXaxis();
  TAxis *xaxis_pt2_corrected = resolution_pt2_corrected->GetXaxis();

  for(unsigned int i=1; i<=6; i++){
    resolution_pt1_corrected->Fill(xaxis_pt1_corrected->GetBinCenter(i), mean_pt_corrected[i]);
    resolution_pt1_corrected->SetBinError(i, mean_error_pt_corrected[i]);
    resolution_pt2_corrected->Fill(xaxis_pt2_corrected->GetBinCenter(i), rms_pt_corrected[i]);
    resolution_pt2_corrected->SetBinError(i, rms_error_pt_corrected[i]);
  }

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------




  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  resolution1->GetYaxis()->SetRangeUser(-0.2, 0.2);
  resolution1->GetXaxis()->SetTitleOffset(1.4);
  resolution1->GetYaxis()->SetTitleOffset(1.4);
  resolution1->GetXaxis()->SetNdivisions(505);
  resolution1->GetYaxis()->SetNdivisions(505);
  if(mass)resolution1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else{
    if(genbins) resolution1->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
    else resolution1->GetXaxis()->SetTitle("p_{T,jet1}^{rec}");
  }
  resolution1->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution1->SetLineWidth(2);
  resolution1->SetLineColor(kRed+1);

  resolution2->GetYaxis()->SetRangeUser(0.0, 0.6);
  resolution2->GetXaxis()->SetTitleOffset(1.4);
  resolution2->GetYaxis()->SetTitleOffset(1.4);
  resolution2->GetXaxis()->SetNdivisions(505);
  resolution2->GetYaxis()->SetNdivisions(505);
  if(mass)resolution2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else{
    if(genbins) resolution2->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
    else resolution2->GetXaxis()->SetTitle("p_{T,jet1}^{rec}");
  }
  resolution2->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2->GetYaxis()->SetTitleOffset(1.1);
  resolution2->SetLineWidth(2);
  resolution2->SetLineColor(kRed+1);

  resolution2_corrected->GetYaxis()->SetRangeUser(0.0, 0.6);
  resolution2_corrected->GetXaxis()->SetTitleOffset(1.4);
  resolution2_corrected->GetYaxis()->SetTitleOffset(1.4);
  resolution2_corrected->GetXaxis()->SetNdivisions(505);
  resolution2_corrected->GetYaxis()->SetNdivisions(505);
  if(mass)resolution2_corrected->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else{
    if(genbins) resolution2_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
    else resolution2_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{rec}");
  }
  resolution2_corrected->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2_corrected->GetYaxis()->SetTitleOffset(1.1);
  resolution2_corrected->SetLineWidth(2);
  resolution2_corrected->SetLineColor(kRed+1);

  resolution_pt1->GetYaxis()->SetRangeUser(-0.1, 0.1);
  resolution_pt1->GetXaxis()->SetTitleOffset(1.4);
  resolution_pt1->GetYaxis()->SetTitleOffset(1.4);
  resolution_pt1->GetXaxis()->SetNdivisions(505);
  resolution_pt1->GetYaxis()->SetNdivisions(505);
  if(mass)resolution_pt1->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else{
    if(genbins) resolution_pt1->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
    else resolution_pt1->GetXaxis()->SetTitle("p_{T,jet1}^{rec}");
  }
  resolution_pt1->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt1->GetYaxis()->SetTitleOffset(1.2);
  resolution_pt1->SetLineWidth(2);
  resolution_pt1->SetLineColor(kRed+1);

  resolution_pt2->GetYaxis()->SetRangeUser(0.0, 0.1);
  resolution_pt2->GetXaxis()->SetTitleOffset(1.4);
  resolution_pt2->GetYaxis()->SetTitleOffset(1.4);
  resolution_pt2->GetXaxis()->SetNdivisions(505);
  resolution_pt2->GetYaxis()->SetNdivisions(505);
  if(mass)resolution_pt2->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  else{
    if(genbins) resolution_pt2->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
    else resolution_pt2->GetXaxis()->SetTitle("p_{T,jet1}^{rec}");
  }
  resolution_pt2->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2->GetYaxis()->SetTitleOffset(1.2);
  resolution_pt2->SetLineWidth(3);
  resolution_pt2->SetLineColor(kRed+1);

  // resolution1_noJEC->GetYaxis()->SetRangeUser(-0.2, 0.2);
  // if(mass)resolution1_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution1_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution1_noJEC->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution1_noJEC->SetLineWidth(3);
  resolution1_noJEC->SetLineColor(kOrange+1);

  // resolution2_noJEC->GetYaxis()->SetRangeUser(0.0, 0.3);
  // if(mass)resolution2_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution2_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution2_noJEC->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2_noJEC->SetLineWidth(3);
  resolution2_noJEC->SetLineColor(kOrange+1);

  // resolution_pt1_noJEC->GetYaxis()->SetRangeUser(-0.1, 0.1);
  // if(mass)resolution_pt1_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution_pt1_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution_pt1_noJEC->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt1_noJEC->SetLineWidth(3);
  resolution_pt1_noJEC->SetLineColor(kOrange+1);

  // resolution_pt2_noJEC->GetYaxis()->SetRangeUser(0.0, 0.1);
  // if(mass)resolution_pt2_noJEC->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution_pt2_noJEC->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution_pt2_noJEC->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2_noJEC->SetLineWidth(3);
  resolution_pt2_noJEC->SetLineColor(kOrange+1);

  // resolution1_corrected->GetYaxis()->SetRangeUser(-0.2, 0.2);
  // if(mass)resolution1_corrected->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution1_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution1_corrected->GetYaxis()->SetTitle("Mean((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution1_corrected->SetLineWidth(3);
  resolution1_corrected->SetLineColor(kAzure+7);

  // resolution2_corrected->GetYaxis()->SetRangeUser(0.0, 0.3);
  // if(mass)resolution2_corrected->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution2_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution2_corrected->GetYaxis()->SetTitle("RMS((M_{jet1}^{rec} - M_{jet1}^{gen}) /M_{jet1}^{gen})");
  resolution2_corrected->SetLineWidth(3);
  resolution2_corrected->SetLineColor(kAzure+7);

  // resolution_pt1_corrected->GetYaxis()->SetRangeUser(-0.1, 0.1);
  // if(mass)resolution_pt1_corrected->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution_pt1_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution_pt1_corrected->GetYaxis()->SetTitle("Mean((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt1_corrected->SetLineWidth(3);
  resolution_pt1_corrected->SetLineColor(kAzure+7);

  // resolution_pt2_corrected->GetYaxis()->SetRangeUser(0.0, 0.1);
  // if(mass)resolution_pt2_corrected->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  // else resolution_pt2_corrected->GetXaxis()->SetTitle("p_{T,jet1}^{gen}");
  // resolution_pt2_corrected->GetYaxis()->SetTitle("RMS((p_{T,jet1}^{rec} - p_{T,jet1}^{gen}) /p_{T,jet1}^{gen})");
  resolution_pt2_corrected->SetLineWidth(3);
  resolution_pt2_corrected->SetLineColor(kAzure+7);


  zero_line = new TLine(0, 0, 2000, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);

  //-----------------------------------------------------------------------
  //------------------- DRAW ----------------------------------------------
  //-----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas();
  Float_t left = 0.1;
  Float_t right = 0.08;
  Float_t bottom = 0.12;
  Float_t top = 0.1;
  gPad->SetMargin (left, right, bottom, top);
  resolution1->Draw("E");
  zero_line->Draw("SAME");
  resolution1_noJEC->Draw("SAME E");
  resolution1_corrected->Draw("SAME E");
  leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(resolution1,"JEC applied","l");
  leg->AddEntry(resolution1_noJEC,"no JEC applied","l");
  leg->AddEntry(resolution1_corrected,"corrected","l");
  leg->Draw();
  gPad->RedrawAxis();
  if(mass)c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_mean_massbin.pdf"); 
  else{
    if(genbins) c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_mean_ptgenbin.pdf"); 
    else c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_mean_ptrecbin.pdf");
  }
  TCanvas *c2 = new TCanvas();
  gPad->SetMargin (left, right, bottom, top);
  resolution_pt1->Draw("E");
  resolution_pt1_noJEC->Draw("SAME E");
  resolution_pt1_corrected->Draw("SAME E");
  leg2 = new TLegend(0.65,0.85,0.85,0.65);
  leg2->AddEntry(resolution_pt1,"JEC applied","l");
  leg2->AddEntry(resolution_pt1_noJEC,"no JEC applied","l");
  leg2->AddEntry(resolution_pt1_corrected,"corrected","l");
  leg2->Draw();
  gPad->RedrawAxis();
  if(mass)c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_mean_massbin.pdf"); 
  else{
    if(genbins) c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_mean_ptgenbin.pdf"); 
    else c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_mean_ptrecbin.pdf");
  }
  TCanvas *c3 = new TCanvas();
  gPad->SetMargin (left, right, bottom, top);
  // resolution2->Draw("E");
  // resolution2_noJEC->Draw("SAME E");
  resolution2_corrected->Draw("SAME E");
  leg3 = new TLegend(0.65,0.85,0.85,0.65);
  // leg3->AddEntry(resolution2,"JEC applied","l");
  // leg3->AddEntry(resolution2_noJEC,"no JEC applied","l");
  leg3->AddEntry(resolution2_corrected," corrected","l");
  leg3->Draw();
  gPad->RedrawAxis();
  if(mass)c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_rms_massbin.pdf"); 
  else{
    if(genbins) c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_rms_ptgenbin.pdf"); 
    else c3->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/mass_rms_ptrecbin.pdf");
  }
  TCanvas *c4 = new TCanvas();
  gPad->SetMargin (left, right, bottom, top);
  resolution_pt2->Draw("E");
  resolution_pt2_noJEC->Draw("SAME E");
  resolution_pt2_corrected->Draw("SAME E");
  leg4 = new TLegend(0.65,0.85,0.85,0.65);
  leg4->AddEntry(resolution_pt2,"JEC applied","l");
  leg4->AddEntry(resolution_pt2_noJEC,"no JEC applied","l");
  leg4->AddEntry(resolution_pt2_corrected,"corrected","l");
  leg4->Draw();
  gPad->RedrawAxis();
  if(mass)c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_rms_massbin.pdf"); 
  else{
    if(genbins) c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_rms_ptgenbin.pdf"); 
    else c4->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Resolution_XCone_RecSel/pt_rms_ptrecbin.pdf");
  }

}

