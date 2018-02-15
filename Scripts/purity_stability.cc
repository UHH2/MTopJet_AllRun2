void purity_stability()
{
  // declare files
  TFile *All = new TFile("uhh2.AnalysisModuleRunner.MC.Background.root");

  // get hist
  TH1F * rec = (TH1F*)All->Get("XCone_GEN_RecOnly/Mass_HadJet33");
  TH1F * gen = (TH1F*)All->Get("XCone_GEN_GenOnly/Mass_HadJet33");
  TH1F * both = (TH1F*)All->Get("XCone_GEN_Both/Mass_HadJet33");

  // new hist
  Float_t bins[] = {0, 100, 150, 200, 250, 300, 500};

  TH1F * rec_binned = new TH1F("rec_binned", "rec_binned",6, bins);
  TH1F * gen_binned = new TH1F("gen_binned", "gen_binned",6, bins);
  TH1F * both_binned = new TH1F("both_binned", "both_binned",6, bins);
  TH1F * purity = new TH1F("purity", "purity",6, bins);
  TH1F * stability = new TH1F("stability", "stability",6, bins);

  // get bin boundries for rebinning
  int bin[] = {rec->GetXaxis()->FindBin(0.0), rec->GetXaxis()->FindBin(100), rec->GetXaxis()->FindBin(150), rec->GetXaxis()->FindBin(200), rec->GetXaxis()->FindBin(250), rec->GetXaxis()->FindBin(300), rec->GetXaxis()->FindBin(500)};

  double rec_val[6];
  double rec_err[6];
  double gen_val[6];
  double gen_err[6];
  double both_val[6];
  double both_err[6];

  for(unsigned int i = 0; i < 6; ++i){
    double val_r = 0;
    double val_g = 0;
    double val_b = 0;
    double err_r = 0;
    double err_g = 0;
    double err_b = 0;
    for(unsigned int j=bin[i]; j < bin[i+1]; ++j){
      val_r += rec->GetBinContent(j);
      err_r += (rec->GetBinError(j))*(rec->GetBinError(j));
      val_g += gen->GetBinContent(j);
      err_g += (gen->GetBinError(j))*(gen->GetBinError(j));
      val_b += both->GetBinContent(j);
      err_b += (both->GetBinError(j))*(both->GetBinError(j));
    }
    rec_val[i] = val_r;
    rec_err[i] = TMath::Sqrt(err_r);
    gen_val[i] = val_g;
    gen_err[i] = TMath::Sqrt(err_g);
    both_val[i] = val_b;
    both_err[i] = TMath::Sqrt(err_b);

  }

  for(unsigned int i = 0; i < 6; ++i){
    rec_binned->Fill(rec_binned->GetXaxis()->GetBinCenter(i+1), rec_val[i]);
    rec_binned->SetBinError(i+1, rec_err[i]);
    gen_binned->Fill(gen_binned->GetXaxis()->GetBinCenter(i+1), gen_val[i]);
    gen_binned->SetBinError(i+1, gen_err[i]);
    both_binned->Fill(both_binned->GetXaxis()->GetBinCenter(i+1), both_val[i]);
    both_binned->SetBinError(i+1, both_err[i]);
  }

  purity->Divide(both_binned, rec_binned, 1., 1., "B");
  stability->Divide(both_binned, gen_binned, 1, 1, "B");


  //
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kTRUE);

  purity->GetYaxis()->SetRangeUser(0, 1);
  purity->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  purity->GetYaxis()->SetTitle("N_{rec+gen} / N_{rec}");
  purity->GetYaxis()->SetTitleOffset(1.4);
  purity->SetLineWidth(2);

  stability->GetYaxis()->SetRangeUser(0, 1);
  stability->GetXaxis()->SetTitle("M_{jet1}^{gen}");
  stability->GetYaxis()->SetTitle("N_{rec+gen} / N_{gen}");
  stability->GetYaxis()->SetTitleOffset(1.4);
  stability->SetLineWidth(2);

  TCanvas *C = new TCanvas();
  C->Divide(3,2);
  C->cd(1);
  rec_binned->Draw();
  C->cd(2);
  both_binned->Draw();
  C->cd(3);
  gen_binned->Draw();
  C->cd(4);
  rec->Draw();
  C->cd(5);
  both->Draw();
  C->cd(6);
  gen->Draw();

  TCanvas *A = new TCanvas();
  purity->Draw();

  TCanvas *B = new TCanvas();
  stability->Draw();


}
