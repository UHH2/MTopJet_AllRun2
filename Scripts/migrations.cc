void migrations()
{

  // declare files
  TFile *TTbar = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_noBTag.root");
 
  TH1F * pt = (TH1F*)TTbar->Get("h_XCone_cor_migration_pt/number_hadjet");
  TH1F * pt350 = (TH1F*)TTbar->Get("h_XCone_cor_migration_pt350/number_hadjet");
  TH1F * mass = (TH1F*)TTbar->Get("h_XCone_cor_migration_mass/number_hadjet");
  TH1F * btag = (TH1F*)TTbar->Get("h_XCone_cor_migration_btag/number_hadjet");

  double pt_number = pt->Integral();
  double pt350_number = pt350->Integral();
  double mass_number = mass->Integral();
  double btag_number = btag->Integral();

  TH1F *h1 = new TH1F("h1", " ", 4, 1, 5);
  h1->Fill(1, pt_number);
  h1->Fill(2, pt350_number);
  h1->Fill(3, mass_number);
  h1->Fill(4, btag_number);

  h1->GetXaxis()->SetLabelSize(.06);
  h1->GetXaxis()->SetBinLabel(1,"no pt cut");
  h1->GetXaxis()->SetBinLabel(2,"pt > 350");
  h1->GetXaxis()->SetBinLabel(3,"no mass cut");
  h1->GetXaxis()->SetBinLabel(4,"no b-tag");

  h1->GetYaxis()->SetTitle("events");
  h1->GetYaxis()->SetTitleOffset(1.7);
  h1->GetYaxis()->SetNdivisions(505);


  h1->SetLineColor(kAzure+7);
  h1->SetLineWidth(3);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  h1->Draw("HIST");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Migrations_Sidebands.pdf");


}

