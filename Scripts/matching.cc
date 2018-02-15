void matching()
{
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TTbar_top = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get mass plots -------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TH1F * JetMass = (TH1F*)TTbar_top->Get("XCone_cor/M_jet1");
  TH1F * JetMass_matched = (TH1F*)TTbar_top->Get("XCone_cor_matched/M_jet1");
  TH1F * JetMass_unmatched = (TH1F*)TTbar_top->Get("XCone_cor_unmatched/M_jet1");
  TH1F * JetMass_matched_fat = (TH1F*)TTbar_top->Get("XCone_cor_matched_fat/M_jet1");
  TH1F * JetMass_unmatched_fat = (TH1F*)TTbar_top->Get("XCone_cor_unmatched_fat/M_jet1");

  TH1F *h1 = JetMass;
  TH1F *h1_m = JetMass_matched;
  TH1F *h1_u = JetMass_unmatched;
  TH1F *h1_m_fat = JetMass_matched_fat;
  TH1F *h1_u_fat = JetMass_unmatched_fat;

  float events_tot, events_matched, fraction, events_matched_fat, fraction_fat;
  events_tot = h1->Integral();
  events_matched = h1_m->Integral();
  events_matched_fat = h1_m_fat->Integral();
  fraction = events_matched/events_tot;
  fraction_fat = events_matched_fat/events_tot;
  cout << endl;
  cout << "total events:                 " << events_tot << endl;
  cout << "events matched to subjet:     " << events_matched << endl;
  cout << "events matched to fatjet:     " << events_matched_fat << endl;
  cout << "------------------------------------------" << endl;
  cout << "fraction matched to subjet:   " << fraction << endl;
  cout << "fraction matched to fatjet:   " << fraction_fat << endl;
  cout << endl;
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  int xlow, xhigh;
  int ylow, yhigh;
  xlow = 0;
  xhigh = 450;
  ylow = 0;
  yhigh = 260;

  h1->SetTitle("XCone matching");
  h1->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  h1->GetYaxis()->SetTitle("events");
  h1->GetYaxis()->SetTitleOffset(1.7);
  h1->GetYaxis()->SetRangeUser(ylow, (h1->GetMaximum())*1.1);
  h1->GetXaxis()->SetNdivisions(505);
  h1->GetYaxis()->SetNdivisions(505);
  // h1->GetYaxis()->SetRangeUser(ylow, yhigh);
  h1->GetXaxis()->SetRangeUser(xlow, xhigh);

  h1->SetLineColor(1);
  h1->SetFillColor(kGray);
  h1_m->SetLineColor(798);
  h1_m->SetLineWidth(4);
  // h1_m->SetLineStyle(9);
  // h1_u->SetLineColor(4);
  // h1_u->SetLineWidth(4);
  h1_m_fat->SetLineColor(kAzure+7);
  h1_m_fat->SetLineWidth(3);
  h1_u_fat->SetLineColor(2);
  h1_u_fat->SetLineWidth(3);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
 

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  leg = new TLegend(0.55,0.65,0.85,0.8);
  leg->AddEntry(h1,"all","f");
  // leg->AddEntry(h1_u,"not matched to subjet","l");
  leg->AddEntry(h1_m_fat,"matched","l");
  leg->AddEntry(h1_m,"matched to subjet","l");
  leg->AddEntry(h1_u_fat,"not matched","l");
  h1->Draw("HIST");
  h1_m_fat->Draw("HIST SAME");
  h1_m->Draw("HIST SAME");
  h1_u_fat->Draw("HIST SAME");
  leg->Draw("");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/XCone_RECO_matching.pdf");
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
