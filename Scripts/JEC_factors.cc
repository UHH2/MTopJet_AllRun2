void JEC_factors()
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
  TH1F * L1 = (TH1F*)file->Get("XCone_subjets/JEC_L1_all_subjets");
  TH1F * L2L3 = (TH1F*)file->Get("XCone_subjets/JEC_L2L3_all_subjets");
  TH1F * all = (TH1F*)file->Get("XCone_subjets/JEC_all_subjets");
  TH1F * L1_ak4 = (TH1F*)file->Get("XCone_subjets/JEC_L1_ak4");
  TH1F * L2L3_ak4 = (TH1F*)file->Get("XCone_subjets/JEC_L2L3_ak4");
  TH1F * all_ak4 = (TH1F*)file->Get("XCone_subjets/JEC_ak4");

  TH1F * pt = (TH1F*)file->Get("XCone_subjets/pt_all_subjets");
  TH1F * eta = (TH1F*)file->Get("XCone_subjets/eta_all_subjets");
  TH1F * pt_had = (TH1F*)file->Get("XCone_subjets/pt_had_subjets");
  TH1F * eta_had = (TH1F*)file->Get("XCone_subjets/eta_had_subjets");
  TH1F * pt_lep = (TH1F*)file->Get("XCone_subjets/pt_lep_subjets");
  TH1F * eta_lep = (TH1F*)file->Get("XCone_subjets/eta_lep_subjets");
  TH1F * pt_ak4 = (TH1F*)file->Get("XCone_subjets/pt_ak4");
  TH1F * eta_ak4 = (TH1F*)file->Get("XCone_subjets/eta_ak4");

  TH1F * area = (TH1F*)file->Get("XCone_subjets/area_all_subjets");
  TH1F * area_ak4 = (TH1F*)file->Get("XCone_subjets/area_ak4");

  // normalize hists
  Double_t norm = 1;
  L1->Scale(norm/(L1->Integral()));
  L2L3->Scale(norm/(L2L3->Integral()));
  all->Scale(norm/(all->Integral()));
  pt->Scale(norm/(pt->Integral()));
  eta->Scale(norm/(eta->Integral()));
  pt_had->Scale(norm/(pt_had->Integral()));
  eta_had->Scale(norm/(eta_had->Integral()));
  pt_lep->Scale(norm/(pt_lep->Integral()));
  eta_lep->Scale(norm/(eta_lep->Integral()));
  area->Scale(norm/(area->Integral()));
  L1_ak4->Scale(norm/(L1_ak4->Integral()));
  L2L3_ak4->Scale(norm/(L2L3_ak4->Integral()));
  all_ak4->Scale(norm/(all_ak4->Integral()));
  pt_ak4->Scale(norm/(pt_ak4->Integral()));
  eta_ak4->Scale(norm/(eta_ak4->Integral()));
  area_ak4->Scale(norm/(area_ak4->Integral()));

  // get means
  cout << endl << "===================================================" <<  endl <<  endl;
  cout << "L1     Mean XCone = " << L1->GetMean() << endl;
  cout << "L1     Mean   AK4 = " << L1_ak4->GetMean() << endl;
  cout << "L2L3   Mean XCone = " << L2L3->GetMean() << endl;
  cout << "L2L3   Mean   AK4 = " << L2L3_ak4->GetMean() << endl;
  cout << "L1L2L3 Mean XCone = " << all->GetMean() << endl;
  cout << "L1L2L3 Mean   AK4 = " << all_ak4->GetMean() << endl;
  cout << "area   Mean XCone = " << area->GetMean() << endl;
  cout << "area   Mean   AK4 = " << area_ak4->GetMean() << endl;
  cout << endl << "===================================================" <<  endl <<  endl;

  // closure tests
  cout << "L1     Integral XCone = " << L1->Integral() << endl;
  cout << "L1     Integral   AK4 = " << L1_ak4->Integral() << endl;
  cout << "L2L3   Integral XCone = " << L2L3->Integral() << endl;
  cout << "L2L3   Integral   AK4 = " << L2L3_ak4->Integral() << endl;
  cout << "L1L2L3 Integral XCone = " << all->Integral() << endl;
  cout << "L1L2L3 Integral   AK4 = " << all_ak4->Integral() << endl;
  cout << "pt     Integral XCone = " << pt->Integral() << endl;
  cout << "pt     Integral   AK4 = " << pt_ak4->Integral() << endl;
  cout << "eta    Integral XCone = " << eta->Integral() << endl;
  cout << "eta    Integral   AK4 = " << eta_ak4->Integral() << endl;
  cout << "area   Integral XCone = " << area->Integral() << endl;
  cout << "area   Integral   AK4 = " << area_ak4->Integral() << endl;
  cout << endl << "===================================================" <<  endl <<  endl;

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  int ylow, yhigh;
  ylow = 0;
  yhigh = 0;
  float xlow, xhigh;
  xlow = 0.6;
  xhigh = 1.4;

  L1->SetTitle("JEC L1");
  L1->GetXaxis()->SetTitle("L1 factor");
  L1->GetYaxis()->SetTitle("#Delta events/events");
  L1->GetYaxis()->SetTitleOffset(1.2);
  L1->GetXaxis()->SetRangeUser(xlow, xhigh);
  L1->GetYaxis()->SetRangeUser(ylow, 0.6);
  L1->SetLineColor(2);
  L1->SetLineWidth(4);
  L1_mean = new TLine(L1->GetMean(),0.,L1->GetMean(),0.6);
  L1_mean->SetLineColor(2);
  L1_mean->SetLineStyle(2);
  L1_mean->SetLineWidth(2);

  L1_ak4->SetLineColor(4);
  L1_ak4->SetLineWidth(4);
  L1_ak4_mean = new TLine(L1_ak4->GetMean(),0.,L1_ak4->GetMean(),0.6);
  L1_ak4_mean->SetLineColor(4);
  L1_ak4_mean->SetLineStyle(2);
  L1_ak4_mean->SetLineWidth(2);

  L2L3->SetTitle("JEC L2L3");
  L2L3->GetXaxis()->SetTitle("L2L3 factor");
  L2L3->GetYaxis()->SetTitle("#Delta events/events");
  L2L3->GetYaxis()->SetTitleOffset(1.2);
  L2L3->GetXaxis()->SetRangeUser(xlow, xhigh);
  L2L3->GetYaxis()->SetRangeUser(ylow, 0.4);
  L2L3->SetLineColor(2);
  L2L3->SetLineWidth(4);
  L2L3_mean = new TLine(L2L3->GetMean(),0.,L2L3->GetMean(),0.4);
  L2L3_mean->SetLineColor(2);
  L2L3_mean->SetLineStyle(2);
  L2L3_mean->SetLineWidth(2);

  L2L3_ak4->SetLineColor(4);
  L2L3_ak4->SetLineWidth(4);
  L2L3_ak4_mean = new TLine(L2L3_ak4->GetMean(),0.,L2L3_ak4->GetMean(),0.4);
  L2L3_ak4_mean->SetLineColor(4);
  L2L3_ak4_mean->SetLineStyle(2);
  L2L3_ak4_mean->SetLineWidth(2);

  all->SetTitle("JEC L1 + L2L3");
  all->GetXaxis()->SetTitle("tot JEC factor");
  all->GetYaxis()->SetTitle("#Delta events/events");
  all->GetYaxis()->SetTitleOffset(1.2);
  all->GetXaxis()->SetRangeUser(xlow, xhigh);
  all->GetYaxis()->SetRangeUser(ylow, 0.4);
  all->SetLineColor(2);
  all->SetLineWidth(4);
  all_mean = new TLine(all->GetMean(),0.,all->GetMean(),0.4);
  all_mean->SetLineColor(2);
  all_mean->SetLineStyle(2);
  all_mean->SetLineWidth(2);

  all_ak4->SetLineColor(4);
  all_ak4->SetLineWidth(4);
  all_ak4_mean = new TLine(all_ak4->GetMean(),0.,all_ak4->GetMean(),0.4);
  all_ak4_mean->SetLineColor(4);
  all_ak4_mean->SetLineStyle(2);
  all_ak4_mean->SetLineWidth(2);

  pt->SetTitle("p_{T} shape");
  pt->GetXaxis()->SetTitle("p_{T}");
  pt->GetYaxis()->SetTitle("#Delta events/events");
  pt->GetYaxis()->SetTitleOffset(1.2);
  pt->GetYaxis()->SetRangeUser(ylow, 0.07);
  pt->SetLineColor(2);
  pt->SetLineWidth(4);

  pt_had->SetTitle("p_{T} shape");
  pt_had->GetXaxis()->SetTitle("p_{T}");
  pt_had->GetYaxis()->SetTitle("#Delta events/events");
  pt_had->GetYaxis()->SetTitleOffset(1.2);
  pt_had->GetYaxis()->SetRangeUser(ylow, 0.07);
  pt_had->SetLineColor(2);
  pt_had->SetLineWidth(4);

  pt_lep->SetTitle("p_{T} shape");
  pt_lep->GetXaxis()->SetTitle("p_{T}");
  pt_lep->GetYaxis()->SetTitle("#Delta events/events");
  pt_lep->GetYaxis()->SetTitleOffset(1.2);
  pt_lep->GetYaxis()->SetRangeUser(ylow, 0.07);
  pt_lep->SetLineColor(2);
  pt_lep->SetLineWidth(4);

  pt_ak4->SetLineColor(4);
  pt_ak4->SetLineWidth(4);

  eta->SetTitle("eta shape");
  eta->GetXaxis()->SetTitle("#eta");
  eta->GetYaxis()->SetTitle("#Delta events/events");
  eta->GetYaxis()->SetTitleOffset(1.2);
  eta->GetYaxis()->SetRangeUser(ylow, 0.05);
  eta->SetLineColor(2);
  eta->SetLineWidth(4);

  eta_had->SetTitle("eta");
  eta_had->GetXaxis()->SetTitle("#eta");
  eta_had->GetYaxis()->SetTitle("#Delta events/events");
  eta_had->GetYaxis()->SetTitleOffset(1.2);
  eta_had->GetYaxis()->SetRangeUser(ylow, 0.05);
  eta_had->SetLineColor(2);
  eta_had->SetLineWidth(4);

  eta_lep->SetTitle("eta shape");
  eta_lep->GetXaxis()->SetTitle("#eta");
  eta_lep->GetYaxis()->SetTitle("#Delta events/events");
  eta_lep->GetYaxis()->SetTitleOffset(1.2);
  eta_lep->GetYaxis()->SetRangeUser(ylow, 0.05);
  eta_lep->SetLineColor(2);
  eta_lep->SetLineWidth(4);

  eta_ak4->SetLineColor(4);
  eta_ak4->SetLineWidth(4);

  area->SetTitle("area shape");
  area->GetXaxis()->SetTitle("area");
  area->GetYaxis()->SetTitle("#Delta events/events");
  area->GetYaxis()->SetTitleOffset(1.2);
  area->GetYaxis()->SetRangeUser(ylow, 0.25);
  area->SetLineColor(2);
  area->SetLineWidth(4);
  area_mean = new TLine(area->GetMean(),0.,area->GetMean(),0.25);
  area_mean->SetLineColor(2);
  area_mean->SetLineStyle(2);
  area_mean->SetLineWidth(2);

  area_ak4->SetLineColor(4);
  area_ak4->SetLineWidth(4);
  area_ak4_mean = new TLine(area_ak4->GetMean(),0.,area_ak4->GetMean(),0.25);
  area_ak4_mean->SetLineColor(4);
  area_ak4_mean->SetLineStyle(2);
  area_ak4_mean->SetLineWidth(2);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
 

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- draw Hists -----------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TCanvas *A = new TCanvas();
  leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(L1,"XCone","l");
  leg->AddEntry(L1_ak4,"AK4","l");
  leg->AddEntry(L1_mean,"XCone mean","l");
  leg->AddEntry(L1_ak4_mean,"AK4 mean","l");
  L1->Draw("HIST");
  L1_ak4->Draw("HIST SAME");
  leg->Draw("");
  L1_mean->Draw("");
  L1_ak4_mean->Draw("");
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/L1.pdf");
  TCanvas *B = new TCanvas();
  leg2 = new TLegend(0.65,0.85,0.85,0.65);
  leg2->AddEntry(L2L3,"XCone","l");
  leg2->AddEntry(L2L3_ak4,"AK4","l");
  leg2->AddEntry(L2L3_mean,"XCone mean","l");
  leg2->AddEntry(L2L3_ak4_mean,"AK4 mean","l");
  L2L3->Draw("HIST");
  L2L3_ak4->Draw("HIST SAME");
  leg2->Draw("");
  L2L3_mean->Draw("");
  L2L3_ak4_mean->Draw("");
  gPad->RedrawAxis();
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/L2L3.pdf");
  TCanvas *C = new TCanvas();
  leg3 = new TLegend(0.65,0.85,0.85,0.65);
  leg3->AddEntry(all,"XCone","l");
  leg3->AddEntry(all_ak4,"AK4","l");
  leg3->AddEntry(all_mean,"XCone mean","l");
  leg3->AddEntry(all_ak4_mean,"AK4 mean","l");
  all->Draw("HIST");
  all_ak4->Draw("HIST SAME");
  leg3->Draw("");
  all_mean->Draw("");
  all_ak4_mean->Draw("");
  gPad->RedrawAxis();
  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/L1L2L3.pdf");
  TCanvas *D = new TCanvas();
  leg4 = new TLegend(0.65,0.85,0.85,0.65);
  leg4->AddEntry(pt,"XCone","l");
  leg4->AddEntry(pt_ak4,"AK4","l");
  pt->Draw("HIST");
  pt_ak4->Draw("HIST SAME");
  leg4->Draw("");
  gPad->RedrawAxis();
  D->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/pt.pdf");
  TCanvas *E = new TCanvas();
  leg5 = new TLegend(0.65,0.85,0.85,0.65);
  leg5->AddEntry(eta,"XCone","l");
  leg5->AddEntry(eta_ak4,"AK4","l");
  eta->Draw("HIST");
  eta_ak4->Draw("HIST SAME");
  leg5->Draw("");
  gPad->RedrawAxis();
  E->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/eta.pdf");
  TCanvas *F = new TCanvas();
  leg4 = new TLegend(0.65,0.85,0.85,0.65);
  leg4->AddEntry(pt_had,"XCone had","l");
  leg4->AddEntry(pt_ak4,"AK4","l");
  pt_had->Draw("HIST");
  pt_ak4->Draw("HIST SAME");
  leg4->Draw("");
  gPad->RedrawAxis();
  F->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/pthad.pdf");
  TCanvas *G = new TCanvas();
  leg5 = new TLegend(0.65,0.85,0.85,0.65);
  leg5->AddEntry(eta_had,"XCone had","l");
  leg5->AddEntry(eta_ak4,"AK4","l");
  eta_had->Draw("HIST");
  eta_ak4->Draw("HIST SAME");
  leg5->Draw("");
  gPad->RedrawAxis();
  G->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/etahad.pdf");
  TCanvas *H = new TCanvas();
  leg6 = new TLegend(0.65,0.85,0.85,0.65);
  leg6->AddEntry(area,"XCone","l");
  leg6->AddEntry(area_ak4,"AK4","l");
  leg6->AddEntry(area_mean,"XCone mean","l");
  leg6->AddEntry(area_ak4_mean,"AK4 mean","l");
  area->Draw("HIST");
  area_ak4->Draw("HIST SAME");
  area_mean->Draw();
  area_ak4_mean->Draw();
  leg6->Draw("");
  gPad->RedrawAxis();
  H->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JEC_factor/area.pdf");
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
