void fit_Wmass(TString mode)
{
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- declare files --------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  TFile *TTbar = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *WJets = new TFile("uhh2.AnalysisModuleRunner.MC.WJets.root");
  TFile *ST = new TFile("uhh2.AnalysisModuleRunner.MC.SingleTop.root");
  TFile *other = new TFile("uhh2.AnalysisModuleRunner.MC.other.root");
  TFile *all = new TFile("uhh2.AnalysisModuleRunner.MC.Background.root");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- get mass plots -------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  std::stringstream hist_title;
  if(mode == "raw")  hist_title << "XCone_raw_subjets_SF/min_mass_Wjet_zoom";
  if(mode == "jec")  hist_title << "XCone_jec_subjets_SF/min_mass_Wjet_zoom";
  if(mode == "cor")  hist_title << "XCone_cor_subjets_SF/min_mass_Wjet_zoom";

  // TH1F * TTbar_ = (TH1F*)TTbar->Get(hist_title.str().c_str());
  TH1F * WJets_ = (TH1F*)WJets->Get(hist_title.str().c_str());
  TH1F * ST_ = (TH1F*)ST->Get(hist_title.str().c_str());
  TH1F * other_ = (TH1F*)other->Get(hist_title.str().c_str());
  TH1F * all_ = (TH1F*)all->Get(hist_title.str().c_str());

  other_->Add(WJets_);
  other_->Add(ST_);




  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- Fit function to peak -- ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  bool do_fit = false;
  TF1* fit;
  if(mode == "jec")fit = new TF1("fit", "gaus", 75, 100); 
  else fit = new TF1("fit", "gaus", 70, 95); 

  all_->Fit(fit, "R");
  fit->SetLineWidth(3);
  fit->SetLineColor(2);
  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------- set up lines and titles ----------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------
  int xlow, xhigh;
  int ylow, yhigh;
  xlow = 65;
  xhigh = 120;
  ylow = 0;
  yhigh = (all_->GetMaximum())*1.1;
  all_->SetTitle(" ");
  all_->GetXaxis()->SetTitle("min(m_{ij}) [GeV]");
  all_->GetYaxis()->SetTitle("events");
  all_->GetYaxis()->SetRangeUser(ylow, yhigh);
  all_->GetXaxis()->SetRangeUser(xlow, xhigh);
  all_->GetXaxis()->SetNdivisions(505);
  all_->GetYaxis()->SetNdivisions(505);
  all_->GetXaxis()->SetTitleSize(0.05);
  all_->GetXaxis()->SetTitleOffset(0.9);
  all_->GetYaxis()->SetTitleSize(0.06);
  all_->GetYaxis()->SetTitleOffset(1.1);

  all_->SetFillColor(kGray);
  all_->SetLineColor(kBlack);

  other_->SetFillColor(kGray+2);
  other_->SetLineColor(kBlack);

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
  all_->Draw("HIST");
  other_->Draw("HIST SAME");
  fit->Draw("SAME");


  leg = new TLegend(0.60,0.65,0.85,0.75);
  leg->AddEntry(all_,"t#bar{t}","f");
  leg->AddEntry(other_,"backgrounds","f");
  leg->AddEntry(fit, "fit","l");

 
  leg->Draw("");

  std::stringstream title;
  title << "Fit mean = " << std::setprecision(4) << fit->GetParameter(1)<<" GeV";
  TLatex *t = new TLatex(.60,.6,title.str().c_str()); 
  t->SetTextSize(0.03);
  t->SetNDC(kTRUE);
  t->Draw("");
  std::stringstream title2;
  title2 << "Fit rms    = " << std::setprecision(4) << fit->GetParameter(2)<<" GeV";
  TLatex *t2 = new TLatex(.60,.55,title2.str().c_str()); 
  t2->SetTextSize(0.03);
  t2->SetNDC(kTRUE);
  t2->Draw("");
  std::stringstream title3;
  title3 << "#chi^{2}/ndf = "  << fit->GetChisquare()/fit->GetNDF();
  TLatex *t3 = new TLatex(.60,.5,title3.str().c_str()); 
  t3->SetTextSize(0.03);
  t3->SetNDC(kTRUE);
  t3->Draw("");

  if(mode == "raw") A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/WMassFit_raw.pdf");
  if(mode == "jec") A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/WMassFit_jec.pdf");
  if(mode == "cor") A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/WMassFit_cor.pdf");



  

  gPad->RedrawAxis(); 
 

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
