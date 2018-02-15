void ttbar_fraction()
{

  // declare files
  TFile *DATA = new TFile("uhh2.AnalysisModuleRunner.DATA.DATA.root");
  TFile *DY = new TFile("uhh2.AnalysisModuleRunner.MC.DY.root");
  TFile *ST = new TFile("uhh2.AnalysisModuleRunner.MC.SingleTop.root");
  TFile *QCD = new TFile("uhh2.AnalysisModuleRunner.MC.QCD.root");
  TFile *WJets = new TFile("uhh2.AnalysisModuleRunner.MC.WJets.root");
  TFile *TTbar = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *DiBoson = new TFile("uhh2.AnalysisModuleRunner.MC.DiBoson.root");

  // calculate integral in every Selection Step
  TH1F * DATA1_ = (TH1F*)DATA->Get("XCone_cor/number_hadjet");
  double DATA_int = DATA1_->Integral();

  TH1F * DY1_ = (TH1F*)DY->Get("XCone_cor/number_hadjet");
  double DY_int = DY1_->Integral();

  TH1F * ST1_ = (TH1F*)ST->Get("XCone_cor/number_hadjet");
  double ST_int = ST1_->Integral();

  TH1F * QCD1_ = (TH1F*)QCD->Get("XCone_cor/number_hadjet");
  double QCD_int = QCD1_->Integral();

  TH1F * WJets1_ = (TH1F*)WJets->Get("XCone_cor/number_hadjet");
  double WJets_int = WJets1_->Integral();
 
  TH1F * TTbar1_ = (TH1F*)TTbar->Get("XCone_cor/number_hadjet");
  double TTbar_int = TTbar1_->Integral();

  TH1F * TTbar2_ = (TH1F*)TTbar->Get("750_xcone/number_hadjet");
  double TTbar_750xcone_int = TTbar2_->Integral();

  TH1F * TTbar3_ = (TH1F*)TTbar->Get("750_ak8/number_jets");
  double TTbar_750ak8_int = TTbar3_->Integral();

  TH1F * DiBoson1_ = (TH1F*)DiBoson->Get("XCone_cor/number_hadjet");
  double DiBoson_int = DiBoson1_->Integral();

  // add all backgrounds which are not TTbar for every selection step
  double back_all = DY_int + ST_int + QCD_int + WJets_int + DiBoson_int + TTbar_int;
  double back_without_ttbar = DY_int + ST_int + QCD_int + WJets_int + DiBoson_int;
  double fraction = TTbar_int / back_all;
  double percent = fraction * 100;
  double DATA_ttbar = DATA_int - back_without_ttbar;
  double SF_ttbar = DATA_ttbar / TTbar_int;
  double ratio_before_SF = DATA_int/back_all;
  double ratio_after_SF = DATA_int/(back_without_ttbar + SF_ttbar*TTbar_int);
  double fraction_ = (SF_ttbar*TTbar_int) / (back_without_ttbar + SF_ttbar*TTbar_int);
  double percent_ = fraction_ * 100;


  cout << "Sample                  | #events  " << endl;
  cout << "----------------------------------------" << endl;
  cout << "DY                      |  " << DY_int << endl;
  cout << "Single Top              |  " << ST_int << endl;
  cout << "WJets                   |  " << WJets_int << endl;
  cout << "QCD                     |  " << QCD_int << endl;
  cout << "DiBoson                 |  " << DiBoson_int << endl;
  cout << "ttbar:                  |  " << TTbar_int << endl;
  cout << "ttbar ak8 750GeV:       |  " << TTbar_750ak8_int << endl;
  cout << "ttbar xcone 750GeV:     |  " << TTbar_750xcone_int << endl;
  cout << "non ttbar background:   |  " << back_without_ttbar << endl;
  cout << "total background:       |  " << back_all << endl;
  cout << "data:                   |  " << DATA_int << endl;
  cout << "data - non ttbar:       |  " << DATA_ttbar << endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << "scale factor TTbar      |  " << SF_ttbar << endl;
  cout << "data/MC ratio before SF |  " << ratio_before_SF << endl;
  cout << "data/MC ratio after SF  |  " << ratio_after_SF << endl;
  cout << "fraction of TTbar events in respect to total background (before SF): " << percent << "%"<<endl;
  cout << "fraction of TTbar events in respect to total background (after SF):  " << percent_ << "%"<<endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << "---------------------------------------------------------------------------" << endl;

}

