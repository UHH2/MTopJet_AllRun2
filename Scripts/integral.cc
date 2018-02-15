void integral()
{

  // declare files
  TFile *TTbar = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *TTbar2 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_2gamma.root");
  TFile *TTbar4 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_4gamma.root");
  TFile *TTbar8 = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_8gamma.root");

  // calculate integral in every Selection Step
  TString hist = "XCone_GEN_GenOnly/number_HadJet33";

  TH1F * h1 = (TH1F*)TTbar->Get(hist);
  double h1_int = h1->Integral();

  TH1F * h2 = (TH1F*)TTbar2->Get(hist);
  double h2_int = h2->Integral();

  TH1F * h4 = (TH1F*)TTbar4->Get(hist);
  double h4_int = h4->Integral();

  TH1F * h8 = (TH1F*)TTbar8->Get(hist);
  double h8_int = h8->Integral();

  cout << "SM: " << h1_int << endl;
  cout << "2G: " << h2_int << endl;
  cout << "4G: " << h4_int << endl;
  cout << "8G: " << h8_int << endl;
}

