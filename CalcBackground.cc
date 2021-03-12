void CalcBackground(){

  vector<TString> years = {"2016v3", "2017v2", "2018"};
  // vector<TString> years = {"2016v3"};

  vector<TString> processes = {"TTbar", "WJets", "SingleTop", "DYJets", "DiBoson", "QCD"};
  vector<TString> channels = {"muon", "elec"};

  TString histname = "XCone_cor/M_jet1_";


  for(auto year: years){
    double mcsum = 0;
    vector<double> integrals;
    for(auto process: processes){
      double integral = 0;
      for(auto channel: channels){
        TFile *file = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/"+channel+"/uhh2.AnalysisModuleRunner.MC."+process+"_"+year+".root");
        TH1F* hist = (TH1F*) file->Get(histname);
        integral += hist->Integral();
        mcsum += hist->Integral();
      }
      integrals.push_back(integral);
    }
    cout << "--------------------" << endl;
    cout << year << endl;
    for(int i=0; i<processes.size();i++) cout << "   - " << processes[i] << " = " << 100*integrals[i]/mcsum << "%"<< endl;
  }

  return;
}
