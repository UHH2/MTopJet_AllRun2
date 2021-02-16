void Normalization(){

  // vector<TString> years = {"2016v3", "2017v2", "2018"};
  vector<TString> years = {"2016v3"};

  vector<TString> bkg = {"WJets", "SingleTop", "other"};

  TString histname = "XCone_cor/M_jet1_";
  TString histname_pre1 = "PreSel01_XCone/M_jet1_";
  TString histname_pre2 = "PreSel02_XCone/M_jet1_";
  TString histname_pre3 = "PreSel03_XCone/M_jet1_";
  TString histname_pre3b = "PreSel03b_XCone/M_jet1_";
  TString histname_pre4 = "PreSel04_XCone/M_jet1_";
  bool domuon = false;

  for(unsigned int s = 0; s<=1; s++){
    if(s == 1) domuon = true;
    for(auto year: years){

      TFile *f_data, *f_ttbar;
      if(domuon){
        f_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.DATA.DATA_"+year+".root");
        f_ttbar = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_"+year+".root");
      }
      else{
        f_data = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.DATA.DATA_"+year+".root");
        f_ttbar = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC.TTbar_"+year+".root");
      }

      TH1F* h_data = (TH1F*) f_data->Get(histname);
      TH1F* h_ttbar = (TH1F*) f_ttbar->Get(histname);

      TH1F* h_data_pre1 = (TH1F*) f_data->Get(histname_pre1);
      TH1F* h_ttbar_pre1 = (TH1F*) f_ttbar->Get(histname_pre1);

      TH1F* h_data_pre2 = (TH1F*) f_data->Get(histname_pre2);
      TH1F* h_ttbar_pre2 = (TH1F*) f_ttbar->Get(histname_pre2);

      TH1F* h_data_pre3 = (TH1F*) f_data->Get(histname_pre3);
      TH1F* h_ttbar_pre3 = (TH1F*) f_ttbar->Get(histname_pre3);

      TH1F* h_data_pre3b = (TH1F*) f_data->Get(histname_pre3b);
      TH1F* h_ttbar_pre3b = (TH1F*) f_ttbar->Get(histname_pre3b);

      TH1F* h_data_pre4 = (TH1F*) f_data->Get(histname_pre4);
      TH1F* h_ttbar_pre4 = (TH1F*) f_ttbar->Get(histname_pre4);

      TH1F* h_bkg;
      TH1F *h_bkg_pre1, *h_bkg_pre2, *h_bkg_pre3, *h_bkg_pre3b, *h_bkg_pre4;

      for(unsigned int i=0; i<bkg.size(); i++){
        TFile *f_bkg;
        if(domuon) f_bkg = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC."+bkg[i]+"_"+year+".root");
        else       f_bkg = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/PostSel/elec/uhh2.AnalysisModuleRunner.MC."+bkg[i]+"_"+year+".root");
        TH1F* hist = (TH1F*) f_bkg->Get(histname);
        TH1F* hist_pre1 = (TH1F*) f_bkg->Get(histname_pre1);
        TH1F* hist_pre2 = (TH1F*) f_bkg->Get(histname_pre2);
        TH1F* hist_pre3 = (TH1F*) f_bkg->Get(histname_pre3);
        TH1F* hist_pre3b = (TH1F*) f_bkg->Get(histname_pre3b);
        TH1F* hist_pre4 = (TH1F*) f_bkg->Get(histname_pre4);

        if(i==0){
          h_bkg = hist;
          h_bkg_pre1 = hist_pre1;
          h_bkg_pre2 = hist_pre2;
          h_bkg_pre3 = hist_pre3;
          h_bkg_pre3b = hist_pre3b;
          h_bkg_pre4 = hist_pre4;
        }
        else{
          h_bkg->Add(hist);
          h_bkg_pre1->Add(hist_pre1);
          h_bkg_pre2->Add(hist_pre2);
          h_bkg_pre3->Add(hist_pre3);
          h_bkg_pre3b->Add(hist_pre3b);
          h_bkg_pre4->Add(hist_pre4);
        }
      }

      cout << "------------------------------------------------" << endl;
      cout << "Year = " << year;
      TString channel = "elec";
      if(domuon) channel = "muon";
      cout << " ("+channel+")" << endl;

      cout << "Events in data (Pre01) = " << h_data_pre1->Integral()<< endl;
      cout << "Events in data (Pre02) = " << h_data_pre2->Integral()<< endl;
      cout << "Events in data (Pre03) = " << h_data_pre3->Integral()<< endl;
      cout << "Events in data (Pre03b) = " << h_data_pre3b->Integral()<< endl;
      cout << "Events in data (Pre04) = " << h_data_pre4->Integral()<< endl;
      cout << "Events in data (Full)  = " << h_data->Integral()<< endl;
      cout << "SF (Pre01) = " << (h_data_pre1->Integral() - h_bkg_pre1->Integral()) / h_ttbar_pre1->Integral() << endl;
      cout << "SF (Pre02) = " << (h_data_pre2->Integral() - h_bkg_pre2->Integral()) / h_ttbar_pre2->Integral() << endl;
      cout << "SF (Pre03) = " << (h_data_pre3->Integral() - h_bkg_pre3->Integral()) / h_ttbar_pre3->Integral() << endl;
      cout << "SF (Pre03b) = " << (h_data_pre3b->Integral() - h_bkg_pre3b->Integral()) / h_ttbar_pre3b->Integral() << endl;
      cout << "SF (Pre04) = " << (h_data_pre4->Integral() - h_bkg_pre4->Integral()) / h_ttbar_pre4->Integral() << endl;
      cout << "SF (Full)  = " << (h_data->Integral() - h_bkg->Integral()) / h_ttbar->Integral() << endl;

      TH1F* h_cutflow_data = new TH1F("cutflow", "Selection steps", 6, 0.5, 6.5);
      TH1F* h_cutflow_ttbar = new TH1F("cutflow", "Selection steps", 6, 0.5, 6.5);

      h_cutflow_data->SetBinContent(1, h_data_pre1->Integral()/h_data_pre1->Integral());
      h_cutflow_data->SetBinContent(2, h_data_pre2->Integral()/h_data_pre1->Integral());
      h_cutflow_data->SetBinContent(3, h_data_pre3->Integral()/h_data_pre1->Integral());
      h_cutflow_data->SetBinContent(3, h_data_pre3b->Integral()/h_data_pre1->Integral());
      h_cutflow_data->SetBinContent(4, h_data_pre4->Integral()/h_data_pre1->Integral());
      h_cutflow_data->SetBinContent(5, h_data->Integral()/h_data_pre1->Integral());

      h_cutflow_ttbar->SetBinContent(1, h_ttbar_pre1->Integral()/h_ttbar_pre1->Integral());
      h_cutflow_ttbar->SetBinContent(2, h_ttbar_pre2->Integral()/h_ttbar_pre1->Integral());
      h_cutflow_ttbar->SetBinContent(3, h_ttbar_pre3->Integral()/h_ttbar_pre1->Integral());
      h_cutflow_ttbar->SetBinContent(3, h_ttbar_pre3b->Integral()/h_ttbar_pre1->Integral());
      h_cutflow_ttbar->SetBinContent(4, h_ttbar_pre4->Integral()/h_ttbar_pre1->Integral());
      h_cutflow_ttbar->SetBinContent(5, h_ttbar->Integral()/h_ttbar_pre1->Integral());

      // TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
      // h_cutflow_data->Draw("HIST");
      // c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2/Cutflow/Cutflow_data_"+year+"_"+channel+".pdf");
      // delete h_cutflow_data;
      // delete c1;
      //
      // TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
      // h_cutflow_ttbar->Draw("HIST");
      // c2->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2/Cutflow/Cutflow_ttbar_"+year+"_"+channel+".pdf");
      // delete h_cutflow_ttbar;
      // delete c2;
    }
  }

  return;
}
