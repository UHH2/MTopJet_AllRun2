void ComparisonToMuon(){

  vector<TString> Elec2016 = {
    "SingleElecB_2016v3",
    "SingleElecC_2016v3",
    "SingleElecD_2016v3",
    "SingleElecE_2016v3",
    "SingleElecF_2016v3",
    "SingleElecG_2016v3",
    "SingleElecH_2016v3",
    "SinglePhotonB_ver2_2016v3",
    "SinglePhotonC_2016v3",
    "SinglePhotonD_2016v3",
    "SinglePhotonE_2016v3",
    "SinglePhotonF_2016v3",
    "SinglePhotonG_2016v3",
    "SinglePhotonH_2016v3"
  };

  vector<TString> Elec2017 = {
    "SingleElecB_2017v2",
    "SingleElecC_2017v2",
    "SingleElecD_2017v2",
    "SingleElecE_2017v2",
    "SingleElecF_2017v2",
    "SinglePhotonB_2017v2",
    "SinglePhotonC_2017v2",
    "SinglePhotonD_2017v2",
    "SinglePhotonE_2017v2",
    "SinglePhotonF_2017v2"
  };

  vector<TString> Elec2018 = {
    "EGammaA_2018",
    "EGammaB_2018",
    "EGammaC_2018",
    "EGammaD_2018_1",
    "EGammaD_2018_2",
    "EGammaD_2018_3"
  };

  vector<TString> Muon2016 = {
    "SingleMuonB_2016v3",
    "SingleMuonC_2016v3",
    "SingleMuonD_2016v3",
    "SingleMuonE_2016v3",
    "SingleMuonF_2016v3",
    "SingleMuonG_2016v3",
    "SingleMuonH_2016v3"
  };

  vector<TString> Muon2017 = {
    "SingleMuonB_2017v2",
    "SingleMuonC_2017v2",
    "SingleMuonD_2017v2",
    "SingleMuonE_2017v2",
    "SingleMuonF_2017v2"
  };

  vector<TString> Muon2018 = {
    "SingleMuonA_2018",
    "SingleMuonB_2018",
    "SingleMuonC_2018",
    "SingleMuonD_2018"
  };

  vector<TString> years = {"2016", "2017", "2018"};

  vector<vector<TString>> Elec;
  vector<vector<TString>> Muon;

  Elec.push_back(Elec2016);
  Elec.push_back(Elec2017);
  Elec.push_back(Elec2018);
  Muon.push_back(Muon2016);
  Muon.push_back(Muon2017);
  Muon.push_back(Muon2018);


  if(Elec.size() != Muon.size()){
    cout << "Muon and Elec don't have the same lenght!" << endl;
    return;
  }

  for(unsigned int i=0; i<Elec.size(); i++){
    vector<vector<TString>> histnames;
    histnames.push_back({"01_Cleaner_Elec/pt","ptlep_BeforeTrigger","Lepton p_{T}"});
    histnames.push_back({"02_Trigger_Elec/pt","ptlep_AfterTrigger","Lepton p_{T}"});
    histnames.push_back({"06_MET_Elec/pt","ptlep","Lepton p_{T}"});
    histnames.push_back({"06_MET_Jets/pt_jet","ptjet","AK4 jet p_{T}"});
    histnames.push_back({"06_MET_Event/MET","MET","p_{T}^{miss}"});
    histnames.push_back({"06_MET_lumi/luminosity","Lumi"," "});
    histnames.push_back({"08_bTag_Elec/pt","ptlep_AfterHEM","Lepton p_{T}"});
    histnames.push_back({"08_bTag_jets/pt_jet","ptjet_AfterHEM","AK4 jet p_{T}"});
    histnames.push_back({"08_bTag_Event/MET","MET_AfterHEM","p_{T}^{miss}"});
    histnames.push_back({"08_bTag_lumi/luminosity","Lumi_AfterHEM"," "});

    // Cutflow only for 2016
    if(i==0){
      // vector<TString> cutflow_names ={"01_Cleaner_Event/MET","02_Trigger_Event/MET","03_Lepton_Event/MET","05_TwoD_Event/MET","06_MET_Event/MET","08_bTag_Event/MET"};
      vector<TString> cutflow_names ={"02_Trigger_Event/MET","03_Lepton_Event/MET","05_TwoD_Event/MET","06_MET_Event/MET","08_bTag_Event/MET"};
      TH1F* cutflow_elec = new TH1F("cutflow", "selection steps", 6, -0.5, 5.5);
      TH1F* cutflow_muon = new TH1F("cutflow", "selection steps", 6, -0.5, 5.5);
      for(unsigned int k=0; k<cutflow_names.size(); k++){
        cout << k << endl;
        double Count_el = 0;
        for(unsigned int j=0; j<Elec[i].size(); j++){
          TFile * file = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/Selection/elec/uhh2.AnalysisModuleRunner.DATA.DATA_"+Elec[i][j]+".root");
          TH1F* hist = (TH1F*) file->Get(cutflow_names[k]);
          Count_el += hist->Integral();
        }
        double Count_mu = 0;
        for(unsigned int j=0; j<Muon[i].size(); j++){
          TFile * file = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/Selection/elec/uhh2.AnalysisModuleRunner.DATA.DATA_"+Elec[i][j]+".root");
          TH1F* hist = (TH1F*) file->Get(cutflow_names[k]);
          Count_mu += hist->Integral();
        }
        cout << Count_el << endl;
        cout << Count_mu << endl;

        int norm_el, norm_mu;
        if(k==0){
          norm_el = Count_el;
          norm_mu = Count_mu;
        }
        cutflow_elec->SetBinContent(k+1, Count_el/norm_el);
        cutflow_muon->SetBinContent(k+1, Count_mu/norm_mu);
      }
      TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
      cutflow_elec->SetLineColor(kBlue);
      cutflow_elec->SetLineWidth(3);
      cutflow_muon->SetLineColor(kRed);
      cutflow_muon->SetLineWidth(3);
      cutflow_elec->Draw("HIST");
      cutflow_muon->Draw("HIST SAME");
      TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.85);
      leg->AddEntry(cutflow_elec, "Electron channel", "l");
      leg->AddEntry(cutflow_muon, "Muon channel", "l");
      leg->Draw();
      c1->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2/Cutflow/Cutflow_data_2016.pdf");
    }
    ////



    for(unsigned int ihist = 0; ihist<histnames.size(); ihist++){

      TH1F* hist_elec;
      for(unsigned int j=0; j<Elec[i].size(); j++){
        TFile * file = new TFile("/nfs/dust/cms/user/schwarzd/MTopJet_Run2/Selection/elec/uhh2.AnalysisModuleRunner.DATA.DATA_"+Elec[i][j]+".root");
        TH1F* hist = (TH1F*) file->Get(histnames[ihist][0]);
        if(histnames[ihist][0].Contains("luminosity")){
          if(years[i] == "2018") hist->Rebin(10);
          else                   hist->Rebin(10);
        }
        if(j==0) hist_elec = hist;
        else     hist_elec->Add(hist);
      }
      TH1F* hist_elec_norm = (TH1F*)hist_elec->Clone();
      hist_elec_norm->Scale(1/hist_elec->Integral());



      TH1F* hist_muon;
      for(unsigned int j=0; j<Muon[i].size(); j++){
        TFile * file = new TFile("/nfs/dust/cms/user/paaschal/MTopJet/Selection/"+years[i]+"/muon/uhh2.AnalysisModuleRunner.DATA.DATA_"+Muon[i][j]+".root");
        TString histname = histnames[ihist][0];
        if(histname.Contains("Elec")) histname = histname.ReplaceAll("Elec", "Muon");
        // cout << histname << endl;
        TH1F* hist = (TH1F*) file->Get(histname);
        if(histnames[ihist][0].Contains("luminosity")){
          if(years[i] == "2018") hist->Rebin(10);
          else                   hist->Rebin(10);
        }
        if(j==0) hist_muon = hist;
        else     hist_muon->Add(hist);
      }
      TH1F* hist_muon_norm = (TH1F*)hist_muon->Clone();
      hist_muon_norm->Scale(1/hist_muon->Integral());

      gStyle->SetOptStat(0);
      gStyle -> SetPadTickX(1);
      gStyle -> SetPadTickY(1);

      TCanvas *c = new TCanvas("c", "c", 600, 600);
      hist_elec->SetTitle("");
      hist_elec->GetXaxis()->SetTitle(histnames[ihist][2]);
      hist_elec->GetYaxis()->SetTitle("Events");
      hist_elec->GetYaxis()->SetRangeUser(0, 1.4*hist_elec->GetMaximum());
      hist_elec->SetMarkerStyle(8);
      hist_elec->SetMarkerColor(kRed);
      hist_elec->SetLineColor(kRed);
      hist_muon->SetMarkerStyle(8);
      hist_muon->SetMarkerColor(kAzure+7);
      hist_muon->SetLineColor(kAzure+7);
      hist_elec->Draw("E1");
      hist_muon->Draw("E1 SAME");
      TLegend * leg = new TLegend(0.5, 0.7, 0.9, 0.9);
      leg->AddEntry(hist_elec, "Electron channel", "pl");
      leg->AddEntry(hist_muon, "Muon channel", "pl");
      leg->Draw();
      c->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2_channels/"+histnames[ihist][1]+"_"+years[i]+".pdf");
      delete c;

      TCanvas *d = new TCanvas("d", "d", 600, 600);
      hist_elec_norm->SetTitle("");
      hist_elec_norm->GetXaxis()->SetTitle(histnames[ihist][2]);
      hist_elec_norm->GetYaxis()->SetTitle("#Delta N / N");
      hist_elec_norm->GetYaxis()->SetRangeUser(0, 1.4*hist_elec_norm->GetMaximum());
      hist_elec_norm->SetMarkerStyle(8);
      hist_elec_norm->SetMarkerColor(kRed);
      hist_elec_norm->SetLineColor(kRed);
      hist_muon_norm->SetMarkerStyle(8);
      hist_muon_norm->SetMarkerColor(kAzure+7);
      hist_muon_norm->SetLineColor(kAzure+7);
      hist_elec_norm->Draw("E1");
      hist_muon_norm->Draw("E1 SAME");
      TLegend * leg_norm = new TLegend(0.5, 0.7, 0.9, 0.9);
      leg_norm->AddEntry(hist_elec_norm, "Electron channel", "pl");
      leg_norm->AddEntry(hist_muon_norm, "Muon channel", "pl");
      leg_norm->Draw();
      d->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MTopRun2_channels/"+histnames[ihist][1]+"_"+years[i]+"_norm.pdf");
      delete d;
    }

  }
  return;
}
