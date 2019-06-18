#include "UHH2/MTopJet/include/ElecTriggerSF.h"

ElecTriggerSF::ElecTriggerSF(uhh2::Context & ctx, std::string var, TString pe){
  pteta = pe;
  if(pteta != "pt" && pteta != "eta" && pteta != "eta_ptbins"){
    cout << "Warning: You should select 'pt', 'eta' or 'eta_ptbins' for ElecTriggerSF class" << endl;
    return;
  }

  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  if (!isMC) {
    cout << "Warning: MCElecScaleFactor will not have an effect on "
         <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }


  TString dir = "/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/ElecTriggerSF/";
  TString filename = dir + "ElecTriggerSF.root";
  TFile *file = new TFile(filename);

  TString histname = "Central";
  if(var == "up")        histname = "Up";
  else if(var == "down") histname = "Down";

  // for the options pt and eta the low and high histograms are exactly the same
  // its just to make it easier to code
  if(pteta == "pt"){
    cout << "ElecTriggerSF: You selected the pt dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_pt");
    h_sf_me = (TH1F*)file->Get(histname+"_pt");
    h_sf_hi = (TH1F*)file->Get(histname+"_pt");
  }
  else if(pteta == "eta"){
    cout << "ElecTriggerSF: You selected the eta dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_eta");
    h_sf_me = (TH1F*)file->Get(histname+"_eta");
    h_sf_hi = (TH1F*)file->Get(histname+"_eta");
  }
  else{
    cout << "ElecTriggerSF: You selected the eta (pt bins) dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_eta_lowpt");
    h_sf_me = (TH1F*)file->Get(histname+"_eta_midpt");
    h_sf_hi = (TH1F*)file->Get(histname+"_eta_highpt");
  }
}

bool ElecTriggerSF::process(uhh2::Event & event){
  if(!isMC) return true;
  if(pteta != "pt" && pteta != "eta" && pteta != "eta_ptbins") return false;

  if(event.electrons->size() < 1) return false;
  double weight = event.weight;
  double pt = event.electrons->at(0).pt();
  double eta = event.electrons->at(0).eta();
  double UsedVariable;
  if(pteta == "pt") UsedVariable = pt;
  else              UsedVariable = eta;

  int bin = 0;
  double sf = 1.0;
  if(pt < 120){
    bin = h_sf_lo->GetXaxis()->FindBin(UsedVariable);
    sf  = h_sf_lo->GetBinContent(bin);
  }
  else if(pt > 120 && pt < 200){
    bin = h_sf_me->GetXaxis()->FindBin(UsedVariable);
    sf  = h_sf_me->GetBinContent(bin);
  }
  else if(pt > 200){
    bin = h_sf_hi->GetXaxis()->FindBin(UsedVariable);
    sf  = h_sf_hi->GetBinContent(bin);
  }
  event.weight = weight * sf;

  return true;
}
