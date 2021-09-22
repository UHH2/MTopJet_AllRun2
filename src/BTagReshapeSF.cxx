#include <UHH2/MTopJet/include/BTagReshapeSF.h>

BTagReshapeSF::BTagReshapeSF(uhh2::Context & ctx, TString year, TString channel){

  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  if (!isMC) {
    cout << "Warning: MCBTagReshapeSFs will not have an effect on "
    <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  TFile *file = new TFile("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/ScaleFactors/BTag/BTagReshapeWeights"+year+".root");
  h_sf = (TH1F*)file->Get("weight_"+channel);

}

bool BTagReshapeSF::process(uhh2::Event & event){
  if(!isMC) return true;

  double njets = event.jets->size();
  int bin = h_sf->GetXaxis()->FindBin(njets);
  double sf = h_sf->GetBinContent(bin);

  double weight = event.weight;
  event.weight = weight * sf;

  return true;
}
