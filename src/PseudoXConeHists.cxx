#include <UHH2/MTopJet/include/PseudoXConeHists.h>


PseudoXConeHists::PseudoXConeHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // book all histograms here
  pt_diff_fat = book<TH1F>("pt_diff_fat", "( p_{T}^{XCone} - p_{T}^{PseudoXCone} ) / p_{T}^{XCone}", 1000, -1, 1);
  eta_diff_fat = book<TH1F>("eta_diff_fat", "( #eta^{XCone} - #eta^{PseudoXCone} ) / #eta^{XCone}", 1000, -1, 1);
  dR_fat = book<TH1F>("dR_fat", "#Delta R (XCone, PseudoXCone)", 350, 0, 3.5);

  pt_diff_sub = book<TH1F>("pt_diff_sub", "( p_{T}^{XCone} - p_{T}^{PseudoXCone} ) / p_{T}^{XCone} in subjets", 1000, -1, 1);
  eta_diff_sub = book<TH1F>("eta_diff_sub", "( #eta^{XCone} - #eta^{PseudoXCone} ) / #eta^{XCone} in subjets", 1000, -1, 1);
  dR_sub = book<TH1F>("dR_sub", "#Delta R (XCone, PseudoXCone) in subjets", 350, 0, 3.5);



  // handle for jets
  h_xcone=ctx.get_handle<std::vector<TopJet>>("xconeCHS");
  h_pseudoxcone=ctx.get_handle<std::vector<TopJet>>("pseudoxconeCHS");
}



void PseudoXConeHists::fill(const Event & event){

  // get weight
  double weight = event.weight;

  //---------------------------------------------------------------------------------------
  //--------------------------------- define needed objects-----------------------------------
  //---------------------------------------------------------------------------------------
  // define all objects needed
  std::vector<TopJet> xconejets = event.get(h_xcone);
  if(xconejets.size() != 2) return;
  if(xconejets.at(0).subjets().size() != 3 || xconejets.at(1).subjets().size() != 3) return;

  std::vector<TopJet> pseudoxconejets = event.get(h_pseudoxcone);
  if(pseudoxconejets.size() < 2) return;
  if(pseudoxconejets.at(0).subjets().size() != 3 || pseudoxconejets.at(1).subjets().size() != 3) return;

  //---------------------------------------------------------------------------------------
  //------------------------ add subjets without and with  pt cut--------------------------
  //---------------------------------------------------------------------------------------
  for(unsigned int i=0; i<xconejets.size(); i++){
    double dR = deltaR(xconejets.at(i), pseudoxconejets.at(i));
    double DeltaPT = xconejets.at(i).v4().Pt() - pseudoxconejets.at(i).v4().Pt();
    double DeltaEta = xconejets.at(i).v4().Eta() - pseudoxconejets.at(i).v4().Eta();
    dR_fat->Fill(dR, weight);
    pt_diff_fat->Fill(DeltaPT/xconejets.at(i).v4().Pt() , weight);
    eta_diff_fat->Fill(DeltaEta/xconejets.at(i).v4().Eta(), weight);
    for(unsigned int j=0; j<xconejets.at(0).subjets().size(); j++){
      double dR_ = deltaR(xconejets.at(i).subjets().at(j), pseudoxconejets.at(i).subjets().at(j));
      double DeltaPT_ = xconejets.at(i).subjets().at(j).v4().Pt() - pseudoxconejets.at(i).subjets().at(j).v4().Pt();
      double DeltaEta_ = xconejets.at(i).subjets().at(j).v4().Eta() - pseudoxconejets.at(i).subjets().at(j).v4().Eta();
      dR_sub->Fill(dR_, weight);
      pt_diff_sub->Fill(DeltaPT_/xconejets.at(i).subjets().at(j).v4().Pt(), weight);
      eta_diff_sub->Fill(DeltaEta_/xconejets.at(i).subjets().at(j).v4().Eta(), weight);
    }
  }
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
  return;
}
