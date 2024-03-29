#include <UHH2/MTopJet/include/CutHists.h>


CutHists::CutHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // book all histograms here
  h_cut = book<TH1F>("cut", "a.u.", 1, -0.5, 0.5);
  // h_cut_postsel = book<TH1F>("cut_postsel", "a.u.", 4, -0.5, 3.5);

  // h_hadjets=ctx.get_handle<std::vector<Jet>>("XCone33_had_Combined");
  // h_lepjets=ctx.get_handle<std::vector<Jet>>("XCone33_lep_Combined");

}



void CutHists::fill(const Event & event){

  /* "cut" variable gives the hist which should be filled

     0 = events after all cuts
     1 = events after all cuts but 2D
     2 = events after all cuts but met
     3 = events after all cuts but btag

     fill N-1 hist for each cut and hist with all cuts
  */


  //---------------------------------------------------------------------------------------
  //--------------------------------- get jets and lepton ---------------------------------
  //---------------------------------------------------------------------------------------
  // std::vector<Jet> hadjets = event.get(h_hadjets);
  // std::vector<Jet> lepjets = event.get(h_lepjets);

  //---------------------------------------------------------------------------------------
  //-------- set Lorentz Vectors of Combined Jets ---------- ------------------------------
  //---------------------------------------------------------------------------------------
  // TLorentzVector hadjet_v4, lepjet_v4;
  // double pxlep, pylep, pzlep, Elep;
  // pxlep = lepjets.at(0).v4().Px();
  // pylep = lepjets.at(0).v4().Py();
  // pzlep = lepjets.at(0).v4().Pz();
  // Elep = lepjets.at(0).v4().E();
  // lepjet_v4.SetPxPyPzE(pxlep, pylep, pzlep, Elep);

  // double pxhad, pyhad, pzhad, Ehad;
  // pxhad = hadjets.at(0).v4().Px();
  // pyhad = hadjets.at(0).v4().Py();
  // pzhad = hadjets.at(0).v4().Pz();
  // Ehad = hadjets.at(0).v4().E();
  // hadjet_v4.SetPxPyPzE(pxhad, pyhad, pzhad, Ehad);
  //---------------------------------------------------------------------------------------
  //---------------------------------- set up Post Sel ------------------------------------
  //---------------------------------------------------------------------------------------
  // bool passed_post = false;
  // if(hadjet_v4.Pt() > 400 && hadjet_v4.M() > lepjet_v4.M()) passed_post = true;

  //---------------------------------------------------------------------------------------
  //--------------------------------- Fill Hists here -------------------------------------
  //---------------------------------------------------------------------------------------

  // get weight
  double weight = event.weight;
  int bin = 0;
  h_cut->Fill(bin, weight);

  //---------------------------------------------------------------------------------------
  //--------------------------------- Clear all used objects ------------------------------
  //---------------------------------------------------------------------------------------
  // hadjet_v4.Delete();
  // lepjet_v4.Delete();
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------

}
