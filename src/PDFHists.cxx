#include "../include/PDFHists.h"


using namespace uhh2;

PDFHists::PDFHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){


  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700_2016v3"  ||
     ctx.get("dataset_version") == "TTbar_Mtt0700to1000_2016v3"  ||
     ctx.get("dataset_version") == "TTbar_Mtt1000toInft_2016v3"  ||
     ctx.get("dataset_version") == "TTbar_Mtt0000to0700_2L2Nu_2017" ||
     ctx.get("dataset_version") == "TTbar_Mtt0000to0700_SemiLep_2017" ||
     ctx.get("dataset_version") == "TTbar_Mtt0000to0700_Hadronic_2017") isTTbar = true;
  else isTTbar = false;

  for(int i=0; i<100; i++){
    std::stringstream ss_name1, ss_name2, ss_name3;
    ss_name1 << "M_jet1_A_PDF_" << i+1;
    ss_name2 << "M_jet1_B_PDF_" << i+1;
    ss_name3 << "M_jet1_C_PDF_" << i+1;

    std::string s_name1 = ss_name1.str();
    std::string s_name2 = ss_name2.str();
    std::string s_name3 = ss_name3.str();

    const char* name1 = s_name1.c_str();
    const char* name2 = s_name2.c_str();
    const char* name3 = s_name3.c_str();

    hist_names1.push_back(s_name1);
    hist_names2.push_back(s_name2);
    hist_names3.push_back(s_name3);


    book<TH1F>(name1, "Leading-jet m_{jet} [GeV]", 25, 0, 500);
    book<TH1F>(name2, "Leading-jet m_{jet} [GeV]", 50, 0, 500);
    book<TH1F>(name3, "Leading-jet m_{jet} [GeV]", 100, 0, 500);

    h_hadjets=ctx.get_handle<std::vector<TopJet>>("XCone33_had_Combined_Corrected");

  }
}

void PDFHists::fill(const Event & event){
  double weight = event.weight;
  std::vector<TopJet> hadjets = event.get(h_hadjets);
  TLorentzVector hadjet_v4;
  double pxhad, pyhad, pzhad, Ehad;
  pxhad = hadjets.at(0).v4().Px();
  pyhad = hadjets.at(0).v4().Py();
  pzhad = hadjets.at(0).v4().Pz();
  Ehad = hadjets.at(0).v4().E();
  hadjet_v4.SetPxPyPzE(pxhad, pyhad, pzhad, Ehad);


  if(isTTbar){    // only fill pdf hists for ttbar MC
    if(event.genInfo->systweights().size()){
      for(int i=0; i<100; i++){
        double pdf_weight = event.genInfo->systweights().at(i+9);
        double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
        const char* name1 = hist_names1[i].c_str();
        const char* name2 = hist_names2[i].c_str();
        const char* name3 = hist_names3[i].c_str();

        hist(name1)->Fill(hadjet_v4.M(),fillweight);
        hist(name2)->Fill(hadjet_v4.M(),fillweight);
        hist(name3)->Fill(hadjet_v4.M(),fillweight);
      }
    }
  }

  return;
}
