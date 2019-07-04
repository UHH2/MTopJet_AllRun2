#include <iostream>
#include <fstream>
#include <TFile.h>
#include "TUnfoldBinningXML.h"
#include <TF2.h>

using namespace std;

int main(int argc, char* argv[])
{
  cout << "create root file" << endl;

  // create root file to store binning schemes

  TString binning_root;
  TString binning_xml;

  binning_root = "Binning.root";
  binning_xml = "Binning.xml";

  TFile *binningSchemes=new TFile(binning_root,"recreate");


  cout << "set up binning" << endl;

  /******************* RECO BINNING ***********************************/
  Double_t BINS_REC_MASS[] = {123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 343, 400};
  Double_t BINS_REC_PT[] = {400, 450, 530, 10000};
  Double_t BINS_REC_MASS_PTMIGRATION[] = {100, 138, 153, 161, 167, 173, 179, 185, 192, 203, 224, 263, 400};
  Double_t BINS_REC_PT_PTMIGRATION[] = {350, 365, 380, 400};
  Double_t BINS_REC_MASS_SUBPTMIGRATION[] = {45, 90, 100, 110, 118, 125, 130, 136, 141, 147, 153, 157, 161, 164, 170, 176, 182, 188, 197, 213, 240, 293, 400};
  Double_t BINS_REC_MASS_MASSMIGRATION[] = {40, 55, 70, 85, 95, 105, 115, 125, 135, 145, 155, 160, 165, 170, 175, 180, 185, 195, 205, 215, 235, 263, 400};
  Double_t BINS_REC_MASS_BTAGMIGRATION[] = {100, 123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 400};
  Double_t BINS_REC_MASS_LEPTONPTMIGRATION[] = {45, 157, 180, 230, 400};

  // calculate total number of bins
  int N_BINS_REC_MASS = sizeof(BINS_REC_MASS)/sizeof(BINS_REC_MASS[0]) - 1;
  int N_BINS_REC_PT = sizeof(BINS_REC_PT)/sizeof(BINS_REC_PT[0]) - 1;
  int N_BINS_REC_MASS_PTMIGRATION = sizeof(BINS_REC_MASS_PTMIGRATION)/sizeof(BINS_REC_MASS_PTMIGRATION[0]) - 1;
  int N_BINS_REC_PT_PTMIGRATION = sizeof(BINS_REC_PT_PTMIGRATION)/sizeof(BINS_REC_PT_PTMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_MASSMIGRATION = sizeof(BINS_REC_MASS_MASSMIGRATION)/sizeof(BINS_REC_MASS_MASSMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_SUBPTMIGRATION = sizeof(BINS_REC_MASS_SUBPTMIGRATION)/sizeof(BINS_REC_MASS_SUBPTMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_BTAGMIGRATION = sizeof(BINS_REC_MASS_BTAGMIGRATION)/sizeof(BINS_REC_MASS_BTAGMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_LEPTONPTMIGRATION = sizeof(BINS_REC_MASS_LEPTONPTMIGRATION)/sizeof(BINS_REC_MASS_LEPTONPTMIGRATION[0]) - 1;

  /******************* GEN BINNING ************************************/
  // here some actions with arrays are needed to have different possible array sizes

  // 1. set-up binning in vectors
  std::vector<Double_t> MASS, MASS20;
  MASS20 = {112, 132, 152, 172, 192, 232}; // original binning
  MASS = {112, 132, 142, 152, 162, 172, 192, 232, 312, 392};


  std::vector<Double_t> MASS_MASSMIGRATION;
  MASS_MASSMIGRATION = {50, 75, 85, 95, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 270, 300};

  std::vector<Double_t> MASS_PTMIGRATION;
  MASS_PTMIGRATION = {112, 132, 142, 152, 162, 172, 192, 232};

  std::vector<Double_t> PT = {400, 500, 10000};
  std::vector<Double_t> PT_PTMIGRATION = {350, 380, 400};

  std::vector<Double_t> MASS_SUBPTMIGRATION;
  MASS_SUBPTMIGRATION = {92, 112, 132, 152, 172, 192, 232};

  std::vector<Double_t> MASS_LEPTONPTMIGRATION;
  MASS_LEPTONPTMIGRATION = MASS20;

  // 2. get number of bins from every vector
  int N_BINS_GEN_MASS = MASS.size() - 1;
  int N_BINS_GEN_MASS20 = MASS20.size() - 1;
  int N_BINS_GEN_PT = PT.size() - 1;
  int N_BINS_GEN_PT_PTMIGRATION = PT_PTMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_MASSMIGRATION = MASS_MASSMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_PTMIGRATION = MASS_PTMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_SUBPTMIGRATION = MASS_SUBPTMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_LEPTONPTMIGRATION = MASS_LEPTONPTMIGRATION.size() - 1;

  // 3. create arrays from vectors arrays with correct size
  Double_t *BINS_GEN_MASS = MASS.data();
  Double_t *BINS_GEN_MASS20 = MASS20.data();
  Double_t *BINS_GEN_PT = PT.data();
  Double_t *BINS_GEN_MASS_PTMIGRATION = MASS_PTMIGRATION.data();
  Double_t *BINS_GEN_PT_PTMIGRATION = PT_PTMIGRATION.data();
  Double_t *BINS_GEN_MASS_MASSMIGRATION = MASS_MASSMIGRATION.data();
  Double_t *BINS_GEN_MASS_SUBPTMIGRATION = MASS_SUBPTMIGRATION.data();
  Double_t *BINS_GEN_MASS_LEPTONPTMIGRATION = MASS_LEPTONPTMIGRATION.data();


  // =======================================================================================================
  //
  // REC BINNING
  //
  TUnfoldBinning *binning_rec=new TUnfoldBinning("binning_rec");

  //
  // define measurement phase space distribution
  //
  TUnfoldBinning *measurement_rec = binning_rec->AddBinning("measurement_rec");
  measurement_rec->AddAxis("mass",N_BINS_REC_MASS,BINS_REC_MASS,
                                true, // underflow bin
                                true // overflow bin
                                );
  measurement_rec->AddAxis("pt",N_BINS_REC_PT,BINS_REC_PT,
                                false, // no underflow bin
                                false // overflow bin
                                );

  //
  // define sideband: migration from lower pt region
  //
  TUnfoldBinning *ptmigration_rec = binning_rec->AddBinning("ptmigration_rec");
  ptmigration_rec->AddAxis("mass",N_BINS_REC_MASS_PTMIGRATION,BINS_REC_MASS_PTMIGRATION,
                                true, // underflow bin
                                true // overflow bin
                                );
  ptmigration_rec->AddAxis("pt",N_BINS_REC_PT_PTMIGRATION,BINS_REC_PT_PTMIGRATION,
                                false, // no underflow bin
                                false // overflow bin
                                );
  //
  // define sideband: migration from lower sub-pt region
  //
  TUnfoldBinning *subptmigration_rec = binning_rec->AddBinning("subptmigration_rec");
  subptmigration_rec->AddAxis("mass",N_BINS_REC_MASS_SUBPTMIGRATION,BINS_REC_MASS_SUBPTMIGRATION,
                                true, // underflow bin
                                true // overflow bin
                                );
  //
  // define sideband: migration from mass cut
  //
  TUnfoldBinning *massmigration_rec = binning_rec->AddBinning("massmigration_rec");
  massmigration_rec->AddAxis("mass",N_BINS_REC_MASS_MASSMIGRATION,BINS_REC_MASS_MASSMIGRATION,
                                true, // underflow bin
                                true // overflow bin
                                );

  //
  // define sideband: migration from lower muon-pt region
  //
  TUnfoldBinning *leptonptmigration_rec = binning_rec->AddBinning("leptonptmigration_rec");
  leptonptmigration_rec->AddAxis("mass",N_BINS_REC_MASS_LEPTONPTMIGRATION,BINS_REC_MASS_LEPTONPTMIGRATION,
                                true, // underflow bin
                                true // overflow bin
                                );
  //
  // define sideband: migration from b tagging
  //
  TUnfoldBinning *btagmigration_rec = binning_rec->AddBinning("btagmigration_rec");
  btagmigration_rec->AddAxis("mass",N_BINS_REC_MASS_BTAGMIGRATION,BINS_REC_MASS_BTAGMIGRATION,
                                true, // underflow bin
                                true // overflow bin
                                );


  // =======================================================================================================
  //
  // GEN BINNING
  //
  TUnfoldBinning *binning_gen=new TUnfoldBinning("binning_gen");

  //
  // define measurement phase space distribution
  //
  TUnfoldBinning *measurement_gen = binning_gen->AddBinning("measurement_gen");
  measurement_gen->AddAxis("mass",N_BINS_GEN_MASS,BINS_GEN_MASS,
                                true, // no underflow bin
                                true // overflow bin
                                );
  measurement_gen->AddAxis("pt",N_BINS_GEN_PT,BINS_GEN_PT,
                                false, // no underflow bin
                                false // overflow bin
                                );
  //
  // define sideband: migration from lower pt region
  //
  TUnfoldBinning *ptmigration_gen = binning_gen->AddBinning("ptmigration_gen");
  ptmigration_gen->AddAxis("mass",N_BINS_GEN_MASS_PTMIGRATION,BINS_GEN_MASS_PTMIGRATION,
                                true, // no underflow bin
                                true // overflow bin
                                );
  ptmigration_gen->AddAxis("pt",N_BINS_GEN_PT_PTMIGRATION,BINS_GEN_PT_PTMIGRATION,
                                false, // no underflow bin
                                false // overflow bin
                                );
  //
  // define sideband: migration from lower sub-pt region
  //
  TUnfoldBinning *subptmigration_gen = binning_gen->AddBinning("subptmigration_gen");
  subptmigration_gen->AddAxis("mass",N_BINS_GEN_MASS_SUBPTMIGRATION,BINS_GEN_MASS_SUBPTMIGRATION,
                                true, // no underflow bin
                                true // overflow bin
                                );
  //
  // define sideband: migration from mass cut
  //
  TUnfoldBinning *massmigration_gen = binning_gen->AddBinning("massmigration_gen");
  massmigration_gen->AddAxis("mass",N_BINS_GEN_MASS_MASSMIGRATION,BINS_GEN_MASS_MASSMIGRATION,
                                true, // no underflow bin
                                true // overflow bin
                                );

  //
  // define sideband: migration from lower muon-pt region
  //
  TUnfoldBinning *leptonptmigration_gen = binning_gen->AddBinning("leptonptmigration_gen");
  leptonptmigration_gen->AddAxis("mass",N_BINS_GEN_MASS_LEPTONPTMIGRATION,BINS_GEN_MASS_LEPTONPTMIGRATION,
                                true, // no underflow bin
                                true // overflow bin
                                );

  // =======================================================================================================
  //
  // GEN BINNING
  //
  TUnfoldBinning *binning_gen20=new TUnfoldBinning("binning_gen20");

  TUnfoldBinning *measurement_gen20 = binning_gen20->AddBinning("measurement_gen20");
  measurement_gen20->AddAxis("mass",N_BINS_GEN_MASS20,BINS_GEN_MASS20,
                                true, // no underflow bin
                                true // overflow bin
                                );
  measurement_gen20->AddAxis("pt",N_BINS_GEN_PT,BINS_GEN_PT,
                                false, // no underflow bin
                                false // overflow bin
                                );


  cout << "wirte binning scheme to root file" << endl;
  binning_rec->Write();
  binning_gen->Write();
  binning_gen20->Write();


  cout << "wirte binning scheme to readable txt file" << endl;
  // also write binning scheme as readable to a txt file
  std::ofstream out("Binning.txt");
  auto coutbuf = std::cout.rdbuf(out.rdbuf());
  cout << "BINNING REC" << endl;
  vector<TUnfoldBinning*> regions_rec = {measurement_rec, ptmigration_rec, subptmigration_rec, massmigration_rec, leptonptmigration_rec, btagmigration_rec};
  vector<TString> names_rec = {"measurement_rec", "ptmigration_rec", "subptmigration_rec", "massmigration_rec", "leptonptmigration_rec", "btagmigration_rec"};
  for(unsigned int i=0; i<regions_rec.size(); i++){
    int start = regions_rec[i]->GetStartBin();
    int end = regions_rec[i]->GetEndBin() - 1;
    cout << "  -- " << names_rec[i] << " (bins "<< start << " - " << end << ")" << endl;
  }
  cout << "-------------------------------------------" << endl;
  cout << "BINNING GEN" << endl;
  vector<TUnfoldBinning*> regions_gen = {measurement_gen, ptmigration_gen, subptmigration_gen, massmigration_gen, leptonptmigration_gen};
  vector<TString> names_gen = {"measurement_gen", "ptmigration_gen", "subptmigration_gen", "massmigration_gen", "leptonptmigration_gen"};
  for(unsigned int i=0; i<regions_gen.size(); i++){
    int start = regions_gen[i]->GetStartBin();
    int end = regions_gen[i]->GetEndBin() - 1;
    cout << "  -- " << names_gen[i] << " (bins "<< start << " - " << end << ")" << endl;
  }
  std::cout.rdbuf(coutbuf);
  ////

  cout << "wirte binning scheme to xml file" << endl;
  ofstream xmlOut(binning_xml);
  TUnfoldBinningXML::ExportXML(*binning_rec,xmlOut,kTRUE,kFALSE);
  TUnfoldBinningXML::ExportXML(*binning_gen20,xmlOut,kFALSE,kFALSE);
  TUnfoldBinningXML::ExportXML(*binning_gen,xmlOut,kFALSE,kTRUE);
  TUnfoldBinningXML::WriteDTD();
  xmlOut.close();

  delete binningSchemes;

  cout << "finished" << endl;
  return 0;
}
