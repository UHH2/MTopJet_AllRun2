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
  bool fine;

  if(strcmp(argv[1], "fine") == 0) fine = true;
  else fine = false;

  if(fine){
    cout << "you selected fine binning" << endl;
    binning_root = "Binning_fine.root";
    binning_xml = "Binning_fine.xml";
  }
  else{
    binning_root = "Binning.root";
    binning_xml = "Binning.xml";
  }

  TFile *binningSchemes=new TFile(binning_root,"recreate");


  cout << "set up binning" << endl;

  /******************* RECO BINNING ***********************************/
  Double_t BINS_REC_MASS[] = {123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 343, 400};
  Double_t BINS_REC_PT[] = {400, 450, 530, 10000};
  Double_t BINS_REC_MASS_PTMIGRATION[] = {100, 123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 400};
  Double_t BINS_REC_PT_PTMIGRATION[] = {350, 380, 400};
  Double_t BINS_REC_MASS_MASSMIGRATION[] = {45, 70, 90, 100, 110, 125, 130, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 400};
  Double_t BINS_REC_MASS_SUBPTMIGRATION[] = {100, 123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 400};
  Double_t BINS_REC_MASS_BTAGMIGRATION[] = {100, 123, 138, 147, 153, 157, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 192, 197, 203, 213, 224, 240, 263, 293, 400};

  // calculate total number of bins
  int N_BINS_REC_MASS = sizeof(BINS_REC_MASS)/sizeof(BINS_REC_MASS[0]) - 1;
  int N_BINS_REC_PT = sizeof(BINS_REC_PT)/sizeof(BINS_REC_PT[0]) - 1;
  int N_BINS_REC_MASS_PTMIGRATION = sizeof(BINS_REC_MASS_PTMIGRATION)/sizeof(BINS_REC_MASS_PTMIGRATION[0]) - 1;
  int N_BINS_REC_PT_PTMIGRATION = sizeof(BINS_REC_PT_PTMIGRATION)/sizeof(BINS_REC_PT_PTMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_MASSMIGRATION = sizeof(BINS_REC_MASS_MASSMIGRATION)/sizeof(BINS_REC_MASS_MASSMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_SUBPTMIGRATION = sizeof(BINS_REC_MASS_SUBPTMIGRATION)/sizeof(BINS_REC_MASS_SUBPTMIGRATION[0]) - 1;
  int N_BINS_REC_MASS_BTAGMIGRATION = sizeof(BINS_REC_MASS_BTAGMIGRATION)/sizeof(BINS_REC_MASS_BTAGMIGRATION[0]) - 1;

  /******************* GEN BINNING ************************************/
  // here some actions with arrays are needed to have different possible array sizes

  // 1. set-up binning in vectors
  std::vector<Double_t> MASS;
  if(fine) MASS = {117, 137, 157, 172, 187, 222};
  // if(fine) MASS = {117, 137, 157, 172, 187, 202, 222};
  // if(fine) MASS = {110, 137, 157, 172, 187, 220, 250};
  else     MASS = {120, 150, 180, 210, 240, 270, 300};

  std::vector<Double_t> MASS_MASSMIGRATION;
  if(fine) MASS_MASSMIGRATION = {50, 120, 135, 150, 165, 180, 195, 210, 225, 240, 270, 300};
  else     MASS_MASSMIGRATION = {50, 90, 120, 150, 180, 210, 240, 270, 300};

  std::vector<Double_t> MASS_PTMIGRATION;
  MASS_PTMIGRATION = MASS;

  std::vector<Double_t> PT = {400, 500, 10000};
  std::vector<Double_t> PT_PTMIGRATION = {350, 400};

  std::vector<Double_t> MASS_SUBPTMIGRATION;
  MASS_SUBPTMIGRATION = MASS;

  // 2. get number of bins from every vector
  int N_BINS_GEN_MASS = MASS.size() - 1;
  int N_BINS_GEN_PT = PT.size() - 1;
  int N_BINS_GEN_PT_PTMIGRATION = PT_PTMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_MASSMIGRATION = MASS_MASSMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_PTMIGRATION = MASS_PTMIGRATION.size() - 1;
  int N_BINS_GEN_MASS_SUBPTMIGRATION = MASS_SUBPTMIGRATION.size() - 1;

  // 3. create arrays from vectors arrays with correct size
  Double_t *BINS_GEN_MASS = MASS.data();
  Double_t *BINS_GEN_PT = PT.data();
  Double_t *BINS_GEN_MASS_PTMIGRATION = MASS_PTMIGRATION.data();
  Double_t *BINS_GEN_PT_PTMIGRATION = PT_PTMIGRATION.data();
  Double_t *BINS_GEN_MASS_MASSMIGRATION = MASS_MASSMIGRATION.data();
  Double_t *BINS_GEN_MASS_SUBPTMIGRATION = MASS_SUBPTMIGRATION.data();


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


  cout << "wirte binning scheme to file" << endl;

  binning_rec->Write();
  binning_gen->Write();

  ofstream xmlOut(binning_xml);
  TUnfoldBinningXML::ExportXML(*binning_rec,xmlOut,kTRUE,kFALSE);
  TUnfoldBinningXML::ExportXML(*binning_gen,xmlOut,kFALSE,kTRUE);
  TUnfoldBinningXML::WriteDTD();
  xmlOut.close();

  delete binningSchemes;

  cout << "finished" << endl;
  return 0;
}
