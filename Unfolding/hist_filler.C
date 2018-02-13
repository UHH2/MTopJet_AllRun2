#include "hist_filler.h"  

int main(int argc, char* argv[]){

 
  bool fine;
  if(strcmp(argv[1], "fine") == 0) fine = true;
  else fine = false;
  if(fine) cout << "The fine binning is selected" << endl;
  else     cout << "The large binning is selected" << endl;


 
  if(strcmp(argv[2], "ttw") == 0) ttweight = true;
  else ttweight = false;
  if(ttweight) cout << "ttbar reweighting is applied" << endl;
  else         cout << "ttbar reweighting is NOT applied" << endl;

  // switch on histogram errors
  TH1::SetDefaultSumw2();

  //=======================================================
  // Step 1: open file to save histograms and binning schemes

  std::string filename;
  if(fine) filename = "Histograms_fine.root";
  else     filename = "Histograms.root";
  outputFile=new TFile(filename.c_str(),"recreate");


 //=======================================================
  // Step 2: read binning from XML
  //         and save them to output file


  outputFile->cd();

  // read binning schemes in XML format
  TDOMParser parser;
  TString binning_xml;
  if(fine) binning_xml = "Binning_fine.xml";
  else     binning_xml = "Binning.xml";

  Int_t error=parser.ParseFile(binning_xml);
  if(error) cout<<"error="<<error<<" from TDOMParser\n";
  TXMLDocument const *XMLdocument=parser.GetXMLDocument();
  binning_rec = TUnfoldBinningXML::ImportXML(XMLdocument,"binning_rec");
  binning_gen = TUnfoldBinningXML::ImportXML(XMLdocument,"binning_gen");


  if(!binning_rec) cout<<"could not read 'rec' binning\n";
  if(!binning_gen) cout<<"could not read 'gen' binning\n";


  // get distributions from measurement phase space and sideband regions
  measurement_rec = binning_rec->FindNode("measurement_rec");
  measurement_gen = binning_gen->FindNode("measurement_gen");

  ptmigration_rec = binning_rec->FindNode("ptmigration_rec");
  ptmigration_gen = binning_gen->FindNode("ptmigration_gen");

  massmigration_rec = binning_rec->FindNode("massmigration_rec");
  massmigration_gen = binning_gen->FindNode("massmigration_gen");

  btagmigration_rec = binning_rec->FindNode("btagmigration_rec");


  //=======================================================
  // create hists from binning
  create_hists();

  //=======================================================
  // fill histograms 

  // define directory
  TString prefix = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/muon/uhh2.AnalysisModuleRunner.";

  TFile *data_File=new TFile(prefix+"DATA.DATA.root");
  fill_data((TTree *) data_File->Get("AnalysisTree"));
  delete data_File;

  // TFile *pseudodata_File=new TFile(prefix+"MC.TTbar.root");
  // TTree *pseudodata_Tree=(TTree *) pseudodata_File->Get("AnalysisTree");


  TFile *mc_matrix_File=new TFile(prefix+"MC.TTbar.root");
  fill_matrix((TTree *) mc_matrix_File->Get("AnalysisTree"));
  delete mc_matrix_File;


  // fill templates here
  std::vector<TFile *> mc_truth_File;
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1665.root"));
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1695_ext2.root"));
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1715.root"));
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1735.root"));
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1755.root"));
  mc_truth_File.push_back(new TFile(prefix+"MC.TTbar_mtop1785.root"));

  fill_template((TTree *) mc_truth_File[0]->Get("AnalysisTree"), "1665");
  fill_template((TTree *) mc_truth_File[1]->Get("AnalysisTree"), "1695");
  fill_template((TTree *) mc_truth_File[2]->Get("AnalysisTree"), "1715");
  fill_template((TTree *) mc_truth_File[3]->Get("AnalysisTree"), "1735");
  fill_template((TTree *) mc_truth_File[4]->Get("AnalysisTree"), "1755");
  fill_template((TTree *) mc_truth_File[5]->Get("AnalysisTree"), "1785");

 
  TFile *mc_bgr_File=new TFile(prefix+"MC.Background_only.root");
  fill_background((TTree *) mc_bgr_File->Get("AnalysisTree"));
  delete mc_bgr_File;

  //=======================================================
  // write file

  outputFile->Write();
  delete outputFile;

  return 0;
}


//=======================================================
//================ fill data ============================
//=======================================================

void fill_data(TTree* tree){
  if(!tree) cout << "could not read 'data' tree\n";
  else      cout << "Filling Data Histograms...\n";

  outputFile->cd();

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("Pt_Rec",&ptRec);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("passed_ptmigration_rec",&passed_ptmigration_rec);
  tree->SetBranchAddress("passed_massmigration_rec",&passed_massmigration_rec);
  tree->SetBranchAddress("passed_btagmigration_rec",&passed_btagmigration_rec);

  tree->SetBranchStatus("*",1);

  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
     if(tree->GetEntry(ievent)<=0) break;

     Int_t binNumber = 0;
     if     (passed_measurement_rec)    binNumber = measurement_rec->GetGlobalBinNumber(massRec,ptRec);
     else if(passed_ptmigration_rec)    binNumber = ptmigration_rec->GetGlobalBinNumber(massRec,ptRec);
     else if(passed_massmigration_rec)  binNumber = massmigration_rec->GetGlobalBinNumber(massRec);
     else if(passed_btagmigration_rec)  binNumber = btagmigration_rec->GetGlobalBinNumber(massRec);


     if(passed_measurement_rec || passed_ptmigration_rec || passed_massmigration_rec || passed_btagmigration_rec){
       h_data->Fill(binNumber);
     }
  }
}

//=======================================================
//================ fill template ========================
//=======================================================

void fill_template(TTree* tree, TString mtop){
  if(!tree) cout << "could not read 'data' tree\n";
  else      cout << "Filling Template Histograms...\n";

  // create hist

  outputFile->cd();

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Gen33",&massGen);
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("passed_measurement_gen",&passed_measurement_gen);
  tree->SetBranchAddress("is_TTbar",&is_TTbar);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("gen_weight_ttfactor",&gen_ttfactor);
  tree->SetBranchStatus("*",1);
 
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(!is_TTbar){
      cout << "ERROR: non-ttbar event in signal sample!" << endl;
      break;
    }
 
    if(ttweight) gen_weight *= gen_ttfactor;
    w_rec = rec_weight * gen_weight;

    if(passed_measurement_gen){
      if(mtop == "1665") h_mc_mtop1665_truth->Fill(massGen, gen_weight);
      if(mtop == "1695") h_mc_mtop1695_truth->Fill(massGen, gen_weight);
      if(mtop == "1715") h_mc_mtop1715_truth->Fill(massGen, gen_weight);
      if(mtop == "1735") h_mc_mtop1735_truth->Fill(massGen, gen_weight);
      if(mtop == "1755") h_mc_mtop1755_truth->Fill(massGen, gen_weight);
      if(mtop == "1785") h_mc_mtop1785_truth->Fill(massGen, gen_weight);
      if(mtop == "1665") h_mc_mtop1665_rec->Fill(massRec, w_rec);
      if(mtop == "1695") h_mc_mtop1695_rec->Fill(massRec, w_rec);
      if(mtop == "1715") h_mc_mtop1715_rec->Fill(massRec, w_rec);
      if(mtop == "1735") h_mc_mtop1735_rec->Fill(massRec, w_rec);
      if(mtop == "1755") h_mc_mtop1755_rec->Fill(massRec, w_rec);
      if(mtop == "1785") h_mc_mtop1785_rec->Fill(massRec, w_rec);
    }
  }
}


//=======================================================
//================ fill background ======================
//=======================================================

void fill_background(TTree *tree){
  if(!tree) cout << "could not read 'mc bgr' tree\n";
  else      cout << "Filling Background Histograms...\n";

  outputFile->cd();



  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("Pt_Rec",&ptRec);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("gen_weight_ttfactor",&gen_ttfactor);
  tree->SetBranchAddress("is_TTbar",&is_TTbar);
  tree->SetBranchAddress("passed_ptmigration_rec",&passed_ptmigration_rec);
  tree->SetBranchAddress("passed_massmigration_rec",&passed_massmigration_rec);
  tree->SetBranchAddress("passed_btagmigration_rec",&passed_btagmigration_rec);
  tree->SetBranchStatus("*",1);



  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(is_TTbar){
      cout << "ERROR: ttbar event in background sample!" << endl;
      break;
    }

    if(ttweight) gen_weight *= gen_ttfactor;
    w_bgr_rec = rec_weight * gen_weight;

    Int_t recBin = 0;
    if     (passed_measurement_rec)    recBin = measurement_rec->GetGlobalBinNumber(massRec,ptRec);
    else if(passed_ptmigration_rec)    recBin = ptmigration_rec->GetGlobalBinNumber(massRec,ptRec);
    else if(passed_massmigration_rec)  recBin = massmigration_rec->GetGlobalBinNumber(massRec);
    else if(passed_btagmigration_rec)  recBin = btagmigration_rec->GetGlobalBinNumber(massRec);

    if(passed_measurement_rec || passed_ptmigration_rec || passed_massmigration_rec || passed_btagmigration_rec){
      h_mc_bgr->Fill(recBin, w_bgr_rec);
      h_mc_rec->Fill(recBin, w_bgr_rec);
    }
  }
}


//=======================================================
//================ fill signal ==========================
//=======================================================

void fill_matrix(TTree* tree){
  if(!tree) cout<<"could not read 'mc signal' tree\n";
  else      cout << "Filling Matrix Histograms...\n";

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("Mass_Gen33",&massGen);
  tree->SetBranchAddress("Pt_Rec",&ptRec);
  tree->SetBranchAddress("Pt_Gen33",&ptGen);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("passed_measurement_gen",&passed_measurement_gen);
  tree->SetBranchAddress("is_TTbar",&is_TTbar);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("gen_weight_ttfactor",&gen_ttfactor);
  tree->SetBranchAddress("passed_ptmigration_rec",&passed_ptmigration_rec);
  tree->SetBranchAddress("passed_ptmigration_gen",&passed_ptmigration_gen);
  tree->SetBranchAddress("passed_massmigration_rec",&passed_massmigration_rec);
  tree->SetBranchAddress("passed_massmigration_gen",&passed_massmigration_gen);
  tree->SetBranchAddress("passed_btagmigration_rec",&passed_btagmigration_rec);
  tree->SetBranchStatus("*",1);
 

  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(!is_TTbar){
      cout << "ERROR: non-ttbar event in signal sample!" << endl;
      break;
    }

    if(ttweight) gen_weight *= gen_ttfactor;

     // get weights for migration matrix
    w_central = rec_weight * gen_weight;
    w_nogen = rec_weight * gen_weight;
    w_norec = gen_weight;
    w_correction = gen_weight * (1 - rec_weight);

    // get weight for gen and rec hists
    w_sig_rec = rec_weight * gen_weight;
    w_gen = gen_weight;


    // get global bins
    Int_t genBin;
    if     (passed_measurement_gen)   genBin = measurement_gen->GetGlobalBinNumber(massGen,ptGen);
    else if(passed_ptmigration_gen)   genBin = ptmigration_gen->GetGlobalBinNumber(massGen,ptGen);
    else if(passed_massmigration_gen) genBin = massmigration_gen->GetGlobalBinNumber(massGen);
    else                              genBin = 0;

    Int_t recBin;
    if     (passed_measurement_rec)    recBin = measurement_rec->GetGlobalBinNumber(massRec,ptRec);
    else if(passed_ptmigration_rec)    recBin = ptmigration_rec->GetGlobalBinNumber(massRec,ptRec);
    else if(passed_massmigration_rec)  recBin = massmigration_rec->GetGlobalBinNumber(massRec);
    else if(passed_btagmigration_rec)  recBin = btagmigration_rec->GetGlobalBinNumber(massRec);
    else                               recBin = 0;

    bool gen_info = false;
    bool rec_info = false;
    if(passed_ptmigration_gen || passed_measurement_gen || passed_massmigration_gen) gen_info = true;
    if(passed_ptmigration_rec || passed_measurement_rec || passed_massmigration_rec || passed_btagmigration_rec) rec_info = true;


    // Fill 1D Hists (for truth distribution only use gen selection)
    if(passed_measurement_gen) h_mc_truth->Fill(massGen, w_gen);

    if(gen_info) h_mc_gen->Fill(genBin, w_gen);
    if(rec_info){
      h_mc_sig->Fill(recBin, w_sig_rec);
      h_mc_rec->Fill(recBin, w_sig_rec);
    }

    // Fill Matrix
    if( rec_info &&  gen_info) h_mc_matrix->Fill(genBin, recBin, w_central);
    if( rec_info && !gen_info) h_mc_matrix->Fill(genBin, recBin, w_nogen);
    if(!rec_info &&  gen_info) h_mc_matrix->Fill(genBin, recBin, w_norec);
    if( rec_info &&  gen_info) h_mc_matrix->Fill(genBin,     0., w_correction); // this is needed because events that dont pass rec, have no rec_weight. 

    //fill hists for purity
    int genBin_recInfo = 0;
    int genBin_recInfo_pt = 0;
    if(passed_measurement_gen && rec_info){
      genBin_recInfo = measurement_gen->GetGlobalBinNumber(massRec,ptGen);
      genBin_recInfo_pt = measurement_gen->GetGlobalBinNumber(massRec,ptRec);
      
      h_purity_all->Fill(massGen, w_gen);
      if(genBin_recInfo == genBin) h_purity_samebin->Fill(massGen, w_gen);
      if(genBin_recInfo_pt == genBin) h_purity_samebin_pt->Fill(massGen, w_gen);
    }
  }
}




void create_hists(){
  h_purity_all = measurement_gen->CreateHistogram("purity_all",kTRUE,0,0,"pt[C]");
  h_purity_samebin = measurement_gen->CreateHistogram("purity_samebin",kTRUE,0,0,"pt[C]");
  h_purity_samebin_pt = measurement_gen->CreateHistogram("purity_samebin_pt",kTRUE,0,0,"pt[C]");


  h_data = binning_rec->CreateHistogram("data");

  h_mc_rec = binning_rec->CreateHistogram("mc_rec");
  h_mc_bgr = binning_rec->CreateHistogram("mc_bgr");
  h_mc_sig = binning_rec->CreateHistogram("mc_sig");
  h_mc_gen = binning_gen->CreateHistogram("mc_gen");
  h_mc_truth = measurement_gen->CreateHistogram("mc_truth",kTRUE,0,0,"pt[C]");

  h_mc_matrix=TUnfoldBinning::CreateHistogramOfMigrations(binning_gen,binning_rec,"mc_matrix");

  h_mc_mtop1665_truth = measurement_gen->CreateHistogram("mc_mtop1665_truth",kTRUE,0,0,"pt[C]");
  h_mc_mtop1695_truth = measurement_gen->CreateHistogram("mc_mtop1695_truth",kTRUE,0,0,"pt[C]");
  h_mc_mtop1715_truth = measurement_gen->CreateHistogram("mc_mtop1715_truth",kTRUE,0,0,"pt[C]");
  h_mc_mtop1735_truth = measurement_gen->CreateHistogram("mc_mtop1735_truth",kTRUE,0,0,"pt[C]");
  h_mc_mtop1755_truth = measurement_gen->CreateHistogram("mc_mtop1755_truth",kTRUE,0,0,"pt[C]");
  h_mc_mtop1785_truth = measurement_gen->CreateHistogram("mc_mtop1785_truth",kTRUE,0,0,"pt[C]");

  h_mc_mtop1665_rec = measurement_rec->CreateHistogram("mc_mtop1665_rec",kTRUE,0,0,"pt[C]");
  h_mc_mtop1695_rec = measurement_rec->CreateHistogram("mc_mtop1695_rec",kTRUE,0,0,"pt[C]");
  h_mc_mtop1715_rec = measurement_rec->CreateHistogram("mc_mtop1715_rec",kTRUE,0,0,"pt[C]");
  h_mc_mtop1735_rec = measurement_rec->CreateHistogram("mc_mtop1735_rec",kTRUE,0,0,"pt[C]");
  h_mc_mtop1755_rec = measurement_rec->CreateHistogram("mc_mtop1755_rec",kTRUE,0,0,"pt[C]");
  h_mc_mtop1785_rec = measurement_rec->CreateHistogram("mc_mtop1785_rec",kTRUE,0,0,"pt[C]");
}
