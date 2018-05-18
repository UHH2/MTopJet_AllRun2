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

  /*
   ██████ ██████  ███████  █████  ████████ ███████      ██████  ██    ██ ████████ ██████  ██    ██ ████████     ███████ ██ ██      ███████
  ██      ██   ██ ██      ██   ██    ██    ██          ██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██      ██ ██      ██
  ██      ██████  █████   ███████    ██    █████       ██    ██ ██    ██    ██    ██████  ██    ██    ██        █████   ██ ██      █████
  ██      ██   ██ ██      ██   ██    ██    ██          ██    ██ ██    ██    ██    ██      ██    ██    ██        ██      ██ ██      ██
   ██████ ██   ██ ███████ ██   ██    ██    ███████      ██████   ██████     ██    ██       ██████     ██        ██      ██ ███████ ███████
  */

  std::string filename;
  if(fine) filename = "Histograms_fine.root";
  else     filename = "Histograms.root";
  outputFile=new TFile(filename.c_str(),"recreate");
  outputFile->cd();

  /*
  ██████  ███████  █████  ██████      ██████  ██ ███    ██ ███    ██ ██ ███    ██  ██████
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██ ████   ██ ████   ██ ██ ████   ██ ██
  ██████  █████   ███████ ██   ██     ██████  ██ ██ ██  ██ ██ ██  ██ ██ ██ ██  ██ ██   ███
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██ ██  ██ ██ ██  ██ ██ ██ ██  ██ ██ ██    ██
  ██   ██ ███████ ██   ██ ██████      ██████  ██ ██   ████ ██   ████ ██ ██   ████  ██████
  */

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

  /*
  ███████ ██ ██      ██          ██   ██ ██ ███████ ████████  ██████   ██████  ██████   █████  ███    ███ ███████
  ██      ██ ██      ██          ██   ██ ██ ██         ██    ██    ██ ██       ██   ██ ██   ██ ████  ████ ██
  █████   ██ ██      ██          ███████ ██ ███████    ██    ██    ██ ██   ███ ██████  ███████ ██ ████ ██ ███████
  ██      ██ ██      ██          ██   ██ ██      ██    ██    ██    ██ ██    ██ ██   ██ ██   ██ ██  ██  ██      ██
  ██      ██ ███████ ███████     ██   ██ ██ ███████    ██     ██████   ██████  ██   ██ ██   ██ ██      ██ ███████
  */

  // define directory
  TString dir = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/muon/";
  TString prefix = "/uhh2.AnalysisModuleRunner.";

  // fill data
  TFile *data_File=new TFile(dir+prefix+"DATA.DATA.root");
  fill_data((TTree *) data_File->Get("AnalysisTree"));
  // fill ttbar
  TFile *mc_matrix_File=new TFile(dir+prefix+"MC.TTbar.root");
  fill_matrix((TTree *) mc_matrix_File->Get("AnalysisTree"), "mc");
  // fill TTbar_amcatnlo
  TFile *pseudodata_File=new TFile(dir+prefix+"MC.TTbar_amcatnlo-pythia.root");
  fill_matrix((TTree *) pseudodata_File->Get("AnalysisTree"), "pseudo1");
  //fill TTbar_powheg
  TFile *pseudodata2_File=new TFile(dir+prefix+"MC.TTbar_powheg-herwig.root");
  fill_matrix((TTree *) pseudodata2_File->Get("AnalysisTree"), "pseudo2");

  // fill SYS
  vector<TString> sys_name = {"jecup", "jecdown", "jerup", "jerdown"};
  vector<TString> subdir = {"JEC_up", "JEC_down", "JER_up", "JER_down"};

  for(unsigned int i=0; i<sys_name.size(); i++){
    TFile *file = new TFile(dir+subdir[i]+prefix+"MC.TTbar.root");
    fill_matrix((TTree *) file->Get("AnalysisTree"), sys_name[i]);
    delete file;
  }

  // fill mass templates
  std::vector<TFile *> mc_truth_File;
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1665.root"));
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1695_ext2.root"));
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1715.root"));
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1735.root"));
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1755.root"));
  mc_truth_File.push_back(new TFile(dir+prefix+"MC.TTbar_mtop1785.root"));

  fill_template((TTree *) mc_truth_File[0]->Get("AnalysisTree"), "1665");
  fill_template((TTree *) mc_truth_File[1]->Get("AnalysisTree"), "1695");
  fill_template((TTree *) mc_truth_File[2]->Get("AnalysisTree"), "1715");
  fill_template((TTree *) mc_truth_File[3]->Get("AnalysisTree"), "1735");
  fill_template((TTree *) mc_truth_File[4]->Get("AnalysisTree"), "1755");
  fill_template((TTree *) mc_truth_File[5]->Get("AnalysisTree"), "1785");

  // fill background
  vector<TString> bkg_name = {"WJets", "SingleTop", "other"};
  for(unsigned int i=0; i<bkg_name.size(); i++){
    TFile *file = new TFile(dir+prefix+"MC."+bkg_name[i]+".root");
    fill_background((TTree *) file->Get("AnalysisTree"), bkg_name[i]);
    delete file;
  }

  return 0;
}


/*
███████ ██ ██      ██          ██████   █████  ████████  █████
██      ██ ██      ██          ██   ██ ██   ██    ██    ██   ██
█████   ██ ██      ██          ██   ██ ███████    ██    ███████
██      ██ ██      ██          ██   ██ ██   ██    ██    ██   ██
██      ██ ███████ ███████     ██████  ██   ██    ██    ██   ██
*/



void fill_data(TTree* tree){
  if(!tree) cout << "could not read 'data' tree\n";
  else      cout << "Filling Data Histograms...\n";

  TH1* h_data = binning_rec->CreateHistogram("data");

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

  h_data->Write();
  delete h_data;
  return;
}

/*
███████ ██ ██      ██          ███    ███  █████  ███████ ███████     ████████ ███████ ███    ███ ██████
██      ██ ██      ██          ████  ████ ██   ██ ██      ██             ██    ██      ████  ████ ██   ██
█████   ██ ██      ██          ██ ████ ██ ███████ ███████ ███████        ██    █████   ██ ████ ██ ██████
██      ██ ██      ██          ██  ██  ██ ██   ██      ██      ██        ██    ██      ██  ██  ██ ██
██      ██ ███████ ███████     ██      ██ ██   ██ ███████ ███████        ██    ███████ ██      ██ ██
*/

void fill_template(TTree* tree, TString mtop){
  if(!tree) cout << "could not read 'data' tree\n";
  else      cout << "Filling Template Histograms...\n";

  TH1* h_mass_truth = measurement_gen->CreateHistogram("mc_mtop"+mtop+"_truth",kTRUE,0,0,"pt[C]");
  TH1* h_mass_rec = measurement_rec->CreateHistogram("mc_mtop"+mtop+"_rec",kTRUE,0,0,"pt[C]");

  outputFile->cd();

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Gen33",&massGen);
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("passed_ptmigration_rec",&passed_ptmigration_rec);
  tree->SetBranchAddress("passed_massmigration_rec",&passed_massmigration_rec);
  tree->SetBranchAddress("passed_btagmigration_rec",&passed_btagmigration_rec);
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
      h_mass_truth->Fill(massGen, gen_weight);
    }
    if(passed_measurement_rec || passed_ptmigration_rec || passed_massmigration_rec || passed_btagmigration_rec){
      h_mass_rec->Fill(massRec, w_rec);
    }
  }
  h_mass_truth->Write();
  h_mass_rec->Write();
  delete h_mass_truth;
  delete h_mass_rec;
  return;
}

/*
███████ ██ ██      ██          ██████   █████   ██████ ██   ██  ██████  ██████   ██████  ██    ██ ███    ██ ██████
██      ██ ██      ██          ██   ██ ██   ██ ██      ██  ██  ██       ██   ██ ██    ██ ██    ██ ████   ██ ██   ██
█████   ██ ██      ██          ██████  ███████ ██      █████   ██   ███ ██████  ██    ██ ██    ██ ██ ██  ██ ██   ██
██      ██ ██      ██          ██   ██ ██   ██ ██      ██  ██  ██    ██ ██   ██ ██    ██ ██    ██ ██  ██ ██ ██   ██
██      ██ ███████ ███████     ██████  ██   ██  ██████ ██   ██  ██████  ██   ██  ██████   ██████  ██   ████ ██████
*/



void fill_background(TTree *tree, TString prefix){
  if(!tree) cout << "could not read 'mc bgr' tree\n";
  else      cout << "Filling Background Histograms...\n";

  TH1* h_mc_bgr = binning_rec->CreateHistogram("BKG_"+prefix);

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
    }
  }
  h_mc_bgr->Write();
  delete h_mc_bgr;
  return;
}


/*
███████ ██ ██      ██          ███    ███  █████  ████████ ██████  ██ ██   ██
██      ██ ██      ██          ████  ████ ██   ██    ██    ██   ██ ██  ██ ██
█████   ██ ██      ██          ██ ████ ██ ███████    ██    ██████  ██   ███
██      ██ ██      ██          ██  ██  ██ ██   ██    ██    ██   ██ ██  ██ ██
██      ██ ███████ ███████     ██      ██ ██   ██    ██    ██   ██ ██ ██   ██
*/

void fill_matrix(TTree* tree, TString prefix){
  if(!tree) cout<<"could not read 'mc signal' tree\n";
  else      cout << "Filling Matrix Histograms ("+prefix+") ...\n";

  // setup hists
  TH1* h_purity_all = measurement_gen->CreateHistogram(prefix + "_purity_all",kTRUE,0,0,"pt[C]");
  TH1* h_purity_samebin = measurement_gen->CreateHistogram(prefix + "_purity_samebin",kTRUE,0,0,"pt[C]");
  TH1* h_purity_samebin_pt = measurement_gen->CreateHistogram(prefix + "_purity_samebin_pt",kTRUE,0,0,"pt[C]");

  TH1* h_mc_sig = binning_rec->CreateHistogram(prefix + "_sig");
  TH1* h_mc_gen = binning_gen->CreateHistogram(prefix + "_gen");
  TH1* h_mc_truth = measurement_gen->CreateHistogram(prefix + "_truth",kTRUE,0,0,"pt[C]");
  TH2* h_mc_matrix = TUnfoldBinning::CreateHistogramOfMigrations(binning_gen,binning_rec, prefix + "_matrix");

  outputFile->cd();
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
  h_purity_all->Write();
  h_purity_samebin->Write();
  h_purity_samebin_pt->Write();
  h_mc_sig->Write();
  h_mc_gen->Write();
  h_mc_truth->Write();
  h_mc_matrix->Write();

  delete h_purity_all;
  delete h_purity_samebin;
  delete h_purity_samebin_pt;
  delete h_mc_sig;
  delete h_mc_gen;
  delete h_mc_truth;
  delete h_mc_matrix;

  return;
}
