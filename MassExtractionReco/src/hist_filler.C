#include "../include/hist_filler.h"

int main(int argc, char* argv[]){


  // switch on histogram errors
  TH1::SetDefaultSumw2();

  /*
   ██████ ██████  ███████  █████  ████████ ███████      ██████  ██    ██ ████████ ██████  ██    ██ ████████     ███████ ██ ██      ███████
  ██      ██   ██ ██      ██   ██    ██    ██          ██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██      ██ ██      ██
  ██      ██████  █████   ███████    ██    █████       ██    ██ ██    ██    ██    ██████  ██    ██    ██        █████   ██ ██      █████
  ██      ██   ██ ██      ██   ██    ██    ██          ██    ██ ██    ██    ██    ██      ██    ██    ██        ██      ██ ██      ██
   ██████ ██   ██ ███████ ██   ██    ██    ███████      ██████   ██████     ██    ██       ██████     ██        ██      ██ ███████ ███████
  */

  std::string filename = "Histograms.root";

  outputFile=new TFile(filename.c_str(),"recreate");
  outputFile->cd();

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
  fill_hist((TTree *) data_File->Get("AnalysisTree"), "data");
  // fill pseudo data
  TFile *pseudo1_File = new TFile(dir+prefix+"MC.TTbar_amcatnlo-pythia.root");
  fill_hist((TTree *) pseudo1_File->Get("AnalysisTree"), "pseudo1");
  // fill ttbar
  TFile *mc_matrix_File=new TFile(dir+prefix+"MC.TTbar.root");
  fill_hist((TTree *) mc_matrix_File->Get("AnalysisTree"), "tt");
  // fill background
  vector<TString> bkg_name = {"WJets", "SingleTop", "other"};
  for(unsigned int i=0; i<bkg_name.size(); i++){
    TFile *file = new TFile(dir+prefix+"MC."+bkg_name[i]+".root");
    fill_hist((TTree *) file->Get("AnalysisTree"), bkg_name[i]);
    delete file;
  }

  // fill SYS
  vector<TString> sys_name = {"btagbcup", "btagbcdown", "btagudsgup", "btagudsgdown", "jecup", "jecdown", "jerup", "jerdown", "corup", "cordown", "muidup", "muiddown", "mutrup", "mutrdown", "puup", "pudown"};
  vector<TString> subdir = {"BTAG_bc_up", "BTAG_bc_down", "BTAG_udsg_up", "BTAG_udsg_down", "JEC_up", "JEC_down", "JER_up", "JER_down", "COR_up", "COR_down", "MUID_up", "MUID_down", "MUTR_up", "MUTR_down", "PU_up", "PU_down"};

  for(unsigned int i=0; i<sys_name.size(); i++){
    TFile *file = new TFile(dir+subdir[i]+prefix+"MC.TTbar.root");
    fill_hist((TTree *) file->Get("AnalysisTree"), sys_name[i]);
    delete file;
  }

  // fill model SYS
  std::vector<TString> model_name = {"Shower", "Generator", "SCALE_upup", "SCALE_upnone", "SCALE_noneup", "SCALE_downdown", "SCALE_downnone", "SCALE_nonedown"};
  for(unsigned int i=0; i<model_name.size(); i++){
    TFile * file;
    if(model_name[i] == "Shower")          file = new TFile(dir+prefix+"MC.TTbar_powheg-herwig.root");
    else if (model_name[i] == "Generator") file = new TFile(dir+prefix+"MC.TTbar_amcatnlo-pythia.root");
    else                                   file = new TFile(dir+model_name[i]+prefix+"MC.TTbar.root");
    fill_hist((TTree *) file->Get("AnalysisTree"), model_name[i]);
  }

  // fill PDF
  fill_pdf( (TTree *) mc_matrix_File->Get("AnalysisTree") );

  // fill mass templates
  std::vector<TFile *> mc_truth_File;
  std::vector<TString> masses = {"1665", "1695", "1715", "1735", "1755", "1785"};
  for(auto mass: masses){
    TFile *file = new TFile(dir+prefix+"MC.TTbar_mtop"+mass+".root");
    fill_hist((TTree *) file->Get("AnalysisTree"), mass);
    delete file;
  }

  outputFile->Close();

  return 0;
}


/*
███████ ██ ██      ██          ██   ██ ██ ███████ ████████
██      ██ ██      ██          ██   ██ ██ ██         ██
█████   ██ ██      ██          ███████ ██ ███████    ██
██      ██ ██      ██          ██   ██ ██      ██    ██
██      ██ ███████ ███████     ██   ██ ██ ███████    ██
*/

void fill_hist(TTree* tree, TString prefix){
  if(!tree) cout << "could not read tree\n";
  else      cout << "Filling " << prefix << " Histograms...\n";

  int nbins = bins.size()-1;
  TH1F* h_data = new TH1F(prefix, "jet mass", nbins, &bins[0]);

  outputFile->cd();

  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchStatus("*",1);

  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;

    double weight = rec_weight * gen_weight;
    if(passed_measurement_rec) h_data->Fill(massRec, weight);
  }

  h_data->Write();
  delete h_data;
  return;
}



/*
███████ ██ ██      ██          ██████  ██████  ███████
██      ██ ██      ██          ██   ██ ██   ██ ██
█████   ██ ██      ██          ██████  ██   ██ █████
██      ██ ██      ██          ██      ██   ██ ██
██      ██ ███████ ███████     ██      ██████  ██
*/



void fill_pdf(TTree* tree){
  if(!tree) cout<<"could not read 'mc signal' tree\n";
  else      cout << "Filling Model SYS Histograms (PDF) ...\n";


  // setup hist names
  vector<TString> histnames;
  for(unsigned int i=0; i<100; i++){
    TString histname = "PDF_";
    histname += i;
    histnames.push_back(histname);
  }

  // setup histograms
  vector<TH1F*> hists;
  int nbins = bins.size()-1;
  TH1F* hist = new TH1F("pdf", "jet mass", nbins, &bins[0]);
  for(unsigned int i=0; i<100; i++){
    hists.push_back((TH1F*)hist->Clone());
  }

  vector<double> *pdf_weights = new vector<double>(100);

  outputFile->cd();
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("Mass_Rec",&massRec);
  tree->SetBranchAddress("passed_measurement_rec",&passed_measurement_rec);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("pdf_weights",&pdf_weights);
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;

    double weight = rec_weight * gen_weight;

    // Fill 1D Hists (for truth distribution only use gen selection)
    for(unsigned int i=0; i<100; i++){
      if(passed_measurement_rec) hists[i]->Fill(massRec, weight*(*pdf_weights)[i]);
    }
  }

  for(unsigned int i=0; i<100; i++){
    hists[i]->Write(histnames[i]);
  }

  return;
}
