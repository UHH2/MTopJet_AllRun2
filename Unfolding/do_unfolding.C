#include "do_unfolding.h"

using namespace std;

TH1* SetSysError(TH1* data_unfolded, TH2* CovTotal);
int FindLargestVariation(vector<TH1*> variations);
TH1* GetModelDelta(TH1* unfolded, TH1* truth);
TH2* CreateCovFromDelta(TH1* delta, TH2* dummyCov);
TH1* GetLumiDelta(TH1* unfolded, double lumiError);
TH1* ConvertToRelative(TH1* sys, TH1* central);
TH1* AddSys(vector<TH1*> sys);

int main(int argc, char* argv[])
{
  if(argc != 3){
    cout << "ERROR: WRONG USAGE!" << endl;
    cout << endl;
    cout << "CORRECT WAY: ./do_unfolding option1 option2" << endl;
    cout << endl;
    cout << "           option1 = data/pseudo1/pseudo2/same" << endl;
    cout << "           option2 = fine/no" << endl;

    return 0;
  }


  bool data = false;
  bool same = false;
  bool pseudo = false;
  bool pseudo1 = false;
  bool pseudo2 = false;

  bool fine;
  if(strcmp(argv[2], "fine") == 0) fine = true;
  else fine = false;

  if(strcmp(argv[1], "data") == 0){
    data = true;
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Data/";
    output_file = "Results_data.root";
  }
  else if(strcmp(argv[1], "same") == 0){
    same = true;
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Same/";
    output_file = "Results_same.root";
  }
  else if(strcmp(argv[1], "pseudo1") == 0){
    pseudo = true;
    pseudo1 = true;
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Pseudo1/";
    output_file = "Results_pseudo.root";
  }
  else if(strcmp(argv[1], "pseudo2") == 0){
    pseudo = true;
    pseudo2 = true;
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Pseudo2/";
    output_file = "Results_pseudo.root";
  }
  else {
    cout << "select 'data', 'same', 'pseudo1', 'pseudo2'"<< endl;
    return 0;
  }

  cout << "Selected Mode: ";
  if(data) cout << "DATA" << endl;
  else if(same) cout << "SAME" << endl;
  else if(pseudo1) cout << "PSEUDO1" << endl;
  else if(pseudo2) cout << "PSEUDO2" << endl;
  else cout << "ERROR: None selected!" << endl;

  if(fine) cout << "You Selected Fine Binning" << endl;

  TH1::SetDefaultSumw2();

  /*
  ██████ ██████  ███████  █████  ████████ ███████     ███████ ██ ██      ███████ ███████
  ██      ██   ██ ██      ██   ██    ██    ██          ██      ██ ██      ██      ██
  ██      ██████  █████   ███████    ██    █████       █████   ██ ██      █████   ███████
  ██      ██   ██ ██      ██   ██    ██    ██          ██      ██ ██      ██           ██
  ██████ ██   ██ ███████ ██   ██    ██    ███████     ██      ██ ███████ ███████ ███████
  */

  cout << "Open Files" << endl;

  TFile *outputFile=new TFile(output_file,"recreate");

  /*
  ██████  ███████  █████  ██████      ██████  ██ ███    ██ ███    ██ ██ ███    ██  ██████
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██ ████   ██ ████   ██ ██ ████   ██ ██
  ██████  █████   ███████ ██   ██     ██████  ██ ██ ██  ██ ██ ██  ██ ██ ██ ██  ██ ██   ███
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██ ██  ██ ██ ██  ██ ██ ██ ██  ██ ██ ██    ██
  ██   ██ ███████ ██   ██ ██████      ██████  ██ ██   ████ ██   ████ ██ ██   ████  ██████
  */

  if(fine) input_file = "Histograms_fine.root";
  else     input_file = "Histograms.root";
  TFile *inputFile=new TFile(input_file);

  outputFile->cd();

  cout << "Get Binning Schemes" << endl;


  // read binning schemes in XML format
  TUnfoldBinning *binning_rec,*binning_gen;
  TDOMParser parser;
  if(fine) binning_xml = "Binning_fine.xml";
  else     binning_xml = "Binning.xml";

  Int_t error = parser.ParseFile(binning_xml);

  if(error) cout<<"error="<<error<<" from TDOMParser\n";
  TXMLDocument const *XMLdocument=parser.GetXMLDocument();
  binning_rec = TUnfoldBinningXML::ImportXML(XMLdocument,"binning_rec");
  binning_gen = TUnfoldBinningXML::ImportXML(XMLdocument,"binning_gen");

  if((!binning_rec)||(!binning_gen)) {
    cout<<"problem to read binning schemes\n";
  }


  // save binning schemes to output file
  binning_rec->Write();
  binning_gen->Write();

  cout << "Get Histograms" << endl;

  /*
  ██████  ███████  █████  ██████      ██   ██ ██ ███████ ████████  ██████   ██████  ██████   █████  ███    ███ ███████
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██ ██         ██    ██    ██ ██       ██   ██ ██   ██ ████  ████ ██
  ██████  █████   ███████ ██   ██     ███████ ██ ███████    ██    ██    ██ ██   ███ ██████  ███████ ██ ████ ██ ███████
  ██   ██ ██      ██   ██ ██   ██     ██   ██ ██      ██    ██    ██    ██ ██    ██ ██   ██ ██   ██ ██  ██  ██      ██
  ██   ██ ███████ ██   ██ ██████      ██   ██ ██ ███████    ██     ██████   ██████  ██   ██ ██   ██ ██      ██ ███████
  */

  TH1D *hist_data,*hist_mc_gen,*hist_pseudo_gen,*hist_mc_bgr,*hist_mc_sig,*hist_mc_rec, *hist_mc_truth;
  TH2 *histMCGenRec, *histMCGenRec_shower, *histMCGenRec_generator;

  TH1D* h_pseudodata_truth;

  TH1D* h_pur_samebin, *h_pur_samebin_pt, *h_pur_all;

  TH1D* mc_mtop1665_truth;
  TH1D* mc_mtop1695_truth;
  TH1D* mc_mtop1715_truth;
  TH1D* mc_mtop1725_truth;
  TH1D* mc_mtop1735_truth;
  TH1D* mc_mtop1755_truth;
  TH1D* mc_mtop1785_truth;
  std::vector<TH1D*> mc_mtop_templates;

  std::vector<TH1D*> hist_PseudoData , hist_PseudoData_gen, hist_PseudoMC_sig, hist_PseudoData_truth, hist_PseudoMC_truth;


  inputFile->GetObject("data",hist_data);
  inputFile->GetObject("mc_gen",hist_mc_gen);
  inputFile->GetObject("mc_sig",hist_mc_sig);

  inputFile->GetObject("mc_truth",hist_mc_truth);
  if(pseudo1) inputFile->GetObject("pseudo1_truth",h_pseudodata_truth);
  if(pseudo2) inputFile->GetObject("pseudo2_truth",h_pseudodata_truth);
  inputFile->GetObject("mc_mtop1665_truth",mc_mtop1665_truth);
  inputFile->GetObject("mc_mtop1695_truth",mc_mtop1695_truth);
  inputFile->GetObject("mc_mtop1715_truth",mc_mtop1715_truth);
  inputFile->GetObject("mc_mtop1735_truth",mc_mtop1735_truth);
  inputFile->GetObject("mc_mtop1755_truth",mc_mtop1755_truth);
  inputFile->GetObject("mc_mtop1785_truth",mc_mtop1785_truth);
  mc_mtop1725_truth = (TH1D*)hist_mc_truth->Clone("mc_mtop1725_truth");

  mc_mtop_templates.push_back(mc_mtop1665_truth);
  mc_mtop_templates.push_back(mc_mtop1695_truth);
  mc_mtop_templates.push_back(mc_mtop1715_truth);
  mc_mtop_templates.push_back(mc_mtop1725_truth);
  mc_mtop_templates.push_back(mc_mtop1735_truth);
  mc_mtop_templates.push_back(mc_mtop1755_truth);
  mc_mtop_templates.push_back(mc_mtop1785_truth);
  std::vector<double> masses = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
  std::vector<bool>     show = {false,  true, false,  true, false,  true, false}; // decides which masspoint is shown
  // std::vector<bool> show = {true, true, true, true, true, true, true};

  // read backgrounds
  vector<TString> bgr_name = {"WJets", "SingleTop", "other"};
  vector<TH1D*> backgrounds;
  for(unsigned int i=0; i<bgr_name.size(); i++){
    backgrounds.push_back((TH1D*)inputFile->Get("BKG_"+bgr_name[i]));
    if(i==0) hist_mc_bgr = (TH1D*)backgrounds[0]->Clone();
    else hist_mc_bgr->Add(backgrounds[i]);
  }

  if(!data) hist_mc_bgr->Reset();
  hist_mc_rec = (TH1D*)hist_mc_sig->Clone();
  if(data) hist_mc_rec->Add(hist_mc_bgr);

  inputFile->GetObject("mc_matrix",histMCGenRec);
  inputFile->GetObject("pseudo1_matrix",histMCGenRec_generator);
  inputFile->GetObject("pseudo2_matrix",histMCGenRec_shower);
  TH2* MatrixDelta_generator = (TH2*) histMCGenRec->Clone();
  MatrixDelta_generator->Add(histMCGenRec_generator, -1);
  TH2* MatrixDelta_shower = (TH2*) histMCGenRec->Clone();
  MatrixDelta_shower->Add(histMCGenRec_shower, -1);
  // purity
  if(pseudo1){
    inputFile->GetObject("pseudo1_purity_all",h_pur_all);
    inputFile->GetObject("pseudo1_purity_samebin",h_pur_samebin);
    inputFile->GetObject("pseudo1_purity_samebin_pt",h_pur_samebin_pt);
  }
  else if(pseudo2){
    inputFile->GetObject("pseudo2_purity_all",h_pur_all);
    inputFile->GetObject("pseudo2_purity_samebin",h_pur_samebin);
    inputFile->GetObject("pseudo2_purity_samebin_pt",h_pur_samebin_pt);
  }
  else{
    inputFile->GetObject("mc_purity_all",h_pur_all);
    inputFile->GetObject("mc_purity_samebin",h_pur_samebin);
    inputFile->GetObject("mc_purity_samebin_pt",h_pur_samebin_pt);
  }

  // if no data is used, background has to be added to input histogram
  if(pseudo1){
    inputFile->GetObject("pseudo1_sig",hist_data);
    inputFile->GetObject("pseudo1_gen",hist_pseudo_gen);
    hist_data->Add(hist_mc_bgr);
  }
  if(pseudo2){
    inputFile->GetObject("pseudo2_sig",hist_data);
    inputFile->GetObject("pseudo2_gen",hist_pseudo_gen);
    hist_data->Add(hist_mc_bgr);
  }
  if(same){
    inputFile->GetObject("mc_sig",hist_data);
    hist_data->Add(hist_mc_bgr);
  }

  // read migrations from variations (vector of every source containing a vector of every variation e.g. up/down)
  vector<TString> btag_name = {"btagbcup", "btagbcdown", "btagudsgup", "btagudsgdown"};
  vector<TString> jec_name = {"jecup", "jecdown"};
  vector<TString> jer_name = {"jerup", "jerdown"};
  vector<TString> muid_name = {"muidup", "muiddown"};
  vector<TString> mutr_name = {"mutrup", "mutrdown"};
  vector<TString> pu_name = {"puup", "pudown"};
  vector< vector<TString> > sys_name = {btag_name, jec_name, jer_name, muid_name, mutr_name, pu_name};
  vector<TString> sys_rel_name = {"b-tagging", "jec", "jer", "MuID", "MuTrigger", "pile-up"}; // used for comparison plot
  vector< vector<TH2*> > sys_matrix;
  for(unsigned int i=0; i<sys_name.size(); i++){
    vector<TH2*> dummy;
    sys_matrix.push_back(dummy);
    for(unsigned int j=0; j<sys_name[i].size(); j++){
      sys_matrix[i].push_back((TH2*)inputFile->Get(sys_name[i][j] + "_matrix"));
    }
  }

  // now read model variations
  vector<TString> shower = {"Shower"};
  vector<TString> generator = {"Generator"};
  vector<TString> scale = {"SCALE_upup", "SCALE_upnone", "SCALE_noneup", "SCALE_downdown", "SCALE_downnone", "SCALE_nonedown"};
  vector< vector<TString> > model_name = {shower, generator, scale};
  vector<TString> model_rel_name = {"shower", "generator", "scale"};
  vector< vector<TH1D*> > model_input;
  vector< vector<TH1D*> > model_truth;
  for(unsigned int i=0; i<model_name.size(); i++){
    vector<TH1D*> dummy;
    model_input.push_back(dummy);
    model_truth.push_back(dummy);
    for(unsigned int j=0; j<model_name[i].size(); j++){
      TH1D* input = (TH1D*)inputFile->Get(model_name[i][j] + "_sig");
      model_input[i].push_back(input);
      model_truth[i].push_back((TH1D*)inputFile->Get(model_name[i][j] + "_truth"));
    }
  }

  vector<TH1D*> GeneratorVariations;
  TH1D* GeneratorTruth = (TH1D*)inputFile->Get("Generator_truth");
  for(int i=0; i<200; i++){
    TString name = "Generator_sig_";
    name += i;
    GeneratorVariations.push_back((TH1D*)inputFile->Get(name));
  }


  hist_data->Write();
  hist_mc_gen->Write();
  hist_mc_bgr->Write();
  hist_mc_sig->Write();
  hist_mc_rec->Write();
  histMCGenRec->Write();

  if((!hist_data)||(!hist_mc_gen)||(!histMCGenRec)) {
    cout<<"problem to read input histograms\n";
  }

  /*
  ██    ██ ███    ██ ███████  ██████  ██      ██████  ██ ███    ██  ██████
  ██    ██ ████   ██ ██      ██    ██ ██      ██   ██ ██ ████   ██ ██
  ██    ██ ██ ██  ██ █████   ██    ██ ██      ██   ██ ██ ██ ██  ██ ██   ███
  ██    ██ ██  ██ ██ ██      ██    ██ ██      ██   ██ ██ ██  ██ ██ ██    ██
  ██████  ██   ████ ██       ██████  ███████ ██████  ██ ██   ████  ██████
  */
  cout << "**************************************************************" << endl;
  cout << "************     STARTING UNFOLDING PROCEDURE     ************" << endl;
  cout << "**************************************************************" << endl;

  int nscan = 100;
  cout << "    -- YOU SELECTED " << nscan << " STEPS TO DETERMINE TAU" << endl << endl;

  bool do_model = true;
  bool do_genvar = false;

  if(!do_model) cout << "    -- MODEL VARIATIONS ARE NOT UNFOLDED!" << endl  << endl;
  if(!do_genvar) cout << "    -- GENERATOR SMEARINGS ARE NOT UNFOLDED!" << endl  << endl;


  TH2 *CovStat, *CovInputStat, *CovMatrixStat;
  vector<TH2*> CovBgrStat;

  vector< vector<TH2*> > CovSys;
  vector< vector<TH1*> > SYS_DELTA;
  vector< TH1* > SYS_rel;
  TH1 *SYS_rel_total;


  vector<TH2*> CovBgrScale;
  vector<TH1*> BGR_DELTA;

  vector< vector<TH2*> > CovModel;
  vector< vector<TH1*> > MODEL_DELTA;
  vector< vector<TH1*> > MODEL_BIAS;
  vector< vector<TH1*> > MODEL_OUTPUT;
  vector< TH1* > MODEL_rel;
  TH1 *MODEL_rel_total;
  vector<TH1*> Genvar_OUTPUT;
  TH1* Genvar_combination;
  vector<TF1*> Genvar_fits;
  vector<TH1D*> Genvar_variations;

  TH1* DeltaLumi;
  TH2* CovLumi;

  TH2 *CovTotal;
  // TH2 *CovTotal_TUnfold;

  TH2 *CorMatrix;
  TH2 *ProbMatrix;

  TH1 *data_unfolded,*data_unfolded_sys,*data_unfolded_stat,*data_unfolded_all;
  TH1 * pseudodata_bias;

  if(same){
    unfolding unfold(hist_data, backgrounds, bgr_name, hist_mc_sig, histMCGenRec, sys_matrix, sys_name, binning_rec, binning_gen, true, 1, 0);
    data_unfolded = unfold.get_output(true);
    data_unfolded_all = unfold.get_output(false);
  }


  if(data || pseudo){
    cout << "***********************" << endl;
    cout << " UNFOLDING OF DATA" << endl;
    if(pseudo) for(auto i: backgrounds) i->Reset();
    unfolding unfold(hist_data, backgrounds, bgr_name, hist_mc_sig, histMCGenRec, sys_matrix, sys_name, binning_rec, binning_gen, false, nscan, 0);
    data_unfolded = unfold.get_output(true);
    data_unfolded_all = unfold.get_output(false);
    CorMatrix = unfold.get_cor_matrix();
    ProbMatrix = unfold.get_prob_matrix();
    pseudodata_bias = unfold.GetBiasDistribution();
    CovInputStat = unfold.GetInputStatCov();
    CovMatrixStat = unfold.GetMatrixStatCov();
    CovBgrStat = unfold.GetBgrStatCov();
    CovBgrScale = unfold.GetBgrScaleCov();
    CovSys = unfold.GetSysCov();
    SYS_DELTA = unfold.get_sys_delta();
    BGR_DELTA = unfold.get_bgr_delta();
    // CovTotal_TUnfold = unfold.GetTotalCov();

    // now unfold every model variation, get difference to truth and fill cov matrices
    // since tau depends only on the Migration Matrix, one has to find taus only in the first unfolding
    // for every following unfolding the same tau is used
    if(do_model){
      for(unsigned int i=0; i<model_name.size(); i++){
        vector<TH1*> dummy;
        MODEL_OUTPUT.push_back(dummy);
        MODEL_DELTA.push_back(dummy);
        MODEL_BIAS.push_back(dummy);
        vector<TH2*> dummy2;
        CovModel.push_back(dummy2);
        for(unsigned int j=0; j<model_name[i].size(); j++){
          cout << "***********************" << endl;
          cout << " UNFOLDING OF " << model_name[i][j] << endl;
          vector<TH1D*> background_dummy = backgrounds; // create empty vector of backgrounds for model uncertainties
          for(auto i: background_dummy) i->Reset();
          unfolding* unfold_model = new unfolding(model_input[i][j], background_dummy, bgr_name, hist_mc_sig, histMCGenRec, sys_matrix, sys_name, binning_rec, binning_gen, false, nscan, 0);
          TH1* output = unfold_model->get_output(true);
          MODEL_OUTPUT[i].push_back(output);
          TH1* delta = GetModelDelta(output, model_truth[i][j]);
          MODEL_DELTA[i].push_back(delta);
          TH1* bias = unfold_model->GetBiasDistribution();
          MODEL_BIAS[i].push_back(bias);
          TH2* cov = CreateCovFromDelta(delta, CovInputStat);
          CovModel[i].push_back(cov);
          //delete output;
          delete unfold_model;
          cout << "unfolding finished" << endl;
        }
      }
    }

    // Unfold 100 smeared Generator Distribution
    // take mean in every bin as central value and rms as error
    // for every following unfolding the same tau is used
    if(do_genvar){
      double tau = 0;
      for(int i=0; i<GeneratorVariations.size(); i++){
        cout << "********************************************" << endl;
        cout << " UNFOLDING OF SMEARED DISTRIBUTION NUMBER " << i << endl;
        vector<TH1D*> background_dummy = backgrounds; // create empty vector of backgrounds for model uncertainties
        for(auto i: background_dummy) i->Reset();
        unfolding* unfold_genvar;
        if(i==0) unfold_genvar = new unfolding(GeneratorVariations[i], background_dummy, bgr_name, hist_mc_sig, histMCGenRec, sys_matrix, sys_name, binning_rec, binning_gen, false, nscan, 0);
        else     unfold_genvar = new unfolding(GeneratorVariations[i], background_dummy, bgr_name, hist_mc_sig, histMCGenRec, sys_matrix, sys_name, binning_rec, binning_gen, false, 0, tau);
        TH1* output = unfold_genvar->get_output(true);
        Genvar_OUTPUT.push_back(output);
        tau = unfold_genvar->get_tau();
        cout << "Log(tau) = " << TMath::Log10(tau) << endl;
      }
      Smearing * smear = new Smearing(Genvar_OUTPUT);
      Genvar_combination = smear->GetResult();
      Genvar_fits = smear->GetFits();
      Genvar_variations = smear->GetVariations();
    }

    // create a delta hist and cov for lumi uncertainty
    DeltaLumi = GetLumiDelta(data_unfolded, 0.025);
    CovLumi = CreateCovFromDelta(DeltaLumi, CovInputStat);

    cout << "******************************************************" << endl;
    // add up stat Cov
    cout << "sum up stat cov matrices" << endl;
    CovStat = (TH2*) CovInputStat->Clone();
    CovStat->Add(CovMatrixStat);
    for(auto bgrcov: CovBgrStat) CovStat->Add(bgrcov);

    // then add sys cov from backgrounds
    cout << "sum up background sys cov matrices" << endl;
    CovTotal = (TH2*) CovStat->Clone();
    for(auto bgrcov: CovBgrScale) CovTotal->Add(bgrcov);

    // then add sys cov (and convert used uncertainty to relative hist)
    cout << "sum up experimental sys cov matrices" << endl;
    for(unsigned int i=0; i<SYS_DELTA.size(); i++){
      int j = FindLargestVariation(SYS_DELTA[i]);
      CovTotal->Add(CovSys[i][j]);
      SYS_rel.push_back(ConvertToRelative(SYS_DELTA[i][j], data_unfolded));
      cout << "using " << sys_name[i][j] << endl;
    }
    SYS_rel_total = AddSys(SYS_rel);

    // Add Lumi
    CovTotal->Add(CovLumi);

    // then add model cov
    if(!do_model) cout << "!!!! ATTENTION: MODEL UNCERTAINTIES SWITCHED OFF!!!!! " <<endl;
    cout << "sum up model sys cov matrices" << endl;
    for(unsigned int i=0; i<MODEL_DELTA.size(); i++){
      int j = FindLargestVariation(MODEL_DELTA[i]);
      if(do_model) CovTotal->Add(CovModel[i][j]);
      MODEL_rel.push_back(ConvertToRelative(MODEL_DELTA[i][j], data_unfolded));
      cout << "using " << model_name[i][j] << endl;
    }

    if(do_model) MODEL_rel_total = AddSys(MODEL_rel);


    cout << "******************************************************" << endl;
    data_unfolded_sys = SetSysError(data_unfolded, CovTotal);
    data_unfolded_stat = SetSysError(data_unfolded, CovStat);
  }



  // get projections
  int ngen = binning_gen->GetEndBin();
  int nrec = binning_rec->GetEndBin();
  TH1D * gen_proj = histMCGenRec->ProjectionX("gen_proj", 0, nrec+1, "e");
  TH1D * rec_proj = histMCGenRec->ProjectionY("rec_proj", 0, ngen+1, "e");

  for(int i = 1; i<= nrec; i++){
    if(rec_proj->GetBinContent(i) > 200){
      std::cout << "Bin Number " << i << std::endl;
      std::cout <<  binning_rec->GetBinName(i)<< std::endl;
    }
  }


  // ======================================================================================================
  // some counting

  // cout << "total Sum of Weights: " << ProbMatrix->GetSumOfWeights() << endl;
  // cout << "total incl overflow/underflow: " << ProbMatrix->Integral(0,ngen+1,0,nrec+1) << endl;


  // get event counts from various bins
  double count_gen = 0;
  double count_rec = 0;
  for(int rec=0; rec<= nrec+1; rec++) count_gen += histMCGenRec->GetBinContent(0, rec);
  for(int gen=0; gen<= ngen+1; gen++) count_rec += histMCGenRec->GetBinContent(gen, 0);

  cout<< "Events in gen underflow  = " << count_gen << endl;
  cout<< "Events in rec underflow  = " << count_rec << endl;
  cout<< "Events in total Matrix   = " << histMCGenRec->Integral(0,ngen+1,0,nrec+1) << endl;
  cout<< "Events without underflow = " << histMCGenRec->Integral(1,ngen+1,1,nrec+1) << endl;
  cout<< "Events in gen measurement= " << hist_mc_gen->Integral(1, ngen+1) << endl;
  cout<< "percent reconstructed    = " << 100*histMCGenRec->Integral(0,ngen+1,1,nrec+1)/hist_mc_gen->Integral(1, ngen) << "%" << endl;





  cout<< "Events in generated in measurement region = " << histMCGenRec->Integral(1, 16,1,nrec+1) << endl;
  cout<< "Events in generated in measurement region = " << histMCGenRec->Integral(1, 16,0,nrec+1) << endl;


  data_unfolded->Write();
  ProbMatrix->Write();

  /*
  ██████ ██   ██ ██     ██████
  ██      ██   ██ ██          ██
  ██      ███████ ██      █████
  ██      ██   ██ ██     ██
  ██████ ██   ██ ██     ███████
  */

  // parameters for chi2 fit
  double lower = 120;
  double upper = 250;
  bool NormToWidth = true;

  // normalise unfolding output
  Normalise * normData_stat = new Normalise(data_unfolded, CovStat, lower, upper, NormToWidth);
  TH1D* data_unfolded_stat_norm = normData_stat->GetHist(); // this hist is just for plotting

  Normalise * normData = new Normalise(data_unfolded, CovTotal, lower, upper, NormToWidth);
  TH1D* data_unfolded_norm = normData->GetHist();
  TH2D* CovMatrix_norm = normData->GetMatrix();

  // normalise mass samples
  vector<TH1D*> mc_mtop_templates_norm;
  vector<Normalise*> normMass;
  for(unsigned int i=0; i<mc_mtop_templates.size(); i++){
    normMass.push_back(new Normalise(mc_mtop_templates[i], lower, upper, NormToWidth));
    mc_mtop_templates_norm.push_back(normMass[i]->GetHist());
  }
  // normalise mc truth
  Normalise * normMC = new Normalise(hist_mc_truth, lower, upper, NormToWidth);
  TH1D* hist_mc_truth_norm  = normMC->GetHist();

  // normalise pseudo data truth
  TH1D* h_pseudodata_truth_norm;
  if(pseudo){
    Normalise * normPseudo = new Normalise(h_pseudodata_truth, lower, upper, NormToWidth);
    h_pseudodata_truth_norm  = normPseudo->GetHist();
  }

  // write a new vector for masses used in the chi2 fit
  vector<TH1D*> chi2_MassSamples;
  std::vector<double> chi2_masses;
  for(unsigned int i=0; i<mc_mtop_templates_norm.size(); i++){
    if(i != 0 && i != mc_mtop_templates_norm.size()-1){
      chi2_MassSamples.push_back(mc_mtop_templates_norm[i]);
      chi2_masses.push_back(masses[i]);
    }
  }

  // perform chi2 fit
  chi2fit * chi2 = new chi2fit(data_unfolded_norm, CovMatrix_norm, chi2_MassSamples, chi2_masses, lower, upper, NormToWidth);
  chi2->CalculateChi2();
  std::vector<double> chi2values = chi2->GetChi2Values();
  TF1* chi2fit = chi2->GetChi2Fit();
  cout << " MASS = " << chi2->GetMass() << " +- " << chi2->GetUncertainty() << std::endl;

  /*
  ██████  ██       ██████  ████████
  ██   ██ ██      ██    ██    ██
  ██████  ██      ██    ██    ██
  ██      ██      ██    ██    ██
  ██      ███████  ██████     ██
  */

  plotter * plot = new plotter(directory);
  plot->draw_chi2(chi2fit, chi2_masses, chi2values, chi2->GetMass(), chi2->GetUncertainty(), "chi2fit");
  plot->draw_matrix(ProbMatrix, "Prob_Matrix", true);
  plot->draw_matrix(CorMatrix, "Cor_Matrix", false);
  plot->draw_matrix(CovStat, "COV_STAT", false);
  plot->draw_matrix(CovInputStat, "COV_INPUT_STAT", false);
  plot->draw_matrix(CovMatrixStat, "COV_MATRIX_STAT", false);
  plot->draw_matrix(CovTotal, "COV_TOTAL", false);
  // plot->draw_matrix(CovTotal_TUnfold, "COV_TOTAL_check", false);
  plot->draw_matrix(histMCGenRec, "Migration_Matrix", true);
  plot->draw_matrix(histMCGenRec_generator, "Migration_Matrix_Generator", true);
  plot->draw_matrix(histMCGenRec_shower, "Migration_Matrix_Shower", true);
  plot->draw_matrix(MatrixDelta_generator, "Migration_Delta_Generator", false);
  plot->draw_matrix(MatrixDelta_shower, "Migration_Delta_Shower", false);

  for(unsigned int i=0; i<sys_name.size(); i++){
    for(unsigned int j=0; j<sys_name[i].size(); j++){
      plot->draw_matrix(CovSys[i][j], "COV_"+sys_name[i][j], false);
      plot->draw_delta(SYS_DELTA[i][j], "DELTA_"+sys_name[i][j]);
    }
  }
  plot->draw_delta_comparison(SYS_rel_total, SYS_rel, sys_rel_name, "exp", "SYS_EXP_COMPARISION");

  plot->draw_matrix(CovLumi, "COV_Lumi", false);
  plot->draw_delta(DeltaLumi, "DELTA_Lumi");

  if(do_model){
    for(unsigned int i=0; i<model_name.size(); i++){
      for(unsigned int j=0; j<model_name[i].size(); j++){
        plot->draw_matrix(CovModel[i][j], "COV_"+model_name[i][j], false);
        plot->draw_delta(MODEL_DELTA[i][j], "DELTA_"+model_name[i][j]);
        plot->draw_bias(MODEL_OUTPUT[i][j], model_truth[i][j], MODEL_BIAS[i][j], "BIAS_"+model_name[i][j]);
      }
    }
    plot->draw_delta_comparison(MODEL_rel_total, MODEL_rel, model_rel_name, "model", "SYS_MODEL_COMPARISION");
  }

  for(unsigned int i=0; i<bgr_name.size(); i++){
    plot->draw_matrix(CovBgrStat[i], "COV_"+bgr_name[i]+"_stat", false);
    plot->draw_matrix(CovBgrScale[i], "COV_"+bgr_name[i]+"_scale", false);
    plot->draw_delta(BGR_DELTA[i], "DELTA_"+bgr_name[i]);
  }

  if(do_genvar){
    plot->draw_output(Genvar_combination, GeneratorTruth, false, "GENERATOR_SMEARING_COMBINATION");
    if(data || pseudo) plot->draw_output_smear(Genvar_OUTPUT, GeneratorTruth, "GeneratorSmearing");
    for(unsigned int i=0; i<Genvar_fits.size(); i++){
      TString name = "SmearFit_Bin";
      int bin = i+1;
      name += bin;
      plot->draw_smearFit(Genvar_variations[i], Genvar_fits[i], name);
    }
  }
  if(pseudo) plot->draw_output_pseudo(data_unfolded, h_pseudodata_truth, hist_mc_truth, false, "Unfold_pseudo");
  if(pseudo) plot->draw_output_pseudo(data_unfolded_norm, h_pseudodata_truth_norm, hist_mc_truth_norm, true, "Unfold_pseudo_norm");
  if(pseudo) plot->draw_bias(data_unfolded, h_pseudodata_truth, pseudodata_bias, "Unfold_pseudo_bias");
  if(pseudo) plot->draw_output_pseudo(data_unfolded_all, hist_pseudo_gen, hist_mc_gen, false, "Unfold_pseudo_all");
  plot->draw_output(data_unfolded_all, hist_mc_gen, false, "Unfold_all");
  plot->draw_output_stat(data_unfolded_sys, data_unfolded_stat, hist_mc_truth, false, "Unfold");
  plot->draw_output(data_unfolded_sys, hist_mc_truth, false, "Unfold_SYS");
  plot->draw_output_stat(data_unfolded_norm, data_unfolded_stat_norm, hist_mc_truth_norm, true, "Unfold_norm");
  plot->draw_output_mass(data_unfolded_norm, mc_mtop_templates_norm, show, true, "Unfold_masspoints_norm");
  plot->draw_output_mass(data_unfolded, mc_mtop_templates, show, false, "Unfold_masspoints");
  plot->draw_projection(gen_proj, hist_mc_gen, "Projection_gen");
  plot->draw_projection(rec_proj, hist_mc_sig, "Projection_rec");
  plot->draw_1D_hist(hist_mc_gen, "Gen");
  plot->draw_rec(hist_data, hist_mc_sig, hist_mc_bgr, "Rec");
  plot->draw_purity(h_pur_samebin, h_pur_all, "Purity");
  plot->draw_purity(h_pur_samebin_pt, h_pur_all, "Purity_pt");

  cout << "finished" << endl;
  return 0;
}

/*
██   ██ ███████ ██      ██████  ███████ ██████  ███████
██   ██ ██      ██      ██   ██ ██      ██   ██ ██
███████ █████   ██      ██████  █████   ██████  ███████
██   ██ ██      ██      ██      ██      ██   ██      ██
██   ██ ███████ ███████ ██      ███████ ██   ██ ███████
*/

// choose variation
int FindLargestVariation(vector<TH1*> variations){
  if(variations.size() < 1) cout << "Vector of variations has size 0" << endl;
  int nbins = variations[0]->GetXaxis()->GetNbins();
  vector<double> entry;
  for(auto i:variations) entry.push_back(0);
  for(unsigned int j=0; j<variations.size(); j++){
    for(int i=1; i<=nbins; i++){
      entry[j] += abs(variations[j]->GetBinContent(i));
    }
  }
  double max_value = 0;
  int position = -1;
  for(unsigned int i=0; i<entry.size(); i++){
    if(entry[i] > max_value){
      max_value = entry[i];
      position = i;
    }
  }
  if(position == -1) cout << "NO MAX VALUE FOUND FOR VARIATION" << endl;
  return position;
}

// set uncert from cov to 1d hist
TH1* SetSysError(TH1* data_unfolded, TH2* CovTotal){
  TH1* hist = (TH1*) data_unfolded->Clone();
  int nbins = hist->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    double error = sqrt(CovTotal->GetBinContent(i,i));
    hist->SetBinError(i, error);
  }
  return hist;
}

TH1* GetModelDelta(TH1* unfolded, TH1* truth){
  TH1* delta = (TH1*) unfolded->Clone();
  delta->Reset();
  int nbins = delta->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    double diff = unfolded->GetBinContent(i) - truth->GetBinContent(i);
    delta->SetBinContent(i, diff);
  }
  return delta;
}

TH2* CreateCovFromDelta(TH1* delta, TH2* dummyCov){
  int nbins = delta->GetXaxis()->GetNbins();
  TH2* cov = (TH2*) dummyCov->Clone();
  cov->Reset();
  for(int i=1; i<=nbins; i++){
    for(int j=1; j<=nbins; j++){
      double entry = delta->GetBinContent(i) * delta->GetBinContent(j);
      cov->SetBinContent(i,j,entry);
    }
  }
  return cov;
}

TH1* GetLumiDelta(TH1* unfolded, double lumiError){
  int nbins = unfolded->GetXaxis()->GetNbins();
  TH1* delta = (TH1*) unfolded->Clone();
  delta->Reset();
  for(int i=1; i<=nbins; i++){
    double entry = lumiError * unfolded->GetBinContent(i);
    delta->SetBinContent(i, entry);
  }
  return delta;
}

TH1* ConvertToRelative(TH1* sys, TH1* central){
  int nbins = central->GetXaxis()->GetNbins();
  TH1* hist = (TH1*)sys->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double cont_central = central->GetBinContent(bin);
    double cont_sys = abs(sys->GetBinContent(bin));
    double percent;
    if(cont_central == 0) percent = 0;
    else percent = 100*cont_sys/cont_central;
    hist->SetBinContent(bin, percent);
    hist->SetBinError(bin, 0.000000001); // just for plotting reasons
  }
  return hist;
}

TH1* AddSys(vector<TH1*> sys){
  int nbins = sys[0]->GetXaxis()->GetNbins();
  TH1* hist = (TH1*)sys[0]->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double sum = 0;
    for(unsigned int j=0; j<sys.size(); j++){
      double cont_sys = abs(sys[j]->GetBinContent(bin));
      sum += cont_sys * cont_sys;
    }
    double entry = sqrt(sum);
    hist->SetBinContent(bin, entry);
  }
  return hist;
}

// get diagonal Cov Matrix from hist errors
TH2D* Tools::GetDiagonalCovMatrix(TH1D* hist){
  int nbins = hist->GetSize();
  TH2D* cov = new TH2D("cov", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  for(int i=1; i <= nbins; i++){
    double error = hist->GetBinError(i) * hist->GetBinError(i);
    cov->SetBinContent(i,i,error);
  }
  return cov;
}
