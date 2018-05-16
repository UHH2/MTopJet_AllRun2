#include "do_unfolding.h"

using namespace std;

int main(int argc, char* argv[])
{

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
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Pseudo/";
    output_file = "Results_pseudo.root";
  }
  else if(strcmp(argv[1], "pseudo2") == 0){
    pseudo = true;
    pseudo2 = true;
    directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Pseudo/";
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

  //==============================================
  // step 1 : open output file
  cout << "Open Files" << endl;

  TFile *outputFile=new TFile(output_file,"recreate");

  //==============================================
  // step 2 : read binning schemes and input histograms

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

  // read histograms
  TH1D *hist_data,*hist_mc_gen,*hist_mc_bgr,*hist_mc_sig,*hist_mc_rec, *hist_mc_truth;
  TH1D *hist_data_count400, *hist_data_count500;
  TH2 *histMCGenRec;

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
  inputFile->GetObject("mc_truth",hist_mc_truth);
  if(pseudo1) inputFile->GetObject("pseudodata1_truth",h_pseudodata_truth);
  if(pseudo2) inputFile->GetObject("pseudodata2_truth",h_pseudodata_truth);
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
  std::vector<bool>     show = {true,  false, false,  true, false, false, true}; // decides which masspoint is shown
  // std::vector<bool> show = {true, true, true, true, true, true, true};

  inputFile->GetObject("mc_bgr",hist_mc_bgr);
  inputFile->GetObject("mc_sig",hist_mc_sig);
  hist_mc_rec = (TH1D*)hist_mc_sig->Clone();
  hist_mc_rec->Add(hist_mc_bgr);

  inputFile->GetObject("mc_matrix",histMCGenRec);

  inputFile->GetObject("mc_purity_all",h_pur_all);
  inputFile->GetObject("mc_purity_samebin",h_pur_samebin);
  inputFile->GetObject("mc_purity_samebin_pt",h_pur_samebin_pt);



  if(pseudo1){
    inputFile->GetObject("pseudo1_matrix",histMCGenRec);
    inputFile->GetObject("mc_sig",hist_data);
  }
  if(pseudo2){
    inputFile->GetObject("pseudo2_matrix",histMCGenRec);
    inputFile->GetObject("mc_sig",hist_data);
  }
  if(same) inputFile->GetObject("mc_sig",hist_data);
  if(pseudo || same) hist_mc_bgr->Reset();

  hist_data->Write();
  hist_mc_gen->Write();
  hist_mc_bgr->Write();
  hist_mc_sig->Write();
  hist_mc_rec->Write();
  histMCGenRec->Write();

  if((!hist_data)||(!hist_mc_gen)||(!histMCGenRec)) {
     cout<<"problem to read input histograms\n";
  }

  // ======================================================================================================
  // unfolding here

  int nscan = 100;

  TH2 *CovMatrix;
  TH2 *CorMatrix;
  TH2 *ProbMatrix;
  TGraph *lcurve;
  // double lcurveX= 0;
  // double lcurveY= 0;
  TH1 *data_unfolded,*data_unfolded_all;

  if(same){
    unfolding unfold(hist_data, hist_mc_bgr, hist_mc_sig, histMCGenRec, binning_rec, binning_gen, true, 1);
    data_unfolded = unfold.get_output(true);
    data_unfolded_all = unfold.get_output(false);
    CorMatrix = unfold.get_cor_matrix();
    CovMatrix = unfold.get_cov_matrix();
    ProbMatrix = unfold.get_prob_matrix();
  }


  if(data || pseudo){
    unfolding unfold(hist_data, hist_mc_bgr, hist_mc_sig, histMCGenRec, binning_rec, binning_gen, false, nscan);
    data_unfolded = unfold.get_output(true);
    data_unfolded_all = unfold.get_output(false);
    CorMatrix = unfold.get_cor_matrix();
    CovMatrix = unfold.get_cov_matrix();
    ProbMatrix = unfold.get_prob_matrix();
    // lcurve = unfold.get_lcurve();
    // lcurveX= unfold.get_best_point("X");
    // lcurveY= unfold.get_best_point("Y");
  }

  // write_syserror("SYS_TEST.txt", output, hist_PseudoData_truth);

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

  // count events in data input (for proper binning)
  int n_bins = hist_data_count400->GetSize() - 2;
  double n_events = 0;
  int lower_bin = 0;
  int upper_bin = 0;
  int events_min = 150;
  int bin_nr = 0;

  cout << endl;
  while(upper_bin < n_bins){
    while(n_events < events_min && upper_bin < n_bins){
      upper_bin ++;
      n_events = hist_data_count400->Integral(lower_bin, upper_bin);
    }
    bin_nr ++;
    cout << n_events << " between bin " << lower_bin << " and bin " << upper_bin << endl;
    lower_bin = upper_bin + 1;
    upper_bin += 2;
    n_events = 0;
  }
  cout << "total number of bins: " << bin_nr << endl << endl;

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

  // ======================================================================================================
  // parameters for chi2 fit
  double lower = 130;
  double upper = 290;
  bool NormToWidth = true;

  // normalise unfolding output
  Normalise * normData = new Normalise(data_unfolded, CovMatrix, lower, upper, NormToWidth);
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

  // perform chi2 fit
  chi2fit * chi2 = new chi2fit(data_unfolded_norm, CovMatrix_norm, mc_mtop_templates_norm, masses, lower, upper, NormToWidth);
  chi2->CalculateChi2();
  std::vector<double> chi2values = chi2->GetChi2Values();
  TF1* chi2fit = chi2->GetChi2Fit();
  cout << " MASS = " << chi2->GetMass() << " +- " << chi2->GetUncertainty() << std::endl;

  // ======================================================================================================
  // make plots

  plotter * plot = new plotter(directory);
  plot->draw_chi2(chi2fit, masses, chi2values, "chi2fit");
  plot->draw_matrix(ProbMatrix, "Prob_Matrix", true);
  plot->draw_matrix(CorMatrix, "Cor_Matrix", false);
  plot->draw_matrix(CovMatrix, "Cov_Matrix", false);
  plot->draw_matrix(histMCGenRec, "Migration_Matrix", true);
  if(pseudo) plot->draw_output_pseudo(data_unfolded, h_pseudodata_truth, hist_mc_truth, false, "Unfold_pseudo");
  if(pseudo) plot->draw_output_pseudo(data_unfolded_norm, h_pseudodata_truth_norm, hist_mc_truth_norm, true, "Unfold_pseudo_norm");
  plot->draw_output(data_unfolded, hist_mc_truth, false, "Unfold");
  plot->draw_output(data_unfolded_norm, hist_mc_truth_norm, true, "Unfold_norm");
  plot->draw_output(data_unfolded_all, hist_mc_gen, false, "Unfold_all");
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
