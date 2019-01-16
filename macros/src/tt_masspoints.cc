#include "../include/CentralInclude.h"

// compile with:
// g++ -o tt_masspoints tt_masspoints.cc `root-config --cflags --glibs`

using namespace std;

TH2D* GetCovFromStat(TH1F* hist, int firstbin, int lastbin, bool binwidth);
TH2D* NormaliseCov(TH2D* old_cov, TH1F* central, int firstbin, int lastbin, bool binwidth);
TH2D* create_cov(TH1F* errorsUP, TH1F* errorsDOWN, int firstbin, int lastbin);
double ComputeChi2(TH1F* data, TH1F* MC, TH2D* cov_data, TH2D* cov_mc, std::vector<int> bins);
void show_matrix(TH2D* matrix, int firstbin, int lastbin);
TH1F* DataWithCorrectErrors(TH1F* JetMass_data, TH2D* cov, int firstbin, int lastbin);
void WriteValuesToFile(int mass, bool gen, double measured, double error);
TH1F* GetNormErrorFromSYS(TH1F* Central, TH1F* Variation, int firstbin, int lastbin);
TH2D* GetCovFromErrors(TH1F* error, int firstbin, int lastbin);
TH1F* GetSymmetricError(TH1F* up, TH1F* down);

int main(int argc, char* argv[])
{
  bool use_data = false;
  bool use_mass = false;
  bool use_pseudo1 = false;
  bool use_pseudo2 = false;
  int mass;
  bool gen = false;
  bool NoModelSys = false;
  bool StatOnly = false;
  bool debug = false;

  if(debug) cout << "[DEBUG] Start Program" << endl;

  if(StatOnly) cout << "[ATTENTION] systematic uncertainties not included!!!" << endl;
  else if(NoModelSys) cout << "[ATTENTION] model uncertainties not included!!!" << endl;

  if(argc < 3){
    cout << "Invalid input" << endl;
    cout << "Usage: ./tt_masspoints data/pseudo1/pseudo2/mass masspoint use_gen" << endl;
  }

  istringstream ss(argv[2]);
  if (!(ss >> mass)) cerr << "Invalid number " << argv[2] << '\n';
  if(strcmp(argv[1], "mass") == 0){
    if(mass != 166 && mass != 169 && mass != 171 && mass != 172 && mass != 173 && mass != 175 && mass != 178){
      cout << "Select a valid Mass" << endl;
      return 0;
    }
  }

  if(strcmp(argv[1], "data") == 0){
    cout << "Input Distribution: DATA" << endl;
    use_data = true;
  }
  else if(strcmp(argv[1], "mass") == 0){
    cout << "Input Distribution: MASS wit mtop = " << mass << ".5 GeV" << endl;
    use_mass = true;
  }
  else if(strcmp(argv[1], "pseudo1") == 0){
    cout << "Input Distribution: aMCatNLO + Pythia" << endl;
    use_pseudo1 = true;
  }
  else if(strcmp(argv[1], "pseudo2") == 0){
    cout << "Input Distribution: Powheg + Herwig++" << endl;
    use_pseudo2 = true;
  }
  else{
    cout << "select data, mass, pseudo1 or pseudo2" << endl;
    return 0;
  }

  if(strcmp(argv[3], "true") == 0){
    cout << "Using GEN distribution" << endl;
    gen = true;
  }
  else{
    cout << "Using REC distribution" << endl;
    gen = false;
  }

  bool binwidth = true;

  /*
  ██████  ███████ ████████     ███████ ██ ██      ███████ ███████
  ██      ██         ██        ██      ██ ██      ██      ██
  ██  ███ █████      ██        █████   ██ ██      █████   ███████
  ██   ██ ██         ██        ██      ██ ██      ██           ██
  ██████  ███████    ██        ██      ██ ███████ ███████ ███████
  */

  if(debug) cout << "[DEBUG] Get Files" << endl;


  TString data_file_string;
  if(use_data)      data_file_string = dir+"uhh2.AnalysisModuleRunner.DATA.DATA.root";
  else if(use_mass){
    if(mass == 166) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1665.root";
    if(mass == 169) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root";
    if(mass == 171) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1715.root";
    if(mass == 172) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root";
    if(mass == 173) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1735.root";
    if(mass == 175) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root";
    if(mass == 178) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1785.root";
  }
  else if(use_pseudo1) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_amcatnlo-pythia.root";
  //else if(use_pseudo1) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_MuRdown_MuFdown.root";
  else if(use_pseudo2) data_file_string = dir+"uhh2.AnalysisModuleRunner.MC.TTbar_powheg-herwig.root";

  TFile *DATA = new TFile(data_file_string);
  TFile *BGR = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.Background_only.root");

  TFile *TT = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile *TT_1665 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1665.root");
  TFile *TT_1695 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1695.root");
  TFile *TT_1715 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1715.root");
  TFile *TT_1735 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1735.root");
  TFile *TT_1755 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1755.root");
  TFile *TT_1785 = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_mtop1785.root");

  vector<TFile*> SYS_FILES;
  if(!StatOnly){
    // SYS_FILES.push_back(new TFile("SYS_DUMMY.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_BTAG.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_BKG.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_PU.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_MUID.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_MUTR.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_COR.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_JEC.root"));
    SYS_FILES.push_back(new TFile(dir+"SYS_JER.root"));
    if(!NoModelSys){
      SYS_FILES.push_back(new TFile(dir+"SYS_SCALE.root"));
      SYS_FILES.push_back(new TFile(dir+"SYS_GENERATOR.root"));
      SYS_FILES.push_back(new TFile(dir+"SYS_SHOWER.root"));
      SYS_FILES.push_back(new TFile(dir+"SYS_PDF.root"));
    }
  }

  /*
  ██████   ███████ ████████    ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██       ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██       ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██       ██   ██ ██      ██    ██         ██
  ██████   ███████    ██       ██   ██ ██ ███████    ██    ███████
  */
  if(debug) cout << "[DEBUG] Get Hists" << endl;

  vector<TH1F*> UP_Variations;
  vector<TH1F*> DOWN_Variations;
  vector<TH2D*> CovMatrix;

  for(auto file:SYS_FILES){
    UP_Variations.push_back((TH1F*)file->Get("EnvelopeUp/hist"));
    DOWN_Variations.push_back((TH1F*)file->Get("EnvelopeDown/hist"));
    CovMatrix.push_back((TH2D*)file->Get("CovMatrix/cov"));
  }


  TString histdir, histname;
  if(gen){
    histdir = "XCone_GEN_GenOnly/";
    histname = "Mass_HadJet33";
  }
  else{
    histdir = "XCone_cor/";
    histname = "M_jet1_B";
  }

  TH1F * JetMass_data = (TH1F*)DATA->Get(histdir + histname);
  TH1F * Central = (TH1F*)TT->Get(histdir + histname);

  // if pseudo data is selected, adjust stat uncertainty to match data
  if(use_pseudo1 || use_pseudo2){
    int nbins = JetMass_data->GetSize() - 2;
    for(int i=1; i<=nbins; i++){
      double cont = JetMass_data->GetBinContent(i);
      // JetMass_data->SetBinContent(i,cont);
      JetMass_data->SetBinError(i, sqrt(cont));
    }
  }


  if(use_data){
    TH1F * JetMass_bgr = (TH1F*)BGR->Get(histdir + histname);
    JetMass_data->Add(JetMass_bgr, -1); // subtract background
  }

  TH1F * JetMass_1665 = (TH1F*)TT_1665->Get(histdir + histname);
  TH1F * JetMass_1695 = (TH1F*)TT_1695->Get(histdir + histname);
  TH1F * JetMass_1715 = (TH1F*)TT_1715->Get(histdir + histname);
  TH1F * JetMass_1725 = (TH1F*)TT->Get(histdir + histname);
  TH1F * JetMass_1735 = (TH1F*)TT_1735->Get(histdir + histname);
  TH1F * JetMass_1755 = (TH1F*)TT_1755->Get(histdir + histname);
  TH1F * JetMass_1785 = (TH1F*)TT_1785->Get(histdir + histname);

  // cout << "rel. uncertainty in peak bin (data) = " << 100*JetMass_data->GetBinError(18)/JetMass_data->GetBinContent(18) << " %" << endl;
  // cout << "rel. uncertainty in peak bin   (mc) = " << 100*JetMass_1725->GetBinError(18)/JetMass_1725->GetBinContent(18) << " %" << endl;


  std::vector<TH1F*> mass_samples = {JetMass_1665, JetMass_1695, JetMass_1715, JetMass_1725, JetMass_1735, JetMass_1755, JetMass_1785};
  if(use_mass){
    if(mass == 166) mass_samples.erase(mass_samples.begin());
    if(mass == 169) mass_samples.erase(mass_samples.begin()+1);
    if(mass == 171) mass_samples.erase(mass_samples.begin()+2);
    if(mass == 172) mass_samples.erase(mass_samples.begin()+3);
    if(mass == 173) mass_samples.erase(mass_samples.begin()+4);
    if(mass == 175) mass_samples.erase(mass_samples.begin()+5);
    if(mass == 178) mass_samples.erase(mass_samples.begin()+6);
  }

  Double_t all_masses[] = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
  const int N_samples = mass_samples.size();
  Double_t masses[N_samples];
  int k = 0;
  for(int i =0; i<N_samples; i++){
    if(mass == 166 && k == 0) k++;
    if(mass == 169 && k == 1) k++;
    if(mass == 171 && k == 2) k++;
    if(mass == 172 && k == 3) k++;
    if(mass == 173 && k == 4) k++;
    if(mass == 175 && k == 5) k++;
    if(mass == 178 && k == 6) continue;
    masses[i] = all_masses[k];
    k++;
  }

  /*
  ██████   █████  ███    ██  ██████  ███████      ██████  ███████     ███████ ██ ████████
  ██   ██ ██   ██ ████   ██ ██       ██          ██    ██ ██          ██      ██    ██
  ██████  ███████ ██ ██  ██ ██   ███ █████       ██    ██ █████       █████   ██    ██
  ██   ██ ██   ██ ██  ██ ██ ██    ██ ██          ██    ██ ██          ██      ██    ██
  ██   ██ ██   ██ ██   ████  ██████  ███████      ██████  ██          ██      ██    ██
  */

  if(debug) cout << "[DEBUG] Set Range" << endl;


  double from, to;
  if(gen){
    from = 150;
    to = 200;
  }
  else{
    from = 150;
    to = 200;
  }
  int firstbin = JetMass_data->GetXaxis()->FindBin(from);
  int lastbin = JetMass_data->GetXaxis()->FindBin(to) - 1;

  cout << lastbin << endl;

  /*
  ██   ██  █████  ███    ██ ██████  ██      ███████     ██    ██ ███    ██  ██████ ███████ ██████  ████████  █████  ██ ███    ██ ████████ ██ ███████ ███████
  ██   ██ ██   ██ ████   ██ ██   ██ ██      ██          ██    ██ ████   ██ ██      ██      ██   ██    ██    ██   ██ ██ ████   ██    ██    ██ ██      ██
  ███████ ███████ ██ ██  ██ ██   ██ ██      █████       ██    ██ ██ ██  ██ ██      █████   ██████     ██    ███████ ██ ██ ██  ██    ██    ██ █████   ███████
  ██   ██ ██   ██ ██  ██ ██ ██   ██ ██      ██          ██    ██ ██  ██ ██ ██      ██      ██   ██    ██    ██   ██ ██ ██  ██ ██    ██    ██ ██           ██
  ██   ██ ██   ██ ██   ████ ██████  ███████ ███████      ██████  ██   ████  ██████ ███████ ██   ██    ██    ██   ██ ██ ██   ████    ██    ██ ███████ ███████
  */
  if(debug) cout << "[DEBUG] Handle Systematics" << endl;

  // calculate covariance matrix with normalisation applied (for data and templates!)
  // careful: do not use scaled distribution here!
  TH2D*cov = GetCovFromStat(JetMass_data, firstbin, lastbin, binwidth);
  vector<TH2D*> cov_mc;
  for(unsigned int i=0; i<mass_samples.size(); i++){
    cov_mc.push_back(GetCovFromStat(mass_samples[i], firstbin, lastbin, binwidth));
  }

  // now get covariance matrices for systematic uncertainties
  // these can then be added to the cov with stat uncertainties
  // no more error propagation is needed!
  vector<TH2D*> SYS_COV;
  vector<TH1F*> NormErrorUp;
  vector<TH1F*> NormErrorDown;
  vector<TH1F*> NormErrorSym;
  TString CovNames[] = {"BTAG", "PU", "MUID", "MUTR", "SCALE", "BKG", "COR", "JEC", "JER", "Generator", "Shower", "PDF"};
  for(auto upvar:UP_Variations){
    NormErrorUp.push_back(GetNormErrorFromSYS(Central, upvar, firstbin, lastbin));
  }
  for(auto downvar:DOWN_Variations){
    NormErrorDown.push_back(GetNormErrorFromSYS(Central, downvar, firstbin, lastbin));
  }

  if(NormErrorDown.size() != NormErrorUp.size()){
    cout << "Not the same number of up and down variations!" << endl;
    return 0;
  }

  for(unsigned int i=0; i<NormErrorUp.size(); i++){
    NormErrorSym.push_back(GetSymmetricError(NormErrorUp[i], NormErrorDown[i]));
    // here the Cov Matrix would be made from normalized variations (ONLY DO ONE - Should not change result)
    //SYS_COV.push_back(GetCovFromErrors(NormErrorSym[i], firstbin, lastbin));
  }

  for(unsigned int i=0; i<CovMatrix.size(); i++){
    // here the CovMatrix is filled before normalisation and then normalized (ONLY DO ONE - Should not change result)
    SYS_COV.push_back(NormaliseCov(CovMatrix[i], Central, firstbin, lastbin, binwidth));
  }

  TH2D* STAT_COV = (TH2D*) cov->Clone("STAT_COV");
  if(!gen) for(auto cov_i:SYS_COV) cov->Add(cov_i);

  /*
  ███    ██  ██████  ██████  ███    ███  █████  ██      ██ ███████ ███████
  ████   ██ ██    ██ ██   ██ ████  ████ ██   ██ ██      ██ ██      ██
  ██ ██  ██ ██    ██ ██████  ██ ████ ██ ███████ ██      ██ ███████ █████
  ██  ██ ██ ██    ██ ██   ██ ██  ██  ██ ██   ██ ██      ██      ██ ██
  ██   ████  ██████  ██   ██ ██      ██ ██   ██ ███████ ██ ███████ ███████
  */
  if(debug) cout << "[DEBUG] Normalise" << endl;
  TH1F* data_ = (TH1F*)JetMass_data->Clone("data");

  // normalise here
  for(auto &hist: mass_samples){
    if(binwidth) hist->Scale(1/hist->Integral(firstbin, lastbin), "width");
    else         hist->Scale(1/hist->Integral(firstbin, lastbin));

  }
  if(binwidth) data_->Scale(1/data_->Integral(firstbin, lastbin),"width");
  else         data_->Scale(1/data_->Integral(firstbin, lastbin));
  //


  // now set correct error bars (diagonal entries of final cov matrix)
  TH1F* data = DataWithCorrectErrors(data_, cov, firstbin, lastbin);
  TH1F* dataSTAT = DataWithCorrectErrors(data_, STAT_COV, firstbin, lastbin);

  // select bins that enter chi2 calculation
  // leave one bin out (chi2 should not change)
  int j = 0;
  std::vector< std::vector<int> > bins;
  bool skip;
  int rand_bin = rand()%(lastbin-firstbin+1)+firstbin; // select random bin that is not used
  for(auto &masspoint: mass_samples){
    std::vector<int> bins_this;
    for(int i=firstbin; i<=lastbin; i++){
      if(i == rand_bin) skip = true;
      else skip = false;
      if(!skip) bins_this.push_back(i);
    }
    bins.push_back(bins_this);
    j++;
  }

  /*
  ██████ ██   ██ ██ ██████
  ██     ██   ██ ██      ██
  ██     ███████ ██  █████
  ██     ██   ██ ██ ██
  ██████ ██   ██ ██ ███████
  */
  if(debug) cout << "[DEBUG] Calculate Chi2" << endl;

  std::vector<double> chi_vec;
  j = 0;
  for(auto &masspoint: mass_samples){
    double chi2 = ComputeChi2(data, masspoint, cov, cov_mc[j], bins[j]);
    cout << "   chi2 (mtop = " << masses[j] <<") = "<< chi2 << endl;
    chi_vec.push_back(chi2);
    j++;
  }

  Double_t y[N_samples];
  for(int i=0; i<N_samples; i++){
    y[i] = chi_vec[i];
  }

  if(debug) cout << "[DEBUG] Do Chi2 Fit" << endl;

  TGraph* chi_hist = new TGraph(N_samples,masses,y);
  TF1 * fit;
  TF1*f1 = new TF1("f1","pol2",0,500);
  chi_hist->Fit("f1","R");
  fit = chi_hist->GetFunction("f1");
  double minY = fit->GetMinimum();
  double minX = fit->GetX(minY, 0, 350);
  double minX_sigmaup = fit->GetX(minY+1, minX, 350);
  double minX_sigmadown = fit->GetX(minY+1, 0, minX);
  cout << endl << "=========================== " << endl;
  cout << "| value = " << std::setprecision(5) << minX << " +- " << std::setprecision(4) << (minX_sigmaup-minX) << "|" << endl;
  cout << "=========================== " << endl << endl;

  if(debug) cout << "[DEBUG] Write Values to File" << endl;
  if(use_mass) WriteValuesToFile(mass, gen, minX, (minX_sigmaup-minX));


  /*
  ██████  ██       ██████  ████████ ████████ ██ ███    ██  ██████
  ██   ██ ██      ██    ██    ██       ██    ██ ████   ██ ██
  ██████  ██      ██    ██    ██       ██    ██ ██ ██  ██ ██   ███
  ██      ██      ██    ██    ██       ██    ██ ██  ██ ██ ██    ██
  ██      ███████  ██████     ██       ██    ██ ██   ████  ██████
  */

  if(debug) cout << "[DEBUG] Start Plotting" << endl;

  double xlow, xhigh;
  double ylow, yhigh;
  xlow = from;
  xhigh = to;
  // xlow = 0;
  // xhigh = 1000;
  ylow = 0;
  yhigh = data->GetMaximum() * 1.1;
  if(mass != 166 ){
    if(JetMass_1665->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1665->GetMaximum() * 1.1;
  }
  else{
    if(JetMass_1695->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1695->GetMaximum() * 1.1;
  }
  if(mass != 178 ){
    if(JetMass_1785->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1785->GetMaximum() * 1.1;
  }
  else{
    if(JetMass_1755->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1755->GetMaximum() * 1.1;
  }

  for(auto &i: mass_samples){
    i->SetLineWidth(3);
    i->SetTitle(" ");
    i->GetXaxis()->SetTitle("Leading-jet m_{jet} [GeV]");
    i->GetYaxis()->SetTitle("a.u.");
    i->GetYaxis()->SetRangeUser(ylow, yhigh);
    i->GetXaxis()->SetRangeUser(xlow, xhigh);
    i->GetXaxis()->SetNdivisions(505);
    i->GetYaxis()->SetNdivisions(505);
    i->GetYaxis()->SetTitleOffset(1.5);
    i->GetXaxis()->SetTitleOffset(1.3);
    i->SetLineWidth(4);
  }

  JetMass_1665->SetLineColor(kRed-4);
  JetMass_1695->SetLineColor(kRed-4);
  JetMass_1715->SetLineColor(kRed-4);

  JetMass_1725->SetLineColor(13);

  JetMass_1735->SetLineColor(kAzure+7);
  JetMass_1755->SetLineColor(kAzure+7);
  JetMass_1785->SetLineColor(kAzure+7);

  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(8);
  data->SetMarkerSize(1);
  dataSTAT->SetLineColor(kBlack);
  dataSTAT->SetMarkerColor(kBlack);
  dataSTAT->SetMarkerStyle(8);
  dataSTAT->SetMarkerSize(1);
  data->GetYaxis()->SetTitleOffset(1.5);
  data->GetXaxis()->SetNdivisions(505);
  data->GetYaxis()->SetNdivisions(505);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);

  if(mass != 166) JetMass_1665->Draw("HIST SAME");
  else            JetMass_1695->Draw("HIST SAME");
  if(mass != 172) JetMass_1725->Draw("HIST SAME");
  if(mass != 178) JetMass_1785->Draw("HIST SAME");
  else            JetMass_1755->Draw("HIST SAME");

  data->Draw("E1 SAME");
  dataSTAT->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
  leg->SetBorderSize(0);
  if(use_data) leg->AddEntry(data,"data","pl");
  else         leg->AddEntry(data,"pseudo data","pl");
  if(mass != 166) leg->AddEntry(JetMass_1665,"t#bar{t} (m_{top} = 166.5 GeV)","l");
  else            leg->AddEntry(JetMass_1695,"t#bar{t} (m_{top} = 169.5 GeV)","l");
  if(mass != 172) leg->AddEntry(JetMass_1725,"t#bar{t} (m_{top} = 172.5 GeV)","l");
  if(mass != 178) leg->AddEntry(JetMass_1785,"t#bar{t} (m_{top} = 178.5 GeV)","l");
  else            leg->AddEntry(JetMass_1755,"t#bar{t} (m_{top} = 175.5 GeV)","l");
  leg->Draw("");
  if(use_data)       A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_data.pdf");
  else if(use_pseudo1){
    if(gen) A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_aMCatNLO+Pythia_gen.pdf");
    else    A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_aMCatNLO+Pythia_rec.pdf");
  }
  else if(use_pseudo2){
    if(gen) A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_Powheg+Herwig_gen.pdf");
    else    A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_Powheg+Herwig_rec.pdf");
  }
  else if(use_mass){
    TString filename;
    if(gen) filename = "/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_mtop" + std::to_string(mass) + "_gen.pdf";
    else    filename = "/afs/desy.de/user/s/schwarzd/Plots/Masspoints/JetMass_mtop" + std::to_string(mass) + "_rec.pdf";
    A->SaveAs(filename);
  }


  TCanvas *B = new TCanvas("B", "B", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  // set axis limits such that full pol2 fit is visible
  double xfrom, xto;
  double chimax = 0;
  double chimin = 1000000;
  int chimax_i = 0;
  for(unsigned int i=0; i<chi_vec.size(); i++){
    if(chi_vec[i] > chimax){
      chimax = chi_vec[i];
      chimax_i = i;
    }
    if(chi_vec[i] < chimin) chimin = chi_vec[i];
  }
  double chimax_mass = masses[chimax_i];
  double dist_to_min = abs(chimax_mass - minX);
  xfrom = minX - dist_to_min - 1;
  xto = minX + dist_to_min + 1;
  chi_hist->GetXaxis()->SetLimits(xfrom, xto);
  if(chimin < minY) chi_hist->GetHistogram()->SetMinimum(chimin);
  else              chi_hist->GetHistogram()->SetMinimum(minY);
  chi_hist->GetHistogram()->SetMaximum(chimax*1.1);
  //
  chi_hist->SetTitle(" ");
  chi_hist->GetXaxis()->SetTitle("m^{MC}_{top} [GeV]");
  chi_hist->GetYaxis()->SetTitle("#chi^{2}");
  chi_hist->GetXaxis()->SetTitleOffset(1.3);
  chi_hist->GetYaxis()->SetTitleOffset(1.5);
  chi_hist->GetXaxis()->SetNdivisions(505);
  chi_hist->GetYaxis()->SetNdivisions(505);
  chi_hist->SetMarkerStyle(20);
  chi_hist->SetMarkerSize(1.5);
  chi_hist->SetLineColor(1);
  chi_hist->Draw("AP");
  if(use_data)       B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_data.pdf");
  else if(use_pseudo1){
    if(gen) B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_aMCatNLO+Pythia_gen.pdf");
    else    B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_aMCatNLO+Pythia_rec.pdf");
  }
  else if(use_pseudo2){
    if(gen) B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_Powheg+Herwig_gen.pdf");
    else    B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_Powheg+Herwig_rec.pdf");
  }
  else if(use_mass){
    TString filename;
    if(gen) filename = "/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_mtop" + std::to_string(mass) + "_gen.pdf";
    else    filename = "/afs/desy.de/user/s/schwarzd/Plots/Masspoints/chi2_mtop" + std::to_string(mass) + "_rec.pdf";
    B->SaveAs(filename);
  }


  for(unsigned int i=0; i<NormErrorUp.size(); i++){
    TCanvas *C = new TCanvas("C", "C", 600, 600);
    double min, max;
    if(NormErrorUp[i]->GetMinimum() < NormErrorDown[i]->GetMinimum()) min = 1.1*NormErrorUp[i]->GetMinimum();
    else  min = 1.1*NormErrorDown[i]->GetMinimum();
    if(NormErrorUp[i]->GetMaximum() > NormErrorDown[i]->GetMaximum()) max = 1.1*NormErrorUp[i]->GetMaximum();
    else  max = 1.1*NormErrorDown[i]->GetMaximum();
    NormErrorUp[i]->SetTitle(CovNames[i]);
    NormErrorUp[i]->GetXaxis()->SetTitleOffset(1.3);
    NormErrorUp[i]->GetYaxis()->SetTitleOffset(1.5);
    NormErrorUp[i]->GetXaxis()->SetRangeUser(xlow-10, xhigh+40);
    NormErrorUp[i]->GetYaxis()->SetRangeUser(min, max);
    NormErrorUp[i]->GetXaxis()->SetTitle("Leading-jet m_{jet} [GeV]");
    NormErrorUp[i]->GetYaxis()->SetTitle("uncertainty");
    NormErrorUp[i]->SetLineWidth(4);
    NormErrorDown[i]->SetLineWidth(4);
    NormErrorSym[i]->SetLineWidth(4);
    NormErrorSym[i]->SetLineStyle(7);
    NormErrorUp[i]->SetLineColor(kRed-4);
    NormErrorDown[i]->SetLineColor(kAzure+7);
    NormErrorSym[i]->SetLineColor(13);
    NormErrorUp[i]->Draw("HIST");
    NormErrorDown[i]->Draw("HIST SAME");
    NormErrorSym[i]->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.63,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(NormErrorUp[i],CovNames[i] + " up","l");
    leg->AddEntry(NormErrorDown[i],CovNames[i] + " down","l");
    leg->AddEntry(NormErrorSym[i],CovNames[i] + " symmetric","l");
    leg->Draw();
    C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/COV/Error_"+CovNames[i]+".pdf");
    delete C;
    delete leg;
  }

  for(unsigned int i=0; i<SYS_COV.size(); i++){
    TCanvas *C = new TCanvas("C", "C", 600, 600);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
    //gPad->SetLogz();
    SYS_COV[i]->SetTitle(CovNames[i]);
    SYS_COV[i]->GetXaxis()->SetTitleOffset(1.1);
    SYS_COV[i]->GetYaxis()->SetTitleOffset(1.5);
    SYS_COV[i]->GetXaxis()->SetTitle("measurement bin");
    SYS_COV[i]->GetYaxis()->SetTitle("measurement bin");
    SYS_COV[i]->GetXaxis()->SetRangeUser(firstbin, lastbin+1);
    SYS_COV[i]->GetYaxis()->SetRangeUser(firstbin, lastbin+1);
    SYS_COV[i]->Draw("COLZ");
    SYS_COV[i]->Draw("BOX SAME");
    C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/COV/"+CovNames[i]+".pdf");
    delete C;
  }

  TCanvas *C = new TCanvas("C", "C", 600, 600);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogz();
  STAT_COV->SetTitle("Statistics");
  STAT_COV->GetXaxis()->SetTitleOffset(1.1);
  STAT_COV->GetYaxis()->SetTitleOffset(1.5);
  STAT_COV->GetXaxis()->SetTitle("measurement bin");
  STAT_COV->GetYaxis()->SetTitle("measurement bin");
  STAT_COV->GetXaxis()->SetRangeUser(firstbin, lastbin+1);
  STAT_COV->GetYaxis()->SetRangeUser(firstbin, lastbin+1);
  STAT_COV->Draw("COLZ");
  STAT_COV->Draw("BOX SAME");
  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Masspoints/COV/STAT.pdf");
  delete C;

  cout << "=====================================" << endl;
  cout << "=========== finished ================" << endl;
  cout << "=====================================" << endl << endl ;

  return 0;
}

/*
██   ██ ███████ ██      ██████  ███████ ██████  ███████
██   ██ ██      ██      ██   ██ ██      ██   ██ ██
███████ █████   ██      ██████  █████   ██████  ███████
██   ██ ██      ██      ██      ██      ██   ██      ██
██   ██ ███████ ███████ ██      ███████ ██   ██ ███████
*/

double ComputeChi2(TH1F* data, TH1F* MC, TH2D* cov_data, TH2D* cov_mc, std::vector<int> bins){
  //double chi2 = 0;
  TH1F* data_sc = (TH1F*)data->Clone("data_sc");
  TH1F* MC_sc = (TH1F*)MC->Clone("MC_sc");

  // add cov from input and cov from template for chi2 value
  TH2D* cov_sc = (TH2D*)cov_data->Clone("cov_sc");
  cov_sc->Add(cov_mc, 1.0);


  // bins
  int N =  bins.size();

  //vectors and matrices
  TVectorD vdata(N);
  TVectorD vMC(N);
  TMatrixDSym mat(N);

  for (unsigned int i=0; i < bins.size(); ++i){
    vdata[i] = data_sc->GetBinContent(bins.at(i));
    vMC[i] = MC_sc->GetBinContent(bins.at(i));
    for (unsigned int j=0; j < bins.size(); ++j){
      mat[i][j] = cov_sc->GetBinContent(bins.at(i),bins.at(j));
      // cout << "cov[" << i << "," << j << "] = " << cov->GetBinContent(bins.at(i),bins.at(j)) << endl;
    }
  }

  //invert the matrix
  // TDecompSVD lu(mat);
  TDecompLU lu(mat);

  TMatrixD imat = TMatrixD(mat);
  lu.Invert(imat);

  TMatrixD test(mat, TMatrixD::kMult, imat);
  TMatrixD unit(TMatrixD::kUnit, mat);
  Bool_t ok = VerifyMatrixIdentity(unit, test, false, 10e-10);
  if (!ok){
    cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << endl;
    printf("\n\n------------------------------- cov matrix ----------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("%3d:  ",i+1);
      for (int j=0; j<N; ++j)  printf("% 5.2f ", mat[i][j] );
      printf("\n");
    }
    printf("\n");
    printf("\n\n------------------------- diag. elements in cov matrix --------------------------- \n");
    for (int i=0; i<N; ++i){
      printf("(%3d,%3d)   %9.3f",i,i,mat[i][i]);
      if ( mat[i][i] < 1.e-10 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // calculate chi2
  double chi2 = 0;

  TVectorD diff = vdata - vMC;
  TVectorD right = imat * diff;
  chi2 = diff * right;

  return chi2;
}


TH2D* GetCovFromStat(TH1F* hist, int firstbin, int lastbin, bool binwidth){
  int nbins = lastbin;
  TH2D* old_cov = new TH2D("old_cov", "old_cov", nbins, 1, lastbin+1, nbins, 1, lastbin+1);
  TH2D* new_cov = new TH2D("new_cov", "new_cov", nbins, 1, lastbin+1, nbins, 1, lastbin+1);
  for(int i=1; i <= lastbin; i++){
    double error = hist->GetBinError(i);
    if(i < firstbin ) old_cov->Fill(i,i,0);
    else  old_cov->Fill(i,i,pow(error,2));
  }
  double integral = hist->Integral(firstbin, lastbin);
  for(int i=1; i <= lastbin; i++){
    for(int j=1; j <= lastbin; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= lastbin; k++){
        for(int l=1; l <= lastbin; l++){
          old_entry = old_cov->GetBinContent(k,l);
          double binwidth_i = hist->GetBinWidth(i);
          double binwidth_j = hist->GetBinWidth(j);
          if(!binwidth){
            binwidth_i = 1.0;
            binwidth_j = 1.0;
          }
          if(i==k) derivation_i = (integral - hist->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (hist->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - hist->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (hist->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      new_cov->Fill(i, j, sum);
    }
  }

  // show_matrix(old_cov, firstbin, lastbin);
  // show_matrix(new_cov, firstbin, lastbin);

  return new_cov;
}

TH2D* NormaliseCov(TH2D* old_cov, TH1F* central, int firstbin, int lastbin, bool binwidth){
  int nbins = lastbin;
  TH2D* new_cov = new TH2D("new_cov", "new_cov", nbins, 1, lastbin+1, nbins, 1, lastbin+1);
  double integral = central->Integral(firstbin, lastbin);
  for(int i=1; i <= lastbin; i++){
    for(int j=1; j <= lastbin; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= lastbin; k++){
        for(int l=1; l <= lastbin; l++){
          // set old entries to 0 if outside bin range
          if(k<firstbin)      old_entry = 0;
          else if(l<firstbin) old_entry = 0;
          else                old_entry = old_cov->GetBinContent(k,l);
          ////
          double binwidth_i = central->GetBinWidth(i);
          double binwidth_j = central->GetBinWidth(j);
          if(!binwidth){
            binwidth_i = 1.0;
            binwidth_j = 1.0;
          }
          if(i==k) derivation_i = (integral - central->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (central->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - central->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (central->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      new_cov->Fill(i, j, sum);
    }
  }

  // show_matrix(old_cov, firstbin, lastbin);
  // show_matrix(new_cov, firstbin, lastbin);

  return new_cov;
}


void show_matrix(TH2D* matrix, int firstbin, int lastbin){

  for(int i=firstbin; i <= lastbin; i++){
    printf("%3d:  ",i);
    for(int j=firstbin; j <= lastbin; j++){
      printf("% 5.2f ",  matrix->GetBinContent(j,i));
    }
    printf("\n");
  }
  cout << "================================================" << endl;

  return;
}



TH1F* DataWithCorrectErrors(TH1F* JetMass_data, TH2D* cov, int firstbin, int lastbin){
  int nbins = lastbin;
  TH1F* data = (TH1F*)JetMass_data->Clone("data");
  for(int i=1; i<=nbins; i++){
    if(i<firstbin){
      data->SetBinContent(i, 0);
      data->SetBinError(i, 0);
    }
    else{
      double error = sqrt(cov->GetBinContent(i,i));
      data->SetBinError(i, error);
    }
  }
  return data;
}


void WriteValuesToFile(int mass, bool gen, double measured, double error){
  string x;
  ifstream file;
  if(gen) file.open("MassValues_GEN.txt");
  else    file.open("MassValues_REC.txt");
  ofstream ofile;
  ofile.open("temp.txt");
  string mass_string;
  if(mass == 166) mass_string = "166.5";
  if(mass == 169) mass_string = "169.5";
  if(mass == 171) mass_string = "171.5";
  if(mass == 172) mass_string = "172.5";
  if(mass == 173) mass_string = "173.5";
  if(mass == 175) mass_string = "175.5";
  if(mass == 178) mass_string = "178.5";
  while (!file.eof()){
    getline(file,x);
    if(x.find(mass_string) != std::string::npos){
      ofile << fixed << setprecision(10) << mass_string << " " << measured << " " << error  << endl;
    }
    else ofile << x << endl;
  }
  file.close();
  if(gen) remove("MassValues_GEN.txt");
  else    remove("MassValues_REC.txt");
  ofile.close();
  if(gen) rename("temp.txt","MassValues_GEN.txt");
  else    rename("temp.txt","MassValues_REC.txt");
  return;
}

TH2D* create_cov(TH1F* errorsUP, TH1F* errorsDOWN, int firstbin, int lastbin){
  int nbins = lastbin-1;
  TH2D* cov = new TH2D("cov", "cov", nbins, 1, lastbin, nbins, 1, lastbin);
  for(int i=1; i <= lastbin; i++){
    double errorup = abs(errorsUP->GetBinContent(i));
    double errordown = abs(errorsDOWN->GetBinContent(i));
    double error;
    if(errorup > errordown) error = errorup;
    else error = errordown;

    if(i < firstbin ) cov->Fill(i,i,0);
    else  cov->Fill(i,i,pow(error,2));
  }
  return cov;
}

TH1F* GetNormErrorFromSYS(TH1F* Central, TH1F* Variation, int firstbin, int lastbin){
  int nbins = lastbin;
  TH1F* Error = (TH1F*)Central->Clone("Error");
  Error->Reset();
  double VariationIntegral = Variation->Integral(firstbin, lastbin);
  double CentralIntegral = Central->Integral(firstbin, lastbin);

  std::vector<double> errors;
  for(unsigned i=1; i<=nbins; i++){
    double relerror_i = ((Variation->GetBinContent(i)/VariationIntegral ) / (Central->GetBinContent(i)/CentralIntegral)) - 1;
    double error_i = relerror_i * (Central->GetBinContent(i)/CentralIntegral) * (1/Central->GetBinWidth(i));
    double xvalue = Central->GetBinCenter(i);
    if(i<firstbin) Error->Fill(xvalue, 0);
    else        Error->Fill(xvalue, error_i);
  }
  return Error;
}

TH2D* GetCovFromErrors(TH1F* error, int firstbin, int lastbin){
  int nbins=lastbin;
  TH2D* cov = new TH2D("cov", "cov", nbins, 1, nbins+1, nbins, 1, nbins+1);
  for(int i=1; i <= nbins; i++){
    for(int j=1; j <= nbins; j++){
      double CovEntry = error->GetBinContent(i) * error->GetBinContent(j);
      if(i<firstbin || j<firstbin) cov->Fill(i,j,0);
      else                   cov->Fill(i,j,CovEntry);
    }
  }
  return cov;
}

TH1F* GetSymmetricError(TH1F* up, TH1F* down){
  // calculate absolute values from up and down Variations
  // but the signe matters because of anti-correlations
  // thus, take the largest absolute value, with the sign of the up variation
  // one could also take the sign of the down variation
  int nbins = up->GetSize() - 2;
  TH1F* sym = (TH1F*)up->Clone("sym");
  sym->Reset();
  for(int i=1; i<=nbins; i++){
    double a = up->GetBinContent(i);
    double b = down->GetBinContent(i);
    double sign;
    if(a != 0) sign = a/abs(a);
    else sign = 1.;
    double err;
    if(abs(a) > abs(b)) err = abs(a);
    else                err = abs(b);
    double xvalue = up->GetBinCenter(i);
    double value = err*sign;
    sym->Fill(xvalue, err*sign);
  }
  return sym;
}
