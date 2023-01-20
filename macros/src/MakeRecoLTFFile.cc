#include "../include/CentralInclude.h"

using namespace std;

// CHI2 copied from unfolding
class chi2fit{
 public:
  chi2fit(TH1F* data, TH2F* cov_matrix, std::vector<TH1F*> masspoints, std::vector<double> masses, double lower, double upper, bool width);

  void CalculateChi2();

  TF1* GetChi2Fit();
  std::vector<double> GetChi2Values();
  double GetMass();
  double GetMin();
  double GetUncertainty();
  TGraph* GetChi2Graph();

 private:
  void SetMasses(std::vector<double> masses);
  double ComputeChi2(TH1F* MC, std::vector<int> bins);
  void SetBins(std::vector<TH1F*> masspoints);
  void SetRange(double lower, double upper);
  void SetMasspoints(std::vector<TH1F*> masspoints);
  void SetData(TH1F* data);
  void SetCovMatrix(TH2F* cov_matrix);


  std::vector<TH1F*> masspoints_;
  TH1F* data_;
  int lower_, upper_;
  bool width_;
  TH2F* cov_;
  std::vector< std::vector<int> > bins_;
  TVectorD chi2_;
  TVectorD masses_;
  TF1 * fit_;

};
////

TH1F* SetError(TH1F* data, TH2F* cov);
void Plot_hist(TH1F* h_data_XS_tot, vector<TH1F*> truth, vector<TString> legnames, TString name);
TH1F* get_hist(TFile* file, TH1F* dummy, TString sel_name, TString weightname);
TH2F* NormalizeMatrix(TH2F* old_cov, TH1F* hist_);
TH2F* CreateMatrix(TH1F* hist);
TH1F* NormalizeHist(TH1F* old_hist);
void Plot_chi2(TF1 * fit_, std::vector<double> masses_, std::vector<double> chi2_, double mass, double uncert, TString file_name);
void PlotMatrix(TH2F* cor, TString name);

int Nbins = 30;
double xmin = 100;
double xmax = 300;

int main(int argc, char* argv[]){

  TFile * outfile = new TFile("RECO.root","RECREATE");

  TString region = "passed_measurement_rec";
  
  TH1F* dummy = new TH1F("mjet", "mjet", Nbins, xmin, xmax);

  TString directory = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/";
  TString prefix_mc = "uhh2.AnalysisModuleRunner.MC.";
  TString prefix_data = "uhh2.AnalysisModuleRunner.DATA.";

  vector<TString> channels = {"elec", "muon"};
  vector<TString> processes = {"DATA", "TTbar", "SingleTop", "WJets", "other", "TTbar_mtop1695", "TTbar_mtop1715", "TTbar_mtop1735", "TTbar_mtop1755"};
  vector<TString> years = {"2016v3", "2017v2", "2018"};
  
  cout << "Fill nominal hists" << endl;
  TH1F *h_tt, *h_st, *h_wj, *h_ot, *h_da;
  TH1F *h_tt_1695, *h_tt_1715, *h_tt_1735, *h_tt_1755;
  
  for(auto process: processes){
    bool firsthist=true;
    for(auto year: years){
      TString prefix = prefix_mc;
      if(process == "DATA") prefix = prefix_data;
      TFile* file_el = new TFile(directory+"/elec/"+prefix+process+"_"+year+".root");
      TFile* file_mu = new TFile(directory+"/muon/"+prefix+process+"_"+year+".root");
      cout << " - fill " << process << " " << year << endl;
      TH1F* hist_el = get_hist(file_el, dummy, region, "none");
      TH1F* hist_mu = get_hist(file_mu, dummy, region, "none");
      // cout << process << ", " << year << ", elec:" << hist_el->Integral() << endl;
      // cout << process << ", " << year << ", muon:" << hist_mu->Integral() << endl;
      
      if(firsthist){
        firsthist = false;
        if(process == "TTbar"){
          h_tt = (TH1F*) hist_el->Clone();
          h_tt->Add(hist_mu);
        }
        else if(process == "SingleTop"){
          h_st = (TH1F*) hist_el->Clone();
          h_st->Add(hist_mu);
        }
        else if(process == "WJets"){
          h_wj = (TH1F*) hist_el->Clone();
          h_wj->Add(hist_mu);
        }
        else if(process == "other"){
          h_ot = (TH1F*) hist_el->Clone();
          h_ot->Add(hist_mu);
        }
        else if(process == "DATA"){
          h_da = (TH1F*) hist_el->Clone();
          h_da->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1695"){
          h_tt_1695 = (TH1F*) hist_el->Clone();
          h_tt_1695->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1715"){
          h_tt_1715 = (TH1F*) hist_el->Clone();
          h_tt_1715->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1735"){
          h_tt_1735 = (TH1F*) hist_el->Clone();
          h_tt_1735->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1755"){
          h_tt_1755 = (TH1F*) hist_el->Clone();
          h_tt_1755->Add(hist_mu);
        }
        else cout << "process " << process << " not known!!" << endl;
      }
      else{
        if(process == "TTbar"){
          h_tt->Add(hist_el);
          h_tt->Add(hist_mu);
        }
        else if(process == "SingleTop"){
          h_st->Add(hist_el);
          h_st->Add(hist_mu);
        }
        else if(process == "WJets"){
          h_wj->Add(hist_el);
          h_wj->Add(hist_mu);
        }
        else if(process == "other"){
          h_ot->Add(hist_el);
          h_ot->Add(hist_mu);
        }
        else if(process == "DATA"){
          h_da->Add(hist_el);
          h_da->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1695"){
          h_tt_1695->Add(hist_el);
          h_tt_1695->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1715"){
          h_tt_1715->Add(hist_el);
          h_tt_1715->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1735"){
          h_tt_1735->Add(hist_el);
          h_tt_1735->Add(hist_mu);
        }
        else if(process == "TTbar_mtop1755"){
          h_tt_1755->Add(hist_el);
          h_tt_1755->Add(hist_mu);
        }
        else cout << "process " << process << " not known!!" << endl;
      }
    }
  }
  
  h_da->Add(h_st, -1);
  h_da->Add(h_wj, -1);
  h_da->Add(h_ot, -1);
  
  TH2F* cov = CreateMatrix(h_da);
  
  TH1F* h_da_norm = NormalizeHist(h_da);
  TH1F* h_tt_norm = NormalizeHist(h_tt);
  TH1F* h_tt1695_norm = NormalizeHist(h_tt_1695);
  TH1F* h_tt1715_norm = NormalizeHist(h_tt_1715);
  TH1F* h_tt1735_norm = NormalizeHist(h_tt_1735);
  TH1F* h_tt1755_norm = NormalizeHist(h_tt_1755);
  
  TH2F* cov_norm = NormalizeMatrix(cov, h_da);
  
  TH1F* h_da_norm_tot = SetError(h_da_norm, cov_norm);
  
  
  outfile->cd();
  h_da_norm_tot->Write("data");
  h_tt_norm->Write("mt1725");
  h_tt1695_norm->Write("mt1695");
  h_tt1715_norm->Write("mt1715");
  h_tt1735_norm->Write("mt1735");
  h_tt1755_norm->Write("mt1755");
  cov_norm->Write("cov");
  outfile->Close();
  
  PlotMatrix(cov, "MtopReco_cov");
  PlotMatrix(cov_norm, "MtopReco_cov_norm");
  
  vector<TString> legnames = {"#it{m}_{t} = 169.5 GeV", "#it{m}_{t} = 172.5 GeV", "#it{m}_{t} = 175.5 GeV"};
  vector<TH1F*> h_mtop_forplot = {h_tt_1695, h_tt, h_tt_1755};
  vector<TH1F*> h_mtop_norm_forplot = {h_tt1695_norm, h_tt_norm, h_tt1755_norm};
  Plot_hist(h_da, h_mtop_forplot, legnames, "MtopReco_hist");
  Plot_hist(h_da_norm_tot, h_mtop_norm_forplot, legnames, "MtopReco_hist_norm");
  
  
  vector<TH1F*> h_mtop_norm = {h_tt1695_norm, h_tt1715_norm, h_tt_norm, h_tt1735_norm, h_tt1755_norm};
  vector<double> masses = {169.5, 171.5, 172.5, 173.5, 175.5};
  double lower = xmin;
  double upper = xmax;
  bool NormToWidth = true;
  
  
  
  chi2fit* chi2;
  std::vector<double> chi2values;
  TF1* chi2_fitfunction;
  chi2 = new chi2fit(h_da_norm_tot, cov_norm, h_mtop_norm, masses, lower, upper, NormToWidth);
  chi2->CalculateChi2();
  chi2values = chi2->GetChi2Values();
  chi2_fitfunction = chi2->GetChi2Fit();
  TGraph* chi2_graph = chi2->GetChi2Graph();
  cout << " MASS = " << chi2->GetMass() << " +- " << chi2->GetUncertainty() << std::endl;
  double totuncert = chi2->GetUncertainty();
  Plot_chi2(chi2_fitfunction, masses, chi2values, chi2->GetMass(), chi2->GetUncertainty(), "chi2fit");
  
  
  return 0;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


TH1F* get_hist(TFile* file, TH1F* dummy, TString sel_name, TString weightname){

  TTree* tree = (TTree *) file->Get("AnalysisTree");

  TH1F* hist = (TH1F*) dummy->Clone();

  Double_t mjet;

  Bool_t passed_selection;
  Double_t weight, rec_weight, gen_weight;
  Float_t additional_factor;


  tree->ResetBranchAddresses();
  tree->SetBranchAddress(sel_name, &passed_selection);
  tree->SetBranchAddress("rec_weight",&rec_weight);
  tree->SetBranchAddress("gen_weight",&gen_weight);
  tree->SetBranchAddress("Mass_Rec", &mjet);

  if(weightname != "none"){
    cout << "  - using weight" << endl;
    tree->SetBranchAddress(weightname,&additional_factor);
  }
  tree->SetBranchStatus("*",1);


  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(weightname != "none") rec_weight *= additional_factor;
    weight = rec_weight * gen_weight;
    if(passed_selection){
      hist->Fill(mjet, weight);
    }
  }
  return hist;
}

TH2F* NormalizeMatrix(TH2F*old_cov, TH1F* hist_){
  int nbins = hist_->GetXaxis()->GetNbins();
  TH2F* new_cov = (TH2F*) old_cov->Clone();
  new_cov->Reset();
  double integral = hist_->Integral();
  cout << "Nbins = " << nbins << ", integral = " << integral << endl;
  // std::cout << "lower bin: "<< lower_ << std::endl;
  // std::cout << "upper bin: "<< upper_ << std::endl;
  for(int i=1; i <= nbins; i++){
    for(int j=1; j <= nbins; j++){
      double sum = 0;
      double derivation_i = 0;
      double derivation_j = 0;
      double old_entry = 0;
      for(int k=1; k <= nbins; k++){
        for(int l=1; l <= nbins; l++){
          old_entry = old_cov->GetBinContent(k,l);
          // double binwidth_i = hist_->GetBinWidth(i);
          // double binwidth_j = hist_->GetBinWidth(j);
          double binwidth_i = 1;
          double binwidth_j = 1;
          if(i==k) derivation_i = (integral - hist_->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (hist_->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - hist_->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (hist_->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      // std::cout << sum << std::endl;
      new_cov->SetBinContent(i, j, sum);
      // if(i==20 || j==20)new_cov->SetBinContent(i, j, 0.0);
    }
  }
  return new_cov;
}


TH1F* NormalizeHist(TH1F* old_hist){
  TH1F* new_hist = (TH1F*) old_hist->Clone();
  // new_hist->Scale(1/old_hist->Integral(),"width");
  new_hist->Scale(1/old_hist->Integral());
  
  return new_hist;
}

TH2F* CreateMatrix(TH1F* hist){
  TH2F* cov = new TH2F("mjet", "mjet", Nbins, xmin, xmax, Nbins, xmin, xmax);
  cov->Reset();
  for(unsigned int i=1; i<=Nbins; i++){
    for(unsigned int j=1; j<=Nbins; j++){
      if(i==j){
        double error = hist->GetBinError(i);
        // cout << i << ", " << j << ": " << error << endl;
        cov->SetBinContent(i,j, error*error);
      }
      else{
        cov->SetBinContent(i,j,0.0);
      }
    }
  }
  return cov;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Chi2 Stuff

chi2fit::chi2fit(TH1F* data, TH2F* cov_matrix, std::vector<TH1F*> masspoints, std::vector<double> masses, double lower, double upper, bool width){
  width_ = width;
  SetData(data);
  SetCovMatrix(cov_matrix);
  SetMasspoints(masspoints);
  SetMasses(masses);
  SetRange(lower, upper);
  SetBins(masspoints);
}


void chi2fit::SetMasspoints(std::vector<TH1F*> masspoints){
  for(unsigned int i = 0; i<masspoints.size(); i++) {
      TH1F *mass = (TH1F*)masspoints[i]->Clone();
      masspoints_.push_back(mass);
  }
  return;
}

void chi2fit::SetData(TH1F* data){
  data_ = (TH1F*)data->Clone("data_");
  return;
}

void chi2fit::SetCovMatrix(TH2F* cov_matrix){
  cov_ = (TH2F*)cov_matrix->Clone("cov_");
  return;
}


void chi2fit::SetRange(double lower, double upper){
  lower_ = data_->GetXaxis()->FindBin(lower);
  upper_ = data_->GetXaxis()->FindBin(upper-1);
  cout << "lower bin = " << lower_ << ", upper bin = " << upper_ << endl;
  return;
}

void chi2fit::SetBins(std::vector<TH1F*> masspoints){
  // when calculating the chi2, one degree of freedom is lost by normalisation
  // therefore covariance matrix is only invertable if one excludes one bin from chi2
  // this bin is selected randomly (chi2 should not change depending on the bin)
  // for every masspoint, the bins that should be used are written into a vector
  bool skip;
  int rand_bin = rand()%(upper_-lower_+1)+lower_;
  // int rand_bin = 20;
  std::cout << "Exclude bin " << rand_bin << " from chi2" << std::endl;
  for(unsigned int j=0; j<masspoints.size(); j++){
    std::vector<int> bins_this;
    for(int i=lower_; i<=upper_; i++){
      if(i == rand_bin) skip = true;
      else skip = false;
      if(!skip) bins_this.push_back(i);
    }
    bins_.push_back(bins_this);
  }
  return;
}


void chi2fit::SetMasses(std::vector<double> masses){
  masses_.ResizeTo(masses.size());
  for(unsigned int i=0; i<masses.size(); i++) masses_[i] = masses[i];
  return;
}

void chi2fit::CalculateChi2(){
  int j = 0;
  const int NMasses = masspoints_.size();
  chi2_.ResizeTo(NMasses);
  for(auto &masspoint: masspoints_){
    double chi2 = ComputeChi2(masspoint, bins_[j]);
    chi2_[j] = chi2;
    j++;
  }

  TGraph* chi2Hist = new TGraph(masses_, chi2_);
  TF1*f1 = new TF1("f1","pol2",0,500);
  chi2Hist->Fit("f1","R");
  fit_ = chi2Hist->GetFunction("f1");

  return;
}

double chi2fit::GetMass(){
  double minY = fit_->GetMinimum();
  double mass = fit_->GetX(minY, 140, 200);
  return mass;
}

double chi2fit::GetMin(){
  double minY = fit_->GetMinimum();
  return minY;
}

double chi2fit::GetUncertainty(){
  double minY = fit_->GetMinimum();
  double minX = fit_->GetX(minY, 140, 200);
  double minX_sigmaup = fit_->GetX(minY+1, minX, 200);
  // double minX_sigmadown = fit_->GetX(minY+1, 140, minX);
  double uncert = minX_sigmaup-minX;
  return uncert;
}


TGraph * chi2fit::GetChi2Graph(){
  TGraph* chi_hist = new TGraph(masses_,chi2_);
  return chi_hist;
}

std::vector<double> chi2fit::GetChi2Values(){
  std::vector<double> chi2;
  int N = chi2_.GetNoElements();
  for(int i=0; i<N; i++) chi2.push_back(chi2_[i]);
  return chi2;
}


TF1* chi2fit::GetChi2Fit(){
  return fit_;
}


//
// compute chi^2
//
double chi2fit::ComputeChi2(TH1F* MC, std::vector<int> bins){
  //double chi2 = 0;
  TH1F* data_sc = (TH1F*)data_->Clone("data_sc");
  TH2F* cov_sc = (TH2F*)cov_->Clone("cov_sc");
  TH1F* MC_sc = (TH1F*)MC->Clone("MC_sc");

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
      // cout << "cov[" << i << "," << j << "] = " << cov_sc->GetBinContent(bins.at(i),bins.at(j)) << endl;
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
    std::cerr << "FcnData::RecalculateInvCovMa. Error inverting the mapped covariance matrix. Aborting." << std::endl;
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
      if ( mat[i][i] < 1.e-20 ) printf("   <<<<<<<< \n");
      else printf("\n");
    }
    exit(3);
  }

  // calculate chi2
  double chi2 = 0;

  TVectorD diff = vdata - vMC;
  TVectorD right = imat * diff;
  chi2 = diff * right;
  std::cout << "Chi2 = " << chi2 << std::endl;
  return chi2;
}


void Plot_chi2(TF1 * fit_, std::vector<double> masses_, std::vector<double> chi2_, double mass, double uncert, TString file_name){
  TString save_path = get_save_path();

  TF1 * fit = (TF1*)fit_->Clone("fit");
  TVectorD masses(masses_.size());
  TVectorD chi2(chi2_.size());
  for(int i=0; i<masses_.size(); i++) masses[i] = masses_[i];
  for(int i=0; i<chi2_.size(); i++) chi2[i] = chi2_[i];

  TGraph* chi_hist = new TGraph(masses,chi2);
  TCanvas *c = new TCanvas("Chi2", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  chi_hist->SetTitle(" ");
  chi_hist->GetXaxis()->SetTitle("m_{t} [GeV]");
  chi_hist->GetYaxis()->SetTitle("#chi^{2}");
  chi_hist->GetYaxis()->SetTitleOffset(1.1);
  chi_hist->GetXaxis()->SetTitleOffset(0.9);
  chi_hist->GetYaxis()->SetTitleSize(0.05);
  chi_hist->GetXaxis()->SetTitleSize(0.05);
  chi_hist->GetXaxis()->SetNdivisions(505);
  chi_hist->GetYaxis()->SetNdivisions(505);
  chi_hist->SetMarkerStyle(20);
  chi_hist->SetMarkerSize(1.5);
  chi_hist->SetLineColor(1);
  chi_hist->Draw("AP");
  fit->Draw("SAME");

  // write extracted mass value into plot
  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextFont(43);
  text.SetTextSize(18);
  char mass_text[32];
  sprintf(mass_text, "%.5g", mass);
  char uncert_text[32];
  if(uncert < 1) sprintf(uncert_text, "%.3g", uncert);
  else           sprintf(uncert_text, "%.4g", uncert);
  TString masstext = "m_{t} = ";
  masstext += mass_text;
  masstext += " #pm ";
  masstext += uncert_text;
  text.DrawLatex(.4,.6, masstext);
  c->SaveAs(save_path+"/Plots/MtopReco_chi2.pdf");
  delete c;
  return;
}

void PlotMatrix(TH2F* cor, TString name){
  TString save_path = get_save_path();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(7);
  gStyle->SetPalette(kSunset);

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetBottomMargin(0.13);
  gPad->SetLeftMargin(0.135);
  gPad->SetRightMargin(0.2);

  cor->SetTitle(" ");

  cor->GetXaxis()->SetTitle("#it{m}_{ jet} [GeV]");
  cor->GetXaxis()->SetTitleOffset(1.1);
  cor->GetXaxis()->SetTitleSize(0.05);
  cor->GetXaxis()->SetLabelSize(0.04);
  cor->GetXaxis()->SetLabelOffset(0.012);
  cor->GetXaxis()->SetNdivisions(506);

  cor->GetYaxis()->SetTitle("#it{m}_{ jet} [GeV]");
  cor->GetYaxis()->SetTitleOffset(1.3);
  cor->GetYaxis()->SetTitleSize(0.05);
  cor->GetYaxis()->SetLabelSize(0.04);
  cor->GetYaxis()->SetLabelOffset(0.012);
  cor->GetYaxis()->SetNdivisions(506);

  // cor->GetZaxis()->SetRangeUser(-1,1);
  cor->GetZaxis()->SetLabelOffset(0.012);

  cor->SetLineColor(1);
  cor->SetLineWidth(2);
  

  cor->Draw("COLZ");
  cor->Draw("BOX SAME");


  gPad->RedrawAxis();
  c->SaveAs(save_path+"/Plots/"+name+".pdf");
}


void Plot_hist(TH1F* h_data_XS_tot, vector<TH1F*> truth, vector<TString> legnames, TString name){
  TString save_path = get_save_path();

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(7);

  bool doTicks = false;

  TString o_tot = "E1";

  if(!doTicks){
    o_tot = "E";
  }

  int titlefont = 43;
  double ymax = h_data_XS_tot->GetMaximum();
  for(auto h: truth){
    if(h->GetMaximum() > ymax) ymax = h->GetMaximum();
  }

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetBottomMargin(0.13); 
  gPad->SetLeftMargin(0.21);
  TGaxis::SetMaxDigits(4);


  h_data_XS_tot->SetTitle(" ");
  h_data_XS_tot->GetYaxis()->SetRangeUser(0, 1.4*ymax);
  h_data_XS_tot->GetYaxis()->SetTitle("a.u.");
  h_data_XS_tot->GetYaxis()->SetTitleOffset(1.32);
  h_data_XS_tot->GetYaxis()->SetTitleSize(0.07);
  h_data_XS_tot->GetYaxis()->SetNdivisions(505);
  h_data_XS_tot->GetXaxis()->SetNdivisions(505);
  h_data_XS_tot->GetXaxis()->SetTitle("#it{m}_{ jet} [GeV]");
  h_data_XS_tot->GetXaxis()->SetTitleOffset(1.1);
  h_data_XS_tot->GetXaxis()->SetTitleSize(0.05);
  h_data_XS_tot->SetLineColor(kBlack);
  h_data_XS_tot->SetMarkerColor(kBlack);
  h_data_XS_tot->SetMarkerStyle(8);
  h_data_XS_tot->SetMarkerSize(1);
  h_data_XS_tot->Draw(o_tot);
  
  
  Int_t color[] = {kRed-4, TColor::GetColor("#0059b3"), 14};
  Int_t style[] = {3, 1, 2};

  for(unsigned int i=0; i<truth.size(); i++){
    truth[i]->SetLineWidth(3);
    truth[i]->SetLineColor(color[i]);
    truth[i]->SetLineStyle(style[i]);
    truth[i]->Draw("HIST SAME");
  }

  h_data_XS_tot->Draw(o_tot + " SAME");

  TLegend *l=new TLegend(0.53,0.55,0.87,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(h_data_XS_tot,"Data","ple");
  for(unsigned int i=0; i<truth.size(); i++){
    l->AddEntry(truth[i],legnames[i],"fl");
  }
  l->SetTextSize(0.04);
  l->Draw();
  gPad->RedrawAxis();


  c->SaveAs(save_path+"/Plots/"+name+".pdf");
  delete c;
  return;
}


TH1F* SetError(TH1F* data, TH2F* cov){
  TH1F* hist = (TH1F*) data->Clone();
  int nbins = hist->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    double error = sqrt(cov->GetBinContent(i,i));
    hist->SetBinError(i, error);
  }
  return hist;
}
