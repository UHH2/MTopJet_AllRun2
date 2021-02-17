#include "../include/CentralInclude.h"

using namespace std;
TH1F* GetRatio(TH1F* h1, TH1F* h2);
TH1F* GetRatioUncert(TH1F* mc);

TH1F* ConvertToCrossSection(TH1F* hist);
TH2F* ConvertToCrossSection(TH1F* hist, TH2F* cov);
TH1F* SetError(TH1F* data, TH2F* cov);
TH2F* NormalizeMatrix(TH2F* old_cov, TH1F* hist_);
TH1F* NormalizeHist(TH1F* old_hist);
void Plot(TH1F* h_data_XS_tot, TH1F* h_data_XS_stat, vector<TH1F*> truth, vector<TH1F*> trut2, vector<TString> legnames, TH1F* h_ratio_data, TH1F* h_ratio_data_stat, TH1F* h_ratio_powheg);
void Plot_norm(TH1F* h_data_XS_tot, TH1F* h_data_XS_stat, vector<TH1F*> truth, vector<TH1F*> trut2, vector<TString> legnames, TH1F* h_ratio_data, TH1F* h_ratio_data_stat, vector<TH1F*> ratios_mtop);
void Plot_chi2(TF1 * fit_, std::vector<double> masses_, std::vector<double> chi2_, double mass, double uncert, TString file_name);

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

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

TString mode = "data";

int main(int argc, char* argv[]){

  if(argc > 1){
    if(strcmp(argv[1], "pseudo1")==0) mode = "pseudo1";
    else if(strcmp(argv[1], "pseudo2")==0) mode = "pseudo2";
  }

  cout << "You selected mode = " << mode <<endl;

  TFile* file16 = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_"+mode+"_2016_combine.root");
  TFile* file17 = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_"+mode+"_2017_combine.root");
  TFile* file18 = new TFile("/nfs/dust/cms/user/schwarzd/TUnfold/Results_"+mode+"_2018_combine.root");


  TH1F * h_data16 = (TH1F*)file16->Get("Unfold_events");
  TH1F * h_data17 = (TH1F*)file17->Get("Unfold_events");
  TH1F * h_data18 = (TH1F*)file18->Get("Unfold_events");
  TH1F * h_data = (TH1F*) h_data16->Clone();
  h_data->Add(h_data17);
  h_data->Add(h_data18);
  TH1F * h_data_XS = ConvertToCrossSection(h_data);

  TH1F * h_truth16 = (TH1F*)file16->Get("Truth_error");
  TH1F * h_truth17 = (TH1F*)file17->Get("Truth_error");
  TH1F * h_truth18 = (TH1F*)file18->Get("Truth_error");
  TH1F * h_truth = (TH1F*) h_truth16->Clone();
  h_truth->Add(h_truth17);
  h_truth->Add(h_truth18);
  TH1F * h_truth_XS = ConvertToCrossSection(h_truth);

  TH2F * cov_tot16 = (TH2F*)file16->Get("CovTotal");
  TH2F * cov_tot17 = (TH2F*)file17->Get("CovTotal");
  TH2F * cov_tot18 = (TH2F*)file18->Get("CovTotal");
  TH2F * cov_tot = (TH2F*) cov_tot16->Clone();
  cov_tot->Add(cov_tot17);
  cov_tot->Add(cov_tot18);
  TH2F * cov_tot_XS = ConvertToCrossSection(h_data_XS, cov_tot);

  TH2F * cov_theo16 = (TH2F*)file16->Get("CovTheo");
  TH2F * cov_theo17 = (TH2F*)file17->Get("CovTheo");
  TH2F * cov_theo18 = (TH2F*)file18->Get("CovTheo");
  TH2F * cov_theo = (TH2F*) cov_theo16->Clone();
  cov_theo->Add(cov_theo17);
  cov_theo->Add(cov_theo18);
  TH2F * cov_theo_XS = ConvertToCrossSection(h_data_XS, cov_theo);
  TH2F * cov_theo_norm = NormalizeMatrix(cov_theo,h_data);

  // Now add Total and Theo for chi2
  TH2F * cov_total_theo = (TH2F*) cov_tot->Clone();
  cov_total_theo->Add(cov_theo);
  TH2F * cov_total_theo_norm = NormalizeMatrix(cov_total_theo, h_data);


  TH2F * cov_stat16 = (TH2F*)file16->Get("CovStat");
  TH2F * cov_stat17 = (TH2F*)file17->Get("CovStat");
  TH2F * cov_stat18 = (TH2F*)file18->Get("CovStat");
  TH2F * cov_stat = (TH2F*) cov_stat16->Clone();
  cov_stat->Add(cov_stat17);
  cov_stat->Add(cov_stat18);
  TH2F * cov_stat_XS = ConvertToCrossSection(h_data_XS, cov_stat);

  TH1F* h_data_XS_tot = SetError(h_data_XS, cov_tot_XS);
  TH1F* h_data_XS_stat = SetError(h_data_XS, cov_stat_XS);

  TH1F * h_ratio_data = GetRatioUncert(h_data_XS_tot);
  TH1F * h_ratio_data_stat = GetRatioUncert(h_data_XS_stat);
  TH1F * h_ratio_powheg = GetRatio(h_truth_XS, h_data_XS_tot);

  // get mtop samples (norm)
  vector<double> masses = {169.5, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5, 175.0, 175.5};
  vector<TString> mass_names = {"1665", "1700", "1705", "1710", "1715", "1720", "1725", "1730", "1735", "1740", "1745", "1750", "1755"};
  vector<TH1F*>  h_mtop_norm;
  for(auto m: mass_names){
    TString histname = "mc_mtop"+m+"_truth";
    if(m=="1725") histname = "mc_truth";
    TH1F * h_16 = (TH1F*)file16->Get(histname);
    TH1F * h_17 = (TH1F*)file17->Get(histname);
    TH1F * h_18 = (TH1F*)file18->Get(histname);
    TH1F * h_all = (TH1F*) h_16->Clone();
    h_all->Add(h_17);
    h_all->Add(h_18);
    TH1F* h_all_norm = NormalizeHist(h_all);
    TH1F* h_all_norm_sys = SetError(h_all_norm, cov_theo_norm);
    h_mtop_norm.push_back(h_all_norm_sys);
  }

  // Normalize
  TH2F* cov_tot_norm = NormalizeMatrix(cov_tot,h_data);
  TH2F* cov_stat_norm = NormalizeMatrix(cov_stat,h_data);
  TH1F* h_data_norm = NormalizeHist(h_data);
  TH1F* h_data_norm_tot = SetError(h_data_norm, cov_tot_norm);
  TH1F* h_data_norm_stat = SetError(h_data_norm, cov_stat_norm);
  TH1F* h_truth_norm = NormalizeHist(h_truth);
  TH1F * h_ratio_data_norm = GetRatioUncert(h_data_norm_tot);
  TH1F * h_ratio_data_norm_stat = GetRatioUncert(h_data_norm_stat);
  TH1F * h_ratio_powheg_norm = GetRatio(h_truth_norm, h_data_norm_tot);
  TH1F * h_ratio_mtop1695_norm = GetRatio(h_mtop_norm[0], h_data_norm_tot);
  TH1F * h_ratio_mtop1725_norm = GetRatio(h_mtop_norm[6], h_data_norm_tot);
  TH1F * h_ratio_mtop1755_norm = GetRatio(h_mtop_norm[12], h_data_norm_tot);



  std::vector<TH1F*> truth, truth2;
  for(auto t: {h_truth_XS}){
    truth.push_back( (TH1F*) t->Clone() );
    truth2.push_back( (TH1F*) t->Clone() );
  }
  vector<TString> legnames = {"POWHEG"};

  vector<TH1F*> mtop, mtop2, ratios_mtop;
  mtop.push_back( (TH1F*) h_mtop_norm[0]->Clone());
  mtop.push_back( (TH1F*) h_mtop_norm[6]->Clone());
  mtop.push_back( (TH1F*) h_mtop_norm[12]->Clone());
  mtop2.push_back( (TH1F*) h_mtop_norm[0]->Clone());
  mtop2.push_back( (TH1F*) h_mtop_norm[6]->Clone());
  mtop2.push_back( (TH1F*) h_mtop_norm[12]->Clone());
  ratios_mtop.push_back( (TH1F*) h_ratio_mtop1695_norm->Clone());
  ratios_mtop.push_back( (TH1F*) h_ratio_mtop1725_norm->Clone());
  ratios_mtop.push_back( (TH1F*) h_ratio_mtop1755_norm->Clone());



  vector<TString> legnames_norm = {"#it{m}_{t} = 169.5 GeV", "#it{m}_{t} = 172.5 GeV", "#it{m}_{t} = 175.5 GeV"};

  Plot(h_data_XS_tot, h_data_XS_stat, truth, truth2, legnames, h_ratio_data, h_ratio_data_stat, h_ratio_powheg);
  Plot_norm(h_data_norm_tot, h_data_norm_stat, mtop, mtop2, legnames_norm, h_ratio_data, h_ratio_data_stat, ratios_mtop);

  cout << "*******************************" << endl;
  cout << "************ chi 2 ************" << endl;
  cout << "*******************************" << endl;
  double lower = 112;
  double upper = 230;
  bool NormToWidth = true;
  chi2fit* chi2;
  std::vector<double> chi2values;
  TF1* chi2_fitfunction;
  chi2 = new chi2fit(h_data_norm_tot, cov_total_theo_norm, h_mtop_norm, masses, lower, upper, NormToWidth);
  chi2->CalculateChi2();
  chi2values = chi2->GetChi2Values();
  chi2_fitfunction = chi2->GetChi2Fit();
  cout << " MASS = " << chi2->GetMass() << " +- " << chi2->GetUncertainty() << std::endl;

  Plot_chi2(chi2_fitfunction, masses, chi2values, chi2->GetMass(), chi2->GetUncertainty(), "chi2fit");



  return 0;
}


TH1F* GetRatio(TH1F* h1, TH1F* h2){
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, 1);
      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      // double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      double error = E1/N2; // only consider uncert from MC
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}

TH1F* GetRatioUncert(TH1F* mc){
  TH1F* mc_uncert = (TH1F*) mc->Clone();
  int Nbins = mc->GetSize() - 2;
  for(int i=1; i<=Nbins;i++){
    mc_uncert->SetBinContent(i, 1.0);
    double central = mc->GetBinContent(i);
    double error = mc->GetBinError(i);
    double error_ratio = error/central;
    mc_uncert->SetBinError(i, error_ratio);
  }
  return mc_uncert;
}

TH1F* ConvertToCrossSection(TH1F* hist){
  double lumi = 137.1;
  TH1F* newhist = (TH1F*) hist->Clone();
  int nbins = hist->GetXaxis()->GetNbins();
  for(int bin=1; bin<=nbins; bin++){
    double events = hist->GetBinContent(bin);
    double error = hist->GetBinError(bin);
    double binwidth = hist->GetBinWidth(bin);
    double cs = events/(binwidth*lumi);
    double er = error/(binwidth*lumi);
    newhist->SetBinContent(bin, cs);
    newhist->SetBinError(bin, er);
  }
  return newhist;
}

TH2F* ConvertToCrossSection(TH1F* hist, TH2F* cov){
  double lumi = 137.1;
  TH1F* newhist = (TH1F*) hist->Clone();
  TH2F* newcov = (TH2F*) cov->Clone();
  newcov->Reset();
  int nbins = hist->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    for(int j=1; j<=nbins; j++){
      double binwidth_i = hist->GetBinWidth(i);
      double binwidth_j = hist->GetBinWidth(j);
      double old_entry = cov->GetBinContent(i,j);
      double new_entry = old_entry/(binwidth_i*binwidth_j*lumi*lumi);
      newcov->SetBinContent(i,j, new_entry);
    }
  }
  return newcov;
}

TH1F* SetError(TH1F* data, TH2F* cov){
  TH1F* hist = (TH1F*) data->Clone();
  int nbins = hist->GetXaxis()->GetNbins();
  for(int i=1; i<=nbins; i++){
    double error = sqrt(cov->GetBinContent(i,i));
    if(error < 0.00001) error = 0.00001;
    hist->SetBinError(i, error);
  }
  return hist;
}

TH2F* NormalizeMatrix(TH2F*old_cov, TH1F* hist_){
  int nbins = hist_->GetXaxis()->GetNbins();
  TH2F* new_cov = (TH2F*) old_cov->Clone();
  new_cov->Reset();
  double integral = hist_->Integral();
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
          double binwidth_i = hist_->GetBinWidth(i);
          double binwidth_j = hist_->GetBinWidth(j);
          if(i==k) derivation_i = (integral - hist_->GetBinContent(i)) / (pow(integral,2)) * (1/binwidth_i);
          else     derivation_i = - (hist_->GetBinContent(i))          / (pow(integral,2)) * (1/binwidth_i);
          if(j==l) derivation_j = (integral - hist_->GetBinContent(j)) / (pow(integral,2)) * (1/binwidth_j);
          else     derivation_j = - (hist_->GetBinContent(j))          / (pow(integral,2)) * (1/binwidth_j);
          sum += derivation_i * derivation_j * old_entry;
        }
      }
      // std::cout << sum << std::endl;
      new_cov->SetBinContent(i, j, sum);
    }
  }
  return new_cov;
}


TH1F* NormalizeHist(TH1F* old_hist){
  TH1F* new_hist = (TH1F*) old_hist->Clone();
  new_hist->Scale(1/old_hist->Integral(),"width");
  return new_hist;
}

void Plot(TH1F* h_data_XS_tot, TH1F* h_data_XS_stat, vector<TH1F*> truth, vector<TH1F*> truth2, vector<TString> legnames, TH1F* h_ratio_data, TH1F* h_ratio_data_stat, TH1F* h_ratio_powheg){
  TString save_path = get_save_path();

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(5);

  bool doTicks = true;

  TString o_tot = "E1";
  TString o_stat = "E1";

  if(!doTicks){
    o_tot = "E";
    o_stat = "E1 X0";
  }

  int titlefont = 43;
  double ymax = 22.0;

  TCanvas *c = new TCanvas("c","",600,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();
  TGaxis::SetMaxDigits(3);


  h_data_XS_tot->SetTitle(" ");
  h_data_XS_tot->GetYaxis()->SetRangeUser(0, ymax);
  h_data_XS_tot->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{m}_{jet}} #left[#frac{fb}{GeV}#right]");
  h_data_XS_tot->GetYaxis()->SetTitleOffset(1.3);
  h_data_XS_tot->GetYaxis()->SetTitleSize(0.06);
  h_data_XS_tot->GetYaxis()->SetNdivisions(505);
  h_data_XS_tot->SetLineColor(kBlack);
  h_data_XS_tot->SetMarkerColor(kBlack);
  h_data_XS_tot->SetMarkerStyle(8);
  h_data_XS_tot->SetMarkerSize(1);
  h_data_XS_tot->Draw(o_tot);
  h_data_XS_stat->SetLineColor(kBlack);
  h_data_XS_stat->SetMarkerColor(kBlack);
  h_data_XS_stat->SetMarkerStyle(8);
  h_data_XS_stat->SetMarkerSize(1);
  Int_t color[] = {TColor::GetColor("#0059b3"), TColor::GetColor("#e67300")};
  Int_t fillcolor[] = {TColor::GetColor("#99ccff"), TColor::GetColor("#ff8000")};
  Int_t style[] = {1, 2};
  Int_t fillstyle[] = {3144, 3153};

  for(unsigned int i=0; i<truth.size(); i++){
    if(i != 0) continue; //skip aMCatNLO
    truth[i]->SetLineWidth(3);
    truth[i]->SetLineColor(color[i]);
    truth[i]->SetLineStyle(style[i]);
    truth[i]->SetMarkerStyle(0);
    truth[i]->SetFillColor(fillcolor[i]);
    if(i==0) truth[i]->Draw("E2 SAME");
    truth2[i]->SetLineWidth(3);
    truth2[i]->SetLineColor(color[i]);
    truth2[i]->SetLineStyle(style[i]);
    truth2[i]->Draw("HIST SAME");
  }
  h_data_XS_stat->Draw(o_stat + " SAME");
  h_data_XS_tot->Draw(o_tot + " SAME");

  TLegend *l=new TLegend(0.55,0.62,0.82,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(h_data_XS_tot,"Data","ple");
  for(unsigned int i=0; i<truth.size(); i++){
    if(i != 0) continue; //skip aMCatNLO
    if(i==0) l->AddEntry(truth[i],legnames[i],"fl");
    else     l->AddEntry(truth[i],legnames[i],"l");
  }
  l->SetTextSize(0.052);
  l->Draw();
  gPad->RedrawAxis();

  TLatex *cmstext = new TLatex(3.5, 24, "CMS");
  cmstext->SetNDC();
  cmstext->SetTextAlign(13);
  cmstext->SetTextFont(62);
  cmstext->SetTextSize(0.08);
  cmstext->SetX(0.24);
  cmstext->SetY(0.84);
  cmstext->Draw();

  TString simtext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, simtext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.24);
  text3->SetTextFont(52);
  text3->SetTextSize(0.06);
  text3->SetY(0.77);
  text3->Draw();

  TString infotext = "137.1 fb^{-1} (13 TeV)";
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(0.9);
  text1->SetY(0.961);
  text1->SetTextSize(0.055);
  text1->Draw();

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  h_data_XS_tot->GetXaxis()->SetLabelSize(0.);
  h_data_XS_tot->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( 112, 0, 112, ymax, 0, ymax, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(titlefont);
  axis->SetLabelSize(21);
  axis->Draw();

  // Ratio Plot kommt in unteres pad
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.38);
  pad2->Draw();
  pad2->cd();

  h_ratio_data->GetXaxis()->SetTickLength(0.07);
  h_ratio_data->GetXaxis()->SetTitleSize(25);
  h_ratio_data->GetXaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetXaxis()->SetTitleOffset(4.0);
  h_ratio_data->GetXaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetXaxis()->SetLabelSize(21);
  h_ratio_data->GetXaxis()->SetLabelOffset(0.035);
  h_ratio_data->GetYaxis()->SetTitle("#frac{Theory}{Data}");
  h_ratio_data->GetYaxis()->CenterTitle();
  h_ratio_data->GetYaxis()->SetTitleSize(22);
  h_ratio_data->GetYaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetYaxis()->SetTitleOffset(2.2);
  h_ratio_data->GetYaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetYaxis()->SetLabelSize(19);
  h_ratio_data->GetYaxis()->SetLabelOffset(0.009);
  h_ratio_data->GetYaxis()->SetNdivisions(505);
  h_ratio_data->SetTitle(" ");
  h_ratio_data->GetYaxis()->SetRangeUser(0.2, 1.8);
  h_ratio_data->GetXaxis()->SetTitle("#it{m}_{jet} #left[GeV#right]");
  h_ratio_data->SetLineColor(kBlack);
  h_ratio_data->SetMarkerColor(kBlack);
  h_ratio_data->SetMarkerStyle(8);
  h_ratio_data->SetMarkerSize(1);
  h_ratio_data->Draw(o_tot);
  h_ratio_data_stat->SetLineColor(kBlack);
  h_ratio_data_stat->SetMarkerColor(kBlack);
  h_ratio_data_stat->SetMarkerStyle(8);
  h_ratio_data_stat->SetMarkerSize(1);

  TH1F* h_ratio_powheg_line = (TH1F*) h_ratio_powheg->Clone();
  h_ratio_powheg->SetMarkerSize(0);
  h_ratio_powheg->SetFillColor(fillcolor[0]);
  h_ratio_powheg->Draw("E2 SAME");
  // h_ratio_amc->SetLineWidth(3);
  // h_ratio_amc->SetLineColor(color[1]);
  // h_ratio_amc->SetLineStyle(style[1]);
  // h_ratio_amc->Draw("HIST SAME");
  h_ratio_data->Draw(o_tot + " SAME");

  gStyle->SetEndErrorSize(5);

  h_ratio_powheg_line->SetLineColor(color[0]);
  h_ratio_powheg_line->SetLineWidth(3);
  h_ratio_powheg_line->SetLineStyle(style[0]);


  // noch einmal alles in der richtigen Reihenfolge zeichnen
  h_ratio_powheg->Draw("E2 SAME");
  h_ratio_powheg_line->Draw("HIST SAME");
  // h_ratio_amc->Draw("HIST SAME");
  h_ratio_data_stat->Draw(o_stat + " SAME");
  h_ratio_data->Draw(o_tot + " SAME");
  gPad->RedrawAxis();

  c->SaveAs(save_path+"/Plots/Unfolding_Run2/Unfold_combination_"+mode+".pdf");
  delete c;
  return;
}

void Plot_norm(TH1F* h_data_XS_tot, TH1F* h_data_XS_stat, vector<TH1F*> truth, vector<TH1F*> truth2, vector<TString> legnames, TH1F* h_ratio_data, TH1F* h_ratio_data_stat, vector<TH1F*> ratios_mtop){
  TString save_path = get_save_path();

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(5);

  bool doTicks = true;

  TString o_tot = "E1";
  TString o_stat = "E1";

  if(!doTicks){
    o_tot = "E";
    o_stat = "E1 X0";
  }

  int titlefont = 43;
  double ymax = 0.04;

  TCanvas *c = new TCanvas("c","",600,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.31, 1, 1.0);
  pad1->SetBottomMargin(0.02); // 0.015
  pad1->SetLeftMargin(0.19);
  pad1->Draw();
  pad1->cd();
  TGaxis::SetMaxDigits(4);


  h_data_XS_tot->SetTitle(" ");
  h_data_XS_tot->GetYaxis()->SetRangeUser(0, ymax);
  h_data_XS_tot->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#it{m}_{jet}} #left[#frac{1}{GeV}#right]");
  h_data_XS_tot->GetYaxis()->SetTitleOffset(1.3);
  h_data_XS_tot->GetYaxis()->SetTitleSize(0.06);
  h_data_XS_tot->GetYaxis()->SetNdivisions(505);
  h_data_XS_tot->SetLineColor(kBlack);
  h_data_XS_tot->SetMarkerColor(kBlack);
  h_data_XS_tot->SetMarkerStyle(8);
  h_data_XS_tot->SetMarkerSize(1);
  h_data_XS_tot->Draw(o_tot);
  h_data_XS_stat->SetLineColor(kBlack);
  h_data_XS_stat->SetMarkerColor(kBlack);
  h_data_XS_stat->SetMarkerStyle(8);
  h_data_XS_stat->SetMarkerSize(1);
  Int_t color[] = {kRed-4, TColor::GetColor("#0059b3"), 14};
  Int_t fillcolor[] = {kRed-10, TColor::GetColor("#99ccff"), 16};
  Int_t style[] = {3, 1, 2};

  for(unsigned int i=0; i<truth.size(); i++){
    truth[i]->SetLineWidth(3);
    truth[i]->SetLineColor(color[i]);
    truth[i]->SetLineStyle(style[i]);
    truth[i]->SetMarkerStyle(0);
    truth[i]->SetFillColor(fillcolor[i]);
    truth[i]->Draw("E2 SAME");
    truth2[i]->SetLineWidth(3);
    truth2[i]->SetLineColor(color[i]);
    truth2[i]->SetLineStyle(style[i]);
    truth2[i]->Draw("HIST SAME");
  }
  truth[1]->Draw("E2 SAME");
  truth2[1]->Draw("HIST SAME");

  h_data_XS_stat->Draw(o_stat + " SAME");
  h_data_XS_tot->Draw(o_tot + " SAME");

  TLegend *l=new TLegend(0.55,0.62,0.82,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(h_data_XS_tot,"Data","ple");
  for(unsigned int i=0; i<truth.size(); i++){
    l->AddEntry(truth[i],legnames[i],"fl");
  }
  l->SetTextSize(0.052);
  l->Draw();
  gPad->RedrawAxis();

  TLatex *cmstext = new TLatex(3.5, 24, "CMS");
  cmstext->SetNDC();
  cmstext->SetTextAlign(13);
  cmstext->SetTextFont(62);
  cmstext->SetTextSize(0.08);
  cmstext->SetX(0.24);
  cmstext->SetY(0.84);
  cmstext->Draw();

  TString simtext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, simtext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.24);
  text3->SetTextFont(52);
  text3->SetTextSize(0.06);
  text3->SetY(0.77);
  text3->Draw();

  TString infotext = "137.1 fb^{-1} (13 TeV)";
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetX(0.9);
  text1->SetY(0.961);
  text1->SetTextSize(0.055);
  text1->Draw();

  // das hier wird nur gemacht, damit die Achsenbeschriftung da ist
  h_data_XS_tot->GetXaxis()->SetLabelSize(0.);
  h_data_XS_tot->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( 112, 0, 112, ymax, 0, ymax, 505,"");
  axis->SetLabelOffset(0.01);
  axis->SetLabelFont(titlefont);
  axis->SetLabelSize(21);
  axis->Draw();

  // Ratio Plot kommt in unteres pad
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.19);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.38);
  pad2->Draw();
  pad2->cd();

  h_ratio_data->GetXaxis()->SetTickLength(0.07);
  h_ratio_data->GetXaxis()->SetTitleSize(25);
  h_ratio_data->GetXaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetXaxis()->SetTitleOffset(4.0);
  h_ratio_data->GetXaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetXaxis()->SetLabelSize(21);
  h_ratio_data->GetXaxis()->SetLabelOffset(0.035);
  h_ratio_data->GetYaxis()->SetTitle("#frac{Theory}{Data}");
  h_ratio_data->GetYaxis()->CenterTitle();
  h_ratio_data->GetYaxis()->SetTitleSize(22);
  h_ratio_data->GetYaxis()->SetTitleFont(titlefont);
  h_ratio_data->GetYaxis()->SetTitleOffset(2.2);
  h_ratio_data->GetYaxis()->SetLabelFont(titlefont);
  h_ratio_data->GetYaxis()->SetLabelSize(19);
  h_ratio_data->GetYaxis()->SetLabelOffset(0.009);
  h_ratio_data->GetYaxis()->SetNdivisions(505);
  h_ratio_data->SetTitle(" ");
  h_ratio_data->GetYaxis()->SetRangeUser(0.2, 1.8);
  h_ratio_data->GetXaxis()->SetTitle("#it{m}_{jet} #left[GeV#right]");
  h_ratio_data->SetLineColor(kBlack);
  h_ratio_data->SetMarkerColor(kBlack);
  h_ratio_data->SetMarkerStyle(8);
  h_ratio_data->SetMarkerSize(1);
  h_ratio_data->Draw(o_tot);
  h_ratio_data_stat->SetLineColor(kBlack);
  h_ratio_data_stat->SetMarkerColor(kBlack);
  h_ratio_data_stat->SetMarkerStyle(8);
  h_ratio_data_stat->SetMarkerSize(1);

  TH1F* h_ratio_powheg_line = (TH1F*) ratios_mtop[1]->Clone();
  ratios_mtop[1]->SetMarkerSize(0);
  ratios_mtop[1]->SetFillColor(fillcolor[1]);
  ratios_mtop[1]->Draw("E2 SAME");
  ratios_mtop[0]->SetLineWidth(3);
  ratios_mtop[0]->SetLineColor(color[0]);
  ratios_mtop[0]->SetLineStyle(style[0]);
  ratios_mtop[0]->Draw("HIST SAME");
  ratios_mtop[2]->SetLineWidth(3);
  ratios_mtop[2]->SetLineColor(color[2]);
  ratios_mtop[2]->SetLineStyle(style[2]);
  ratios_mtop[2]->Draw("HIST SAME");
  h_ratio_data->Draw(o_tot+" SAME");

  gStyle->SetEndErrorSize(5);

  h_ratio_powheg_line->SetLineColor(color[1]);
  h_ratio_powheg_line->SetLineWidth(3);
  h_ratio_powheg_line->SetLineStyle(style[1]);


  // noch einmal alles in der richtigen Reihenfolge zeichnen
  ratios_mtop[1]->Draw("E2 SAME");
  h_ratio_powheg_line->Draw("HIST SAME");
  ratios_mtop[0]->Draw("HIST SAME");
  h_ratio_data_stat->Draw(o_stat + " SAME");
  h_ratio_data->Draw(o_tot + " SAME");
  gPad->RedrawAxis();

  c->SaveAs(save_path+"/Plots/Unfolding_Run2/Unfold_combination_norm_"+mode+".pdf");
  delete c;
  return;
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
  upper_ = data_->GetXaxis()->FindBin(upper);
  return;
}

void chi2fit::SetBins(std::vector<TH1F*> masspoints){
  // when calculating the chi2, one degree of freedom is lost by normalisation
  // therefore covariance matrix is only invertable if one excludes one bin from chi2
  // this bin is selected randomly (chi2 should not change depending on the bin)
  // for every masspoint, the bins that should be used are written into a vector
  bool skip;
  int rand_bin = rand()%(upper_-lower_+1)+lower_;
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
  c->SaveAs(save_path+"/Plots/Unfolding_Run2/Unfold_combination_chi2_"+mode+".pdf");
  delete c;
  return;
}
