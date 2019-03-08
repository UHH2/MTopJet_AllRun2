#include "../include/mass.h"

int main(int argc, char* argv[])
{

  if(argc != 2){
    cout << "ERROR: WRONG USAGE!" << endl;
    cout << endl;
    cout << "CORRECT WAY: ./mass option" << endl;
    cout << endl;
    cout << "           option = data/pseudo1" << endl;
    return 0;
  }

  bool use_data = false;
  bool use_pseudo1 = false;
  bool use_pseudo1695 = false;
  bool use_pseudo1715 = false;
  bool use_pseudo1735 = false;
  bool use_pseudo1755 = false;

  if(strcmp(argv[1], "data") == 0)            use_data = true;
  else if(strcmp(argv[1], "pseudo1") == 0)    use_pseudo1 = true;
  else if(strcmp(argv[1], "pseudo1695") == 0) use_pseudo1695 = true;
  else if(strcmp(argv[1], "pseudo1715") == 0) use_pseudo1715 = true;
  else if(strcmp(argv[1], "pseudo1735") == 0) use_pseudo1735 = true;
  else if(strcmp(argv[1], "pseudo1755") == 0) use_pseudo1755 = true;
  else{
    cout << "[ERROR] select a valid option!" << endl;
    return 0;
  }

  TString dir = argv[1];

  /*
  ██████  ███████ ████████     ███████ ██ ██      ███████ ███████
  ██      ██         ██        ██      ██ ██      ██      ██
  ██  ███ █████      ██        █████   ██ ██      █████   ███████
  ██   ██ ██         ██        ██      ██ ██      ██           ██
  ██████  ███████    ██        ██      ██ ███████ ███████ ███████
  */

  TFile* file = new TFile("Histograms.root");

  vector<TFile*> SYS_FILES;
  // SYS_FILES.push_back(new TFile("SYS_DUMMY.root"));
  SYS_FILES.push_back(new TFile("SYS_BTAG.root"));
  SYS_FILES.push_back(new TFile("SYS_BKG.root"));
  SYS_FILES.push_back(new TFile("SYS_PU.root"));
  SYS_FILES.push_back(new TFile("SYS_MUID.root"));
  SYS_FILES.push_back(new TFile("SYS_MUTR.root"));
  SYS_FILES.push_back(new TFile("SYS_COR.root"));
  SYS_FILES.push_back(new TFile("SYS_JEC.root"));
  SYS_FILES.push_back(new TFile("SYS_JER.root"));
  SYS_FILES.push_back(new TFile("SYS_SCALE.root"));
  SYS_FILES.push_back(new TFile("SYS_GENERATOR.root"));
  SYS_FILES.push_back(new TFile("SYS_SHOWER.root"));
  // SYS_FILES.push_back(new TFile("SYS_PDF.root"));


  /*
  ██████   ███████ ████████    ██   ██ ██ ███████ ████████ ███████
  ██       ██         ██       ██   ██ ██ ██         ██    ██
  ██   ███ █████      ██       ███████ ██ ███████    ██    ███████
  ██    ██ ██         ██       ██   ██ ██      ██    ██         ██
  ██████   ███████    ██       ██   ██ ██ ███████    ██    ███████
  */

  vector<TH1F*> UP_Variations;
  vector<TH1F*> DOWN_Variations;
  vector<TH2D*> CovMatrix;

  for(auto f:SYS_FILES){
    UP_Variations.push_back((TH1F*)f->Get("EnvelopeUp/hist"));
    DOWN_Variations.push_back((TH1F*)f->Get("EnvelopeDown/hist"));
    CovMatrix.push_back((TH2D*)f->Get("CovMatrix/cov"));
  }


  // set data distribution
  TString dataname;
  if(use_data)            dataname = "data";
  else if(use_pseudo1)    dataname = "pseudo1";
  else if(use_pseudo1695) dataname = "1695";
  else if(use_pseudo1715) dataname = "1715";
  else if(use_pseudo1735) dataname = "1735";
  else if(use_pseudo1755) dataname = "1755";
  TH1F * JetMass_data = (TH1F*)file->Get(dataname);

  // set central tt distribution
  TH1F * Central = (TH1F*)file->Get("tt");

  // subtract backgrounds from data
  vector<TString> bkg_name = {"WJets", "SingleTop", "other"};
  for(auto bkg: bkg_name){
    TH1F * JetMass_bgr = (TH1F*)file->Get(bkg);
    JetMass_data->Add(JetMass_bgr, -1);
  }

  TH1F * JetMass_1665 = (TH1F*)file->Get("1665");
  TH1F * JetMass_1695 = (TH1F*)file->Get("1695");
  TH1F * JetMass_1715 = (TH1F*)file->Get("1715");
  TH1F * JetMass_1725 = (TH1F*)file->Get("tt");
  TH1F * JetMass_1735 = (TH1F*)file->Get("1735");
  TH1F * JetMass_1755 = (TH1F*)file->Get("1755");
  TH1F * JetMass_1785 = (TH1F*)file->Get("1785");

  std::vector<TH1F*> mass_samples_temp = {JetMass_1665, JetMass_1695, JetMass_1715, JetMass_1725, JetMass_1735, JetMass_1755, JetMass_1785};
  std::vector<Double_t> masses_temp = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
  std::vector<TH1F*> mass_samples;
  std::vector<Double_t> masses;
  // if mass is extracted from mass sample, it has to be removed from mass templates
  for(int i=0; i<masses_temp.size(); i++){
    if(use_data || use_pseudo1){
      if(i == 0 || i==6) continue;
    }
    else if(use_pseudo1695){
      if(i == 1 || i==6) continue;
    }
    else if(use_pseudo1715){
      if(i == 2 || i==6) continue;
    }
    else if(use_pseudo1735){
      if(i == 4 || i==0) continue;
    }
    else if(use_pseudo1755){
      if(i == 5 || i==0) continue;
    }
    mass_samples.push_back(mass_samples_temp[i]);
    masses.push_back(masses_temp[i]);
  }
  int N_samples = mass_samples.size();

  /*
  ██████   █████  ███    ██  ██████  ███████      ██████  ███████     ███████ ██ ████████
  ██   ██ ██   ██ ████   ██ ██       ██          ██    ██ ██          ██      ██    ██
  ██████  ███████ ██ ██  ██ ██   ███ █████       ██    ██ █████       █████   ██    ██
  ██   ██ ██   ██ ██  ██ ██ ██    ██ ██          ██    ██ ██          ██      ██    ██
  ██   ██ ██   ██ ██   ████  ██████  ███████      ██████  ██          ██      ██    ██
  */


  double from = 82;
  double to = 322;

  int firstbin = JetMass_data->GetXaxis()->FindBin(from);
  int lastbin = JetMass_data->GetXaxis()->FindBin(to) - 1;


  /*
  ██   ██  █████  ███    ██ ██████  ██      ███████     ██    ██ ███    ██  ██████ ███████ ██████  ████████  █████  ██ ███    ██ ████████ ██ ███████ ███████
  ██   ██ ██   ██ ████   ██ ██   ██ ██      ██          ██    ██ ████   ██ ██      ██      ██   ██    ██    ██   ██ ██ ████   ██    ██    ██ ██      ██
  ███████ ███████ ██ ██  ██ ██   ██ ██      █████       ██    ██ ██ ██  ██ ██      █████   ██████     ██    ███████ ██ ██ ██  ██    ██    ██ █████   ███████
  ██   ██ ██   ██ ██  ██ ██ ██   ██ ██      ██          ██    ██ ██  ██ ██ ██      ██      ██   ██    ██    ██   ██ ██ ██  ██ ██    ██    ██ ██           ██
  ██   ██ ██   ██ ██   ████ ██████  ███████ ███████      ██████  ██   ████  ██████ ███████ ██   ██    ██    ██   ██ ██ ██   ████    ██    ██ ███████ ███████
  */

  // calculate covariance matrix with normalisation applied (for data and templates!)
  // careful: do not use scaled distribution here!
  TH2D*cov = GetCovFromStat(JetMass_data, firstbin, lastbin, true);
  vector<TH2D*> cov_mc;
  for(unsigned int i=0; i<mass_samples.size(); i++){
    cov_mc.push_back(GetCovFromStat(mass_samples[i], firstbin, lastbin, true));
  }

  // now get covariance matrices for systematic uncertainties
  // these can then be added to the cov with stat uncertainties
  // no more error propagation is needed!
  vector<TH2D*> SYS_COV;
  vector<TH1F*> NormErrorUp;
  vector<TH1F*> NormErrorDown;
  vector<TH1F*> NormErrorSym;
  TString CovNames[] = {"BTAG", "PU", "MUID", "MUTR", "SCALE", "BKG", "COR", "JEC", "JER", "Generator", "Shower"};
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
    SYS_COV.push_back(NormaliseCov(CovMatrix[i], Central, firstbin, lastbin, true));
  }

  TH2D* STAT_COV = (TH2D*) cov->Clone("STAT_COV");
  if(use_data) for(auto cov_i:SYS_COV) cov->Add(cov_i);

  /*
  ███    ██  ██████  ██████  ███    ███  █████  ██      ██ ███████ ███████
  ████   ██ ██    ██ ██   ██ ████  ████ ██   ██ ██      ██ ██      ██
  ██ ██  ██ ██    ██ ██████  ██ ████ ██ ███████ ██      ██ ███████ █████
  ██  ██ ██ ██    ██ ██   ██ ██  ██  ██ ██   ██ ██      ██      ██ ██
  ██   ████  ██████  ██   ██ ██      ██ ██   ██ ███████ ██ ███████ ███████
  */
  TH1F* data_ = (TH1F*)JetMass_data->Clone("data");

  // normalise here
  for(auto &hist: mass_samples){
    hist->Scale(1/hist->Integral(firstbin, lastbin), "width");
  }
  data_->Scale(1/data_->Integral(firstbin, lastbin),"width");
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

  TGraph* chi_hist = new TGraph(N_samples,&masses[0],y);
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

  std::ofstream out("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/Mass.txt");
  auto coutbuf = std::cout.rdbuf(out.rdbuf());
  cout << minX << " " << (minX_sigmaup-minX) << " ";
  std::cout.rdbuf(coutbuf);

  /*
  ██████  ██       ██████  ████████ ████████ ██ ███    ██  ██████
  ██   ██ ██      ██    ██    ██       ██    ██ ████   ██ ██
  ██████  ██      ██    ██    ██       ██    ██ ██ ██  ██ ██   ███
  ██      ██      ██    ██    ██       ██    ██ ██  ██ ██ ██    ██
  ██      ███████  ██████     ██       ██    ██ ██   ████  ██████
  */

  double xlow, xhigh;
  double ylow, yhigh;
  xlow = from;
  xhigh = to;
  // xlow = 0;
  // xhigh = 1000;
  ylow = 0;
  yhigh = data->GetMaximum() * 1.1;
  if(JetMass_1695->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1695->GetMaximum() * 1.1;
  if(JetMass_1755->GetMaximum() * 1.1 > yhigh)  yhigh = JetMass_1755->GetMaximum() * 1.1;

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

  JetMass_1695->SetLineColor(kRed-4);
  JetMass_1715->SetLineColor(kRed-4);

  JetMass_1725->SetLineColor(13);

  JetMass_1735->SetLineColor(kAzure+7);
  JetMass_1755->SetLineColor(kAzure+7);

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

  JetMass_1695->Draw("HIST SAME");
  JetMass_1725->Draw("HIST SAME");
  JetMass_1755->Draw("HIST SAME");

  data->Draw("E1 SAME");
  dataSTAT->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.3,0.15,0.7,0.35);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(data,"data","pl");
  leg->AddEntry(JetMass_1695,"t#bar{t} (m_{top} = 169.5 GeV)","l");
  leg->AddEntry(JetMass_1725,"t#bar{t} (m_{top} = 172.5 GeV)","l");
  leg->AddEntry(JetMass_1755,"t#bar{t} (m_{top} = 175.5 GeV)","l");
  leg->Draw("");
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/JetMass.pdf");

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
  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextFont(43);
  text.SetTextSize(18);
  char mass_text[32];
  sprintf(mass_text, "%.5g", minX);
  char uncert_text[32];
  if((minX_sigmaup-minX) < 1) sprintf(uncert_text, "%.3g", (minX_sigmaup-minX));
  else           sprintf(uncert_text, "%.4g", (minX_sigmaup-minX));
  TString masstext = "m_{top} = ";
  masstext += mass_text;
  masstext += " #pm ";
  masstext += uncert_text;
  text.DrawLatex(.4,.6, masstext);
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/chi2.pdf");

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
    C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/COV/Error_"+CovNames[i]+".pdf");
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
    C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/COV/"+CovNames[i]+".pdf");
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
  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/MassExtractionReco/"+dir+"/COV/STAT.pdf");
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
