#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

using namespace std;

TH1F* Normalize(TH1F* hist);
void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown);
void PlotFit(TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin);
void PlotChi2(TF1* chi2function, vector<double> FSRvalues);
vector<double> ExtractFSRValues(TF1* chi2function);
void PlotResults(vector<vector<double>> Values, vector<TString> masswindows);
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname);
void PlotError(TH1F* hist, TString sysname);

TString year;
TString save_path;

int main(int argc, char* argv[]){

  bool debug = false;
  if(argc != 2){
    cout << "Usage: ./FSRuncertainty <year>" << endl;
    cout << "  year = 2016/2017/2018/combine" << endl;
    return 0;
  }
  year = argv[1];
  SetupGlobalStyle();

  vector<TString> masswindows = {"140toInf", "140to160", "160to170", "170to180", "180to190", "190to200", "200toInf"};
  // vector<TString> masswindows = {"140to200"};
  vector<vector<double>> f_fsr;

  for(auto mass_string: masswindows){

    save_path = get_save_path()+"/Plots/FSRuncertainty/";
    save_path = creat_folder_and_path(save_path, year);
    save_path = creat_folder_and_path(save_path, mass_string);

    TFile* infile = new TFile("FSR_hists_mjet"+mass_string+".root");

    // #################################################################################################
    // Get Hists #######################################################################################
    TH1F *data = (TH1F*) infile->Get("data_"+year);
    TH1F *ttbar = (TH1F*) infile->Get("nominal_"+year);
    vector<TH1F*> bkg;
    bkg.push_back((TH1F*) infile->Get("bgr_st_"+year));
    bkg.push_back((TH1F*) infile->Get("bgr_wj_"+year));
    bkg.push_back((TH1F*) infile->Get("bgr_ot_"+year));
    vector<double> bkgsys = {0.23, 0.19, 1.0};

    TH1F *jecup = (TH1F*) infile->Get("jecup_"+year);
    TH1F *jecdown = (TH1F*) infile->Get("jecdown_"+year);
    TH1F *corup = (TH1F*) infile->Get("corup_"+year);
    TH1F *cordown = (TH1F*) infile->Get("cordown_"+year);
    TH1F *jmsup = (TH1F*) infile->Get("jmsup_"+year);
    TH1F *jmsdown = (TH1F*) infile->Get("jmsdown_"+year);

    vector<TH1F*> FSRup, FSRdown;
    vector<double> FSRvalues;
    if(year == "2016"){
      // FSRvalues = {1./2, 1., 2.};
      FSRvalues = {0.25, 1., 4.}; // squared
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2_"+year));
    }
    else{
      // FSRvalues = {1./4, 1./2, 1./sqrt(2), 1., sqrt(2), 2., 4.};
      FSRvalues = {0.0625, 0.25, 0.5, 1., 2., 4., 16.}; // squared
      FSRdown.push_back((TH1F*) infile->Get("fsrdown4_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdown2_"+year));
      FSRdown.push_back((TH1F*) infile->Get("fsrdownsqrt2_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrupsqrt2_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup2_"+year));
      FSRup.push_back((TH1F*) infile->Get("fsrup4_"+year));
    }

    // #################################################################################################
    // Subtract backgournd #############################################################################
    TH1F* data_sub = SubtractBackgrounds(data, bkg, bkgsys);

    // #################################################################################################
    // Rebin ###########################################################################################
    int rebin = 5;
    if(rebin > 1){
      data_sub->Rebin(rebin);
      ttbar->Rebin(rebin);
      jecup->Rebin(rebin);
      jecdown->Rebin(rebin);
      corup->Rebin(rebin);
      cordown->Rebin(rebin);
      jmsup->Rebin(rebin);
      jmsdown->Rebin(rebin);
      for(auto fsr: FSRup) fsr->Rebin(rebin);
      for(auto fsr: FSRdown) fsr->Rebin(rebin);
      cout << "Rebin with a factor of " << rebin << " - resulting in a bin width of " <<  data_sub->GetBinWidth(1) << endl;
    }

    // #################################################################################################
    // Normalize #######################################################################################
    TH1F *data_norm = Normalize(data_sub);
    TH1F *ttbar_norm = Normalize(ttbar);
    TH1F *jecup_norm = Normalize(jecup);
    TH1F *jecdown_norm = Normalize(jecdown);
    TH1F *corup_norm = Normalize(corup);
    TH1F *cordown_norm = Normalize(cordown);
    TH1F *jmsup_norm = Normalize(jmsup);
    TH1F *jmsdown_norm = Normalize(jmsdown);
    vector<TH1F*> FSRup_norm, FSRdown_norm;
    for(auto fsr: FSRup) FSRup_norm.push_back(Normalize(fsr));
    for(auto fsr: FSRdown) FSRdown_norm.push_back(Normalize(fsr));
    PlotTau32(data_norm, ttbar_norm, FSRup_norm, FSRdown_norm);
    PlotError(data_norm, "DATA");
    PlotError(ttbar_norm, "TTbar");
    // #################################################################################################
    // Get additional uncertainties ####################################################################
    vector<TH1F*> sys;
    sys.push_back(GetSYS(ttbar_norm, jecup_norm, jecdown_norm, "jec"));
    sys.push_back(GetSYS(ttbar_norm, corup_norm, cordown_norm, "cor"));
    sys.push_back(GetSYS(ttbar_norm, jmsup_norm, jmsdown_norm, "jms"));


    // #################################################################################################
    // Fit per bin #####################################################################################
    vector<vector<double>> fitparameters;
    vector<int> validbins;
    int nbins = ttbar_norm->GetXaxis()->GetNbins();
    TString fitformula = "[0] + [1]*log(x) + [2]*x";
    if(year == "2016") fitformula = "[0] + [1]*log(x)";
    int Npar = 1;
    if(fitformula.Contains("[1]")) Npar = 2;
    if(fitformula.Contains("[2]")) Npar = 3;
    if(fitformula.Contains("[3]")) Npar = 4;
    if(fitformula.Contains("[4]")) Npar = 5;
    if(fitformula.Contains("[5]")) Npar = 6;

    vector<TGraphErrors*> graphs, bands1, bands2;
    vector<TF1*> fits;
    for(int bin=1; bin<=nbins; bin++){
      // Ignore empty bins
      double mincontent = 0.0001;
      bool ignorebin = false;
      if(ttbar_norm->GetBinContent(bin) < mincontent) ignorebin = true;
      if(data_norm->GetBinContent(bin) < mincontent) ignorebin = true;
      for(auto fsr: FSRup_norm){
        if(fsr->GetBinContent(bin) < mincontent) ignorebin = true;
      }
      for(auto fsr: FSRdown_norm){
        if(fsr->GetBinContent(bin) < mincontent) ignorebin = true;
      }
      if(ignorebin){
        cout << "Too few entries, ignore bin " << bin << "!" << endl;
        continue;
      }

      validbins.push_back(bin);
      // Write bin contents and error in vectors in correct order
      vector<double> bincontents;
      vector<double> binerrors;
      vector<double> binerrors_large;
      vector<double> xerrors;

      for(auto fsr: FSRdown_norm){
        bincontents.push_back(fsr->GetBinContent(bin));
        binerrors.push_back(fsr->GetBinError(bin));
        binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
        xerrors.push_back(0.0);
      }
      bincontents.push_back(ttbar_norm->GetBinContent(bin));
      binerrors.push_back(ttbar_norm->GetBinError(bin));
      binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
      xerrors.push_back(0.0);
      for(auto fsr: FSRup_norm){
        bincontents.push_back(fsr->GetBinContent(bin));
        binerrors.push_back(fsr->GetBinError(bin));
        binerrors_large.push_back(FSRup_norm[FSRup_norm.size()-1]->GetBinError(bin)); // use always largest error
        xerrors.push_back(0.0);
      }

      // Create TGraph
      int npoints = FSRvalues.size();
      TGraphErrors* graph = new TGraphErrors(npoints, &FSRvalues[0], &bincontents[0], &xerrors[0], &binerrors[0]);
      TF1 *fit = new TF1("fit", fitformula);
      graph->Fit("fit", "Q");
      vector<double> params;
      for(int i=0; i<Npar; i++) params.push_back(fit->GetParameter(i));

      // Set Fit uncertainty
      TGraphErrors* band1 = (TGraphErrors*) graph->Clone();
      TFitter* fitter = (TFitter*) TVirtualFitter::GetFitter();
      fitter->GetConfidenceIntervals(band1, 0.68);
      TGraphErrors* band2 = new TGraphErrors(npoints, &FSRvalues[0], &bincontents[0], &xerrors[0], &binerrors_large[0]);


      fitparameters.push_back(params);
      graphs.push_back(graph);
      fits.push_back(fit);
      bands1.push_back(band1);
      bands2.push_back(band2);
    }
    for(int i=0; i<validbins.size();i++){
      PlotFit(graphs[i], fits[i], bands1[i], bands2[i], validbins[i], validbins[0], validbins[validbins.size()-1]);
    }

    // #################################################################################################
    // Now calculate the Chi2 ##########################################################################
    TString chi2formula;
    for(int i=0; i<validbins.size();i++){
      double dataentry = data_norm->GetBinContent(validbins[i]);
      double dataerror = data_norm->GetBinError(validbins[i]);
      // double fiterror = ttbar_norm->GetBinError(validbins[i]);
      double fiterror = FSRup_norm[FSRup_norm.size()-1]->GetBinError(validbins[i]);
      double error2 = pow(dataerror,2) + pow(fiterror,2);
      if(debug){
        cout << "----------------------------" << endl;
        cout << "Bin " << validbins[i] << endl;
        cout << "Data Error " << dataerror << endl;
        cout << "Fit  Error " << fiterror << endl;
      }

      for(auto s: sys) error2 += pow(s->GetBinContent(validbins[i]), 2);

      double error = sqrt(error2);
      TString modifiedformula = fitformula;
      for(int par=0; par<Npar; par++){
        TString paramstring = "["+to_string(par)+"]";
        modifiedformula.ReplaceAll(paramstring, to_string(fitparameters[i][par]));
      }

      TString formula = "(("+to_string(dataentry)+"-(" + modifiedformula +"))/"+to_string(error)+")^2";
      if(debug) cout << formula << endl;
      if(i==0) chi2formula = formula;
      else{
        chi2formula += "+";
        chi2formula += formula;
      }
    }
    double min = 0.01;
    double max = 100.;
    TF1* chi2function = new TF1("chi2", chi2formula, min, max);
    vector<double> Values = ExtractFSRValues(chi2function);
    f_fsr.push_back(Values);
    PlotChi2(chi2function, Values);

    cout << "--------------------------------------------------------" << endl;
    cout << "f(FSR) = " << Values[0] << " + " << Values[1] << " - " << Values[2] << ", chi2min = " << Values[5] << endl;
    cout << "--------------------------------------------------------" << endl;


    fstream fsr_txt;
    fsr_txt.open(save_path+"/fsr_factor.txt", ios::out);
    fsr_txt << "f_FSR  = " << Values[0] << " + " << Values[1] << " - " << Values[2] << endl;
    fsr_txt << "f_up   = " << Values[3] << endl;
    fsr_txt << "f_down = " << Values[4] << endl;
    fsr_txt.close();

  }

  PlotResults(f_fsr, masswindows);
  return 0;

}

// #################################################################################################
// #################################################################################################
// #################################################################################################
TH1F* Normalize(TH1F* hist){
  TH1F* norm = (TH1F*) hist->Clone();
  int nbins = hist->GetXaxis()->GetNbins();
  double integral = hist->Integral();
  norm->Scale(1/integral);
  // Now do error propagation:
  // Bin content of norm hist: b_i = N_i / (N_1 + N_2 + N_3 + ...)
  // Two derivations:
  // 1. db_i / dN_i = (N_2 + N_3 + ...) / (N_1 + N_2 + N_3 + ...)^2 = (Integral - N_i) / Integral^2
  // 2. db_i / dN_j = - N_1 / (N_1 + N_2 + N_3 + ...)^2

  for(int i=1; i<=nbins; i++){
    double error2 = 0;
    for(int j=1; j<=nbins; j++){
      double additionalterm;
      if(i==j) additionalterm = (integral - hist->GetBinContent(j)) / (integral * integral);
      else     additionalterm = - hist->GetBinContent(j) / (integral * integral);
      error2 += additionalterm*additionalterm * hist->GetBinError(j)*hist->GetBinError(j);
    }
    double error = sqrt(error2);
    norm->SetBinError(i, error);
  }
  return norm;
}


void PlotFit(TGraphErrors* graph, TF1* fit, TGraphErrors* band1, TGraphErrors* band2, int bin, int firstbin, int lastbin){
  TString binnr = to_string(bin);
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graph->SetMarkerStyle(8);
  graph->SetTitle("Bin "+binnr);
  graph->GetXaxis()->SetTitle("(f^{FSR})^{2}");
  graph->GetYaxis()->SetTitle("a.u.");
  graph->Draw("AP");
  //fit->SetLineColor(kBlack);

  band2->SetLineColorAlpha(kOrange-4,0.5);
  band2->SetFillColorAlpha(kOrange-4,0.5);
  band1->SetLineColorAlpha(kAzure-4,0.8);
  band1->SetFillColorAlpha(kAzure-4,0.8);
  band2->Draw("L3 SAME");
  // band1->Draw("L3 SAME");

  graph->Draw("P SAME");
  fit->Draw("SAME");

  c->SaveAs(save_path+"/Bin"+binnr+".pdf");
  if(bin == firstbin)     c->Print(save_path+"/AllBins.pdf(","pdf");
  else if(bin == lastbin) c->Print(save_path+"/AllBins.pdf)","pdf");
  else                    c->Print(save_path+"/AllBins.pdf","pdf");
  delete c;
}

void PlotChi2(TF1* chi2function, vector<double> FSRvalues){
  TCanvas* c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  // gPad->SetLogx();
  chi2function->SetTitle(" ");
  chi2function->GetXaxis()->SetRangeUser(0.1, 25);
  chi2function->GetYaxis()->SetRangeUser(0, 100);
  chi2function->GetXaxis()->SetTitle("(f^{FSR})^{2}");
  chi2function->GetYaxis()->SetTitle("#chi^{2}");
  chi2function->GetYaxis()->SetTitleOffset(1.2);
  chi2function->Draw();
  c->SaveAs(save_path+"/chi2.pdf");
  delete c;
}

vector<double> ExtractFSRValues(TF1* chi2function){
  double minChi2      = chi2function->GetMinimum();
  double f_FSR2       = chi2function->GetX(minChi2, 0.001, 1000);
  double f_FSR        = sqrt(f_FSR2);
  double f_up2        = chi2function->GetX(minChi2+1, f_FSR2, 1000);
  double f_down2      = chi2function->GetX(minChi2+1, 0.001, f_FSR2);
  double f_up          = sqrt(f_up2);
  double f_down        = sqrt(f_down2);

  double sigmaup      = f_up - f_FSR;
  double sigmadown    = f_FSR - f_down;
  vector<double> values = {f_FSR, sigmaup, sigmadown, f_up, f_down, minChi2};
  return values;
}

void PlotResults(vector<vector<double>> f_fsr, vector<TString> masswindows){
  vector<double> xvalues, xup, xdown, central, up, down;
  for(unsigned int i=0; i<f_fsr.size(); i++){
    central.push_back(f_fsr[i][0]);
    up.push_back(f_fsr[i][1]);
    down.push_back(f_fsr[i][2]);
    xvalues.push_back(i+1);
    xup.push_back(0.);
    xdown.push_back(0.);
  }

  double ymin = 0.0;
  double ymax = 5.0;

  TH1F* dummy = new TH1F("dummy", "dummy", f_fsr.size(), 0.5, 0.5+f_fsr.size());
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  for(unsigned int i=0; i<masswindows.size(); i++) dummy->GetXaxis()->SetBinLabel(i+1, masswindows[i]);
  dummy->SetTitle("");
  dummy->GetYaxis()->SetTitle("#it{f}^{FSR}");
  dummy->GetXaxis()->SetTitle("#it{m}_{jet}");
  dummy->GetYaxis()->SetRangeUser(ymin, ymax);
  dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->Draw();
  TGraphAsymmErrors * results = new TGraphAsymmErrors(f_fsr.size(), &xvalues[0], &central[0], &xdown[0], &xup[0], &down[0], &up[0]);
  results->SetMarkerStyle(8);
  results->Draw("P");
  TLine *errbandup   = new TLine(0.5, central[0]+up[0], 0.5+f_fsr.size(), central[0]+up[0]);
  TLine *errbanddown = new TLine(0.5, central[0]-down[0], 0.5+f_fsr.size(), central[0]-down[0]);
  errbandup->SetLineColor(13);
  errbanddown->SetLineColor(13);
  errbandup->SetLineStyle(2);
  errbanddown->SetLineStyle(2);
  errbandup->Draw("SAME");
  errbanddown->Draw("SAME");
  results->Draw("P SAME");

  double x = .7;
  for(int i=0; i<f_fsr.size();i++){
    // round min chi2 and store in TString
    std::ostringstream out;
    out.precision(2);
    out << std::fixed << f_fsr[i][5];
    TString label = "#chi_{min}^{2} = " + out.str();
    ////
    TLatex latex;
    latex.SetTextFont(62);
    latex.SetTextSize(0.02);
    latex.DrawLatex(x, .5, label);
    x += 1.0;
  }

  c->SaveAs(save_path+"/Plots/FSRuncertainty/Results_"+year+".pdf");
  delete c;
}

void PlotTau32(TH1F* data, TH1F* ttbar, vector<TH1F*> FSRup, vector<TH1F*> FSRdown){
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(.2);
  data->SetMarkerStyle(8);
  data->SetMarkerColor(1);
  data->SetLineColor(1);
  ttbar->SetLineWidth(2);
  ttbar->SetLineColor(14);
  ttbar->SetTitle(" ");
  ttbar->GetXaxis()->SetTitle("#tau_{32}");
  ttbar->GetYaxis()->SetTitle("a.u.");
  ttbar->GetYaxis()->SetRangeUser(0., 1.4*ttbar->GetMaximum());

  ttbar->Draw("HIST");
  Int_t color[] = {kAzure+7, kGreen-4, 810};
  int j = FSRdown.size()-1;
  for(int i=0; i<FSRup.size(); i++){
    FSRup[i]->SetLineWidth(2);
    FSRdown[j]->SetLineWidth(2);
    FSRup[i]->SetLineStyle(1);
    FSRdown[j]->SetLineStyle(2);
    FSRup[i]->SetLineColor(color[i]);
    FSRdown[j]->SetLineColor(color[i]);
    if(FSRup.size() == 1){
      FSRup[i]->Draw("HIST SAME");
      FSRdown[j]->Draw("HIST SAME");
    }
    if(i == 2){
      FSRup[i]->Draw("HIST SAME");
      FSRdown[j]->Draw("HIST SAME");
    }
    j--;
  }
  data->Draw("EP SAME");
  TLegend * leg = new TLegend(.25, .65, .50, .85);
  leg->AddEntry(data, "Data - background", "pe");
  leg->AddEntry(ttbar, "t#bar{t} nominal", "l");
  vector<TString> legnamesdown = {"t#bar{t} f^{FSR} = #frac{1}{4}", "t#bar{t} f^{FSR} = #frac{1}{2}", "t#bar{t} f^{FSR} = #frac{1}{#sqrt{2}}"};
  vector<TString> legnamesup = {"t#bar{t} f^{FSR} = #sqrt{2}", "t#bar{t} f^{FSR} = 2", "t#bar{t} f^{FSR} = 4"};
  if( year == "2016" ){
    legnamesdown = {"t#bar{t} f^{FSR} = #frac{1}{2}"};
    legnamesup = {"t#bar{t} f^{FSR} = 2"};
  }
  for(int i=0; i<FSRdown.size(); i++){
    if(i == 0) leg->AddEntry(FSRdown[i], legnamesdown[i], "l");
  }
  for(int i=0; i<FSRup.size(); i++){
    if(i == FSRup.size()-1) leg->AddEntry(FSRup[i], legnamesup[i], "l");
  }
  leg->Draw();
  gPad->RedrawAxis();
  c->SaveAs(save_path+"/Tau32.pdf");
  delete c;
}

TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize){
  int nbins = data->GetXaxis()->GetNbins();
  TH1F* result = (TH1F*) data->Clone();
  for(unsigned int i=0; i<bgr.size(); i++){
    result->Add(bgr[i], -1);
  }
  for(int bin=1; bin<=nbins; bin++){
    double syserror2 = 0;
    for(unsigned int i=0; i<bgr.size(); i++){
      syserror2 += pow(bgr[i]->GetBinContent(bin) * syssize[i], 2);
    }
    double olderror2 = pow(result->GetBinError(bin), 2);
    result->SetBinError(bin, sqrt(syserror2+olderror2));
  }
  return result;
}

TH1F* GetSYS(TH1F* hist, TH1F* up, TH1F* down, TString sysname){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double diff1 = fabs( hist->GetBinContent(bin) - up->GetBinContent(bin) );
    double diff2 = fabs( hist->GetBinContent(bin) - down->GetBinContent(bin) );
    if(diff1 > diff2) sys->SetBinContent(bin, diff1);
    else              sys->SetBinContent(bin, diff2);
  }
  TCanvas *c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("#tau_{32}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(0., 0.005);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST");
  c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  delete c;
  return sys;
}


void PlotError(TH1F* hist, TString sysname){
  int nbins = hist->GetXaxis()->GetNbins();
  TH1F* sys = (TH1F*) hist->Clone();
  sys->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double error = hist->GetBinError(bin);
    sys->SetBinContent(bin, error);
  }
  TCanvas *c = new TCanvas("", "", 600, 600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  sys->SetTitle(sysname);
  sys->GetXaxis()->SetTitle("#tau_{32}");
  sys->GetYaxis()->SetTitle("Uncertainty");
  sys->GetYaxis()->SetRangeUser(0., 0.005);
  sys->SetLineColor(kRed);
  sys->SetFillColor(kRed);
  sys->Draw("HIST");
  c->SaveAs(save_path+"/SYS_"+sysname+".pdf");
  delete c;
}
