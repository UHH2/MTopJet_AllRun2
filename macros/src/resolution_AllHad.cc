#include "../include/CentralInclude.h"

using namespace std;

// compile with:
// g++ -o resolution_subjets_flavorJEC resolution_subjets_flavorJEC.cc `root-config --cflags --glibs`

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown);
vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error);

int main(int argc, char* argv[]){
  TString save_path = get_save_path();
  bool use_median = false;
  if(argc >1 && strcmp(argv[1], "median") == 0) use_median = true;

  //= declare files ============================================================
  bool year_16 = false;
  bool year_17 = false;
  bool year_18 = false;
  TString filename, year;
  TFile *f_tt;
  if(argc > 1 && strcmp(argv[1], "2016") == 0) {
    year = "2016";
    f_tt = new TFile(dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2016v3.root");
  }
  if(argc > 1 && strcmp(argv[1], "2017") == 0){
    year = "2017";
    f_tt = new TFile(dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2017v2.root");
  }
  if(argc > 1 && strcmp(argv[1], "2018") == 0){
    year = "2018";
    f_tt = new TFile(dir+"AllHad/uhh2.AnalysisModuleRunner.MC.TTbar_allHad_2018.root");
  }

  string year_s = (string) year;

  // set binning
  int n_ptbin = 10;
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TString bin_strings[] = {"0", "50", "80", "120", "170", "220", "270", "320", "370", "420", "600"};

  // this is for naming of pdf files
  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  // now get plots from file
  vector<TH1F*> reso, reso_noJEC, reso_cor;
  vector<TH1F*> reso_jet1, reso_jet2, reso_jet3;
  vector<TH1F*> ptrec;

  std::string dir_jec = "RecGenHists_allHad/";
  std::string dir_jec1 = "RecGenHists_allHad_jet1/";
  std::string dir_jec2 = "RecGenHists_allHad_jet2/";
  std::string dir_jec3 = "RecGenHists_allHad_jet3/";
  std::string dir_raw = "RecGenHists_allHad_noJEC/";
  std::string dir_cor = "RecGenHists_allHad_corrected/";

  std::string name_jec;
  std::string name_jec1;
  std::string name_jec2;
  std::string name_jec3;
  std::string name_raw;
  std::string name_cor;
  std::string name_ptrec;


  for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){

    name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
    name_jec1 = dir_jec1 + "PtResolution_" + std::to_string(ptbin);
    name_jec2 = dir_jec2 + "PtResolution_" + std::to_string(ptbin);
    name_jec3 = dir_jec3 + "PtResolution_" + std::to_string(ptbin);
    name_raw = dir_raw + "PtResolution_" + std::to_string(ptbin);
    name_cor = dir_cor + "PtResolution_" + std::to_string(ptbin);
    name_ptrec = dir_jec + "PtRec_" + std::to_string(ptbin);

    reso.push_back( (TH1F*)f_tt->Get(name_jec.c_str()) );
    reso_jet1.push_back( (TH1F*)f_tt->Get(name_jec1.c_str()) );
    reso_jet2.push_back( (TH1F*)f_tt->Get(name_jec2.c_str()) );
    reso_jet3.push_back( (TH1F*)f_tt->Get(name_jec3.c_str()) );
    reso_noJEC.push_back( (TH1F*)f_tt->Get(name_raw.c_str()) );
    reso_cor.push_back( (TH1F*)f_tt->Get(name_cor.c_str()) );
    ptrec.push_back( (TH1F*)f_tt->Get(name_ptrec.c_str()) );
  }

  // fit histograms and extract mean/median value
  TH1F * resolution  = GetResoPlot(reso, use_median);
  TH1F * resolution_jet1  = GetResoPlot(reso_jet1, use_median);
  TH1F * resolution_jet2  = GetResoPlot(reso_jet2, use_median);
  TH1F * resolution_jet3  = GetResoPlot(reso_jet3, use_median);
  TH1F * resolution_noJEC = GetResoPlot(reso_noJEC, use_median);
  TH1F * resolution_cor = GetResoPlot(reso_cor, use_median);
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // calculate non-closure as function of ptrec
  // arguments: (hists, use median?, do ptrec?, mean or error)
  vector<double> pt = GetMeans(ptrec, use_median, true, "mean");
  vector<double> pt_err = GetMeans(ptrec, use_median, true, "error");
  vector<double> MeanForUncert = GetMeans(reso_cor, use_median, false, "mean");
  vector<double> MeanForUncert_err = GetMeans(reso_cor, use_median, false, "error");
  ////
  TGraphErrors* non_closure = new TGraphErrors(n_ptbin, &pt[0], &MeanForUncert[0], &pt_err[0], &MeanForUncert_err[0]);
  // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
  // do same to reach ptrec = 600
  non_closure->SetPoint(non_closure->GetN(), 0, MeanForUncert[0]);
  int n_last = MeanForUncert.size()-1;
  non_closure->SetPoint(non_closure->GetN(), 600, MeanForUncert[n_last]);

  TGraph* area = AreaFromEnvelope(non_closure, "area");
  TGraph* upvar = AreaFromEnvelope(non_closure, "up");
  TGraph* downvar = AreaFromEnvelope(non_closure, "down");

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // non-closure in root file
  TString rootname = "../CorrectionFile/Correction_SysFromResolution_"+year+".root";
  TFile * outfile = new TFile(rootname,"RECREATE");;
  upvar->Write("Up");
  downvar->Write("Down");
  outfile->Close();

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // start plotting

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);


  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  resolution->GetYaxis()->SetRangeUser(-0.15, 0.15);
  resolution->GetXaxis()->SetNdivisions(505);
  resolution->GetYaxis()->SetNdivisions(505);
  resolution->GetXaxis()->SetTitleSize(0.06);
  resolution->GetYaxis()->SetTitleSize(0.04);
  resolution->GetXaxis()->SetTitleOffset(0.9);
  resolution->GetYaxis()->SetTitleOffset(1.7);
  resolution->GetXaxis()->SetLabelSize(0.05);
  resolution->GetYaxis()->SetLabelSize(0.05);
  resolution->GetXaxis()->SetTitle("p_{T,}^{gen} [GeV]");
  resolution->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution->SetLineWidth(4);
  resolution->SetLineColor(kAzure+7);
  resolution_noJEC->SetLineWidth(4);
  resolution_noJEC->SetLineColor(kOrange+1);
  resolution_cor->SetLineWidth(4);
  resolution_cor->SetLineColor(kRed+1);


  resolution_jet1->GetYaxis()->SetRangeUser(-0.15, 0.15);
  resolution_jet1->GetXaxis()->SetNdivisions(505);
  resolution_jet1->GetYaxis()->SetNdivisions(505);
  resolution_jet1->GetXaxis()->SetTitleSize(0.06);
  resolution_jet1->GetYaxis()->SetTitleSize(0.04);
  resolution_jet1->GetXaxis()->SetTitleOffset(0.9);
  resolution_jet1->GetYaxis()->SetTitleOffset(1.7);
  resolution_jet1->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
  resolution_jet1->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
  resolution_jet1->SetLineWidth(4);
  resolution_jet2->SetLineWidth(4);
  resolution_jet3->SetLineWidth(4);

  resolution_jet1->SetLineColor(kOrange+1);
  resolution_jet2->SetLineColor(kAzure+7);
  resolution_jet3->SetLineColor(kRed+1);

  non_closure->GetXaxis()->SetRangeUser(0, 600);
  non_closure->GetYaxis()->SetRangeUser(-10, 10);
  non_closure->GetXaxis()->SetNdivisions(505);
  non_closure->GetYaxis()->SetNdivisions(505);
  non_closure->GetXaxis()->SetTitleSize(0.05);
  non_closure->GetYaxis()->SetTitleSize(0.04);
  non_closure->GetXaxis()->SetTitleOffset(0.9);
  non_closure->GetYaxis()->SetTitleOffset(1.5);
  non_closure->GetXaxis()->SetTitle("p_{T}^{rec}");
  non_closure->GetYaxis()->SetTitle("non-closure uncertainty [%]");
  non_closure->SetTitle(" ");
  non_closure->SetMarkerStyle(20);
  non_closure->SetMarkerSize(0.8);


  TCanvas *c1 = new TCanvas();
  gPad->SetTopMargin(0.02);
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.13);
  resolution->Draw("E1");
  zero_line->Draw("SAME");
  resolution->Draw("SAME E1");
  resolution_noJEC->Draw("SAME E1");
  TLegend *leg1 = new TLegend(0.5,0.7,0.80,0.95);
  leg1->AddEntry(resolution_noJEC,"no correction","le");
  leg1->AddEntry(resolution,"AK4 correction","le");
  leg1->SetTextSize(0.05);
  leg1->Draw();
  CMSSimLabel(true, 0.22, 0.93);
  TLatex text1;
  text1.SetNDC(kTRUE);
  text1.SetTextFont(43);
  text1.SetTextSize(18);
  text1.DrawLatex(.22,.2, "all hadronic t#bar{t}");
  gPad->RedrawAxis();
  // first save without additional correction
  c1->SaveAs(save_path+"/Resolution_Subjets/muon/AllHadronic/"+year+"/pt_"+mean_median+"_noAdditional_"+year+".pdf");
  // and once again with the additional correction
  resolution_cor->Draw("SAME E1");
  leg1->AddEntry(resolution_cor,"AK4 + additional correction","le");
  // leg1->AddEntry((TObject*)0, "all hadronic t#bar{t}", "");
  c1->SaveAs(save_path+"/Resolution_Subjets/muon/AllHadronic/"+year+"/pt_"+mean_median+"_"+year+".pdf");


  TCanvas *c2 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  resolution_jet1->Draw("E1");
  zero_line->Draw("SAME");
  resolution_jet1->Draw("SAME E1");
  resolution_jet2->Draw("SAME E1");
  resolution_jet3->Draw("SAME E1");
  TLegend *leg2 = new TLegend(0.45,0.85,0.80,0.65);
  leg2->AddEntry(resolution_jet1,"only leading subjet","l");
  leg2->AddEntry(resolution_jet2,"only second subjet","l");
  leg2->AddEntry(resolution_jet3,"only third subjet","l");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  TLatex text2;
  text2.SetNDC(kTRUE);
  text2.SetTextFont(43);
  text2.SetTextSize(18);
  text2.DrawLatex(.2,.2, "all hadronic");
  gPad->RedrawAxis();
  c2->SaveAs(save_path+"/Resolution_Subjets/muon/AllHadronic/"+year+"/pt_"+mean_median+"_seperateJets_"+year+".pdf");

  TCanvas *c3 = new TCanvas();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  TGraphErrors* non_closure_percent = (TGraphErrors*) non_closure->Clone();
  for (int i=0;i<non_closure_percent->GetN();i++) non_closure_percent->GetY()[i] *= 100;
  TGraphErrors* area_percent = (TGraphErrors*) area->Clone();
  for (int i=0;i<area_percent->GetN();i++) area_percent->GetY()[i] *= 100;
  non_closure_percent->RemovePoint(n_last+2); // remove points at x=600 for drawing
  non_closure_percent->RemovePoint(n_last+1); // remove points at x=0 for drawing
  non_closure_percent->GetYaxis()->SetRangeUser(-5, 5);
  non_closure_percent->Draw("AP");
  area_percent->SetFillColor(16);
  area_percent->Draw("f SAME");
  zero_line->SetLineColor(kRed);
  zero_line->Draw("SAME");
  non_closure_percent->Draw("P SAME");
  TLegend *leg3 = new TLegend(0.45,0.85,0.80,0.65);
  leg3->AddEntry(non_closure_percent, mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]","pl");
  leg3->AddEntry(area_percent, "non-closure", "f");
  leg3->Draw();
  gPad->RedrawAxis();
  c3->SaveAs(save_path+"/Resolution_Subjets/muon/AllHadronic/"+year+"/nonClosure"+mean_median+"_"+year+".pdf");

  return 0;
}
//------------------------------------------------------------------------------
/*
.██████  ███████ ████████     ██████  ███████ ███████  ██████
██       ██         ██        ██   ██ ██      ██      ██    ██
██   ███ █████      ██        ██████  █████   ███████ ██    ██
██    ██ ██         ██        ██   ██ ██           ██ ██    ██
.██████  ███████    ██        ██   ██ ███████ ███████  ██████
*/

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median){
  int n_bins = hists.size();
  Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
  TH1F * reso = new TH1F("reso"," ",n_bins, bins);
  vector<double> mean, error;
  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1 = new TF1("f1","gaus",-0.2,0.2);
    hists[i]->Fit("f1", "QR");
    double m = f1->GetParameter(1);
    double e = f1->GetParError(1);
    if(!use_median) mean.push_back(m);
    error.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  TAxis *xaxis = reso->GetXaxis();
  for(unsigned int i=0; i<n_bins; i++){
    reso->Fill(xaxis->GetBinCenter(i+1), mean[i]);
    reso->SetBinError(i+1, error[i]);
  }
  return reso;
}
//------------------------------------------------------------------------------

/*
.██████  ███████ ████████     ███    ███ ███████  █████  ███    ██ ███████
██       ██         ██        ████  ████ ██      ██   ██ ████   ██ ██
██   ███ █████      ██        ██ ████ ██ █████   ███████ ██ ██  ██ ███████
██    ██ ██         ██        ██  ██  ██ ██      ██   ██ ██  ██ ██      ██
.██████  ███████    ██        ██      ██ ███████ ██   ██ ██   ████ ███████
*/

vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error){
  int n_bins = hists.size();
  vector<double> mean, err;

  for(unsigned int i=0; i<n_bins; i++){
    TF1* f1;
    if(do_ptrec) f1 = new TF1("f1","gaus");
    else         f1 = new TF1("f1","gaus", -0.2, 0.2);
    TString options = "QR";
    if(do_ptrec) options = "Q";
    hists[i]->Fit("f1", options);
    double m = f1->GetParameter(1);
    if(!use_median) mean.push_back(m);
    double e = f1->GetParError(1);
    err.push_back(e);

    // GET MEDIAN
    if(use_median){
      double values[] = {0};     // this array is filled by the GetQuantiles function
      double quantile[] = {0.5}; // this has to be an array with all the quantiles (only 0.5 for median)
      double Nquantiles = 1;     // length of array
      hists[i]->GetQuantiles(1, values, quantile);
      double median = values[0];
      mean.push_back(median);
    }
    ////
  }
  if(error == "error") return err;
  else return mean;
}

//------------------------------------------------------------------------------
/*
███    ██  ██████  ███    ██        ██████ ██       ██████  ███████ ██    ██ ██████  ███████
████   ██ ██    ██ ████   ██       ██      ██      ██    ██ ██      ██    ██ ██   ██ ██
██ ██  ██ ██    ██ ██ ██  ██ █████ ██      ██      ██    ██ ███████ ██    ██ ██████  █████
██  ██ ██ ██    ██ ██  ██ ██       ██      ██      ██    ██      ██ ██    ██ ██   ██ ██
██   ████  ██████  ██   ████        ██████ ███████  ██████  ███████  ██████  ██   ██ ███████
*/
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown){
  int steps = 150;
  double stepsize = 700/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  // first generate new TGraph with all positive values points
  const int Npoints = uncert->GetN();
  double xPoint[Npoints];
  double yPoint[Npoints];
  for(unsigned int i=0; i<Npoints; i++){
    double y;
    uncert->GetPoint(i, xPoint[i], y);
    yPoint[i] = sqrt(y*y); // take abs for y-value
  }
  TGraph* uncert_abs = new TGraph( Npoints, xPoint, yPoint );

  //

  for(int i=0; i<steps; i++){
    x = stepsize * i;
    up[i] = sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    down[i] = - sqrt(uncert_abs->Eval(x) * uncert_abs->Eval(x));
    xvalues[i] = x;
  }
  TGraph *uncertfit;
  if(updown == "area"){
    uncertfit = new TGraph(2*steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
      uncertfit->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
    }
  }
  else if(updown == "up"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],up[i]);
    }
  }
  else if(updown == "down"){
    uncertfit = new TGraph(steps);
    for (int i=0;i<steps;i++) {
      uncertfit->SetPoint(i,xvalues[i],down[i]);
    }
  }
  return uncertfit;
}
