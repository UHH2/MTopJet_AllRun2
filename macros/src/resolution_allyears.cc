#include "../include/CentralInclude.h"

using namespace std;

// compile with:
// g++ -o resolution_subjets_flavorJEC resolution_subjets_flavorJEC.cc `root-config --cflags --glibs`

TH1F* GetResoPlot(vector<TH1F*> hists, bool use_median);
TGraph *AreaFromEnvelope(TGraph* uncert, TString updown);
vector<double> GetMeans(vector<TH1F*> hists, bool use_median, bool do_ptrec, TString error);

int main(int argc, char* argv[]){
  TString save_path = get_save_path();
  TString channel = "muon";

  // this is for naming of pdf files
  bool use_median = false;
  // if(argc > 1 && strcmp(argv[1], "median") == 0) use_median = true;
  TString mean_median;
  if(use_median) mean_median = "median";
  else           mean_median = "mean";

  vector<TString> years = {"2016", "2017", "2018"};
  vector<TH1F*>   hists_raw, hists_jec, hists_xcone;

  for(auto year: years){
    // declare files
    TString filedir;
    if(channel == "muon") filedir = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/"+year+"/muon/";
    else filedir = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/"+year+"/elec/";
    TFile *f_tt = new TFile(filedir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
    // set binning
    int n_ptbin = 10;
    Float_t bins[] = {0, 50, 80, 120, 170, 220, 270, 320, 370, 420, 600};
    TString bin_strings[] = {"0", "50", "80", "120", "170", "220", "270", "320", "370", "420", "600"};

    // now get plots from file
    vector<TH1F*> reso, reso_noJEC, reso_cor;
    vector<TH1F*> ptrec;

    std::string dir_jec = "RecGenHists_subjets/";
    std::string dir_raw = "RecGenHists_subjets_noJEC/";
    std::string dir_cor = "RecGenHists_subjets_corrected/";

    std::string name_jec;
    std::string name_raw;
    std::string name_cor;
    std::string name_ptrec;

    for(int ptbin = 1; ptbin <= n_ptbin; ptbin++){

      name_jec = dir_jec + "PtResolution_" + std::to_string(ptbin);
      name_raw = dir_raw + "PtResolution_" + std::to_string(ptbin);
      name_cor = dir_cor + "PtResolution_" + std::to_string(ptbin);
      name_ptrec = dir_jec + "PtRec_" + std::to_string(ptbin);

      reso.push_back( (TH1F*)f_tt->Get(name_jec.c_str()) );
      reso_noJEC.push_back( (TH1F*)f_tt->Get(name_raw.c_str()) );
      reso_cor.push_back( (TH1F*)f_tt->Get(name_cor.c_str()) );
      ptrec.push_back( (TH1F*)f_tt->Get(name_ptrec.c_str()) );
    }

    // fit histograms and extract mean/median value
    TH1F * resolution       = GetResoPlot(reso, use_median);
    TH1F * resolution_noJEC = GetResoPlot(reso_noJEC, use_median);
    TH1F * resolution_cor   = GetResoPlot(reso_cor, use_median);

    hists_raw.push_back(resolution_noJEC);
    hists_jec.push_back(resolution);
    hists_xcone.push_back(resolution_cor);
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // calculate non-closure as function of ptrec
    // arguments: (hists, use median?, do ptrec?, mean or error)
    vector<double> pt                = GetMeans(ptrec, use_median, true, "mean");
    vector<double> pt_err            = GetMeans(ptrec, use_median, true, "error");
    vector<double> MeanForUncert     = GetMeans(reso_cor, use_median, false, "mean");
    vector<double> MeanForUncert_err = GetMeans(reso_cor, use_median, false, "error");
    ////
    TGraphErrors* non_closure = new TGraphErrors(n_ptbin, &pt[0], &MeanForUncert[0], &pt_err[0], &MeanForUncert_err[0]);
    // now insert a point at ptrec = 0 to get a const uncertainty from 0 to first point of ptrec
    // do same to reach ptrec = 600
    non_closure->SetPoint(non_closure->GetN(), 0, MeanForUncert[0]);
    int n_last = MeanForUncert.size()-1;
    non_closure->SetPoint(non_closure->GetN(), 600, MeanForUncert[n_last]);

    TGraph* area    = AreaFromEnvelope(non_closure, "area");
    TGraph* upvar   = AreaFromEnvelope(non_closure, "up");
    TGraph* downvar = AreaFromEnvelope(non_closure, "down");

    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // non-closure in root file
    // TString rootname = "files/Correction_SysFromResolution.root";
    // TFile * outfile = new TFile(rootname,"RECREATE");;
    // upvar->Write("Up");
    // downvar->Write("Down");
    // outfile->Close();
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // start plotting
  }

  vector<vector<TH1F*>> hists_all = {hists_raw, hists_jec, hists_xcone};

  TLine *zero_line = new TLine(0, 0, 600, 0);
  zero_line->SetLineColor(15);
  zero_line->SetLineWidth(3);
  zero_line->SetLineStyle(7);

  // Canvas properties
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetLegendBorderSize(0);

  int index=0;
  for(auto hists: hists_all){
    hists[0]->GetYaxis()->SetRangeUser(-0.1, 0.1);
    hists[0]->GetXaxis()->SetNdivisions(505);
    hists[0]->GetYaxis()->SetNdivisions(505);
    hists[0]->GetXaxis()->SetTitleSize(0.05);
    hists[0]->GetYaxis()->SetTitleSize(0.04);
    hists[0]->GetXaxis()->SetTitleOffset(0.9);
    hists[0]->GetYaxis()->SetTitleOffset(1.5);
    hists[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
    hists[0]->GetYaxis()->SetTitle(mean_median+" #left[ #frac{p_{T}^{rec} - p_{T}^{gen}}{p_{T}^{gen}} #right]");
    hists[0]->SetLineWidth(4);
    hists[0]->SetLineColor(kAzure+7);
    hists[1]->SetLineWidth(4);
    hists[1]->SetLineColor(kOrange+1);
    hists[2]->SetLineWidth(4);
    hists[2]->SetLineColor(kRed+1);

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    hists[0]->Draw("E1");
    zero_line->Draw("SAME");
    hists[0]->Draw("SAME E1");
    hists[1]->Draw("SAME E1");
    hists[2]->Draw("SAME E1");
    TLegend *leg1 = new TLegend(0.45,0.85,0.80,0.65);
    leg1->AddEntry(hists[0],"2016","l");
    leg1->AddEntry(hists[1],"2017","l");
    leg1->AddEntry(hists[2],"2018","l");
    leg1->SetTextSize(0.05);
    leg1->Draw();
    TLatex text1;
    text1.SetNDC(kTRUE);
    text1.SetTextFont(43);
    text1.SetTextSize(18);
    text1.DrawLatex(.2,.2, "lepton+jets");
    gPad->RedrawAxis();
    // first save without additional correction
    if(index==0) c1->SaveAs(save_path+"/Plots/Resolution_Subjets/"+channel+"/pt_"+mean_median+"_all_years_raw.pdf");
    if(index==1) c1->SaveAs(save_path+"/Plots/Resolution_Subjets/"+channel+"/pt_"+mean_median+"_all_years_jec.pdf");
    if(index==2) c1->SaveAs(save_path+"/Plots/Resolution_Subjets/"+channel+"/pt_"+mean_median+"_all_years_xcone.pdf");

    index++;
  }
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
