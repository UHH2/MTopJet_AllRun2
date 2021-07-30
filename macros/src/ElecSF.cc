#include "../include/CentralInclude.h"


using namespace std;

// variables from root tree
Double_t pt;
Double_t eta;
Double_t weight;
Double_t weight_sfpt;
Double_t weight_sfeta;
Double_t weight_sfetapt;
Double_t weight_sfetaptUP;
Double_t weight_sfetaptDOWN;
Int_t run;
Int_t lumi;
Int_t eventnr;
bool passed;
bool debug = false;
TString weighttag;
TString year, year_v;

// filling hists
void fill_pteta(TTree* tree, vector<TH1F*> h_pt, vector<TH1F*> h_eta);
// Calculate SF
vector<TH1F*> GetSF(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TH1F* hist);
// Plotter
void PlotEfficiency(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TString xaxis, TString histname);
void PlotHist(TH1F* h_data, TString xaxis, TString histname);
void PlotSF(vector<TH1F*> h_SF, TString xaxis, TString histname);


int main(int argc, char* argv[]){

  for(auto s:argv) cout << s << endl;

  year = argv[1];
  if(year.EqualTo("2016")) year_v = "_2016v3";
  else if(year.EqualTo("2017")) year_v = "_2017v2";
  else if(year.EqualTo("2018")) year_v = "_2018";
  else throw runtime_error("I need the correct year; 2016, 2017 or 2018");

  if(argc == 2){
    cout << "Calculate ElecTrigger SF for "+year+" dependent on pt and eta..." << endl << endl;
    weighttag = "weight";
  }
  else if(argc == 3){
    if(strcmp(argv[2], "pt") == 0){
      cout << "Cross check: apply pt SF and calculate again..." << endl << endl;
      weighttag = "weight_sfpt";
    }
    else if(strcmp(argv[2], "eta") == 0){
      cout << "Cross check: apply eta SF and calculate again..." << endl << endl;
      weighttag = "weight_sfeta";
    }
    else if(strcmp(argv[2], "etapt") == 0){
      cout << "Cross check: apply eta (pt binned) SF and calculate again..." << endl << endl;
      weighttag = "weight_sfetapt";
    }
    else if(strcmp(argv[2], "etaptup") == 0){
      cout << "Cross check: apply eta (pt binned) SF (UP variation) and calculate again..." << endl << endl;
      weighttag = "weight_sfetaptUP";
    }
    else if(strcmp(argv[2], "etaptdown") == 0){
      cout << "Cross check: apply eta (pt binned) SF (DOWN variation) and calculate again..." << endl << endl;
      weighttag = "weight_sfetaptDOWN";
    }
    else{
      cout << "[ERROR] no valid option selected!" << endl;
      return 0;
    }
  }

  vector<double> pt_bins = {55, 75, 95, 115, 135, 155, 175, 200, 300, 1500};
  vector<TH1F*> h_pt_data, h_pt_mc;
  h_pt_data.push_back(new TH1F("h_pt_all_data","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_data.push_back(new TH1F("h_pt_pass_data","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_data.push_back(new TH1F("h_pt_all_data_weirdbin","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_data.push_back(new TH1F("h_pt_pass_data_weirdbin","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_data.push_back(new TH1F("h_pt_all_data_barrel","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_data.push_back(new TH1F("h_pt_pass_data_barrel","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_data.push_back(new TH1F("h_pt_all_data_endcap","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_data.push_back(new TH1F("h_pt_pass_data_endcap","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_mc.push_back(new TH1F("h_pt_all_mc","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_mc.push_back(new TH1F("h_pt_pass_mc","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_mc.push_back(new TH1F("h_pt_all_mc_weirdbin","pt", pt_bins.size()-1, &pt_bins[0]));
  h_pt_mc.push_back(new TH1F("h_pt_pass_mc_weirdbin","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_mc.push_back(new TH1F("h_pt_all_mc_barrel","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_mc.push_back(new TH1F("h_pt_pass_mc_barrel","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_mc.push_back(new TH1F("h_pt_all_mc_endcap","pt", pt_bins.size()-1, &pt_bins[0]));
  // h_pt_mc.push_back(new TH1F("h_pt_pass_mc_endcap","pt", pt_bins.size()-1, &pt_bins[0]));

  vector<double> eta_bins = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
  vector<TH1F*> h_eta_data, h_eta_mc;
  h_eta_data.push_back(new TH1F("h_eta_all_data","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_pass_data","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_all_data_lowpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_pass_data_lowpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_all_data_midpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_pass_data_midpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_all_data_highpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_pass_data_highpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_data.push_back(new TH1F("h_eta_fail_data_highpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_all_mc","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_pass_mc","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_all_mc_lowpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_pass_mc_lowpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_all_mc_midpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_pass_mc_midpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_all_mc_highpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_pass_mc_highpt","eta", eta_bins.size()-1, &eta_bins[0]));
  h_eta_mc.push_back(new TH1F("h_eta_fail_mc_highpt","eta", eta_bins.size()-1, &eta_bins[0]));

  TFile *f_data=new TFile("/nfs/dust/cms/user/paaschal/MTopJet_Run2/ElecSF/uhh2.AnalysisModuleRunner.DATA.DATA"+year_v+".root");
  fill_pteta((TTree *) f_data->Get("AnalysisTree"), h_pt_data, h_eta_data);

  TFile *f_tt=new TFile("/nfs/dust/cms/user/paaschal/MTopJet_Run2/ElecSF/uhh2.AnalysisModuleRunner.MC.TTbar"+year_v+".root");
  fill_pteta((TTree *) f_tt->Get("AnalysisTree"), h_pt_mc, h_eta_mc);

  TFile *f_st=new TFile("/nfs/dust/cms/user/paaschal/MTopJet_Run2/ElecSF/uhh2.AnalysisModuleRunner.MC.SingleTop"+year_v+".root");
  fill_pteta((TTree *) f_st->Get("AnalysisTree"), h_pt_mc, h_eta_mc);

  cout << h_pt_mc[1]->GetEntries() << "\t" << h_pt_mc[0]->GetEntries() << endl;

  if(debug) cout << "Plots Hists ... " << endl;
  if(debug) cout << "\n\t ... data" << endl;
  PlotHist(h_pt_data[0], "p_{T}", "data_pt_all");
  PlotHist(h_pt_data[1], "p_{T}", "data_pt_pass");
  PlotHist(h_pt_data[2], "p_{T}", "data_pt_fail_highpt");
  PlotHist(h_pt_data[3], "p_{T}", "data_pt_pass_highpt");
  PlotHist(h_eta_data[0], "#eta", "data_eta_all");
  PlotHist(h_eta_data[1], "#eta", "data_eta_pass");
  PlotHist(h_eta_data[2], "#eta", "data_eta_lowpt_all");
  PlotHist(h_eta_data[3], "#eta", "data_eta_lowpt_pass");
  PlotHist(h_eta_data[4], "#eta", "data_eta_midpt_all");
  PlotHist(h_eta_data[5], "#eta", "data_eta_midpt_pass");
  PlotHist(h_eta_data[6], "#eta", "data_eta_highpt_all");
  PlotHist(h_eta_data[7], "#eta", "data_eta_highpt_pass");
  PlotHist(h_eta_data[8], "#eta", "data_eta_highpt_fail");

  if(debug) cout << "\n\t ... MC" << endl;
  PlotHist(h_pt_mc[0], "p_{T}", "mc_pt_all");
  PlotHist(h_pt_mc[1], "p_{T}", "mc_pt_pass");
  PlotHist(h_pt_mc[2], "p_{T}", "mc_pt_fail_highpt");
  PlotHist(h_pt_mc[3], "p_{T}", "mc_pt_pass_highpt");
  PlotHist(h_eta_mc[0], "#eta", "mc_eta_all");
  PlotHist(h_eta_mc[1], "#eta", "mc_eta_pass");
  PlotHist(h_eta_mc[2], "#eta", "mc_eta_lowpt_all");
  PlotHist(h_eta_mc[3], "#eta", "mc_eta_lowpt_pass");
  PlotHist(h_eta_mc[4], "#eta", "mc_eta_midpt_all");
  PlotHist(h_eta_mc[5], "#eta", "mc_eta_midpt_pass");
  PlotHist(h_eta_mc[6], "#eta", "mc_eta_highpt_all");
  PlotHist(h_eta_mc[7], "#eta", "mc_eta_highpt_pass");
  PlotHist(h_eta_mc[8], "#eta", "mc_eta_highpt_fail");

  if(debug) cout << "Create efficiency graphs ... " << endl;
  TGraphAsymmErrors* h_effi_pt_data = new TGraphAsymmErrors(h_pt_data[1], h_pt_data[0],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_data = new TGraphAsymmErrors(h_eta_data[1], h_eta_data[0],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_lowpt_data = new TGraphAsymmErrors(h_eta_data[3], h_eta_data[2],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_midpt_data = new TGraphAsymmErrors(h_eta_data[5], h_eta_data[4],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_highpt_data = new TGraphAsymmErrors(h_eta_data[7], h_eta_data[6],"cl=0.683 b(1,1) mode");

  cout << h_pt_mc[1]->GetEntries() << "\t" << h_pt_mc[0]->GetEntries() << endl;
  for(int i=1; i<=h_pt_mc[0]->GetNbinsX();i++) cout << "\t" << h_pt_mc[0]->GetBinContent(i);
  cout << endl;
  for(int i=1; i<=h_pt_mc[0]->GetNbinsX();i++) cout << "\t" << h_pt_mc[0]->GetBinCenter(i) << " (" << h_pt_mc[0]->GetBinWidth(i) << ")";
  cout << endl;
  cout << endl;
  for(int i=1; i<=h_pt_mc[1]->GetNbinsX();i++) cout << "\t" << h_pt_mc[1]->GetBinContent(i);
  cout << endl;
  for(int i=1; i<=h_pt_mc[1]->GetNbinsX();i++) cout << "\t" << h_pt_mc[1]->GetBinCenter(i) << " (" << h_pt_mc[1]->GetBinWidth(i) << ")";
  cout << endl;
  TGraphAsymmErrors* h_effi_pt_mc = new TGraphAsymmErrors(h_pt_mc[1], h_pt_mc[0],"cl=0.683 b(1,1) mode");

  cout << h_eta_mc[1]->GetNbinsX() << "\t" << h_eta_mc[0]->GetNbinsX() << endl;
  for(int i=1; i<=h_eta_mc[0]->GetNbinsX();i++) cout << h_eta_mc[1]->GetBinContent(i) << "\t" << h_eta_mc[0]->GetBinContent(i) << endl;
  TGraphAsymmErrors* h_effi_eta_mc = new TGraphAsymmErrors(h_eta_mc[1], h_eta_mc[0],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_lowpt_mc = new TGraphAsymmErrors(h_eta_mc[3], h_eta_mc[2],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_midpt_mc = new TGraphAsymmErrors(h_eta_mc[5], h_eta_mc[4],"cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* h_effi_eta_highpt_mc = new TGraphAsymmErrors(h_eta_mc[7], h_eta_mc[6],"cl=0.683 b(1,1) mode");

  if(debug) cout << "Plot efficiency graphs ... " << endl;
  PlotEfficiency(h_effi_pt_data, h_effi_pt_mc, "p_{T}", "effi_pt");
  PlotEfficiency(h_effi_eta_data, h_effi_eta_mc, "#eta", "effi_eta");
  PlotEfficiency(h_effi_eta_lowpt_data, h_effi_eta_lowpt_mc, "#eta", "effi_eta_lowpt");
  PlotEfficiency(h_effi_eta_midpt_data, h_effi_eta_midpt_mc, "#eta", "effi_eta_midpt");
  PlotEfficiency(h_effi_eta_highpt_data, h_effi_eta_highpt_mc, "#eta", "effi_eta_highpt");

  if(debug) cout << "\t ... pt" << endl;
  vector<TH1F*> h_SF_pt = GetSF(h_effi_pt_data, h_effi_pt_mc, h_pt_data[0]);
  if(debug) cout << "\t ... eta" << endl;
  vector<TH1F*> h_SF_eta = GetSF(h_effi_eta_data, h_effi_eta_mc, h_eta_data[0]);
  if(debug) cout << "\t ... lowpt" << endl;
  vector<TH1F*> h_SF_eta_lowpt = GetSF(h_effi_eta_lowpt_data, h_effi_eta_lowpt_mc, h_eta_data[2]);
  if(debug) cout << "\t ... midpt" << endl;
  vector<TH1F*> h_SF_eta_midpt = GetSF(h_effi_eta_midpt_data, h_effi_eta_midpt_mc, h_eta_data[4]);
  if(debug) cout << "\t ... highpt" << endl;
  vector<TH1F*> h_SF_eta_highpt = GetSF(h_effi_eta_highpt_data, h_effi_eta_highpt_mc, h_eta_data[4]);


  if(debug) cout << "Get SF ... " << endl;
  int bin1 = h_SF_eta[0]->GetXaxis()->FindBin(0.25);
  int bin2 = h_SF_eta[0]->GetXaxis()->FindBin(0.1);
  cout <<h_SF_eta[0]->GetBinContent(bin1) << endl;
  cout <<h_SF_eta[0]->GetBinContent(bin2) << endl;


  if(debug) cout << "Plot SF ... " << endl;
  PlotSF(h_SF_pt, "p_{T}", "SF_pt");
  PlotSF(h_SF_eta, "#eta", "SF_eta");
  PlotSF(h_SF_eta_lowpt, "#eta", "SF_eta_lowpt");
  PlotSF(h_SF_eta_midpt, "#eta", "SF_eta_midpt");
  PlotSF(h_SF_eta_highpt, "#eta", "SF_eta_highpt");

  if(argc == 2){
    cout << "Creating file..." << endl;
    TFile * f_out = new TFile("files/ElecTriggerSF"+year+".root","RECREATE");;
    h_SF_pt[0]->Write("Central_pt");
    h_SF_pt[1]->Write("Up_pt");
    h_SF_pt[2]->Write("Down_pt");
    h_SF_eta[0]->Write("Central_eta");
    h_SF_eta[1]->Write("Up_eta");
    h_SF_eta[2]->Write("Down_eta");
    h_SF_eta_lowpt[0]->Write("Central_eta_lowpt");
    h_SF_eta_lowpt[1]->Write("Up_eta_lowpt");
    h_SF_eta_lowpt[2]->Write("Down_eta_lowpt");
    h_SF_eta_midpt[0]->Write("Central_eta_midpt");
    h_SF_eta_midpt[1]->Write("Up_eta_midpt");
    h_SF_eta_midpt[2]->Write("Down_eta_midpt");
    h_SF_eta_highpt[0]->Write("Central_eta_highpt");
    h_SF_eta_highpt[1]->Write("Up_eta_highpt");
    h_SF_eta_highpt[2]->Write("Down_eta_highpt");
    f_out->Close();
  }

  return 0;
}



void fill_pteta(TTree* tree, vector<TH1F*> h_pt, vector<TH1F*> h_eta){
  if(!tree) cout << "could not read tree\n";
  else      cout << "Filling Histograms...\n";

  // outputFile->cd();
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("eta",&eta);
  tree->SetBranchAddress("passed",&passed);
  tree->SetBranchAddress("weight_sfpt",&weight_sfpt);
  tree->SetBranchAddress("weight_sfeta",&weight_sfeta);
  tree->SetBranchAddress("weight_sfetapt",&weight_sfetapt);
  tree->SetBranchAddress("weight_sfetaptUP",&weight_sfetaptUP);
  tree->SetBranchAddress("weight_sfetaptDOWN",&weight_sfetaptDOWN);
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("eventnr",&eventnr);

  tree->SetBranchStatus("*",1);

  int ipass = 0; int ifull = 0;
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    ifull++;
    h_pt[0]->Fill(pt, weight);
    h_eta[0]->Fill(eta, weight);
    if(pt < 120)             h_eta[2]->Fill(eta, weight);
    else if(pt>120 && pt<200)h_eta[4]->Fill(eta, weight);
    else if(pt>200)          h_eta[6]->Fill(eta, weight);

    if(!passed&&debug){
      if(pt>200){
        h_pt[2]->Fill(pt, weight);
        h_eta[8]->Fill(eta, weight);
        cout << "--------------------" << endl;
        cout << "run:        " << run << endl;
        cout << "lumi block: " << lumi << endl;
        cout << "event nr:   " << eventnr << endl;
      }
    }

    if(passed){
      ipass++;
      // cout << pt << ", " << eta << endl;
      // cout  << weight << ", "<< weight_sfpt << ", " << weight_sfeta << ", " << weight_sfetapt << endl;
      double weight_pass = 0;
      if(weighttag == "weight_sfpt") weight_pass = weight_sfpt;
      else if(weighttag == "weight_sfeta") weight_pass = weight_sfeta;
      else if(weighttag == "weight_sfetapt") weight_pass = weight_sfetapt;
      else if(weighttag == "weight_sfetaptUP") weight_pass = weight_sfetaptUP;
      else if(weighttag == "weight_sfetaptDOWN") weight_pass = weight_sfetaptDOWN;
      else weight_pass = weight;
      // if(year=="2017"&&ievent==10) cout << "weight: " << weight << "\t" << "weight_pass: " << weight_pass << endl;
      if(weight<weight_pass) cout << "UHHHHHHHH" << endl;
      h_pt[1]->Fill(pt, weight_pass);
      h_eta[1]->Fill(eta, weight_pass);
      if(pt < 120)             h_eta[3]->Fill(eta, weight_pass);
      else if(pt>120 && pt<200)h_eta[5]->Fill(eta, weight_pass);
      else if(pt>200)          h_eta[7]->Fill(eta, weight_pass);

      if(pt>200) h_pt[3]->Fill(pt, weight);
    }
  }
  cout << h_pt[1]->GetBinContent(9) << "\t" << h_pt[0]->GetBinContent(9) << endl;
  cout << ifull << "\t" << ipass << endl;
  return;
}

// calculate central value and up/down variation of SF and return in vector
// [0] - central
// [1] - up
// [2] - down
vector<TH1F*> GetSF(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TH1F* hist){
  TH1F* h_central = (TH1F*) hist->Clone();
  TH1F* h_up = (TH1F*) hist->Clone();
  TH1F* h_down = (TH1F*) hist->Clone();
  int nbins = hist->GetSize() -2;

  for(unsigned int i=0; i<nbins; i++){
    double xvalue = hist->GetBinCenter(i+1);
    double value_data = h_data->Eval(xvalue);
    double e_data_hi = h_data->GetErrorYhigh(i);
    double e_data_lo = h_data->GetErrorYlow(i);

    double value_mc = h_mc->Eval(xvalue);
    double e_mc_hi = h_mc->GetErrorYhigh(i);
    double e_mc_lo = h_mc->GetErrorYlow(i);

    double central = 1;
    if(value_data != 0 && value_mc != 0) central = value_data/value_mc;

    double up = 0;
    if(value_mc != 0) up = -value_data/(value_mc*value_mc) * e_mc_hi + e_data_hi/value_mc;

    double down = 0;
    if(value_mc != 0) down = -value_data/(value_mc*value_mc) * e_mc_lo + e_data_lo/value_mc;

    // add 1% uncertainty to cover non-closure
    up += 0.01;
    down += 0.01;

    // cout << i<< ", " << xvalue << ", " << value_mc << ", " << value_data << ", "<< central << endl;
    h_central->SetBinContent(i+1, central);
    h_up->SetBinContent(i+1, central+up);
    h_down->SetBinContent(i+1, central-down);
    // bin errors are just set very small to not annoy in plotting
    h_central->SetBinError(i+1, 0.000000000000001);
    h_up->SetBinError(i+1, 0.000000000000001);
    h_down->SetBinError(i+1, 0.000000000000001);
  }
  vector<TH1F*> h_sf;
  h_sf.push_back(h_central);
  h_sf.push_back(h_up);
  h_sf.push_back(h_down);
  return h_sf;
}


void PlotEfficiency(TGraphAsymmErrors* h_data, TGraphAsymmErrors* h_mc, TString xaxis, TString histname){
  // if(histname == "effi_pt") cout << h_data->Eval(25) << endl;
  // if(histname == "effi_pt") cout << h_mc->Eval(25) << endl;

  h_data->SetTitle(" ");
  h_data->GetXaxis()->SetTitle(xaxis);
  h_data->GetYaxis()->SetTitle("efficiency");
  h_data->GetYaxis()->SetTitleOffset(1.6);
  h_data->GetXaxis()->SetTitleOffset(1.3);
  h_data->GetXaxis()->SetNdivisions(505);
  h_data->GetYaxis()->SetNdivisions(505);
  h_data->GetYaxis()->SetRangeUser(0.0, 1.0);

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);

  h_mc->SetLineColor(kRed);
  h_mc->SetMarkerColor(kRed);
  h_mc->SetMarkerStyle(8);
  h_mc->SetMarkerSize(1);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h_data->Draw("AP");
  h_mc->Draw("P SAME");
  TLegend *leg = new TLegend(0.33,0.20,0.66,0.33);
  leg->AddEntry(h_data,"data","pl");
  leg->AddEntry(h_mc,"simulation","pl");
  leg->Draw("");
  gPad->RedrawAxis();
  if(weighttag == "weight_sfpt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_pt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfeta") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_eta/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetapt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etapt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptUP") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptUP/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptDOWN") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptDOWN/"+histname+"_"+year+".pdf");
  else A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/"+histname+"_"+year+".pdf");
  delete A;
  return;
}

void PlotHist(TH1F* hist, TString xaxis, TString histname){
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle(xaxis);
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.6);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum()*1.2);
  hist->SetLineColor(kBlack);
  hist->SetFillColor(13);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  if(weighttag == "weight_sfpt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_pt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfeta") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_eta/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetapt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etapt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptUP") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptUP/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptDOWN") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptDOWN/"+histname+"_"+year+".pdf");
  else A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/"+histname+"_"+year+".pdf");
  delete A;
  return;
}

void PlotSF(vector<TH1F*> h_SF, TString xaxis, TString histname){
  if(h_SF.size() != 3) cout << "THERE ARE NOT EXACTLY 3 SF hists" << endl;
  h_SF[0]->SetTitle(" ");
  h_SF[0]->GetXaxis()->SetTitle(xaxis);
  h_SF[0]->GetYaxis()->SetTitle("scale factor");
  h_SF[0]->GetYaxis()->SetTitleOffset(1.6);
  h_SF[0]->GetXaxis()->SetTitleOffset(1.3);
  h_SF[0]->GetXaxis()->SetNdivisions(505);
  h_SF[0]->GetYaxis()->SetNdivisions(505);
  h_SF[0]->GetYaxis()->SetRangeUser(0.7, 1.05);
  h_SF[0]->SetLineColor(kBlack);
  h_SF[0]->SetLineWidth(3);
  h_SF[0]->SetFillColor(13);
  h_SF[0]->SetMarkerStyle(0);

  Color_t uncert_col = kOrange+9;

  h_SF[1]->SetLineColor(uncert_col);
  h_SF[1]->SetFillStyle(3144);
  h_SF[1]->SetFillColor(uncert_col);
  h_SF[2]->SetLineColor(uncert_col);
  h_SF[2]->SetFillColor(10);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.1);
  h_SF[0]->Draw("E1");
  h_SF[1]->Draw("HIST SAME");
  h_SF[2]->Draw("HIST SAME");
  h_SF[0]->Draw("E1 SAME");

  TLegend *leg = new TLegend(0.33,0.20,0.66,0.33);
  leg->AddEntry(h_SF[0],"scale factor","l");
  leg->AddEntry(h_SF[1],"uncertainty","f");
  leg->Draw();

  gPad->RedrawAxis();
  if(weighttag == "weight_sfpt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_pt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfeta") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_eta/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetapt") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etapt/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptUP") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptUP/"+histname+"_"+year+".pdf");
  else if(weighttag == "weight_sfetaptDOWN") A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/CLOSURE_etaptDOWN/"+histname+"_"+year+".pdf");
  else A->SaveAs("/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/ElecTriggerSF/"+histname+"_"+year+".pdf");
  delete A;
  return;
}
