#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

int main(int argc, char* argv[]){
  TString path_save = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS";
  TString path_plot = "/nfs/dust/cms/user/paaschal/MTopJet/PostSel/";
  vector<TString> years = {"2016", "2017", "2018"};
  vector<TH1F*> hists;
  TString hist_pt = "comparison_topjet_xcone_pass_rec/wjet_pt_match_S1divW";

  for(int year=0; year<3; year++){
    path_plot += years[year];
    TFile *file = new TFile(path_plot+"/muon/uhh2.AnalysisModuleRunner.MC.TTbar.root");
    hists.push_back((TH1F*)file->Get(hist_pt));
  }

  int n=15;
  TH1F *hist=AddHists(hists, 1);
  cout << hist->GetNbinsX() << endl;
  cout << hist->Integral() << endl;
  cout << hist->Integral(1, n) << endl;
  cout << hist->Integral(n+1, 51) << endl;

  double integral1 = hist->Integral(0, n) ;
  double integral2 = hist->Integral(n+1, 50);
  cout << integral1+integral2<<endl;


  // double sum=0; double sum1=0; double sum2=0;
  // for(int bin=0; bin<hist->GetNbinsX(); bin++){
  //   cout << "-------------------------------------- " << bin << endl;
  //   cout << hist->Integral(bin+1, bin+1) << endl;
  //   cout << hist->GetBinContent(bin+1) << endl;
  // }
}
