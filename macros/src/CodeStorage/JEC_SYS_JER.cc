#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/Plotting.h"
#include "../include/tdrstyle_all.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TFitter.h"

#include "TDecompLU.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>

using namespace std;

MapHH GetHistograms(TString process, TString h_name);
void SubtractBackgroundsMap(MapHH& m_data, MapHH& m_st, MapHH& m_wjets, MapHH& m_other);
TH1F* SubtractBackgrounds(TH1F* data, vector<TH1F*> bgr, vector<double> syssize);
TMatrixD* ProduceCovMatrix(TH1F& map, const TString& save, const TString& process);
MapM GetCovMatrixMap(MapH& map, const TString& save, const TString& process);

typedef map<TString, TMatrixD> MapM;

TString year, channel, ichannel, var;
TString save, save_path, save_general, save_nfs;

bool debug = false;

int factor = 1e+6;
int precision = 5; // For double to string function - dtos()
int binTOremove = -1;

TString c_year = "combine";
TString c_channel = "combine";

TCanvas *c_sys = new TCanvas("sys", "sys", 600, 600); // to Plot all systematics in on pdf

int main(int argc, char* argv[]){

  gErrorIgnoreLevel = kError;
  gStyle->SetOptTitle(0);

  gErrorIgnoreLevel = kWarning;
  SetupGlobalStyle();

  save_path = get_save_path()+"/FSRuncertainty/"+year+"/"+ichannel;
  // if(!weights) save_path += "/noWeights";
  save_general = save_path; // For mass/pt windows

  save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/JER/";
  save = save_nfs;
  // if(!weights) save_nfs += "noWeights/";
  CreateSavePath((string) (save_nfs+"/"+year));

  TString hist = "JetMassScaleHists/wjet_jms_mass";
  MapHH m_data = GetHistograms("data", hist);
  MapHH m_ttbar = GetHistograms("ttbar", hist);
  MapHH m_st = GetHistograms("st", hist);
  MapHH m_wjets = GetHistograms("wjets", hist);
  MapHH m_other = GetHistograms("other", hist);
  MapHH m_JERup = GetHistograms("JERup", hist);
  MapHH m_JERdown = GetHistograms("JERdown", hist);

  MapHH m_bkg;
  for(auto channel: m_bkg){
    TString c = channel.first;
    for(auto year: channel.second){
      TString y = year.first;
      m_bkg; m_bkg[c];
      m_bkg[c][y] = (TH1F*) m_other[c][y]->Clone();
    }
  }

  cout << "Add Bkg ... " << endl;
  for(auto m: {m_wjets, m_st}){
    for(auto channel: m_bkg){
      for(auto year: m_bkg[channel.first]){
        m_bkg[channel.first][year.first]->Add(m[channel.first][year.first], 1);
      }
    }
  }


  cout << "Subtract Bkg from data ... " << endl;
  SubtractBackgroundsMap(m_data, m_st, m_wjets, m_other); // Subtract Bkg from data

  TH1F* h_data = Normalize(m_data[c_year][c_channel]);
  TH1F* h_ttbar = Normalize(m_ttbar[c_year][c_channel]);
  TH1F* h_JERup = Normalize(m_JERup[c_year][c_channel]);
  TH1F* h_JERdown = Normalize(m_JERdown[c_year][c_channel]);

  // Draw
  TH1F* ratio_ttbar   = GetRatio(h_ttbar,   h_ttbar, true);
  TH1F* ratio_JERup   = GetRatio(h_JERup,   h_ttbar, false);
  TH1F* ratio_JERdown = GetRatio(h_JERdown, h_ttbar, false);
  TH1F* ratio_data    = GetRatio(h_data,    h_ttbar, false);

  cout << "\t ... hh: draw plots" << endl;
  SetupCanvas(true);
  // void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
  Cosmetics(h_ttbar,      "#it{m}_{W} [GeV]", kBlack,  kSolid,  0, 2, true);
  Cosmetics(h_JERup,        "#it{m}_{W} [GeV]", kAzure+7, kSolid,  0, 2, true);
  Cosmetics(h_JERdown,        "#it{m}_{W} [GeV]", kAzure+7, kDashed,  0, 2, true);
  Cosmetics(h_data,       "#it{m}_{W} [GeV]", kBlack,   kSolid,  8, 2, true);

  // void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
  RatioCosmetics(ratio_ttbar, "#it{m}_{W} [GeV]", kBlack,  kSolid,  2, 0);
  RatioCosmetics(ratio_JERup,   "#it{m}_{W} [GeV]", kAzure+7, kSolid,  2, 0);
  RatioCosmetics(ratio_JERdown,   "#it{m}_{W} [GeV]", kAzure+7, kDashed,  2, 0);
  RatioCosmetics(ratio_data,  "#it{m}_{W} [GeV]", kBlack,   kSolid,  2, 8);

  m_rp1_top->cd();
  h_ttbar->GetXaxis()->SetRangeUser(30, 150);
  h_ttbar->GetYaxis()->SetRangeUser(0.001, h_ttbar->GetMaximum()*1.2);
  h_ttbar->Draw("hist ][");
  h_JERup->Draw("hist same ][");
  h_JERdown->Draw("hist same ][");
  h_data->Draw("pe same ][");
  gPad->RedrawAxis();
  DrawLegend({h_data, h_ttbar, h_JERup, h_JERdown}, "wmass_jer");

  // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
  DrawLumi(138., true, true, false, "wmass_jer");

  m_rp1->cd();
  ratio_ttbar->GetXaxis()->SetRangeUser(30, 150);
  ratio_ttbar->GetYaxis()->SetRangeUser(0.91, 1.19);
  ratio_ttbar->Draw("hist ][");
  ratio_JERup->Draw("hist same ][");
  ratio_JERdown->Draw("hist same ][");
  ratio_data->Draw("pe same ][");
  gPad->RedrawAxis();

  m_can->Print(save+"wmass_jer.pdf");

  int mDim = m_data[c_year][c_channel]->GetNbinsX();

  auto start = high_resolution_clock::now(); // Calculation time - start
  TString y,c,h;
  c = c_channel;
  y = c_year;
  h = "test";
  TMatrixD mtemp = GetCovMatrix(m_ttbar[c_year][c_channel], debug);
  TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, "ttbar");
  // DrawCov(htemp, save+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
  // bool onlyDiag = process.Contains("data")?false:true;
  // TMatrixD norm = NormCovMatrix(map[h], mtemp, false, onlyDiag);
  // TMatrixD* norm_pointer = new TMatrixD(norm);
  // TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
  // DrawCov(htemp_norm, save+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
  // auto stop = high_resolution_clock::now();  // Calculation time - stop
  // auto duration = duration_cast<seconds>(stop - start);
  // if(!onlyDiag) cout << process << " " << h << " took " << RED << duration.count() << "s" << RESET << endl;
  // int dim = map[h]->GetNbinsX();
  // covs[h].ResizeTo(dim,dim);
  //
  // covs[h] = norm;

}

// #################################################################################################
// #################################################################################################
// #################################################################################################

MapHH GetHistograms(TString process, TString h_name){
  VecTS collection_muon, collection_elec;
  if(process.EqualTo("data")){ collection_muon = data_muon; collection_elec = data_elec;}
  else if(process.EqualTo("ttbar")){ collection_muon = ttbar_muon; collection_elec = ttbar_elec;}
  else if(process.EqualTo("wjets")){ collection_muon = wjets_muon; collection_elec = wjets_elec;}
  else if(process.EqualTo("st")){ collection_muon = st_muon; collection_elec = st_elec;}
  else if(process.EqualTo("other")){ collection_muon = other_muon; collection_elec = other_elec;}
  else if(process.EqualTo("JERup")){ collection_muon = jer_up_muon; collection_elec = jer_up_elec;}
  else if(process.EqualTo("JERdown")){ collection_muon = jer_down_muon; collection_elec = jer_down_elec;}
  else throw runtime_error("Check the process ´"+process+"´ to obtain the histograms");

  MapHH map;
  vector<TH1F*> h_muon = get_all_hists(collection_muon, h_name);
  vector<TH1F*> h_elec = get_all_hists(collection_elec, h_name);
  vector<TH1F*> h_combine = combine_channels(h_muon, h_elec);

  for(auto nhist:h_muon) nhist->SetTitle(h_name);
  for(auto nhist:h_elec) nhist->SetTitle(h_name);
  for(auto nhist:h_combine) nhist->SetTitle(h_name);

  map["muon"]["2016"] = h_muon[0]; if(h_muon[0]->IsZombie()){throw runtime_error("NOOO h_muon[0]");};
  map["muon"]["2017"] = h_muon[1]; if(h_muon[1]->IsZombie()){throw runtime_error("NOOO h_muon[1]");};
  map["muon"]["2018"] = h_muon[2]; if(h_muon[2]->IsZombie()){throw runtime_error("NOOO h_muon[2]");};
  map["muon"]["combine"] = h_muon[3]; if(h_muon[3]->IsZombie()){throw runtime_error("NOOO h_muon[3]");};

  map["elec"]["2016"] = h_elec[0]; if(h_elec[0]->IsZombie()){throw runtime_error("NOOO h_elec[0]");};
  map["elec"]["2017"] = h_elec[1]; if(h_elec[1]->IsZombie()){throw runtime_error("NOOO h_elec[1]");};
  map["elec"]["2018"] = h_elec[2]; if(h_elec[2]->IsZombie()){throw runtime_error("NOOO h_elec[2]");};
  map["elec"]["combine"] = h_elec[3]; if(h_elec[3]->IsZombie()){throw runtime_error("NOOO h_elec[3]");};

  map["combine"]["2016"] = h_combine[0]; if(h_combine[0]->IsZombie()){throw runtime_error("NOOO ombine[0]");};
  map["combine"]["2017"] = h_combine[1]; if(h_combine[1]->IsZombie()){throw runtime_error("NOOO ombine[1]");};
  map["combine"]["2018"] = h_combine[2]; if(h_combine[2]->IsZombie()){throw runtime_error("NOOO ombine[2]");};
  map["combine"]["combine"] = h_combine[3]; if(h_combine[3]->IsZombie()){throw runtime_error("NOOO ombine[3]");};

  return map;
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
  for(int bin=1; bin<=nbins; bin++){
    if(result->GetBinContent(bin)<0){
      result->SetBinContent(bin, 0);
      result->SetBinError(bin, 0);
    }
  }

  return result;
}

void SubtractBackgroundsMap(MapHH& m_data, MapHH& m_st, MapHH& m_wjets, MapHH& m_other){
  // MapHHH map;
  for(auto channel: m_data){
    TString c = channel.first;
    for(auto year: channel.second){
      TString y = year.first;
      vector<TH1F*> bkgs = {m_st[c][y], m_wjets[c][y], m_other[c][y]};
      vector<double> bkgrate = {0.23, 0.19, 1.00};
      TH1F* hist = SubtractBackgrounds(m_data[c][y], bkgs, bkgrate);
      m_data[c][y] = hist;
    }
  }
}

// TMatrixD* ProduceCovMatrix(TH1F& map, const TString& save, const TString& process){
//   TMatrixD* covs = new TMatrixD();
//   auto start = high_resolution_clock::now(); // Calculation time - start
//
//   TString h = hist.first; TString c = channel; TString y = year;
//   TMatrixD mtemp = GetCovMatrix(map[h], debug);
//
//   TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
//   TString wbin = h; wbin.ReplaceAll(cut, "");
//   DrawCov(htemp, save+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
//
//   bool onlyDiag = process.Contains("data")?false:true;
//   TMatrixD norm = NormCovMatrix(map[h], mtemp, false, onlyDiag);
//   TMatrixD* norm_pointer = new TMatrixD(norm);
//   TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
//   DrawCov(htemp_norm, save+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);
//
//   auto stop = high_resolution_clock::now();  // Calculation time - stop
//   auto duration = duration_cast<seconds>(stop - start);
//   if(!onlyDiag) cout << process << " " << h << " took " << RED << duration.count() << "s" << RESET << endl;
//
//   int dim = map[h]->GetNbinsX();
//   covs[h].ResizeTo(dim,dim);
//   covs[h] = norm;
//
//   return covs;
// }

MapM GetCovMatrixMap(MapH& map, const TString& save, const TString& process){
  MapM covs;
  for(auto hist: map){
    auto start = high_resolution_clock::now(); // Calculation time - start

    TString h = hist.first; TString c = channel; TString y = year;
    TMatrixD mtemp = GetCovMatrix(map[h], debug);

    TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
    TString wbin = h; wbin.ReplaceAll(cut, "");
    DrawCov(htemp, save+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    bool onlyDiag = process.Contains("data")?false:true;
    TMatrixD norm = NormCovMatrix(map[h], mtemp, false, onlyDiag);
    TMatrixD* norm_pointer = new TMatrixD(norm);
    TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
    DrawCov(htemp_norm, save+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    auto stop = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<seconds>(stop - start);
    if(!onlyDiag) cout << process << " " << h << " took " << RED << duration.count() << "s" << RESET << endl;

    int dim = map[h]->GetNbinsX();
    covs[h].ResizeTo(dim,dim);
    covs[h] = norm;
  }
  return covs;
}
