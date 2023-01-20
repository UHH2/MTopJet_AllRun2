#include "../include/CentralInclude.h"
#include "../include/CovarianzMatrix.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/tdrstyle_all.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <algorithm>

using namespace std;

typedef map<int, TGraph2DErrors*> MapI2DGe;
typedef map<TString, MapI2DGe> MapSI2DGe;

typedef map<TString, TMatrixD*> MapM;
typedef map<TString, MapM> MapMM;
typedef map<TString, MapMM> MapMMM;

typedef map<TString, TF2*> MapF2;
typedef map<TString, MapF2> MapFF2;
typedef map<int, TF2*> MapIF2;
typedef map<TString, MapIF2> MapSIF2;

// ---------------------------------------------------------------------------
// Declare functions

// MapHHH SetUpMap(VecTS hists);
MapHH GetHistograms(TString process, TString h_name);
void AddTwoMaps(MapHHH map1, MapHHH map2, int option);
MapHHH RebinMap(MapHHH process, int width);
MapMMM GetCovMatrixMap(TFile* file, MapHHH map, TString save, TString process);
TH1F* SubtractBackgrounds(const TH1F* data, const VecH &bgr);

bool debug = true;
TString save, year, channel;
TString schannel = "";
VecTS hists = {
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_hl",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_lh",
  "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_ll"
};
TString cut = "comparison_topjet_xcone_pass_rec/wmass_match_ptdiv_";
int space;

int main(int argc, char* argv[]){

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  // Input ---------------------------------------------------------------------
  year = argv[1];
  channel = argv[2];
  if(channel.EqualTo("muon")||channel.EqualTo("elec")) schannel = "_"+channel;
  // TString sub = argv[3]; bool isSUB = false;
  // if(sub.EqualTo("1")) isSUB = true;
  // cout << isSUB << endl;
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptTitle(0);
  SetupGlobalStyle();

  // year ---------------------------------------------------------------------
  bool is16=false; bool is17=false; bool is18=false; bool isAll=false; bool is1718=false;
  if     (strcmp(year, "2016")==0)     is16  = true;
  else if(strcmp(year, "2017")==0)     is17  = true;
  else if(strcmp(year, "2018")==0)     is18  = true;
  else if(strcmp(year, "combine")==0)  isAll = true;
  // else cerr("Give me the correct year please (2016, 2017, 2018 or combine)");

  // Bins ----------------------------------------------------------------------
  int bin_width = 1;
  int number_bins = 180/bin_width;
  TString str_number_bins = to_string(number_bins); // For creating folders

  // Limit for bin content -----------------------------------------------------
  int limit = (isAll)?100:40;
  TString slimit = to_string(limit);


  // Create Directories --------------------------------------------------------
  if(debug) cout << "Create Directiories ..." << endl;

  save = get_save_path();
  save += "/JetCorrections/fit/"+year+"/"+channel+"/"+slimit+"/"+"BinWidth_"+to_string(bin_width)+"/";
  CreateSavePath((string) save);
  cout << save << endl;
  TString save_root = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/files/";

  cout << "Start collecting Histograms ..." << endl;
  // ---------------------------------------------------------------------------
  // I load all histograms; This code was created after studies ended. Only one
  // option is available. If necessary include other options as well (no bins etc.)

  TString w_mass_hh = hists[0];
  TString w_mass_hl = hists[1];
  TString w_mass_lh = hists[2];
  TString w_mass_ll = hists[3];

  MapHHHH histograms;
  MapHHH m_data, m_ttbar, m_st, m_wjets, m_other, m_JECup, m_JECdown, m_CORup, m_CORdown;

  // MapHHH GetHistograms(channel, muon, elec, combine, nhist)
  if(debug) cout << "\t ... Get Histograms" << endl;
  for(TString hist: hists)
  {
    m_data[hist] = GetHistograms("data", hist);
    m_ttbar[hist] = GetHistograms("ttbar", hist);
    m_st[hist] = GetHistograms("st", hist);
    m_wjets[hist] = GetHistograms("wjets", hist);
    m_other[hist] = GetHistograms("other", hist);
    m_JECup[hist] = GetHistograms("JECup", hist);
    m_JECdown[hist] = GetHistograms("JECdown", hist);
    m_CORup[hist] = GetHistograms("CORup", hist);
    m_CORdown[hist] = GetHistograms("CORdown", hist);
  }

  if(debug) cout << "\t ... Get Background" << endl;
  MapHHH m_bkg = m_other;
  for(MapHHH map: {m_wjets, m_st}) AddTwoMaps(m_bkg, map, 1);

  AddTwoMaps(m_data, m_bkg, -1); // Subtract Bkg from data


  if(debug) cout << "\t ... Rebin" << endl;
  MapHHH m_data_rebin = RebinMap(m_data, bin_width);
  MapHHH m_ttbar_rebin = RebinMap(m_ttbar, bin_width);
  MapHHH m_JECup_rebin = RebinMap(m_JECup, bin_width);
  MapHHH m_JECdown_rebin = RebinMap(m_JECdown, bin_width);
  MapHHH m_CORup_rebin = RebinMap(m_CORup, bin_width);
  MapHHH m_CORdown_rebin = RebinMap(m_CORdown, bin_width);

  ////////////////////////////////////////////////////////////////////////////////
  // Save covariance Matrix in root file. Otherwise it takes to long for debugging
  // In JEC_SYS

  TString file = "CovMatrixJMS_"+channel+"_"+year+".root";
  TFile* f_out = new TFile(file,"update");

  if(debug) cout << "\t ... Get Cov Matrix" << endl;
  MapMMM m_cov_norm_data = GetCovMatrixMap(f_out, m_data_rebin, save, "data");
  MapMMM m_cov_norm_ttbar = GetCovMatrixMap(f_out, m_ttbar_rebin, save, "ttbar");
  MapMMM m_cov_norm_JECup = GetCovMatrixMap(f_out, m_JECup_rebin, save, "JECup");
  MapMMM m_cov_norm_JECdown = GetCovMatrixMap(f_out, m_JECdown_rebin, save, "JECdown");
  MapMMM m_cov_norm_CORup = GetCovMatrixMap(f_out, m_CORup_rebin, save, "CORup");
  MapMMM m_cov_norm_CORdown = GetCovMatrixMap(f_out, m_CORdown_rebin, save, "CORdown");

  f_out->cd();
  f_out->Close();

  return 0;
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// Define functions

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
MapHH GetHistograms(TString process, TString h_name){
  VecTS collection_muon, collection_elec;
  if(process.EqualTo("data")){ collection_muon = data_muon; collection_elec = data_elec;}
  else if(process.EqualTo("ttbar")){ collection_muon = ttbar_muon; collection_elec = ttbar_elec;}
  else if(process.EqualTo("wjets")){ collection_muon = wjets_muon; collection_elec = wjets_elec;}
  else if(process.EqualTo("st")){ collection_muon = st_muon; collection_elec = st_elec;}
  else if(process.EqualTo("other")){ collection_muon = other_muon; collection_elec = other_elec;}
  else if(process.EqualTo("JECup")){ collection_muon = jec_up_muon; collection_elec = jec_up_elec;}
  else if(process.EqualTo("JECdown")){ collection_muon = jec_down_muon; collection_elec = jec_down_elec;}
  else if(process.EqualTo("CORup")){ collection_muon = cor_up_muon; collection_elec = cor_up_elec;}
  else if(process.EqualTo("CORdown")){ collection_muon = cor_down_muon; collection_elec = cor_down_elec;}
  else throw runtime_error("Check the process ´"+process+"´ to obtain the histograms");

  MapHH map;
  vector<TH1F*> h_muon = get_all_hists(collection_muon, h_name);
  vector<TH1F*> h_elec = get_all_hists(collection_elec, h_name);
  vector<TH1F*> h_combine = combine_channels(h_muon, h_elec);

  for(auto nhist:h_muon) nhist->SetTitle(h_name);
  for(auto nhist:h_elec) nhist->SetTitle(h_name);
  for(auto nhist:h_combine) nhist->SetTitle(h_name);

  map["muon"]["2016"] = h_muon[0];
  map["muon"]["2017"] = h_muon[1];
  map["muon"]["2018"] = h_muon[2];
  map["muon"]["combine"] = h_muon[3];

  map["elec"]["2016"] = h_elec[0];
  map["elec"]["2017"] = h_elec[1];
  map["elec"]["2018"] = h_elec[2];
  map["elec"]["combine"] = h_elec[3];

  map["combine"]["2016"] = h_combine[0];
  map["combine"]["2017"] = h_combine[1];
  map["combine"]["2018"] = h_combine[2];
  map["combine"]["combine"] = h_combine[3];

  return map;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================
void AddTwoMaps(MapHHH map1, MapHHH map2, int option){
  // MapHHH map;
  for(auto hist: map1){
    TString h = hist.first;
    for(auto channel: hist.second){
      TString c = channel.first;
      for(auto year: channel.second){
        TString y = year.first;
        // map[h][c][y] = AddHists(year.second, bkg[h][c][y], 1);
        year.second->Add(map2[h][c][y], option);
      }
    }
  }
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapMMM GetCovMatrixMap(TFile* f_out, MapHHH map, TString save, TString process){
  MapMMM covs;
  for(auto hist: map){
    auto start = high_resolution_clock::now(); // Calculation time - start

    TString h = hist.first; TString c = channel; TString y = year;
    covs[h]; covs[h][c];
    TMatrixD mtemp = GetCovMatrix(map[h][c][y], debug);
    covs[h][c][y] = &mtemp;

    TH2D* htemp = TMatrixDtoTH2D(mtemp, mtemp.GetNcols(), 0, 180, h+c+y+process);
    TString wbin = h; wbin.ReplaceAll(cut, "");
    DrawCov(htemp, save+process+"_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    TMatrixD norm = NormCovMatrix(map[h][c][y], mtemp, false);
    TH2D* htemp_norm = TMatrixDtoTH2D(norm, norm.GetNcols(), 0, 180, h+c+y+process+"norm");
    DrawCov(htemp_norm, save+process+"_norm_"+wbin+"_"+c+"_"+y, "#it{m}_{W}", 0.3);

    auto stop = high_resolution_clock::now();  // Calculation time - stop
    auto duration = duration_cast<seconds>(stop - start);
    cout << process << " " << h << " took " << duration.count() << "s" << endl;

    f_out->cd();
    htemp->Write(process+"_"+wbin, TObject::kOverwrite);
    htemp_norm->Write(process+"_norm_"+wbin, TObject::kOverwrite);

  }
  return covs;
}

// ======================================================================================================
// ===                                                                                                ===
// ======================================================================================================

MapHHH RebinMap(MapHHH process, int width){
  MapHHH map;
  for(auto hist: process){
    TString h = hist.first; TString c = channel; TString y = year;
    map[h]; map[h][c];
    map[h][c][y] = (TH1F*) process[h][c][y]->Rebin(width);
  }
  return map;
}
