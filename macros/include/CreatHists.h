#include "CentralInclude.h"
#include "HistogramUtils.h"

using namespace std;

// ----------------------------------------------------------------------------- Functions

TH1F* combine_years(vector<TH1F*> hists){
  TH1F* new_hist = AddHists(hists, 1);
  return new_hist;
}

vector<TH1F*> combine_channels(vector<TH1F*> h_muon, vector<TH1F*> h_elec){
  vector<TH1F*> new_hists = AddHists(h_muon, h_elec, 1);
  return new_hists;
}

// Creating Hists
vector<TH1F*> get_all_hists(vector<TString> names, TString hist){
  vector<TFile*> files;
  vector<TH1F*> hists;
  for(TString file: names){files.push_back(new TFile(file));}
  for(TFile* file: files){hists.push_back((TH1F*)file->Get(hist));}
  TH1F* h_combined = combine_years(hists);
  hists.push_back(h_combined);
  return hists;
}

TH1F* get_hist(TString name, TString hist){
    TFile* f = new TFile(name);
    TH1F* h = (TH1F*) f->Get(hist);
    return h;
}

// TH1F* get_single_hists(vector<TString> names, TString hist_name, int index){
//   TFile* file = new TFile(names[index]);
//   TH1F* hist = (TH1F*)file->Get(hist_name);
//   return hist;
// }

// ----------------------------------------------------------------------------- Nominal
// TString dir = get_save_path();

// =====================================================================================
// === Constructor                                                                   ===
// =====================================================================================
// In macros only use get_hist() and construct the histrogram names with the given TStrigns
// Old definitions are kept to not rewrite older macros.

// channels
TString muon = "muon";
TString elec = "elec";
TString muon_elecJMS = "muon_elecJMS";
TString elec_muonJMS = "elec_muonJMS";

// Years
TString y16      = "_2016v3";  int index_16  = 0;
TString y17      = "_2017v2";  int index_17  = 1;
TString y18      = "_2018";    int index_18  = 2;
TString yCombine = "_combine"; int index_Com = 3;

VecTS years = {"2016", "2017", "2018", "combine"};
VecTS channels = {"muon", "elec", "combine"};

// Hist files
TString data_f   = "uhh2.AnalysisModuleRunner.DATA.DATA";
TString ttbar_f  = "uhh2.AnalysisModuleRunner.MC.TTbar";
TString wjets_f  = "uhh2.AnalysisModuleRunner.MC.WJets";
TString st_f     = "uhh2.AnalysisModuleRunner.MC.SingleTop";
TString other_f  = "uhh2.AnalysisModuleRunner.MC.other";
TString uhh_mc   = "uhh2.AnalysisModuleRunner.MC.";

// systematics
TString fsr_up_16 = "uhh2.AnalysisModuleRunner.MC.TTbar_fsrup_2016v3.root";
TString fsr_down_16 = "uhh2.AnalysisModuleRunner.MC.TTbar_fsrdown_2016v3.root";

TString fsr_up_sqrt2 = "FSRup_sqrt2";
TString fsr_up_2 = "FSRup_2";
TString fsr_up_4 = "FSRup_4";
TString fsr_down_sqrt2 = "FSRdown_sqrt2";
TString fsr_down_2 = "FSRdown_2";
TString fsr_down_4 = "FSRdown_4";
TString jms_uu   = "JMS_upup";
TString jms_ud   = "JMS_updown";
TString jms_du   = "JMS_downup";
TString jms_dd   = "JMS_downdown";
TString jms_u    = "JMS_up";
TString jms_d    = "JMS_down";
TString jms_u_corr = "JMS_up_acorr";
TString jms_d_corr = "JMS_down_acorr";
TString jec_up   = "JEC_up";
TString jec_down = "JEC_down";
TString cor_up   = "COR_up";
TString cor_down = "COR_down";

// =====================================================================================
// === Commonly used                                                                 ===
// =====================================================================================

// Hist paths
TString data_muon_16 = dir+muon+"/"+data_f+y16+".root";
TString data_muon_17 = dir+muon+"/"+data_f+y17+".root";
TString data_muon_18 = dir+muon+"/"+data_f+y18+".root";
TString data_elec_16 = dir+elec+"/"+data_f+y16+".root";
TString data_elec_17 = dir+elec+"/"+data_f+y17+".root";
TString data_elec_18 = dir+elec+"/"+data_f+y18+".root";
vector<TString> data_muon = {data_muon_16, data_muon_17, data_muon_18};
vector<TString> data_elec = {data_elec_16, data_elec_17, data_elec_18};

TString ttbar_muon_16 = dir+muon+"/"+ttbar_f+y16+".root";
TString ttbar_muon_17 = dir+muon+"/"+ttbar_f+y17+".root";
TString ttbar_muon_18 = dir+muon+"/"+ttbar_f+y18+".root";
TString ttbar_elec_16 = dir+elec+"/"+ttbar_f+y16+".root";
TString ttbar_elec_17 = dir+elec+"/"+ttbar_f+y17+".root";
TString ttbar_elec_18 = dir+elec+"/"+ttbar_f+y18+".root";
vector<TString> ttbar_muon = {ttbar_muon_16, ttbar_muon_17, ttbar_muon_18};
vector<TString> ttbar_elec = {ttbar_elec_16, ttbar_elec_17, ttbar_elec_18};

TString wjets_muon_16 = dir+muon+"/"+wjets_f+y16+".root";
TString wjets_muon_17 = dir+muon+"/"+wjets_f+y17+".root";
TString wjets_muon_18 = dir+muon+"/"+wjets_f+y18+".root";
TString wjets_elec_16 = dir+elec+"/"+wjets_f+y16+".root";
TString wjets_elec_17 = dir+elec+"/"+wjets_f+y17+".root";
TString wjets_elec_18 = dir+elec+"/"+wjets_f+y18+".root";
vector<TString> wjets_muon = {wjets_muon_16, wjets_muon_17, wjets_muon_18};
vector<TString> wjets_elec = {wjets_elec_16, wjets_elec_17, wjets_elec_18};

TString st_muon_16 = dir+muon+"/"+st_f+y16+".root";
TString st_muon_17 = dir+muon+"/"+st_f+y17+".root";
TString st_muon_18 = dir+muon+"/"+st_f+y18+".root";
TString st_elec_16 = dir+elec+"/"+st_f+y16+".root";
TString st_elec_17 = dir+elec+"/"+st_f+y17+".root";
TString st_elec_18 = dir+elec+"/"+st_f+y18+".root";
vector<TString> st_muon = {st_muon_16, st_muon_17, st_muon_18};
vector<TString> st_elec = {st_elec_16, st_elec_17, st_elec_18};

TString other_muon_16 = dir+muon+"/"+other_f+y16+".root";
TString other_muon_17 = dir+muon+"/"+other_f+y17+".root";
TString other_muon_18 = dir+muon+"/"+other_f+y18+".root";
TString other_elec_16 = dir+elec+"/"+other_f+y16+".root";
TString other_elec_17 = dir+elec+"/"+other_f+y17+".root";
TString other_elec_18 = dir+elec+"/"+other_f+y18+".root";
vector<TString> other_muon = {other_muon_16, other_muon_17, other_muon_18};
vector<TString> other_elec = {other_elec_16, other_elec_17, other_elec_18};

// ----------------------------------------------------------------------------- mTop

TString mtop_muon_1665_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1665_2016v3.root";
TString mtop_muon_1665_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1665_2017v2.root";
TString mtop_muon_1665_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1665_2018.root";
TString mtop_elec_1665_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1665_2016v3.root";
TString mtop_elec_1665_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1665_2017v2.root";
TString mtop_elec_1665_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1665_2018.root";
vector<TString> mtop_1665_muon = {mtop_muon_1665_16, mtop_muon_1665_17, mtop_muon_1665_18};
vector<TString> mtop_1665_elec = {mtop_elec_1665_16, mtop_elec_1665_17, mtop_elec_1665_18};

TString mtop_muon_1695_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1695_2016v3.root";
TString mtop_muon_1695_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1695_2017v2.root";
TString mtop_muon_1695_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1695_2018.root";
TString mtop_elec_1695_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1695_2016v3.root";
TString mtop_elec_1695_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1695_2017v2.root";
TString mtop_elec_1695_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1695_2018.root";
vector<TString> mtop_1695_muon = {mtop_muon_1695_16, mtop_muon_1695_17, mtop_muon_1695_18};
vector<TString> mtop_1695_elec = {mtop_elec_1695_16, mtop_elec_1695_17, mtop_elec_1695_18};

TString mtop_muon_1715_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1715_2016v3.root";
TString mtop_muon_1715_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1715_2017v2.root";
TString mtop_muon_1715_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1715_2018.root";
TString mtop_elec_1715_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1715_2016v3.root";
TString mtop_elec_1715_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1715_2017v2.root";
TString mtop_elec_1715_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1715_2018.root";
vector<TString> mtop_1715_muon = {mtop_muon_1715_16, mtop_muon_1715_17, mtop_muon_1715_18};
vector<TString> mtop_1715_elec = {mtop_elec_1715_16, mtop_elec_1715_17, mtop_elec_1715_18};

TString mtop_muon_1735_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1735_2016v3.root";
TString mtop_muon_1735_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1735_2017v2.root";
TString mtop_muon_1735_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1735_2018.root";
TString mtop_elec_1735_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1735_2016v3.root";
TString mtop_elec_1735_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1735_2017v2.root";
TString mtop_elec_1735_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1735_2018.root";
vector<TString> mtop_1735_muon = {mtop_muon_1735_16, mtop_muon_1735_17, mtop_muon_1735_18};
vector<TString> mtop_1735_elec = {mtop_elec_1735_16, mtop_elec_1735_17, mtop_elec_1735_18};

TString mtop_muon_1755_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1755_2016v3.root";
TString mtop_muon_1755_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1755_2017v2.root";
TString mtop_muon_1755_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1755_2018.root";
TString mtop_elec_1755_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1755_2016v3.root";
TString mtop_elec_1755_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1755_2017v2.root";
TString mtop_elec_1755_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1755_2018.root";
vector<TString> mtop_1755_muon = {mtop_muon_1755_16, mtop_muon_1755_17, mtop_muon_1755_18};
vector<TString> mtop_1755_elec = {mtop_elec_1755_16, mtop_elec_1755_17, mtop_elec_1755_18};

TString mtop_muon_1785_16 = dir+muon+"/"+uhh_mc+"TTbar_mtop1785_2016v3.root";
TString mtop_muon_1785_17 = dir+muon+"/"+uhh_mc+"TTbar_mtop1785_2017v2.root";
TString mtop_muon_1785_18 = dir+muon+"/"+uhh_mc+"TTbar_mtop1785_2018.root";
TString mtop_elec_1785_16 = dir+elec+"/"+uhh_mc+"TTbar_mtop1785_2016v3.root";
TString mtop_elec_1785_17 = dir+elec+"/"+uhh_mc+"TTbar_mtop1785_2017v2.root";
TString mtop_elec_1785_18 = dir+elec+"/"+uhh_mc+"TTbar_mtop1785_2018.root";
vector<TString> mtop_1785_muon = {mtop_muon_1785_16, mtop_muon_1785_17, mtop_muon_1785_18};
vector<TString> mtop_1785_elec = {mtop_elec_1785_16, mtop_elec_1785_17, mtop_elec_1785_18};

// ----------------------------------------------------------------------------- FSR

TString fsr_up_muon_16 = dir+muon+"/"+fsr_up_16;
TString fsr_up_elec_16 = dir+elec+"/"+fsr_up_16;
TString fsr_down_muon_16 = dir+muon+"/"+fsr_down_16;
TString fsr_down_elec_16 = dir+elec+"/"+fsr_down_16;

TString fsr_up_sqrt2_muon_17 = dir+muon+"/"+fsr_up_sqrt2+"/"+ttbar_f+y17+".root";
TString fsr_up_sqrt2_muon_18 = dir+muon+"/"+fsr_up_sqrt2+"/"+ttbar_f+y18+".root";
TString fsr_up_sqrt2_elec_17 = dir+elec+"/"+fsr_up_sqrt2+"/"+ttbar_f+y17+".root";
TString fsr_up_sqrt2_elec_18 = dir+elec+"/"+fsr_up_sqrt2+"/"+ttbar_f+y18+".root";
vector<TString> fsr_up_sqrt2_muon = {fsr_up_sqrt2_muon_17, fsr_up_sqrt2_muon_18};
vector<TString> fsr_up_sqrt2_elec = {fsr_up_sqrt2_elec_17, fsr_up_sqrt2_elec_18};

TString fsr_down_sqrt2_muon_17 = dir+muon+"/"+fsr_down_sqrt2+"/"+ttbar_f+y17+".root";
TString fsr_down_sqrt2_muon_18 = dir+muon+"/"+fsr_down_sqrt2+"/"+ttbar_f+y18+".root";
TString fsr_down_sqrt2_elec_17 = dir+elec+"/"+fsr_down_sqrt2+"/"+ttbar_f+y17+".root";
TString fsr_down_sqrt2_elec_18 = dir+elec+"/"+fsr_down_sqrt2+"/"+ttbar_f+y18+".root";
vector<TString> fsr_down_sqrt2_muon = {fsr_down_sqrt2_muon_17, fsr_down_sqrt2_muon_18};
vector<TString> fsr_down_sqrt2_elec = {fsr_down_sqrt2_elec_17, fsr_down_sqrt2_elec_18};

TString fsr_up_2_muon_17 = dir+muon+"/"+fsr_up_2+"/"+ttbar_f+y17+".root";
TString fsr_up_2_muon_18 = dir+muon+"/"+fsr_up_2+"/"+ttbar_f+y18+".root";
TString fsr_up_2_elec_17 = dir+elec+"/"+fsr_up_2+"/"+ttbar_f+y17+".root";
TString fsr_up_2_elec_18 = dir+elec+"/"+fsr_up_2+"/"+ttbar_f+y18+".root";
vector<TString> fsr_up_2_muon = {fsr_up_2_muon_17, fsr_up_2_muon_18};
vector<TString> fsr_up_2_elec = {fsr_up_2_elec_17, fsr_up_2_elec_18};

TString fsr_down_2_muon_17 = dir+muon+"/"+fsr_down_2+"/"+ttbar_f+y17+".root";
TString fsr_down_2_muon_18 = dir+muon+"/"+fsr_down_2+"/"+ttbar_f+y18+".root";
TString fsr_down_2_elec_17 = dir+elec+"/"+fsr_down_2+"/"+ttbar_f+y17+".root";
TString fsr_down_2_elec_18 = dir+elec+"/"+fsr_down_2+"/"+ttbar_f+y18+".root";
vector<TString> fsr_down_2_muon = {fsr_down_2_muon_17, fsr_down_2_muon_18};
vector<TString> fsr_down_2_elec = {fsr_down_2_elec_17, fsr_down_2_elec_18};

TString fsr_up_4_muon_17 = dir+muon+"/"+fsr_up_4+"/"+ttbar_f+y17+".root";
TString fsr_up_4_muon_18 = dir+muon+"/"+fsr_up_4+"/"+ttbar_f+y18+".root";
TString fsr_up_4_elec_17 = dir+elec+"/"+fsr_up_4+"/"+ttbar_f+y17+".root";
TString fsr_up_4_elec_18 = dir+elec+"/"+fsr_up_4+"/"+ttbar_f+y18+".root";
vector<TString> fsr_up_4_muon = {fsr_up_4_muon_17, fsr_up_4_muon_18};
vector<TString> fsr_up_4_elec = {fsr_up_4_elec_17, fsr_up_4_elec_18};

TString fsr_down_4_muon_17 = dir+muon+"/"+fsr_down_4+"/"+ttbar_f+y17+".root";
TString fsr_down_4_muon_18 = dir+muon+"/"+fsr_down_4+"/"+ttbar_f+y18+".root";
TString fsr_down_4_elec_17 = dir+elec+"/"+fsr_down_4+"/"+ttbar_f+y17+".root";
TString fsr_down_4_elec_18 = dir+elec+"/"+fsr_down_4+"/"+ttbar_f+y18+".root";
vector<TString> fsr_down_4_muon = {fsr_down_4_muon_17, fsr_down_4_muon_18};
vector<TString> fsr_down_4_elec = {fsr_down_4_elec_17, fsr_down_4_elec_18};

// ----------------------------------------------------------------------------- JMS

// Default
TString jms_uu_muon_16 = dir+muon+"/"+jms_uu+"/"+ttbar_f+y16+".root";
TString jms_uu_muon_17 = dir+muon+"/"+jms_uu+"/"+ttbar_f+y17+".root";
TString jms_uu_muon_18 = dir+muon+"/"+jms_uu+"/"+ttbar_f+y18+".root";
TString jms_uu_elec_16 = dir+elec+"/"+jms_uu+"/"+ttbar_f+y16+".root";
TString jms_uu_elec_17 = dir+elec+"/"+jms_uu+"/"+ttbar_f+y17+".root";
TString jms_uu_elec_18 = dir+elec+"/"+jms_uu+"/"+ttbar_f+y18+".root";
vector<TString> jms_uu_muon = {jms_uu_muon_16, jms_uu_muon_17, jms_uu_muon_18};
vector<TString> jms_uu_elec = {jms_uu_elec_16, jms_uu_elec_17, jms_uu_elec_18};

TString jms_ud_muon_16 = dir+muon+"/"+jms_ud+"/"+ttbar_f+y16+".root";
TString jms_ud_muon_17 = dir+muon+"/"+jms_ud+"/"+ttbar_f+y17+".root";
TString jms_ud_muon_18 = dir+muon+"/"+jms_ud+"/"+ttbar_f+y18+".root";
TString jms_ud_elec_16 = dir+elec+"/"+jms_ud+"/"+ttbar_f+y16+".root";
TString jms_ud_elec_17 = dir+elec+"/"+jms_ud+"/"+ttbar_f+y17+".root";
TString jms_ud_elec_18 = dir+elec+"/"+jms_ud+"/"+ttbar_f+y18+".root";
vector<TString> jms_ud_muon = {jms_ud_muon_16, jms_ud_muon_17, jms_ud_muon_18};
vector<TString> jms_ud_elec = {jms_ud_elec_16, jms_ud_elec_17, jms_ud_elec_18};

TString jms_du_muon_16 = dir+muon+"/"+jms_du+"/"+ttbar_f+y16+".root";
TString jms_du_muon_17 = dir+muon+"/"+jms_du+"/"+ttbar_f+y17+".root";
TString jms_du_muon_18 = dir+muon+"/"+jms_du+"/"+ttbar_f+y18+".root";
TString jms_du_elec_16 = dir+elec+"/"+jms_du+"/"+ttbar_f+y16+".root";
TString jms_du_elec_17 = dir+elec+"/"+jms_du+"/"+ttbar_f+y17+".root";
TString jms_du_elec_18 = dir+elec+"/"+jms_du+"/"+ttbar_f+y18+".root";
vector<TString> jms_du_muon = {jms_du_muon_16, jms_du_muon_17, jms_du_muon_18};
vector<TString> jms_du_elec = {jms_du_elec_16, jms_du_elec_17, jms_du_elec_18};

TString jms_dd_muon_16 = dir+muon+"/"+jms_dd+"/"+ttbar_f+y16+".root";
TString jms_dd_muon_17 = dir+muon+"/"+jms_dd+"/"+ttbar_f+y17+".root";
TString jms_dd_muon_18 = dir+muon+"/"+jms_dd+"/"+ttbar_f+y18+".root";
TString jms_dd_elec_16 = dir+elec+"/"+jms_dd+"/"+ttbar_f+y16+".root";
TString jms_dd_elec_17 = dir+elec+"/"+jms_dd+"/"+ttbar_f+y17+".root";
TString jms_dd_elec_18 = dir+elec+"/"+jms_dd+"/"+ttbar_f+y18+".root";
vector<TString> jms_dd_muon = {jms_dd_muon_16, jms_dd_muon_17, jms_dd_muon_18};
vector<TString> jms_dd_elec = {jms_dd_elec_16, jms_dd_elec_17, jms_dd_elec_18};

// Propagation of Uncertainty; rho = -0.21
TString jms_u_muon_16 = dir+muon+"/"+jms_u+"/"+ttbar_f+y16+".root";
TString jms_u_muon_17 = dir+muon+"/"+jms_u+"/"+ttbar_f+y17+".root";
TString jms_u_muon_18 = dir+muon+"/"+jms_u+"/"+ttbar_f+y18+".root";
TString jms_u_elec_16 = dir+elec+"/"+jms_u+"/"+ttbar_f+y16+".root";
TString jms_u_elec_17 = dir+elec+"/"+jms_u+"/"+ttbar_f+y17+".root";
TString jms_u_elec_18 = dir+elec+"/"+jms_u+"/"+ttbar_f+y18+".root";
vector<TString> jms_u_muon = {jms_u_muon_16, jms_u_muon_17, jms_u_muon_18};
vector<TString> jms_u_elec = {jms_u_elec_16, jms_u_elec_17, jms_u_elec_18};

TString jms_d_muon_16 = dir+muon+"/"+jms_d+"/"+ttbar_f+y16+".root";
TString jms_d_muon_17 = dir+muon+"/"+jms_d+"/"+ttbar_f+y17+".root";
TString jms_d_muon_18 = dir+muon+"/"+jms_d+"/"+ttbar_f+y18+".root";
TString jms_d_elec_16 = dir+elec+"/"+jms_d+"/"+ttbar_f+y16+".root";
TString jms_d_elec_17 = dir+elec+"/"+jms_d+"/"+ttbar_f+y17+".root";
TString jms_d_elec_18 = dir+elec+"/"+jms_d+"/"+ttbar_f+y18+".root";
vector<TString> jms_d_muon = {jms_d_muon_16, jms_d_muon_17, jms_d_muon_18};
vector<TString> jms_d_elec = {jms_d_elec_16, jms_d_elec_17, jms_d_elec_18};

// fully anticorrelated; just for comparison
TString jms_d_muon_16_corr = dir+muon+"/"+jms_d_corr+"/"+ttbar_f+y16+".root";
TString jms_d_muon_17_corr = dir+muon+"/"+jms_d_corr+"/"+ttbar_f+y17+".root";
TString jms_d_muon_18_corr = dir+muon+"/"+jms_d_corr+"/"+ttbar_f+y18+".root";
TString jms_d_elec_16_corr = dir+elec+"/"+jms_d_corr+"/"+ttbar_f+y16+".root";
TString jms_d_elec_17_corr = dir+elec+"/"+jms_d_corr+"/"+ttbar_f+y17+".root";
TString jms_d_elec_18_corr = dir+elec+"/"+jms_d_corr+"/"+ttbar_f+y18+".root";
vector<TString> jms_d_muon_corr = {jms_d_muon_16_corr, jms_d_muon_17_corr, jms_d_muon_18_corr};
vector<TString> jms_d_elec_corr = {jms_d_elec_16_corr, jms_d_elec_17_corr, jms_d_elec_18_corr};

TString jms_u_muon_16_corr = dir+muon+"/"+jms_u_corr+"/"+ttbar_f+y16+".root";
TString jms_u_muon_17_corr = dir+muon+"/"+jms_u_corr+"/"+ttbar_f+y17+".root";
TString jms_u_muon_18_corr = dir+muon+"/"+jms_u_corr+"/"+ttbar_f+y18+".root";
TString jms_u_elec_16_corr = dir+elec+"/"+jms_u_corr+"/"+ttbar_f+y16+".root";
TString jms_u_elec_17_corr = dir+elec+"/"+jms_u_corr+"/"+ttbar_f+y17+".root";
TString jms_u_elec_18_corr = dir+elec+"/"+jms_u_corr+"/"+ttbar_f+y18+".root";
vector<TString> jms_u_muon_corr = {jms_u_muon_16_corr, jms_u_muon_17_corr, jms_u_muon_18_corr};
vector<TString> jms_u_elec_corr = {jms_u_elec_16_corr, jms_u_elec_17_corr, jms_u_elec_18_corr};

// ----------------------------------------------------------------------------- JEC

TString jec_up_muon_16 = dir+muon+"/"+jec_up+"/"+ttbar_f+y16+".root";
TString jec_up_muon_17 = dir+muon+"/"+jec_up+"/"+ttbar_f+y17+".root";
TString jec_up_muon_18 = dir+muon+"/"+jec_up+"/"+ttbar_f+y18+".root";
TString jec_up_elec_16 = dir+elec+"/"+jec_up+"/"+ttbar_f+y16+".root";
TString jec_up_elec_17 = dir+elec+"/"+jec_up+"/"+ttbar_f+y17+".root";
TString jec_up_elec_18 = dir+elec+"/"+jec_up+"/"+ttbar_f+y18+".root";
vector<TString> jec_up_muon = {jec_up_muon_16, jec_up_muon_17, jec_up_muon_18};
vector<TString> jec_up_elec = {jec_up_elec_16, jec_up_elec_17, jec_up_elec_18};

TString jec_down_muon_16 = dir+muon+"/"+jec_down+"/"+ttbar_f+y16+".root";
TString jec_down_muon_17 = dir+muon+"/"+jec_down+"/"+ttbar_f+y17+".root";
TString jec_down_muon_18 = dir+muon+"/"+jec_down+"/"+ttbar_f+y18+".root";
TString jec_down_elec_16 = dir+elec+"/"+jec_down+"/"+ttbar_f+y16+".root";
TString jec_down_elec_17 = dir+elec+"/"+jec_down+"/"+ttbar_f+y17+".root";
TString jec_down_elec_18 = dir+elec+"/"+jec_down+"/"+ttbar_f+y18+".root";
vector<TString> jec_down_muon = {jec_down_muon_16, jec_down_muon_17, jec_down_muon_18};
vector<TString> jec_down_elec = {jec_down_elec_16, jec_down_elec_17, jec_down_elec_18};

// ----------------------------------------------------------------------------- COR

TString cor_up_muon_16 = dir+muon+"/"+cor_up+"/"+ttbar_f+y16+".root";
TString cor_up_muon_17 = dir+muon+"/"+cor_up+"/"+ttbar_f+y17+".root";
TString cor_up_muon_18 = dir+muon+"/"+cor_up+"/"+ttbar_f+y18+".root";
TString cor_up_elec_16 = dir+elec+"/"+cor_up+"/"+ttbar_f+y16+".root";
TString cor_up_elec_17 = dir+elec+"/"+cor_up+"/"+ttbar_f+y17+".root";
TString cor_up_elec_18 = dir+elec+"/"+cor_up+"/"+ttbar_f+y18+".root";
vector<TString> cor_up_muon = {cor_up_muon_16, cor_up_muon_17, cor_up_muon_18};
vector<TString> cor_up_elec = {cor_up_elec_16, cor_up_elec_17, cor_up_elec_18};

TString cor_down_muon_16 = dir+muon+"/"+cor_down+"/"+ttbar_f+y16+".root";
TString cor_down_muon_17 = dir+muon+"/"+cor_down+"/"+ttbar_f+y17+".root";
TString cor_down_muon_18 = dir+muon+"/"+cor_down+"/"+ttbar_f+y18+".root";
TString cor_down_elec_16 = dir+elec+"/"+cor_down+"/"+ttbar_f+y16+".root";
TString cor_down_elec_17 = dir+elec+"/"+cor_down+"/"+ttbar_f+y17+".root";
TString cor_down_elec_18 = dir+elec+"/"+cor_down+"/"+ttbar_f+y18+".root";
vector<TString> cor_down_muon = {cor_down_muon_16, cor_down_muon_17, cor_down_muon_18};
vector<TString> cor_down_elec = {cor_down_elec_16, cor_down_elec_17, cor_down_elec_18};
