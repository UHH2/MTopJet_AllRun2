#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/tdrstyle_all.h"

using namespace std;

TH1F* GetSF(TH1F* h_without, TH1F* h_with);
TH1F* CombineMCs(MapHHH map, TString c, TString s);
TH1F* CombineChannels(MapHHH map, TString s);
void Plot2D(VecH h_wo, VecH h_w, VecH h_check, double xmin, double xmax, double offset, double offset2, TString process, TString hist, bool isNorm);
MapHHH CollectHists(MapFF files, TString hist);
MapHHH NormalizeMap(MapHHH map);
MapHHH GetRatioMap(MapHHH map);

TString indir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/";
TString uhh2 = "/uhh2.AnalysisModuleRunner.MC.";
TString outdir = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/BTagReshape/";
TString btag_w = "BTagReshape_nom/";
TString btag_alt = "BTagReshape_alt/";
TString btag_wo = "BTagReshape_no/";
TString btag_btag = "ttbar/";

TString study = "njets";
TString year;

VecTS process = {"TTbar", "SingleTop", "WJets", "other", "combine"};

MapFF files;
MapHHH h_njets, h_btag;

bool is16, is17, is18;

int main(int argc, char* argv[]){

  // =====================================================================================
  // === Preperations                                                                  ===
  // =====================================================================================

  year = argv[1];
  is16 = (year.EqualTo("2016"))?true:false;
  is17 = (year.EqualTo("2017"))?true:false;
  is18 = (year.EqualTo("2018"))?true:false;

  cout << "Get Files ..." << endl;
  channels = {"muon", "elec"};
  TString version = is16?"2016v3.root":is17?"2017v2.root":is18?"2018.root":"NULL.root";

  for(auto c: channels){
    for(TString p: process){
      if(p.EqualTo("combine")) continue;
      files[c][p] = new TFile(indir+c+uhh2+p+"_"+version);
    }
  }

  TString hist = "Jet/number";
  h_njets = CollectHists(files, hist);

  hist = "XCone/h_deepjet";
  h_btag = CollectHists(files, hist);

  // TH1F *h_wo, *h_w, *h_check;
  cout << "Fill Hists ..." << endl;
  MapHHH place = h_njets; TH1F* h_place;
  h_njets["muon"]["combine"]["without"] = CombineMCs(place, "muon", "without"); h_njets["muon"]["combine"]["with"] = CombineMCs(place, "muon", "with"); h_njets["muon"]["combine"]["check"] = CombineMCs(place, "muon", "check");
  h_njets["elec"]["combine"]["without"] = CombineMCs(place, "elec", "without"); h_njets["elec"]["combine"]["with"] = CombineMCs(place, "elec", "with"); h_njets["elec"]["combine"]["check"] = CombineMCs(place, "elec", "check");
  h_njets["combine"]["combine"]["without"] = (TH1F*) h_njets["muon"]["combine"]["without"]->Clone(); h_njets["combine"]["combine"]["without"]->Add(h_njets["elec"]["combine"]["without"], 1);
  h_njets["combine"]["combine"]["with"] = (TH1F*) h_njets["muon"]["combine"]["with"]->Clone(); h_njets["combine"]["combine"]["with"]->Add(h_njets["elec"]["combine"]["with"], 1);
  h_njets["combine"]["combine"]["check"] = (TH1F*) h_njets["muon"]["combine"]["check"]->Clone(); h_njets["combine"]["combine"]["check"]->Add(h_njets["elec"]["combine"]["check"], 1);

  place = h_btag;
  h_btag["muon"]["combine"]["without"] = CombineMCs(place, "muon", "without"); h_btag["muon"]["combine"]["with"] = CombineMCs(place, "muon", "with"); h_btag["muon"]["combine"]["check"] = CombineMCs(place, "muon", "check");
  h_btag["elec"]["combine"]["without"] = CombineMCs(place, "elec", "without"); h_btag["elec"]["combine"]["with"] = CombineMCs(place, "elec", "with"); h_btag["elec"]["combine"]["check"] = CombineMCs(place, "elec", "check");


  cout << "\t ... normalize hists" << endl;
  MapHHH h_njets_norm = NormalizeMap(h_njets);
  MapHHH h_btag_norm = NormalizeMap(h_btag);

  MapHHH h_njets_ratio = GetRatioMap(h_njets);
  MapHHH h_btag_ratio = GetRatioMap(h_btag);

  MapHHH h_njets_ratio_norm = GetRatioMap(h_njets_norm);
  MapHHH h_btag_ratio_norm = GetRatioMap(h_btag_norm);

  // =====================================================================================
  // === Draw Hists                                                                    ===
  // =====================================================================================

  // void Plot2D(TH1F* h_wo, TH1F* h_w, double xmin, double xmax, double offset, TString process, TString hist, bool isNorm)
  Plot2D({h_njets["muon"]["combine"]["without"],h_njets_ratio["muon"]["combine"]["without"]}, {h_njets["muon"]["combine"]["with"],h_njets_ratio["muon"]["combine"]["with"]}, {h_njets["muon"]["combine"]["check"],h_njets_ratio["muon"]["combine"]["check"]}, 0, 10, -0.05, 1.05, "combine", "njets", false);
  Plot2D({h_njets_norm["muon"]["combine"]["without"],h_njets_ratio_norm["muon"]["combine"]["without"]}, {h_njets_norm["muon"]["combine"]["with"],h_njets_ratio_norm["muon"]["combine"]["with"]}, {h_njets_norm["muon"]["combine"]["check"],h_njets_ratio["muon"]["combine"]["check"]}, 0, 10, -0.05, 1.05, "combine", "njets", true);
  Plot2D({h_njets_norm["muon"]["TTbar"]["without"],h_njets_ratio_norm["muon"]["TTbar"]["without"]}, {h_njets_norm["muon"]["TTbar"]["with"],h_njets_ratio_norm["muon"]["TTbar"]["with"]}, {h_njets_norm["muon"]["TTbar"]["check"],h_njets_ratio["muon"]["TTbar"]["check"]}, 0, 10, -0.05, 1.05, "ttbar", "njets", true);

  Plot2D({h_btag["muon"]["combine"]["without"],h_btag_ratio["muon"]["combine"]["without"]}, {h_btag["muon"]["combine"]["with"],h_btag_ratio["muon"]["combine"]["with"]}, {h_btag["muon"]["combine"]["check"],h_btag_ratio["muon"]["combine"]["check"]}, 0, 1, -0.05, 1.05, "combine", "btag", false);
  Plot2D({h_btag_norm["muon"]["combine"]["without"],h_btag_ratio_norm["muon"]["combine"]["without"]}, {h_btag_norm["muon"]["combine"]["with"],h_btag_ratio_norm["muon"]["combine"]["with"]}, {h_btag_norm["muon"]["combine"]["check"],h_btag_ratio["muon"]["combine"]["check"]}, 0, 1, -0.05, 1.05, "combine", "btag", true);

  // =====================================================================================
  // === Creat ROOT files with SF                                                      ===
  // =====================================================================================
  cout << "Get SFs ..." << endl;
  MapHH m_SFs;
  for(auto c: {"muon", "elec"}){
    for(auto p: process){
      m_SFs[c][p] = GetSF(h_njets_norm[c][p]["without"], h_njets_norm[c][p]["with"]);
    }
  }

  cout << "\t ... Creating file" << endl;
  TFile * f_out = new TFile("files/BTagReshapeWeights"+year+".root","RECREATE");;
  m_SFs["muon"]["combine"]->Write("weight_muon");
  m_SFs["elec"]["combine"]->Write("weight_elec");
  f_out->Close();
}

// =====================================================================================
// =====================================================================================
// =====================================================================================

// =====================================================================================
void Plot2D(VecH h_wo, VecH h_w, VecH h_check, double xmin, double xmax, double offset, double offset2, TString process, TString hist, bool isNorm){
  TLegend* leg = tdrLeg(0.58, 0.7, 0.9, 0.9, 0.04);
  TCanvas* canv;

  gErrorIgnoreLevel = kError;
  TString xaxis = (hist.EqualTo("njets"))?"n_{jets}":hist;
  TString yaxis = (isNorm)?"a.u.":"Events";
  int lumi_period = is16?1:is17?2:is18?3:4;
  cout << lumi_period << endl;
  canv = tdrDiCanvas2("btagreshape", xmin, xmax, 0.01, h_w[0]->GetMaximum()*1.3, 0.94, 1.06, xaxis, yaxis, "#frac{without}{with}", kRectangular, lumi_period, 11, offset, offset2);
  canv->SetTickx(0);
  canv->SetTicky(0);

  canv->cd(1);
  tdrDraw(h_wo[0], "hist same", 1, kBlack, kSolid, kBlue+2, 0);
  tdrDraw(h_w[0], "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(h_check[0], "hist same", 1, kBlack, kDashed, kGreen+2, 0);

  leg->AddEntry(h_wo[0], "Without b-tagging SFs", "l");
  leg->AddEntry(h_w[0], "With b-tagging SFs", "l");
  leg->AddEntry(h_check[0], "Applied weights from n_{jet}", "l");
  leg->Draw("same");

  canv->cd(2);
  tdrDraw(h_w[1], "hist same", 1, kBlack, kSolid, kRed+2, 0);
  tdrDraw(h_wo[1], "hist same", 1, kBlack, kSolid, kBlue+2, 0);
  tdrDraw(h_check[1], "hist same", 1, kBlack, kDashed, kGreen+2, 0);

  canv->Update();
  TString save_info = hist+"_"+process+"_"+year+"_"+((isNorm)?"norm":"");
  canv->Print(outdir+"btagreshape_"+save_info+".pdf","pdf");
  leg->Clear();
  delete canv;
}


// =====================================================================================
TH1F* GetSF(TH1F* h_without, TH1F* h_with){
  TH1F* h_central = (TH1F*) h_without->Clone();
  int nbins = h_central->GetNbinsX();
  for(unsigned int i=1; i<=nbins; i++){

    double value_wo = h_without->GetBinContent(i);
    double value_w = h_with->GetBinContent(i);

    double central = 1;
    if(value_wo != 0 && value_w != 0) central = value_wo/value_w;
    // cout << "bin " << i << "\t" << central << "\t" << value_wo << "\t" << value_w << endl;

    h_central->SetBinContent(i, central);
    h_central->SetBinError(i, 0.000000000000001);
  }
  return h_central;
}


// =====================================================================================
MapHHH CollectHists(MapFF files, TString hist){
  cout << "Collect hists " << hist << endl;
  MapHHH map;
  for(auto c: channels){
    map[c];
    for(auto p: process){
      map[c][p]; // initilize
      if(p.EqualTo("combine")) continue;
      map[c][p]["without"] = (TH1F*) files[c][p]->Get("PreSel03_"+hist);
      map[c][p]["with"] = (TH1F*) files[c][p]->Get("PreSel03b_"+hist);
      map[c][p]["check"] = (TH1F*) files[c][p]->Get("PreSelReshapeSF_"+hist);
    }
  }
  return map;
}

// =====================================================================================
TH1F* CombineMCs(MapHHH map, TString c, TString s){
  TH1F* h_s = (TH1F*) map[c]["TTbar"][s]->Clone();
  for(TString p: process){
    if(p.Contains("combine") || p.Contains("TTbar")) continue;
    h_s->Add(map[c][p][s], 1);
  }
  return h_s;
}

// =====================================================================================
TH1F* CombineChannels(MapHHH map, TString s){
  TH1F* h_s_m = (TH1F*) map["muon"]["TTbar"][s]->Clone();
  TH1F* h_s_e = (TH1F*) map["elec"]["TTbar"][s]->Clone();
  for(TString p: process){
    if(p.Contains("combine") || p.Contains("TTbar")) continue;
    h_s_m->Add(map["muon"][p][s], 1);
    h_s_e->Add(map["elec"][p][s], 1);
  }
  h_s_m->Add(h_s_e, 1);
  return h_s_m;
}

// =====================================================================================
MapHHH NormalizeMap(MapHHH map){
  MapHHH h_norm;
  for(auto c: channels){
    for(auto p: process){
      for(auto s: {"without", "with", "check"}){
        h_norm[c][p][s] = Normalize(map[c][p][s]);
      }
    }
  }
  return h_norm;
}

// =====================================================================================
MapHHH GetRatioMap(MapHHH map){
  MapHHH ratio;
  for(auto c: channels){
    for(auto p: process){
      ratio[c][p]["with"] = GetRatio(map[c][p]["with"], map[c][p]["with"], true);
      ratio[c][p]["without"] = GetRatio(map[c][p]["without"], map[c][p]["with"], false);
      ratio[c][p]["check"] = GetRatio(map[c][p]["check"], map[c][p]["with"], false);
    }
  }
  return ratio;
}
