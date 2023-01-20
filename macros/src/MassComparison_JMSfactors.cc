#include "../include/CentralInclude.h"
#include "../include/CreatHists.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"

void Plot(const MapHH &hists, TString y);

MapHH m_combine;
MapHHH m_ttbar;

int main(int argc, char* argv[]){
  VecTS channels = {"muon", "elec", "elec_muonJMS", "muon_elecJMS"};
//   VecTS jms = {"upup", "downdown", "nominal"};
//   TString hist = "JetMassScaleHists/hadjet_jms_mass";
//   TString hist_old = "XCone_cor/M_jet1";
//
//   for(TString c: channels){
//       cout << 1 << c << endl;
//       TH1F* ttbar_16 = get_hist(dir+c+"/"+ttbar_f+y16+".root", hist);
//       cout << 1 << c << endl;
//       TH1F* ttbar_16_uu = get_hist(dir+c+"/"+jms_uu+"/"+ttbar_f+y16+".root", hist);
//       cout << 1 << c << endl;
//       TH1F* ttbar_16_dd = get_hist(dir+c+"/"+jms_dd+"/"+ttbar_f+y16+".root", hist);
//       m_ttbar["2016"][c]["nominal"] = ttbar_16;
//       m_ttbar["2016"][c]["upup"] = ttbar_16_uu;
//       m_ttbar["2016"][c]["downdown"] = ttbar_16_dd;
//
//       cout << 2 << c << endl;
//       TH1F* ttbar_17 = get_hist(dir+c+"/"+ttbar_f+y17+".root", hist);
//       cout << 2 << c << endl;
//       TH1F* ttbar_17_uu = get_hist(dir+c+"/"+jms_uu+"/"+ttbar_f+y17+".root", hist);
//       cout << 2 << c << endl;
//       TH1F* ttbar_17_dd = get_hist(dir+c+"/"+jms_dd+"/"+ttbar_f+y17+".root", hist);
//       m_ttbar["2017"][c]["nominal"] = ttbar_17;
//       m_ttbar["2017"][c]["upup"] = ttbar_17_uu;
//       m_ttbar["2017"][c]["downdown"] = ttbar_17_dd;
//
//       cout << 3 << c << endl;
//       TH1F* ttbar_18 = get_hist(dir+c+"/"+ttbar_f+y18+".root", hist);
//       cout << 3 << c << endl;
//       TH1F* ttbar_18_uu = get_hist(dir+c+"/"+jms_uu+"/"+ttbar_f+y18+".root", hist);
//       cout << 3 << c << endl;
//       TH1F* ttbar_18_dd = get_hist(dir+c+"/"+jms_dd+"/"+ttbar_f+y18+".root", hist);
//       m_ttbar["2018"][c]["nominal"] = ttbar_18;
//       m_ttbar["2018"][c]["upup"] = ttbar_18_uu;
//       m_ttbar["2018"][c]["downdown"] = ttbar_18_dd;
//   }
//
//   for(TString y: years){
//     if(y.EqualTo("combine")) continue;
//     cout << 4 << y << endl;
//     m_ttbar[y]["combine"]["nominal"] = AddHists(m_ttbar[y]["muon"]["nominal"], m_ttbar[y]["elec"]["nominal"], 1);
//     cout << 4 << y << endl;
//     m_ttbar[y]["combine"]["upup"] = AddHists(m_ttbar[y]["muon"]["upup"], m_ttbar[y]["elec"]["upup"], 1);
//     cout << 4 << y << endl;
//     m_ttbar[y]["combine"]["downdown"] = AddHists(m_ttbar[y]["muon"]["downdown"], m_ttbar[y]["elec"]["downdown"], 1);
//   }
//
//   cout << 5 << endl;
//   m_combine["combine"]["nominal"] = AddHists({m_ttbar["2016"]["combine"]["nominal"], m_ttbar["2017"]["combine"]["nominal"], m_ttbar["2018"]["combine"]["nominal"]}, 1);
//   cout << 5 << endl;
//   m_combine["elec_muonJMS"]["nominal"] = AddHists({m_ttbar["2016"]["elec_muonJMS"]["nominal"], m_ttbar["2017"]["elec_muonJMS"]["nominal"], m_ttbar["2018"]["elec_muonJMS"]["nominal"]}, 1);
//   cout << 5 << endl;
//   m_combine["muon_elecJMS"]["nominal"] = AddHists({m_ttbar["2016"]["muon_elecJMS"]["nominal"], m_ttbar["2017"]["muon_elecJMS"]["nominal"], m_ttbar["2018"]["muon_elecJMS"]["nominal"]}, 1);
//
//   cout << 6 << endl;
//   m_combine["combine"]["upup"] = AddHists({m_ttbar["2016"]["combine"]["upup"], m_ttbar["2017"]["combine"]["upup"], m_ttbar["2018"]["combine"]["upup"]}, 1);
//   cout << 6 << endl;
//   m_combine["elec_muonJMS"]["upup"] = AddHists({m_ttbar["2016"]["elec_muonJMS"]["upup"], m_ttbar["2017"]["elec_muonJMS"]["upup"], m_ttbar["2018"]["elec_muonJMS"]["upup"]}, 1);
//   cout << 6 << endl;
//   m_combine["muon_elecJMS"]["upup"] = AddHists({m_ttbar["2016"]["muon_elecJMS"]["upup"], m_ttbar["2017"]["muon_elecJMS"]["upup"], m_ttbar["2018"]["muon_elecJMS"]["upup"]}, 1);
//
//   cout << 7 << endl;
//   m_combine["combine"]["downdown"] = AddHists({m_ttbar["2016"]["combine"]["downdown"], m_ttbar["2017"]["combine"]["downdown"], m_ttbar["2018"]["combine"]["downdown"]}, 1);
//   cout << 7 << endl;
//   m_combine["elec_muonJMS"]["downdown"] = AddHists({m_ttbar["2016"]["elec_muonJMS"]["downdown"], m_ttbar["2017"]["elec_muonJMS"]["downdown"], m_ttbar["2018"]["elec_muonJMS"]["downdown"]}, 1);
//   cout << 7 << endl;
//   m_combine["muon_elecJMS"]["downdown"] = AddHists({m_ttbar["2016"]["muon_elecJMS"]["downdown"], m_ttbar["2017"]["muon_elecJMS"]["downdown"], m_ttbar["2018"]["muon_elecJMS"]["downdown"]}, 1);
//
//
// }
//
// void Plot(const MapHH &hists, TString y){
//   TH1F* combine_JMSc_nom = hists[y]["nominal"];
//   TH1F* combine_JMSc_uu = hists[y]["upup"];
//   TH1F* combine_JMSc_dd = hists[y]["downdown"];
//
//   TH1F* elec_JMSm_nom = hists[y]["nominal"];
//   TH1F* elec_JMSm_uu = hists[y]["upup"];
//   TH1F* elec_JMSm_dd = hists[y]["downdown"];
//
//   TH1F* muon_JMSe_nom = hists[y]["nominal"];
//   TH1F* muon_JMSe_uu = hists[y]["upup"];
//   TH1F* muon_JMSe_dd = hists[y]["downdown"];
//
//   TH1F* elec_JMSc_nom = hists[y]["nominal"];
//   TH1F* elec_JMSc_uu = hists[y]["upup"];
//   TH1F* elec_JMSc_dd = hists[y]["downdown"];
//
//   TH1F* muon_JMSc_nom = hists[y]["nominal"];
//   TH1F* muon_JMSc_uu = hists[y]["upup"];
//   TH1F* muon_JMSc_dd = hists[y]["downdown"];
}
