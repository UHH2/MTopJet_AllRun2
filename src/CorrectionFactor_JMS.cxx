#include <UHH2/MTopJet/include/CorrectionFactor_JMS.h>
#include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"

using namespace uhh2;
using namespace std;


std::unique_ptr<FactorizedJetCorrector> build_corrector(const std::vector<std::string> & filenames);
JetCorrectionUncertainty* corrector_uncertainty(const std::vector<std::string> & filenames);

double eta_bins[13] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191};
std::vector<TF1*> CentralCorrectionFunctions;
std::vector<TGraph*> UpDownCorrectionGraphs;
bool debug = false;

double x_cor_err = 0.26344;
double x_jec_err = 0.100036; // average

/*
Here only the Up Variation is used for both coreections. The points/factors for the BestFit method
is calculated seperatly and have to be inserted by hand. The up an down factors are symmetric around the nominal value.
Therefore only the differences of f_nominal and f_up is used to calculate the new factor point_i.

JMS = JetMassScale
*/

// #################################################################################################################
// #################################################################################################################
// #################################################################################################################

CorrectionFactor_JMS::CorrectionFactor_JMS(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen, Year year_):
h_oldjets(ctx.get_handle<std::vector<TopJet>>(jet_collection_rec)),
h_genjets(ctx.get_handle<std::vector<GenTopJet>>(jet_collection_gen)),
year(year_)
{
  // const Year & year = extract_year(ctx);
  // JEC -----------------------------------------------------------------------------------------------------------
  isMC = (ctx.get("dataset_type") == "MC");

  string jec_tag_2016 = "Summer16_07Aug2017";
  string jec_ver_2016 = "20";

  string jec_tag_2017 = "Fall17_17Nov2017";
  string jec_ver_2017 = "32";

  string jec_tag_2018 = "Autumn18";
  string jec_ver_2018 = "19";

  string jec_jet_coll = "AK4PFchs";

  if (isMC){
    jet_corrector_MC.reset(new YearSwitcher(ctx));
    jet_corrector_MC->setup2016v3(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2017v2(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2018(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll),jet_collection_rec));

    corrector_MC_2016 = build_corrector(JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll));
    uncertainty_MC_2016 = corrector_uncertainty(JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll));

    corrector_MC_2017 = build_corrector(JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll));
    uncertainty_MC_2017 = corrector_uncertainty(JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll));

    corrector_MC_2018 = build_corrector(JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll));
    uncertainty_MC_2018 = corrector_uncertainty(JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll));

  } else{
    throw runtime_error("JMS should not be applied to data!");
  }

  // JER -----------------------------------------------------------------------------------------------------------
  JERSmearing::SFtype1 JER_sf;

  std::string jer_tag = "";
  std::string jetCollection = "AK4PFchs";
  const Year & year = extract_year(ctx);

  if (year == Year::is2016v3) {
    jer_tag = "Summer16_25nsV1";
  } else if (year == Year::is2017v2) {
    jer_tag = "Fall17_V3";
  } else if (year == Year::is2018) {
    jer_tag = "Autumn18_V7";
  } else {
    throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
  }

  JER_Smearer.reset(new GenericJetResolutionSmearer(ctx,  jet_collection_rec, jet_collection_gen, JERFiles::JERPathStringMC(jer_tag,jetCollection,"SF"), JERFiles::JERPathStringMC(jer_tag,jetCollection,"PtResolution")));

  // XCone ---------------------------------------------------------------------------------------------------------
  if(year == Year::is2016v3)      str_year = "2016";
  else if(year == Year::is2017v2) str_year = "2017";
  else if(year == Year::is2018)   str_year = "2018";

  get_function(str_year); // read file here, not in every event
  get_additionalSYS();
}

// #################################################################################################################
// #################################################################################################################
// #################################################################################################################

bool CorrectionFactor_JMS::process(uhh2::Event & event){ // Dummy
  double weight = event.weight; // avoid warning
  weight += 0;                  // avoid warning
  return false;
}

vector<Jet> CorrectionFactor_JMS::GetSubjetsJMS(uhh2::Event & event, const VecD& points, const int &PoU){
  if(debug) std::cout << "\nBestFit: Start Get_Mass_BestFit";
  std::vector<TopJet> oldjets = event.get(h_oldjets);

  // Get points from input vector ----------------------------------------------------------------------------------
  double point_J = points[0];
  double point_X = points[1];
  vector<float> factors_had, factors_X, factors_J;
  if(debug) cout << point_J  << "     " << point_X << endl;
  // Get new correction factors for XCone and JEC ------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Get oldsubjets";
  std::vector<Jet> oldsubjets = oldjets.at(0).subjets();
  if(debug) std::cout << "\nBestFit: number oldsubjets " << oldsubjets.size() << endl;
  // first get JEC factors
  VecDD corrections_JEC;
  if(debug) std::cout << "\nBestFit: Get J_factor";
  if(year == Year::is2016v3) corrections_JEC = get_corrections_JEC(*corrector_MC_2016, oldsubjets, event, uncertainty_MC_2016);
  else if(year == Year::is2017v2) corrections_JEC = get_corrections_JEC(*corrector_MC_2017, oldsubjets, event, uncertainty_MC_2017);
  else if(year == Year::is2018) corrections_JEC = get_corrections_JEC(*corrector_MC_2018, oldsubjets, event, uncertainty_MC_2018);

  VecD factor_J = get_factor_JEC(corrections_JEC, point_J);

  // Now get XCone factor using corrected jet
  VecD factor_X;
  VecDD corrections_XCone;
  if(debug) std::cout << "\nBestFit: Get X_factor";
  for (unsigned int i=0; i<oldsubjets.size(); i++){
    double pt = oldsubjets[i].pt()*factor_J[i];
    double eta = oldsubjets[i].eta();
    VecD single_corrections = get_corrections_XCone(pt, eta);
    corrections_XCone.push_back(single_corrections);
    factor_X.push_back(get_factor_XCone(pt, single_corrections[0], single_corrections[1], point_X));
  }

  // Apply new correction factors for each subjet ------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Set newsubjets";
  std::vector<TLorentzVector> newsubjets_v4;
  std::vector<Jet> newsubjets;
  LorentzVector xcone_had_jms_noJER_v4;
  VecD sigma_JMS = error_JMS(factor_J, factor_X, corrections_JEC, corrections_XCone);

  for(unsigned int subjet=0; subjet<factor_J.size(); subjet++)
  {
    if(debug) cout << "Factor JEC   = " << factor_J[subjet] << endl;
    if(debug) cout << "Factor XCone = " << factor_X[subjet] << endl;

    newsubjets_v4.push_back(jet_to_tlorentz(oldsubjets[subjet]));
    double f_JMS = factor_J[subjet]*factor_X[subjet];
    if(debug) cout << "factor JMS: " << f_JMS << endl;
    if(debug) cout << "error JMS:  " << sigma_JMS[subjet] << endl;
    if(PoU==1)      f_JMS += sigma_JMS[subjet];
    else if(PoU==2) f_JMS -= sigma_JMS[subjet];
    newsubjets_v4[subjet] *= f_JMS;

    // Set Jets for JER
    LorentzVector newsubjets_lorentz = tlorentz_to_lorentz(newsubjets_v4[subjet]); // .set_v4 does not take TLorentzVector!
    Jet newsubjet;
    newsubjets.push_back(newsubjet);
    newsubjets[subjet].set_v4(newsubjets_lorentz);
    xcone_had_jms_noJER_v4 +=  newsubjets_lorentz;
  }

  if(debug) cout << "\nNo JER mass: " << xcone_had_jms_noJER_v4.M();
  // Apply JER for each subjet -------------------------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: JER";
  std::vector<GenTopJet>* gen_topjets = &event.get(h_genjets);
  // match gentopjet to rectopjet
  double dR0 = deltaR(oldjets.at(0), gen_topjets->at(0));
  double dR1 = deltaR(oldjets.at(0), gen_topjets->at(1));
  unsigned int index_gen=0;
  if(dR0 > dR1) index_gen=1;

  std::vector<GenJet> gen_subjets = gen_topjets->at(index_gen).subjets();
  JER_Smearer->apply_JER_smearing(newsubjets, gen_subjets, 0.4, event.rho);

  return newsubjets;
}

// #################################################################################################################
// #################################################################################################################
// #################################################################################################################

/* Copy frome CorrectionFactor.cc */

VecD CorrectionFactor_JMS::get_corrections_XCone(double pt, double eta){
  // first transform eta to absolute value
  if(eta < 0) eta *= -1;

  // make sure that there is not jet with eta > 2.5
  // (this should not happen because in CombineXCone this should be sorted out)
  if(eta > 2.5) eta = 2.5;

  // get eta bin
  int etabin = 0;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    if(eta > eta_bins[i]) etabin = i;
  }

  // get factor from function
  if(pt > 430) pt = 430;
  if(pt < 30)  pt = 30;

  double f_c  = CentralCorrectionFunctions[etabin]->Eval(pt); // first get central factor
  double f_ud = UpDownCorrectionGraphs[etabin]->Eval(pt);     // get up/down factor
  return {f_c, f_ud};
}

double CorrectionFactor_JMS::get_factor_XCone(double pt, double f_c, double f_ud, double point_x){
  double df = f_ud - f_c;
  // get additional sys (this already is a difference)
  double dg = AdditionalSys->Eval(pt);
  // to get the total factor one has to:
  // 1. add df and dg in quadrature (this gives the total diff to central factor)
  // 2. add this to central factor
  // Dont need CorUp or CorDown.
  // Variations are calculated the same way just with another point_x
  double f_sys = f_c + sqrt(df*df + dg*dg);
  // get calculated factor from chi2 best fit for xcone corrections (point_x)
  double factor = f_c + point_x*fabs(f_c-f_sys);

  return factor;
}


void CorrectionFactor_JMS::get_function(TString & year){

  TString dir = "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_allHad_"+year+".root";
  TFile *file = new TFile(filename);

  TString histname_updown;
  TString histname_central;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    stringstream BinString;
    BinString << i;
    histname_updown = BinString.str();
    histname_central = BinString.str();

    histname_updown += "/Up";
    UpDownCorrectionGraphs.push_back((TGraph*)file->Get(histname_updown));

    histname_central += "/Central";
    CentralCorrectionFunctions.push_back((TF1*)file->Get(histname_central));
  }
}

void CorrectionFactor_JMS::get_additionalSYS(){

  TString dir = "/nfs/dust/cms/user/schwarzd/CMSSW10/CMSSW_10_2_17/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_SysFromResolution_"+str_year+".root";
  TFile *file = new TFile(filename);

  TString histname = "Up";
  AdditionalSys = (TGraph*)file->Get(histname);

  return;
}

// #################################################################################################################
// #################################################################################################################
// #################################################################################################################

/*
Implemented in UHH2/common/src/JetCorrections.cxx due to the corrector.
Additional Code is implemented in YearRunSwitchers (access the get_factor_JEC function in JetCorrection)
and AnalysisModule (dummy function).
==============> Check if all deleted!!!!!
*/

std::unique_ptr<FactorizedJetCorrector> build_corrector(const std::vector<std::string> & filenames){
  std::vector<JetCorrectorParameters> pars;
  for(const auto & filename : filenames){
    pars.emplace_back(locate_file(filename));
  }
  return uhh2::make_unique<FactorizedJetCorrector>(pars);
}

VecDD CorrectionFactor_JMS::get_corrections_JEC(FactorizedJetCorrector & corrector, vector<Jet> subjets, const Event & event, JetCorrectionUncertainty* jec_unc){
  VecDD JEC_corrections;
  for(unsigned int i=0; i<subjets.size(); i++){

    auto factor_raw = subjets[i].JEC_factor_raw();
    corrector.setJetPt(subjets[i].pt() * factor_raw);
    corrector.setJetEta(subjets[i].eta());
    corrector.setJetE(subjets[i].energy() * factor_raw);
    corrector.setJetA(subjets[i].jetArea());
    corrector.setJetPhi(subjets[i].phi());
    corrector.setRho(event.rho);
    auto correctionfactors = corrector.getSubCorrections();
    auto correctionfactor = correctionfactors.back();

    LorentzVector jet_v4_corrected = subjets[i].v4() * (factor_raw *correctionfactor);

    // get factor for up variation
    // ignore jets with very low pt or high eta, avoiding a crash from the JESUncertainty tool
    double pt = jet_v4_corrected.Pt();
    double eta = jet_v4_corrected.Eta();
    if(pt<5. || fabs(eta)>5.){
      JEC_corrections.push_back({0.,0.});
    }
    else{
      jec_unc->setJetEta(eta);
      jec_unc->setJetPt(pt);

      double unc = 0.;
      unc = jec_unc->getUncertainty(true);
      double correctionfactor_up = correctionfactor * (1 + fabs(unc));
      JEC_corrections.push_back({correctionfactor, correctionfactor_up});
    }
  }
  return JEC_corrections;
}

VecD CorrectionFactor_JMS::get_factor_JEC(VecDD JEC_corrections, double point_J){
  VecD JEC_factors;
  for(unsigned int i=0; i<JEC_corrections.size(); i++){
    double correctionfactor_jms = JEC_corrections[i][0] + point_J * fabs(JEC_corrections[i][1]-JEC_corrections[i][0]);
    JEC_factors.push_back(correctionfactor_jms);
  }
  return JEC_factors;
}

JetCorrectionUncertainty* corrector_uncertainty(const std::vector<std::string> & filenames){
  //take name from the L1FastJet correction (0th element of filenames) and replace "L1FastJet" by "UncertaintySources" to get the proper name of the uncertainty file
  TString unc_file = locate_file(filenames[0]);
  if (unc_file.Contains("L1FastJet")) {
    unc_file.ReplaceAll("L1FastJet","UncertaintySources");
  }
  else if (unc_file.Contains("L2Relative")) {
    unc_file.ReplaceAll("L2Relative","UncertaintySources");
  }
  else {
    throw runtime_error("WARNING No JEC Uncertainty File found!");
  }
  JetCorrectionUncertainty* jec_uncertainty = new JetCorrectionUncertainty(*(new JetCorrectorParameters(unc_file.Data(), "Total")));
  return jec_uncertainty;
}

// #################################################################################################################
// #################################################################################################################
// #################################################################################################################

VecD CorrectionFactor_JMS::error_JMS(VecD factor_J, VecD factor_X, VecDD corrections_J, VecDD corrections_X){
  VecD errors;
  double rho = -0.21;
  for(unsigned int i=0; i<factor_J.size(); i++){
    double delta_fJ = x_jec_err*fabs(corrections_J[i][0]-corrections_J[i][1]);
    double delta_fX = x_cor_err*fabs(corrections_X[i][0]-corrections_X[i][1]);
    double PoU_jec  = delta_fJ*factor_X[i];
    double PoU_cor  = delta_fX*factor_J[i];
    double error    = sqrt(PoU_jec*PoU_jec+PoU_cor*PoU_cor+2*rho*PoU_jec*PoU_cor);
    errors.push_back(error);
  }
  return errors;
}

double CorrectionFactor_JMS::get_wmass_BestFit(vector<Jet> newsubjets, const std::vector<int>& WSubIndices){
  LorentzVector Wjet_v4_match = newsubjets[WSubIndices[0]].v4() + newsubjets[WSubIndices[1]].v4();
  if(debug) cout << "mw: " << Wjet_v4_match.M() << endl;
  return Wjet_v4_match.M();
}

double CorrectionFactor_JMS::get_mass_BestFit(vector<Jet> newsubjets){
  // Creat new hadjet ----------------------------------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Set new fatjet";
  if(newsubjets.size() < 3) return -1.0;
  LorentzVector xcone_had_jms_JER_v4 = newsubjets[0].v4() + newsubjets[1].v4() + newsubjets[2].v4();
  if(debug) cout << "\nJER JMS mass: " << xcone_had_jms_JER_v4.M() << endl;
  // return mass ---------------------------------------------------------------------------------------------------
  return xcone_had_jms_JER_v4.M();
}
