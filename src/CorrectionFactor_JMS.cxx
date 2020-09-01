#include <UHH2/MTopJet/include/CorrectionFactor_JMS.h>

/*
Here only the Up Variation is used for both coreections. The points/factors for the BestFit method
is calculated seperatly and have to be inserted by hand. The up an down factors are symmetric around the nominal value.
Therefore only the differences of f_nominal and f_up is used to calculate the new factor point_i.

JMS = JetMassScale
*/

/*
.██   ██  ██████  ██████  ███    ██ ███████
. ██ ██  ██      ██    ██ ████   ██ ██
.  ███   ██      ██    ██ ██ ██  ██ █████
. ██ ██  ██      ██    ██ ██  ██ ██ ██
.██   ██  ██████  ██████  ██   ████ ███████
*/

/* Copy frome CorrectionFactor.cc */

double CorrectionFactor_JMS::get_factor_XCone(double pt, double eta, double point_x){

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

  // if(CorUp || CorDown){
  // first get central factor
  double f_c = CentralCorrectionFunctions[etabin]->Eval(pt);
  // get up/down factor
  double f_ud = UpDownCorrectionGraphs[etabin]->Eval(pt);
  // get difference
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


void CorrectionFactor_JMS::get_function(TString year){

  TString dir = "/nfs/dust/cms/user/paaschal/CMSSW_10_2_X/CMSSW_10_2_16/src/UHH2/MTopJet/CorrectionFile/";
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

  TString dir = "/nfs/dust/cms/user/paaschal/CMSSW_10_2_X/CMSSW_10_2_16/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_SysFromResolution_"+str_year+".root";
  TFile *file = new TFile(filename);

  TString histname = "Up";
  AdditionalSys = (TGraph*)file->Get(histname);

  return;
}

// #################################################################################################################
// #################################################################################################################

/*
.     ██ ███████  ██████
.     ██ ██      ██
.     ██ █████   ██
.██   ██ ██      ██
. █████  ███████  ██████
*/

/*
Implemented in UHH2/common/src/JetCorrections.cxx due to the corrector.
Additional Code is implemented in YearRunSwitchers (access the get_factor_JEC function in JetCorrection)
and AnalysisModule (dummy function).
*/


// #################################################################################################################
// #################################################################################################################

/*
.██████  ██████  ███    ███ ██████  ██ ███    ██ ███████ ██████
██      ██    ██ ████  ████ ██   ██ ██ ████   ██ ██      ██   ██
██      ██    ██ ██ ████ ██ ██████  ██ ██ ██  ██ █████   ██   ██
██      ██    ██ ██  ██  ██ ██   ██ ██ ██  ██ ██ ██      ██   ██
.██████  ██████  ██      ██ ██████  ██ ██   ████ ███████ ██████
*/

CorrectionFactor_JMS::CorrectionFactor_JMS(uhh2::Context & ctx, const std::string& jet_collection_rec, const std::string& jet_collection_gen, Year year_):
h_oldjets(ctx.get_handle<std::vector<TopJet>>(jet_collection_rec)),
h_genjets(ctx.get_handle<std::vector<GenTopJet>>(jet_collection_gen)),
year(year_)
{
  // const Year & year = extract_year(ctx);
  // JEC -----------------------------------------------------------------------------------------------------------
  // setup JEC for XCone -------------------------------------------------------------------------------------------
  isMC = (ctx.get("dataset_type") == "MC");

  string jec_tag_2016 = "Summer16_07Aug2017";
  string jec_ver_2016 = "20";

  string jec_tag_2017 = "Fall17_17Nov2017";
  string jec_ver_2017 = "32";

  string jec_tag_2018 = "Autumn18";
  string jec_ver_2018 = "19";

  string jec_jet_coll = "AK4PFchs";

  if (isMC)
  {
    jet_corrector_MC.reset(new YearSwitcher(ctx));
    jet_corrector_MC->setup2016v3(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2017v2(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll),jet_collection_rec));
    jet_corrector_MC->setup2018(std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll),jet_collection_rec));
  }

  else
  {
    jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      jec_switcher_16->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr),jet_collection_rec));
    }

    jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    for (const auto & runItr : runPeriods2017) {
      jec_switcher_17->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr),jet_collection_rec));
    }

    jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    for (const auto & runItr : runPeriods2018) {
      jec_switcher_18->setupRun(runItr, std::make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr),jet_collection_rec));
    }

    jet_corrector_data.reset(new YearSwitcher(ctx));
    jet_corrector_data->setup2016v3(jec_switcher_16);
    jet_corrector_data->setup2017v2(jec_switcher_17);
    jet_corrector_data->setup2018(jec_switcher_18);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // JER -----------------------------------------------------------------------------------------------------------
  JERSmearing::SFtype1 JER_sf;
  std::string filenameAppend = "AK4PFPuppi.txt";

  std::string resFilename = "";
  if (year == Year::is2016v3) {
    JER_sf = JERSmearing::SF_13TeV_Summer16_25nsV1;
    resFilename = "2016/Summer16_25nsV1_MC_PtResolution_"+filenameAppend;
  } else if (year == Year::is2017v2) {
    JER_sf = JERSmearing::SF_13TeV_Fall17_V3;
    resFilename = "2017/Fall17_V3_MC_PtResolution_"+filenameAppend;
  } else if (year == Year::is2018) {
    JER_sf = JERSmearing::SF_13TeV_Autumn18_RunABCD_V4;
    resFilename = "2018/Autumn18_V4_MC_PtResolution_"+filenameAppend;
  } else {
    throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
  }

  // Use Constructer just to setup the class. The handels do not matter since the
  // function apply_JER_smearing is not using them
  JER_Smearer.reset(new GenericJetResolutionSmearer(ctx, jet_collection_rec, jet_collection_gen, JER_sf, resFilename));

  // ---------------------------------------------------------------------------------------------------------------
  // XCone ---------------------------------------------------------------------------------------------------------
  if(year == Year::is2016v3)      str_year = "2016";
  else if(year == Year::is2017v2) str_year = "2017";
  else if(year == Year::is2018)   str_year = "2018";

  get_function(str_year); // read file here, not in every event
  get_additionalSYS();
}

// #################################################################################################################
// #################################################################################################################
bool CorrectionFactor_JMS::process(uhh2::Event & event){ // Dummy
  double weight = event.weight; // avoid warning
  weight += 0;                  // avoid warning
  return false;
}

// #################################################################################################################
// #################################################################################################################
double CorrectionFactor_JMS::get_mass_BestFit(uhh2::Event & event, vector<double> points){
  double debug = false;
  if(debug) std::cout << "\nBestFit: Start Get_Mass_BestFit";
  std::vector<TopJet> oldjets = event.get(h_oldjets);

  // Get points from input vector ----------------------------------------------------------------------------------
  double point_J = points[0];
  double point_X = points[1];
  vector<float> factors_had, factors_X, factors_J;

  // Get new correction factors for XCone and JEC ------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Get oldsubjets";
  std::vector<Jet> oldsubjets = oldjets.at(0).subjets();
  if(debug) std::cout << "\nBestFit: number oldsubjets " << oldsubjets.size() << endl;
  vector<double> factor_X, factor_J;
  if(debug) std::cout << "\nBestFit: Get X_factor";
  for (auto & jet : oldsubjets)
  {
    factor_X.push_back(get_factor_XCone(jet.pt(), jet.eta(), point_X));
  }
  if(debug) std::cout << "\nBestFit: Get J_factor";
  factor_J = jet_corrector_MC->get_factor_JEC_year(event, point_J);

  // Apply new correction factors for each subjet ------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Set newsubjets";
  std::vector<TLorentzVector> newsubjets_v4;
  std::vector<Jet> newsubjets;

  for(unsigned int subjet=0; subjet<factor_J.size(); subjet++)
  {
    newsubjets_v4.push_back(jet_to_tlorentz(oldsubjets[subjet]));
    newsubjets_v4[subjet] *= factor_J[subjet]*factor_X[subjet];

    // Set Jets for JER
    LorentzVector newsubjets_lorentz = tlorentz_to_lorentz(newsubjets_v4[subjet]); // .set_v4 does not take TLorentzVector!
    Jet newsubjet;
    newsubjets.push_back(newsubjet);
    newsubjets[subjet].set_v4(newsubjets_lorentz);
  }

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

  // Creat new hadjet ----------------------------------------------------------------------------------------------
  if(debug) std::cout << "\nBestFit: Set new fatjet";
  TopJet xcone_had_jms;
  xcone_had_jms.set_subjets(newsubjets);
  LorentzVector xcone_had_jms_v4 = xcone_had_jms.v4();

  // return mass _--------------------------------------------------------------------------------------------------
  return xcone_had_jms_v4.M();
}
