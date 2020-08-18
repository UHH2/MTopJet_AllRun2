#include <UHH2/MTopJet/include/CorrectionFactor_BestFit.h>

/* Here onlyt the Up Variation is used for both coreections. The points/factors for the BestFit method
is calculated seperatlly. The up an down factors are symmetric around the nominal value.
Therefore only the differences of f_nominal and f_up is used to calculate the new factor point_i.*/


// to share some code between JetCorrector and JetLeptonCleaner, provide some methods
// dealing with jet energy corrections here:
std::unique_ptr<FactorizedJetCorrector> build_corrector(const std::vector<std::string> & filenames){
  std::vector<JetCorrectorParameters> pars;
  for(const auto & filename : filenames){
    pars.emplace_back(locate_file(filename));
  }
  return std::make_unique<FactorizedJetCorrector>(pars);
}

JetCorrectionUncertainty* corrector_uncertainty(uhh2::Context & ctx, const std::vector<std::string> & filenames, int &direction){

  auto dir = ctx.get("jecsmear_direction", "nominal");
  if (ctx.get("dataset_type") != "MC") {
    direction = 0;
  }
  else if(dir == "up"){
    direction = 1;
  }
  else if(dir == "down"){
    direction = -1;
  }
  else if(dir != "nominal"){
    // direction = 0 is default
    throw runtime_error("JetCorrector: invalid value jecsmear_direction='" + dir + "' (valid: 'nominal', 'up', 'down')");
  }

  // Get optional source of JEC, defaults the total uncertainty if the user doesn't specify one
  std::string source = ctx.get("jecsmear_source", "Total");

  //initialize JetCorrectionUncertainty if shift direction is not "nominal", else return NULL pointer
  if(direction!=0){
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
    JetCorrectionUncertainty* jec_uncertainty = new JetCorrectionUncertainty(*(new JetCorrectorParameters(unc_file.Data(), source)));
    return jec_uncertainty;
  }
  return NULL;

}

/*
.██   ██  ██████  ██████  ███    ██ ███████
. ██ ██  ██      ██    ██ ████   ██ ██
.  ███   ██      ██    ██ ██ ██  ██ █████
. ██ ██  ██      ██    ██ ██  ██ ██ ██
.██   ██  ██████  ██████  ██   ████ ███████
*/

/* Copy frome CorrectionFactor.cc */

double CorrectionFactor_BestFit::get_factor_XCone(double pt, double eta, double point_x){

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
  double factor = f_c + point_x*abs(f_c-f_sys);


  return factor;
}


void CorrectionFactor_BestFit::get_function(TString year){

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

void CorrectionFactor_BestFit::get_additionalSYS(){

  TString dir = "/nfs/dust/cms/user/paaschal/CMSSW_10_2_X/CMSSW_10_2_16/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_SysFromResolution_"+str_year+".root";
  TFile *file = new TFile(filename);

  TString histname = "Up";
  AdditionalSys = (TGraph*)file->Get(histname);

  return;
}

/*
.     ██ ███████  ██████
.     ██ ██      ██
.     ██ █████   ██
.██   ██ ██      ██
. █████  ███████  ██████
*/

/* Implemented in UHH2/common/src/JetCorrections.cxx due to the corrector. */


// #################################################################################################################
// #################################################################################################################



/*
.██████  ██████  ███    ███ ██████  ██ ███    ██ ███████ ██████
██      ██    ██ ████  ████ ██   ██ ██ ████   ██ ██      ██   ██
██      ██    ██ ██ ████ ██ ██████  ██ ██ ██  ██ █████   ██   ██
██      ██    ██ ██  ██  ██ ██   ██ ██ ██  ██ ██ ██      ██   ██
.██████  ██████  ██      ██ ██████  ██ ██   ████ ███████ ██████
*/

CorrectionFactor_BestFit::CorrectionFactor_BestFit(uhh2::Context & ctx,  const std::string& jet_collection_rec, Year year_):
h_oldjets(ctx.get_handle<std::vector<TopJet>>(jet_collection_rec)),
year(year_)
{
  // JEC -----------------------------------------------------------------------
  // setup JEC for XCone -------------------------------------------------------
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

  // XCone ---------------------------------------------------------------------
  if(year == Year::is2016v3)      str_year = "2016";
  else if(year == Year::is2017v2) str_year = "2017";
  else if(year == Year::is2018)   str_year = "2018";

  get_function(str_year); // read file here, not in every event
  get_additionalSYS();
}

// #################################################################################################################
// #################################################################################################################
bool CorrectionFactor_BestFit::process(uhh2::Event & event){
  return false;
}

// #################################################################################################################
// #################################################################################################################
double CorrectionFactor_BestFit::get_mass_BestFit(uhh2::Event & event, vector<TopJet> jet, vector<double> points){
  cout << "start" << endl;
  std::vector<TopJet> oldjets = event.get(h_oldjets);
  std::vector<TopJet> newjets;
  for(unsigned int i=0; i < oldjets.size(); i++) newjets.push_back(oldjets.at(i));
  cout << "points" << endl;

  double point_X = points[0];
  double point_J = points[1];
  vector<float> factors_had, factors_X, factors_J;

  // now correct had subjets
  std::vector<Jet> oldsubjets = oldjets.at(0).subjets();
  std::vector<Jet> newsubjets;
  LorentzVector newjet_v4, newjet_v4_all;
  Jet newjet;
  vector<double> factor_X, factor_J;
  cout << "start loop subjets" << endl;

  for(unsigned int i=0; i < oldsubjets.size(); i++){
    cout << "subjet " << i << endl;
    double f_X = get_factor_XCone(oldsubjets.at(i).pt(), oldsubjets.at(i).eta(), point_X);
    cout << f_X << endl;
    factor_X.push_back(f_X);
  }
  factor_J = jet_corrector_MC->get_factor_JEC(event, point_J);
  cout << factor_J.size() << "   " << factor_X.size() << endl;

  // newjet_v4 = oldsubjets.at(i).v4() * factor_X * factor_J;
  // newjet_v4_all += newjet_v4;
  // factors_had.push_back(factor_X*factor_J);
  // factors_X.push_back(factor_X);
  // factors_J.push_back(factor_J);

  // newjet.set_v4(newjet_v4);  // Because of this, all JEC factors are set to the defaul value of 1


  double mass =0;
  return mass;
}
