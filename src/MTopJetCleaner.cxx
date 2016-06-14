// #include "UHH2/core/include/Event.h"
// #include "UHH2/MTopJet/include/MTopJetCleaner.h"
// #include "UHH2/common/include/Utils.h"
// #include "UHH2/core/include/Utils.h"
// #include "UHH2/common/include/JetCorrections.h"
// #include <UHH2/core/include/LorentzVector.h>

// #include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h"
// #include "UHH2/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"
// #include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h"
// #include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h"

// #include <UHH2/core/include/LorentzVector.h>

// #include <string>

// using namespace uhh2;
// using namespace std;

// namespace {
    
// // to share some code between JetCorrector and JetLeptonCleaner, provide some methods
// // dealing with jet energy corrections here:
// std::unique_ptr<FactorizedJetCorrector> build_corrector(const std::vector<std::string> & filenames){
//     std::vector<JetCorrectorParameters> pars;
//     for(const auto & filename : filenames){
//         pars.emplace_back(locate_file(filename));
//     }
//     return uhh2::make_unique<FactorizedJetCorrector>(pars);
// }

// void correct_jet(FactorizedJetCorrector & corrector, Jet & jet, const Event & event, JetCorrectionUncertainty* jec_unc = NULL, int jec_unc_direction=0){
//     auto factor_raw = jet.JEC_factor_raw();
//     corrector.setJetPt(jet.pt() * factor_raw);
//     corrector.setJetEta(jet.eta());
//     corrector.setJetE(jet.energy() * factor_raw);
//     corrector.setJetA(jet.jetArea());
//     corrector.setRho(event.rho);
//     auto correctionfactor = corrector.getCorrection();

//     LorentzVector jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);

//     if(jec_unc_direction!=0){
//       if (jec_unc==NULL){
// 	std::cerr << "JEC variation should be applied, but JEC uncertainty object is NULL! Abort." << std::endl;
// 	exit(EXIT_FAILURE);
//       }
//       // ignore jets with very low pt or high eta, avoiding a crash from the JESUncertainty tool
//       double pt = jet_v4_corrected.Pt();
//       double eta = jet_v4_corrected.Eta();
//       if (!(pt<5. || fabs(eta)>5.)) {
      
// 	jec_unc->setJetEta(eta);
// 	jec_unc->setJetPt(pt);
	
// 	double unc = 0.;	  
// 	if (jec_unc_direction == 1){
// 	  unc = jec_unc->getUncertainty(1);
// 	  correctionfactor *= (1 + fabs(unc));
// 	} else if (jec_unc_direction == -1){
// 	  unc = jec_unc->getUncertainty(-1);
// 	  correctionfactor *= (1 - fabs(unc));
// 	}
// 	jet_v4_corrected = jet.v4() * (factor_raw *correctionfactor);
//       }
//     }

//     jet.set_v4(jet_v4_corrected);
//     jet.set_JEC_factor_raw(1. / correctionfactor);
// }

// JetCorrectionUncertainty* corrector_uncertainty(uhh2::Context & ctx, const std::vector<std::string> & filenames, int &direction){
    
//     auto dir = ctx.get("jecsmear_direction", "nominal");
//     if(dir == "up"){
//         direction = 1;
//     }
//     else if(dir == "down"){
//         direction = -1;
//     }
//     else if(dir != "nominal"){
//         // direction = 0 is default
//         throw runtime_error("JetCorrector: invalid value jecsmear_direction='" + dir + "' (valid: 'nominal', 'up', 'down')");
//     }

//     //initialize JetCorrectionUncertainty if shift direction is not "nominal", else return NULL pointer
//     if(direction!=0){
//       //take name from the L1FastJet correction (0th element of filenames) and replace "L1FastJet" by "Uncertainty" to get the proper name of the uncertainty file
//       TString unc_file = locate_file(filenames[0]);
//       if (unc_file.Contains("L1FastJet")) {
//         unc_file.ReplaceAll("L1FastJet","Uncertainty");
//       }
//       else if (unc_file.Contains("L2Relative")) {
//         unc_file.ReplaceAll("L2Relative","Uncertainty");
//       }
//       else {
//         throw runtime_error("WARNING No JEC Uncertainty File found!");
//       }
//       JetCorrectionUncertainty* jec_uncertainty = new JetCorrectionUncertainty(unc_file.Data());
//       return jec_uncertainty;
//     }
//     return NULL;
    
// }

// }



// // ** JetLeptonCleaner
// TopJetLeptonCleaner::TopJetLeptonCleaner(uhh2::Context & ctx, const std::vector<std::string> & filenames){
//     corrector = build_corrector(filenames);
//     direction = 0;
//     jec_uncertainty = corrector_uncertainty(ctx, filenames, direction) ;
// }

// bool TopJetLeptonCleaner::process(uhh2::Event & event){
//     assert(event.jets);
//     if(event.muons){
//         for(const auto & mu : *event.muons){
//             if(mu_id && !(mu_id(mu, event))) continue;
//             for(auto & jet : *event.jets){
//                 if(deltaR(jet, mu) < drmax && jet.muonMultiplicity() > 0){
//                     auto jet_p4_raw = jet.v4() * jet.JEC_factor_raw();
//                     // note that muon energy fraction as stored in the jet refers to the raw jet energy.
//                     double muon_energy_in_jet = jet_p4_raw.E() * jet.muonEnergyFraction();
//                     double new_muon_energy_in_jet = muon_energy_in_jet - mu.energy();
                    
//                     // test compatibility of the hypothesis that the muon has been clustered to the jet with
//                     // the jet information. The hypothesis is rejected if the muon energy in the jet is too small
//                     // (but allow 10% off). Note that in general (for muon multiplicity > 1), the muon energy in
//                     // the jet might be larger than from the single muon; make sure to consider that in the test
//                     // by requiring one direction in the comparison only in case the muon multiplicity is 1.
//                     if(new_muon_energy_in_jet < -0.1 * mu.energy() || (jet.muonMultiplicity() == 1 && new_muon_energy_in_jet > 0.1 * mu.energy())){
//                         continue;
//                     }
//                     jet_p4_raw -= mu.v4();
//                     // if that subtraction flipped the jet direction (angle between new and old > 90 degrees or pi/2), emit a warning and set its momentum to 0.
//                     // Only warn if pt > 5GeV (otherwise, the jet is 0 anyway for all practical purposes).
//                     if(jet_p4_raw.pt() > 5 && deltaR(jet_p4_raw, jet) > M_PI/2){
//                         cout << "Warning: subtracting lepton flipped jet direction" << endl;
//                         jet.set_v4(LorentzVector());
//                         continue;
//                     }
//                     // re-correct jet. First, set p4_raw = p4_corrected such that
//                     // the 'correct_jet' method does what it should do if using JEC_factor_raw ...
//                     jet.set_JEC_factor_raw(1.0);
//                     jet.set_v4(jet_p4_raw);
//                     // set new muon multiplicity and muon energy fraction:
//                     jet.set_muonMultiplicity(jet.muonMultiplicity() - 1);
//                     jet.set_muonEnergyFraction(max(new_muon_energy_in_jet / jet_p4_raw.E(), 0.0));
//                     correct_jet(*corrector, jet, event, jec_uncertainty, direction);
//                 }
//             }
//         }
//     }
//     if(event.electrons){
//         for(const auto & ele : *event.electrons){
//             if(ele_id && !(ele_id(ele, event))) continue;
//             for(auto & jet : *event.jets){
//                 if(deltaR(jet, ele) < drmax && jet.electronMultiplicity() > 0){
//                     auto jet_p4_raw = jet.v4() * jet.JEC_factor_raw();
//                     double electron_energy_in_jet = jet_p4_raw.E() * jet.chargedEmEnergyFraction();
//                     double new_electron_energy_in_jet = electron_energy_in_jet - ele.energy();
                    
//                     if(new_electron_energy_in_jet < -0.1 * ele.energy() || (jet.chargedEmEnergyFraction() == 1 && new_electron_energy_in_jet > 0.1 * ele.energy())){
//                         continue;
//                     }
//                     jet_p4_raw -= ele.v4();
//                     if(jet_p4_raw.pt() > 5 && deltaR(jet_p4_raw, jet) > M_PI/2){
//                         cout << "Warning: subtracting lepton flipped jet direction" << endl;
//                         jet.set_v4(LorentzVector());
//                         continue;
//                     }
//                     // re-correct jet:
//                     jet.set_JEC_factor_raw(1.0);
//                     jet.set_v4(jet_p4_raw);
//                     jet.set_electronMultiplicity(jet.electronMultiplicity() - 1);
//                     jet.set_chargedEmEnergyFraction(max(new_electron_energy_in_jet / jet_p4_raw.E(), 0.0));
//                     correct_jet(*corrector, jet, event, jec_uncertainty, direction);
//                 }
//             }
//         }
//     }
//     return true;
// }

// // see ~JetCorrector
// TopJetLeptonCleaner::~TopJetLeptonCleaner(){}

