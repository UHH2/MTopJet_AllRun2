/* #pragma once */

/* #include "UHH2/core/include/AnalysisModule.h" */
/* #include <UHH2/core/include/Event.h> */
/* #include "UHH2/common/include/ObjectIdUtils.h" */
/* #include <UHH2/core/include/NtupleObjects.h> */
/* #include "UHH2/common/include/Utils.h" */
/* #include "UHH2/common/include/JetCorrections.h" */
/* #include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h" */
/* #include "UHH2/JetMETObjects/interface/FactorizedJetCorrector.h" */
/* #include "UHH2/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h" */
/* #include "UHH2/JetMETObjects/interface/JetCorrectorParameters.h" */

/* class FactorizedJetCorrector; */
/* class FactorizedJetCorrectorCalculator; */


/* class TopJetLeptonCleaner: public uhh2::AnalysisModule { */
/* public: */
/*     // jec_filenames is teh same as for the JetCorrector. */
/*     explicit TopJetLeptonCleaner(uhh2::Context & ctx, const std::vector<std::string> & jec_filenames); */
    
/*     void set_muon_id(const MuonId & mu_id_){ */
/*         mu_id = mu_id_; */
/*     } */
    
/*     void set_electron_id(const ElectronId & ele_id_){ */
/*         ele_id = ele_id_; */
/*     } */
    
/*     void set_drmax(double drmax_){ */
/*         drmax = drmax_; */
/*     } */
    
/*     virtual bool process(uhh2::Event & event) override; */
    
/*     virtual ~TopJetLeptonCleaner(); */
    
/* private: */
/*     std::unique_ptr<FactorizedJetCorrector> corrector; */
/*     MuonId mu_id; */
/*     ElectronId ele_id; */
/*     double drmax = 0.4; */
/*     JetCorrectionUncertainty* jec_uncertainty; */
/*     int direction = 0; // -1 = down, +1 = up, 0 = nominal */
/* }; */
