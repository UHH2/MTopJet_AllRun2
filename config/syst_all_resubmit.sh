# sframe_batch -r $1/MTopJetPostSelection_elec.xml
# sframe_batch -r $1/MTopJetPostSelection_muon.xml
# sframe_batch -r $1/MTopJetPostSelection_JMSmuon_elec.xml
# sframe_batch -r $1/MTopJetPostSelection_JMSelec_muon.xml

source syst_resubmit.sh 2016MTopJetmuon
source syst_resubmit.sh 2016MTopJetelec
source syst_resubmit.sh 2016MTopJetJMSelecmuon
source syst_resubmit.sh 2016MTopJetJMSmuonelec

source syst_resubmit.sh 2017MTopJetmuon
source syst_resubmit.sh 2017MTopJetelec
source syst_resubmit.sh 2017MTopJetJMSelecmuon
source syst_resubmit.sh 2017MTopJetJMSmuonelec

source syst_resubmit.sh 2018MTopJetmuon
source syst_resubmit.sh 2018MTopJetelec
source syst_resubmit.sh 2018MTopJetJMSelecmuon
source syst_resubmit.sh 2018MTopJetJMSmuonelec
