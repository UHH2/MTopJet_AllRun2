#!/bin/bash

sframe_batch.py -s $1/MTopJetPostSelection_muon.xml
sframe_batch.py -s $1/MTopJetPostSelection_elec.xml
sframe_batch.py -s $1/MTopJetPostSelection_JMSelec_muon.xml
sframe_batch.py -s $1/MTopJetPostSelection_JMSmuon_elec.xml
python sframe_syst_batch.py $1/MTopJetPostSelection_JMSelec_SYS_muon.xml
python sframe_syst_batch.py $1/MTopJetPostSelection_JMSmuon_SYS_elec.xml

# sframe_batch.py -s 2017/MTopJetPostSelection_muon.xml
# sframe_batch.py -s 2017/MTopJetPostSelection_elec.xml
# sframe_batch.py -s 2017/MTopJetPostSelection_JMSelec_muon.xml
# sframe_batch.py -s 2017/MTopJetPostSelection_JMSmuon_elec.xml
# python sframe_syst_batch.py 2017/MTopJetPostSelection_JMSelec_SYS_muon.xml
# python sframe_syst_batch.py 2017/MTopJetPostSelection_JMSmuon_SYS_elec.xml
#
# sframe_batch.py -s 2018/MTopJetPostSelection_muon.xml
# sframe_batch.py -s 2018/MTopJetPostSelection_elec.xml
# sframe_batch.py -s 2018/MTopJetPostSelection_JMSelec_muon.xml
# sframe_batch.py -s 2018/MTopJetPostSelection_JMSmuon_elec.xml
# python sframe_syst_batch.py 2017/MTopJetPostSelection_JMSelec_SYS_muon.xml
# python sframe_syst_batch.py 2017/MTopJetPostSelection_JMSmuon_SYS_elec.xml
