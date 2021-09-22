#!/usr/bin/python

from datetime import datetime
import sys, time, subprocess
from ROOT import *
from multiprocessing import Pool, Value
import random
import shutil
import glob
import os
import shutil

year = "2016"
channel = "muon"

def main():

    global year; global channel;
    year = sys.argv[1]
    channel = sys.argv[2]

    v="" # 2018
    if year == "2016": v = "v3"
    if year == "2017": v = "v2"

    print v

    root = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/"+channel+"/uhh2.AnalysisModuleRunner.MC."
    processes = ["DYJets_", "DiBoson_", "QCD_", "SingleTop_", "TTbar_", "WJets_"]

    os.chdir("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/scripts")
    current_path = os.getcwd()
    print current_path

    job_list = []
    command = "./create-btag-effiency-histos.py "
    for p in processes:
        command += root+p+year+v+".root "

    os.system(command)

    os.system("mv BTagMCEfficiencyHists.root ../MTopJet/btag_effi/BTagMCEfficiencyHists_"+year+"_"+channel+".root")



# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
