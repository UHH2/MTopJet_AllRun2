#!/usr/bin/python

from datetime import datetime
import sys, time, subprocess
#from ROOT import *
from multiprocessing import Pool, Value
import random
import shutil
import glob
import os
import shutil
import json

debug = False

global year
channels = ["muon", "elec"]
years = ["2016", "2017", "2018"]

hists = {
    'JMSflavor/' : {
        'mjet_flavor_2016.pdf':     'jec/flavor/mjet_flavor_2016.pdf',
        'mjet_flavor_2017.pdf':     'jec/flavor/mjet_flavor_2017.pdf',
        'mjet_flavor_2018.pdf':     'jec/flavor/mjet_flavor_2018.pdf',
        'mjet_flavor_combine.pdf':  'jec/flavor/mjet_flavor_combine.pdf'
    },
    'FSR/combine/' : {
        'Bin6.pdf':                 'FSRuncertainty/combine/combine/mjet_140toInf/Bin6.pdf',
        'chi2Cor.pdf':              'FSRuncertainty/combine/combine/mjet_140toInf/chi2Cor.pdf',
        'covData_norm_combine.pdf': 'FSRuncertainty/combine/combine/mjet_140toInf/covData_norm_combine.pdf',
    },
    'FSR/2016/' : {
        'Bin6.pdf':                 'FSRuncertainty/2016/combine/mjet_140toInf/Bin6.pdf',
        'chi2Cor.pdf':              'FSRuncertainty/2016/combine/mjet_140toInf/chi2Cor.pdf',
        'covData_norm_2016.pdf':    'FSRuncertainty/2016/combine/mjet_140toInf/covData_norm_2016.pdf',
    },
    'PaperPlots/' : {
        'tau32_2016.pdf':           'FSRuncertainty/2016/combine/mjet_140toInf/Tau32_norm.pdf',
        'tau32_combine.pdf':        'FSRuncertainty/combine/combine/mjet_140toInf/Tau32_norm.pdf',
        'wmass_hh.pdf':             'Wjet_mass_all_hh_combine_norm.pdf',
        'wmass_hl.pdf':             'Wjet_mass_all_hl_combine_norm.pdf',
        'wmass_lh.pdf':             'Wjet_mass_all_lh_combine_norm.pdf',
        'wmass_ll.pdf':             'Wjet_mass_all_ll_combine_norm.pdf',
        'chi2_JMS.pdf':             'jec/chi2_result.pdf',
    },
}


# ------------------------------------------------------------------------------
# Convert To PDF
def MoveToAN(): # 1: year - 2: channel - 3: histclass - 4: histname

    path_in  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/'
    path_out = "/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/AN/notes/FullRun2/figs/"
    tmp_nfs = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/'
    tmp_afs = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/AN/notes/FullRun2/figs/'

    print "Start moving ... "
    cmd  = 'cp'+' '
    cmd += path_in
    for key1 in hists:
        cmd1  = cmd
        cmd1 += key1
        for key2 in hists[key1]:
            cmd2  = cmd1
            cmd2 += key2+' '
            cmd2 += path_out
            cmd2 += hists[key1][key2]
            print_tmp = cmd2.replace(tmp_nfs, '')
            print_cmd = print_tmp.replace(tmp_afs, '')
            print print_cmd
            os.system(cmd2)


# ------------------------------------------------------------------------------
#

if __name__ == "__main__":
    MoveToAN()
