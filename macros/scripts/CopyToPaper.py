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

hists = {
        'tau32_2016.pdf':           'fsr_tau32_2016.pdf',
        'tau32_combine.pdf':        'fsr_tau32_combine.pdf',
        'wmass_hh.pdf':             'Wjet_mass_all_hh_combine_norm.pdf',
        'wmass_hl.pdf':             'Wjet_mass_all_hl_combine_norm.pdf',
        'wmass_lh.pdf':             'Wjet_mass_all_lh_combine_norm.pdf',
        'wmass_ll.pdf':             'Wjet_mass_all_ll_combine_norm.pdf',
        'chi2_JMS.pdf':             'chi2_JMS.pdf',
}


# ------------------------------------------------------------------------------
# Convert To PDF
def MoveToPaper(): # 1: year - 2: channel - 3: histclass - 4: histname

    path_in  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/'
    path_out = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/Papers/TOP-21-012/figs/'
    # Save space for print
    tmp_nfs  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots'
    tmp_afs  = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/Papers'

    print "Start moving ... "
    cmd  = 'cp'+' '
    cmd += path_in
    for key1 in hists:

        cmd1  = cmd
        cmd1 += key1+' '
        cmd1 += path_out
        cmd1 += hists[key1]

        print_tmp = cmd1.replace(tmp_nfs, '')
        print_cmd = print_tmp.replace(tmp_afs, '')
        print print_cmd
        os.system(cmd1)


# ------------------------------------------------------------------------------
#

if __name__ == "__main__":
    MoveToPaper()
