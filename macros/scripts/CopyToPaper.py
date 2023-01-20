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

hists_paper = {
        'tau32_2016.pdf':           'fsr_tau32_2016.pdf',
        'tau32_combine.pdf':        'fsr_tau32_combine.pdf',
        'wmass_hh.pdf':             'Wjet_mass_all_hh_combine_norm.pdf',
        'wmass_hl.pdf':             'Wjet_mass_all_hl_combine_norm.pdf',
        'wmass_lh.pdf':             'Wjet_mass_all_lh_combine_norm.pdf',
        'wmass_ll.pdf':             'Wjet_mass_all_ll_combine_norm.pdf',
        'chi2_JMS.pdf':             'chi2_JMS.pdf',
}

hists_journal = {
        'tau32_2016.pdf':           'Figure_008-a.pdf',
        'tau32_combine.pdf':        'Figure_008-b.pdf',
        'wmass_hh.pdf':             'Figure_004-d.pdf',
        'wmass_hl.pdf':             'Figure_004-c.pdf',
        'wmass_lh.pdf':             'Figure_004-b.pdf',
        'wmass_ll.pdf':             'Figure_004-a.pdf',
        'chi2_JMS.pdf':             'Figure_005.pdf',
}

# ------------------------------------------------------------------------------
# Convert To PDF
def MoveToPaper(isPaper, isPre, isEPJC): # 1: year - 2: channel - 3: histclass - 4: histname

    if isEPJC:
        print 'Set paper folder by default'
        isPaper = 1

    folder = "Papers" if isPaper else "PAS"
    prelim = "preliminary/" if isPre and isPaper else ""
    path_in  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/'+prelim
    path_out_paper = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/'+folder+'/TOP-21-012/figs/'+prelim
    path_out_epjc = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/'+folder+'/TOP-21-012/'
    path_out = path_out_epjc if isEPJC else path_out_paper
    hists = hists_journal if isEPJC else hists_paper
    # Save space for print
    tmp_nfs  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots'
    tmp_afs  = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/'+folder

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
    print("Copy to Paper")
    MoveToPaper(1, 0, 1) # Paper
    # print("Copy preliminary to Paper")
    # MoveToPaper(1, 1) # Paper Prelim
    # print("Copy to PAS")
    # MoveToPaper(0, 0) # PAS
