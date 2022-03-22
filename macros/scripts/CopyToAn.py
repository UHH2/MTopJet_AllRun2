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

FSRsys = ['cor', 'jec', 'jms', 'hdamp', 'isr2', 'tune', 'qcdbased']
hists = {
    'JMSflavor/' : {
        'mjet_flavor_2016.pdf':     'jec/flavor/mjet_flavor_2016.pdf',
        'mjet_flavor_2017.pdf':     'jec/flavor/mjet_flavor_2017.pdf',
        'mjet_flavor_2018.pdf':     'jec/flavor/mjet_flavor_2018.pdf',
        'mjet_flavor_combine.pdf':  'jec/flavor/mjet_flavor_combine.pdf'
    },
    'WMass/' : {
        'wmass_data_hh.pdf':  'jec/wmass/single_years_wmass_data_hh.pdf',
        'wmass_data_hl.pdf':  'jec/wmass/single_years_wmass_data_hl.pdf',
        'wmass_data_lh.pdf':  'jec/wmass/single_years_wmass_data_lh.pdf',
        'wmass_data_ll.pdf':  'jec/wmass/single_years_wmass_data_ll.pdf'
    },
    'FSR/' : {
        'Results_mjet_combine.pdf':               'FSRuncertainty/combine/combine/Results_mjet_combine.pdf',
        'Results_pt_combine.pdf':                 'FSRuncertainty/combine/combine/Results_pt_combine.pdf',
        'Results_mjet_2016.pdf':                  'FSRuncertainty/2016/combine/Results_mjet_2016.pdf',
        'Results_pt_2016.pdf':                    'FSRuncertainty/2016/combine/Results_pt_2016.pdf',
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
    'FSR/SYS/' : {
        'tau32_SYS_PLACEHOLDER_2016_norm.pdf':    'FSRuncertainty/SYS/tau32_SYS_PLACEHOLDER_2016_norm.pdf',
        'tau32_SYS_PLACEHOLDER_combine_norm.pdf': 'FSRuncertainty/SYS/tau32_SYS_PLACEHOLDER_combine_norm.pdf',
    },
    'JetCorrections/Ellipse/' : {
        'Year_comparison_cor.pdf': 'jec/ellipse_year_comparison_cor.pdf',
        'linear_comparison.pdf':       'jec/ellipse_linear_comparison.pdf',
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
    path_in2 = '/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/'
    path_out = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/AN/notes/FullRun2/figs/'

    # Save space for print
    tmp_nfs  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/'
    tmp_afs  = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/AN/notes/FullRun2/figs/'
    tmp_afs2 = '/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/'

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

            cmd_list = [cmd2]
            if 'FSR/SYS' in key1:
                cmd_list = [cmd2.replace('PLACEHOLDER', h) for h in FSRsys]
            if 'JetCorrections' in key1:
                cmd_list = [cmd2.replace(path_in, path_in2)]

            if debug:
                print cmd_list

            for c in cmd_list:
                print_tmp = c.replace(tmp_nfs, '')
                print_cmd = print_tmp.replace(tmp_afs, '')
                print print_cmd
                os.system(c)


# ------------------------------------------------------------------------------
#

if __name__ == "__main__":
    MoveToAN()
