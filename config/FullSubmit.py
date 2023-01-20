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
internal_option = ""
command = "python"
program = "sframe_syst_batch.py"
n_jobs = 6
Ndone = 0 # keep track of finished jobs by each worker

def main():
    print "USAGE: ./steer.py <option> <year> (<jobs>)"

    global year; global command; global program; global internal_option; global n_jobs
    delete = False

    option = sys.argv[1]
    if option == "submit":
        internal_option = "-s"
        command = "python"
        program = "sframe_syst_batch.py"
    elif option == "resubmit":
        internal_option = "-r"
        command = "source"
        program = "syst_resubmit.sh"
    elif option == "add" or option == "merge":
        internal_option = "-f"
        command = "source"
        program = "syst_merge.sh"
    elif option == "list" or option == "status":
        internal_option = "-l"
        command = "source"
        program = "syst_status.sh"
    elif option == "remove" or option == "delete":
        delete = True
    else:
        internal_option = ""

    year = sys.argv[2]

    if len(sys.argv) == 4:
        n_jobs = int(sys.argv[3])

    pool = Pool(processes = n_jobs)


    xmls = []
    xmls.append("MTopJetPostSelection_muon.xml")
    xmls.append("MTopJetPostSelection_elec.xml")
    xmls.append("MTopJetPostSelection_SYS_muon.xml")
    xmls.append("MTopJetPostSelection_SYS_elec.xml")

    xmls.append("MTopJetPostSelection_JMSmuon_elec.xml")
    xmls.append("MTopJetPostSelection_JMSelec_muon.xml")
    xmls.append("MTopJetPostSelection_JMSelec_SYS_muon.xml")
    xmls.append("MTopJetPostSelection_JMSmuon_SYS_elec.xml")

    print "----------------------"
    print "  -- start script... "
    print "  -- submit", year
    print "  -- njobs", n_jobs


    # setup global counter
    Ndone = Value('i', 0)

    # create lists for every worker including the xmls that should be processed
    # also pair it with index of worker (for keeping track of done jobs)
    job_lists = []
    for i in range(n_jobs):
        job_lists.append( create_job_list(i, n_jobs, xmls) )

    print(job_lists)

    if delete:
        result = pool.map_async(delete_workdirs, job_lists, chunksize=1)
    else:
        result = pool.map_async(submit_job, job_lists, chunksize=1)

    pool.close()
    pool.join()

# ------------------------------------------------------------------------------
#
def submit_job(dirlist):
    for arg in dirlist:
        isSYS = False
        if 'SYS' in arg: isSYS = True

        # print 'starting',arg
        if not isSYS:
            # print 'sframe_batch', internal_option, year+'/'+arg
            subprocess.call(['sframe_batch', internal_option, year+'/'+arg], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
        else:
            # print command, program, year+'/'+arg
            subprocess.call([command, program, year+'/'+arg], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
    return 1

def delete_workdirs(dirlist):
    dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/"
    for arg in dirlist:
        isSYS = False
        if 'SYS' in arg: isSYS = True

        Workdir = ""
        channel = ""
        if "MTopJetPostSelection_muon" in arg:
            Workdir = "Workdir_PostSelMu_"+year
            channel = "muon"
        elif "MTopJetPostSelection_elec" in arg:
            channel = "elec"
            Workdir = "Workdir_PostSelEl_"+year
        elif "MTopJetPostSelection_JMSmuon_elec" in arg:
            channel = "elec_muonJMS"
            Workdir = "Workdir_PostSelEl_"+year+"_JMSmuon"
        elif "MTopJetPostSelection_JMSelec_muon" in arg:
            channel = "muon_elecJMS"
            Workdir = "Workdir_PostSelEl_"+year+"_JMSelec"

        elif "MTopJetPostSelection_SYS_muon" in arg:
            channel = "muon"
            Workdir = "SFrameUncerts_"+year+"MTopJetmuon"
        elif "MTopJetPostSelection_SYS_elec" in arg:
            channel = "elec"
            Workdir = "SFrameUncerts_"+year+"MTopJetelec"
        elif "MTopJetPostSelection_JMSelec_SYS_muon" in arg:
            channel = "muon_elecJMS"
            Workdir = "SFrameUncerts_"+year+"MTopJetJMSelecmuon"
        elif "MTopJetPostSelection_JMSmuon_SYS_elec" in arg:
            channel = "elec_muonJMS"
            Workdir = "SFrameUncerts_"+year+"MTopJetJMSmuonelec"

        subprocess.call(['rm', '-rf', Workdir], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
        if not isSYS: subprocess.call("rm -r "+dir+channel+"/"+Workdir);

    return 1



# ------------------------------------------------------------------------------
# for every worker an array containing different dirs is created
# it is taken care of distributing the job as equally as possible
def create_job_list(i, nworkers, dirs):
        integer_div = len(dirs) // nworkers
        rest = len(dirs) - (integer_div * nworkers)
        new_list = []
        # the rest should not just be added to the last worker but distributed equally
        # therefore every worker with i+1 gets one additional job
        if i+1 <= rest:
            min_index = i*integer_div + i
            max_index = (i+1)*integer_div + 1 + i
        else:
            min_index = i*integer_div + rest
            max_index = (i+1)*integer_div + rest

        # loop over total list and only keep element if index is between min and max
        for index, val in enumerate(dirs):
            if index >= min_index and index < max_index:
                new_list.append( dirs[index] )

        return new_list

# ------------------------------------------------------------------------------
# needed to setup global counter
def init(args):
    global Ndone
    Ndone = args
    # global year

# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
