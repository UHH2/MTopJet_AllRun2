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

year = '2016'
internal_option = ""
command = ["python", "sframe_syst_batch.py", ""]
program = "sframe_syst_batch.py"
n_jobs = 8 # can be change via shell
Ndone = 0 # keep track of finished jobs by each worker

# ------------------------------------------------------------------------------
# ---           Main                                                         ---
# ------------------------------------------------------------------------------

def main():
    print "------------------------------------------------"
    print "USAGE: ./SubmitSteer.py <option> <year> (<jobs>)"
    print "  -- change xmls in file"
    print "------------------------------------------------"

    global year; global command; global program; global internal_option; global n_jobs
    delete = False

    option = sys.argv[1]
    if option == "submit":
        internal_option = "-s"
        command = ["python", "sframe_syst_batch.py"]

    elif option == "resubmit":
        internal_option = "-r"
        command = ["./syst_resubmit.sh"]

    elif option == "add" or option == "merge":
        internal_option = "-f"
        command = ["./syst_merge.sh"]

    elif option == "list" or option == "status":
        internal_option = ""
        command = ["./syst_status.sh"]

    elif option == "remove" or option == "delete":
        command = ["./delete_workdir.sh"]
        delete = True

    else:
        internal_option = ""

    year = sys.argv[2]

    if len(sys.argv) == 4:
        n_jobs = int(sys.argv[3])

    pool = Pool(processes = n_jobs)

    print internal_option, year, n_jobs

    xmls = []
    # xmls.append("MTopJetElecSF.xml")

    xmls.append("MTopJetPostSelection_muon.xml")
    xmls.append("MTopJetPostSelection_elec.xml")
    # xmls.append("MTopJetPostSelection_muon_V20.xml")
    # xmls.append("MTopJetPostSelection_elec_V20.xml")
    # xmls.append("MTopJetPostSelection_JMSmuon_elec.xml")
    # xmls.append("MTopJetPostSelection_JMSelec_muon.xml")

    xmls.append("MTopJetPostSelection_SYS_muon.xml")
    xmls.append("MTopJetPostSelection_SYS_elec.xml")
    # xmls.append("MTopJetPostSelection_SYS_muon_V20.xml")
    # xmls.append("MTopJetPostSelection_SYS_elec_V20.xml")
    # xmls.append("MTopJetPostSelection_JMSelec_SYS_muon.xml")
    # xmls.append("MTopJetPostSelection_JMSmuon_SYS_elec.xml")

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
# ---           Submit                                                       ---
# ------------------------------------------------------------------------------
#
def submit_job(dirlist):
    for arg in dirlist:
        isSYS = True if 'SYS' in arg else False

        logname = arg
        logname = logname.replace("MTopJetPostSelection_", "")
        logname = logname.replace(".xml", "")

        cmd = command # to avoid overlapping commands
        print 'starting',arg,logname
        if not isSYS:
            subprocess.call(['sframe_batch.py', internal_option, year+'/'+arg], stdout=open("./logfiles/log_"+year+"_"+logname+".txt","w"), stderr=subprocess.STDOUT)
        else:
            if internal_option == "-s":
                cmd.append(year+'/'+arg)
            else:
                Workdir = ""
                if "MTopJetPostSelection_SYS_muon" in arg: Workdir = year+"MTopJetmuon"
                elif "MTopJetPostSelection_SYS_elec" in arg: Workdir = year+"MTopJetelec"
                elif "MTopJetPostSelection_JMSelec_SYS_muon" in arg: Workdir = year+"MTopJetJMSelecmuon"
                elif "MTopJetPostSelection_JMSmuon_SYS_elec" in arg: Workdir = year+"MTopJetJMSmuonelec"
                cmd.append(Workdir)

            # print cmd
            subprocess.call(cmd, stdout=open("./logfiles/log_"+year+"_"+logname+".txt","w"), stderr=subprocess.STDOUT)

# ------------------------------------------------------------------------------
# ---           Delete                                                       ---
# ------------------------------------------------------------------------------
#
def delete_workdirs(dirlist):
    for arg in dirlist:
        Workdir = ""
        if "MTopJetPostSelection_muon" in arg: Workdir = "Workdir_PostSelMu_"+year
        elif "MTopJetPostSelection_elec" in arg: Workdir = "Workdir_PostSelEl_"+year
        elif "MTopJetPostSelection_JMSmuon_elec" in arg: Workdir = "Workdir_PostSelEl_"+year+"_JMSmuon"
        elif "MTopJetPostSelection_JMSelec_muon" in arg: Workdir = "Workdir_PostSelMu_"+year+"_JMSelec"

        elif "MTopJetPostSelection_SYS_muon" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetmuon"
        elif "MTopJetPostSelection_SYS_elec" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetelec"
        elif "MTopJetPostSelection_JMSelec_SYS_muon" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetJMSelecmuon"
        elif "MTopJetPostSelection_JMSmuon_SYS_elec" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetJMSmuonelec"

        elif "MTopJetElecSF" in arg: Workdir = "Workdir_ElecSF_"+year

        cmd = command
        command.append(Workdir)
        print cmd
        # os.system(cmd)
        subprocess.call(cmd, stdout=open("/test.txt","w"), stderr=subprocess.STDOUT) # Popen runs command in backgorund

# ------------------------------------------------------------------------------
# ---           Creat list                                                   ---
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
# ---           Init                                                         ---
# ------------------------------------------------------------------------------
# needed to setup global counter
def init(args):
    global Ndone
    Ndone = args
    # global year

# ------------------------------------------------------------------------------
# ---           START                                                        ---
# ------------------------------------------------------------------------------
#
if __name__ == "__main__":
    main()
