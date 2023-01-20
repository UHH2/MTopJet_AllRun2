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

config = os.environ["CMSSW_BASE"]+"/src/UHH2/MTopJet/config/"
year = '2016'
internal_option = ""
command = ["python", "sframe_syst_batch.py", ""]
program = "sframe_syst_batch.py"
syst = ["JMS", "JMS_flavor", "JER", "JEC", "COR", "FSR", "ISR"]
syst_uncert = {
    'JMS':        { 'para': 'JetMassScale_direction',  'nominal': 'nominal'},
    'JMS_flavor': { 'para': 'JetMassScale_Flavor',     'nominal': 'nominal'},
    'JEC':        { 'para': 'jecsmear_direction',      'nominal': 'nominal'},
    'JER':        { 'para': 'jersmear_direction',      'nominal': 'nominal'},
    'COR':        { 'para': 'JetCorrection_direction', 'nominal': 'nominal'},
    'FSR':        { 'para': 'PS_variation',            'nominal': 'none'},
    'ISR':        { 'para': 'PS_variation',            'nominal': 'none'},
}
var = ["up", "down"]
ps_weights = ["sqrt2", "2", "4"]
n_jobs = 4 # can be change via shell
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
    create = False

    option = sys.argv[1]
    if option == "submit":
        internal_option = "-s"
        command = ["./syst_submit.sh"]

    elif option == "create" or option == "split":
        internal_option = ""
        command = [""]
        create = True

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
        command = ["rm", "-r"]
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

    # xmls.append("MTopJetPostSelection_tt_muon.xml")
    # xmls.append("MTopJetPostSelection_tt_elec.xml")

    xmls.append("MTopJetPostSelection_SYS_muon.xml")
    xmls.append("MTopJetPostSelection_SYS_elec.xml")

    # xmls.append("MTopJetPostSelection_muon_V20.xml")
    # xmls.append("MTopJetPostSelection_elec_V20.xml")
    # xmls.append("MTopJetPostSelection_JMSmuon_elec.xml")
    # xmls.append("MTopJetPostSelection_JMSelec_muon.xml")

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

    xmls_sys = [s for s in xmls if 'SYS' in s]
    print xmls_sys

    if delete:
        result = pool.map_async(delete_workdirs, job_lists, chunksize=1)
    if create:
        create_workdirs(xmls_sys)
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

        channel = "muon" if "muon" in arg else "elec"

        cmd = command+[year+'_'+channel] # to avoid overlapping commands
        print 'starting',arg,logname
        if not isSYS:
            subprocess.call(['sframe_batch.py', internal_option, year+'/'+arg], stdout=open("./logfiles/log_"+year+"_"+logname+".txt","w"), stderr=subprocess.STDOUT)
        else:
            subprocess.call(cmd, stdout=open("./logfiles/log_"+year+"_"+logname+".txt","w"), stderr=subprocess.STDOUT)

# ------------------------------------------------------------------------------
# ---           Submit                                                       ---
# ------------------------------------------------------------------------------
#
def create_workdirs(dirlist):
    for arg in dirlist:
        print arg
        isSYS = True if 'SYS' in arg else False
        if not isSYS:
            continue

        logname = arg
        logname = logname.replace("MTopJetPostSelection_", "")
        logname = logname.replace(".xml", "")

        cmd = command # to avoid overlapping commands
        print 'starting',arg,logname
        for s in syst:
            isPS = True if 'SR' in s else False
            if year == '2016' and isPS:
                continue
            add = ps_weights if isPS else [""]
            for v in var:
                for a in add:
                    print '\n-ls    ------',s,v,a
                    sys = s+'_'+v
                    direction = v # up down
                    if isPS:
                        sys = s+v+'_'+a
                        direction = sys
                    elif 'JMS' in s and not 'flavor' in s:
                        sys = s+'_'+v+v
                        direction = v+v
                    elif 'JMS' in s and 'flavor' in s:
                        sys = s+'_'+v
                    channel = 'muon' if 'muon' in arg else 'elec'
                    workdir = './SFrameUncerts_'+year+'_'+channel+'/'+sys+'/'
                    cmd = 'mkdir -p '+workdir
                    os.system(cmd)
                    cmd = 'cp '+year+'/MTopJetPostSelection_SYS_'+channel+'.xml '+workdir+'config.xml'
                    os.system(cmd)
                    cmd = 'cp JobConfig.dtd '+workdir
                    os.system(cmd)
                    cmd = "sed -i '/<PARA>/s/<NOM>/<VAR>/g' "+workdir+"config.xml;"
                    cmd += "sed -i 's/<SYST>/"+sys+"/g' "+workdir+"config.xml"
                    cmdtmp = cmd.replace('<PARA>', syst_uncert[s]['para'])
                    cmdtmp2 = cmdtmp.replace('<NOM>', syst_uncert[s]['nominal'])
                    cmd = cmdtmp2.replace('<VAR>', direction)
                    # print cmd
                    os.system(cmd)
                    cmd = "grep -n '"+syst_uncert[s]['para']+"' "+workdir+"config.xml"
                    os.system(cmd)
                    # for test in syst_uncert:
                    #     cmd = "grep -n '"+syst_uncert[test]['para']+"' "+workdir+"config.xml"
                    #     os.system(cmd)
                    cmd = "grep -n '<!ENTITY OUTdir' "+workdir+"config.xml"
                    os.system(cmd)

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

        elif "MTopJetPostSelection_SYS_muon" in arg: Workdir = "SFrameUncerts_"+year+"_muon"
        elif "MTopJetPostSelection_SYS_elec" in arg: Workdir = "SFrameUncerts_"+year+"_elec"
        elif "MTopJetPostSelection_JMSelec_SYS_muon" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetJMSelecmuon"
        elif "MTopJetPostSelection_JMSmuon_SYS_elec" in arg: Workdir = "SFrameUncerts_"+year+"MTopJetJMSmuonelec"

        elif "MTopJetElecSF" in arg: Workdir = "Workdir_ElecSF_"+year

        channel = 'muon' if '_muon' in arg else 'elec'
        dirROOT = '/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/'+channel+'/'
        dirXML = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/config/'

        cmd = command
        if not "SFrameUncerts" in Workdir:
            cmd.append(dirROOT+Workdir)
        else:
            cmd.append(dirROOT+"*/workdir")
        cmd.append(dirXML+Workdir)
        print cmd
        subprocess.call(cmd, stdout=open("logfiles/delete.txt","w"), stderr=subprocess.STDOUT) # Popen runs command in backgorund

# def DeleteWorkdirs(dirlist):
#     for (path, dirs, files) in os.walk(config+"SFrameUncerts "):


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
