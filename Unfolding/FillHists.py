#!/usr/bin/python
from datetime import datetime
import sys, multiprocessing, time, subprocess
from ROOT import *

n_jobs = 2


def main():
    global n_jobs                                                       # Needed to modify global copy of n_jobs
    pool_muon = multiprocessing.Pool(processes = n_jobs)
    job_list = [ 'muon', 'elec' ]
    result = pool_muon.map_async(submit_job, job_list)
    pool_muon.close()
    pool_muon.join()
    print 'now combine...'
    subprocess.call(['./combine.sh'])

def submit_job(arg):
    print 'fill hists', arg, 'submitted'
    subprocess.call(['./hist_filler',arg], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
    print 'fill hists', arg, 'finished!'
    return 1


if __name__ == "__main__":
    main()
