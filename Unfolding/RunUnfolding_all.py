#!/usr/bin/python
from datetime import datetime
import sys, multiprocessing, time, subprocess
from ROOT import *

n_jobs = 6


def main():
    global n_jobs                                                       # Needed to modify global copy of n_jobs
    pool_muon = multiprocessing.Pool(processes = n_jobs)
    job_list = [
        ['data', 'muon'],
        ['pseudo1', 'muon'],
        ['pseudo1695', 'muon'],
        ['pseudo1715', 'muon'],
        ['pseudo1735', 'muon'],
        ['pseudo1755', 'muon'],
        ['data', 'elec'],
        ['pseudo1', 'elec'],
        ['pseudo1695', 'elec'],
        ['pseudo1715', 'elec'],
        ['pseudo1735', 'elec'],
        ['pseudo1755', 'elec'],
        ['data', 'combine'],
        ['pseudo1', 'combine'],
        ['pseudo1695', 'combine'],
        ['pseudo1715', 'combine'],
        ['pseudo1735', 'combine'],
        ['pseudo1755', 'combine']
    ]
    result = pool_muon.map_async(submit_job, job_list)
    pool_muon.close()
    pool_muon.join()

def submit_job(arg):
    print 'Job', arg[0], arg[1], 'submitted'
    subprocess.call(['./do_unfolding',arg[0],arg[1]], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
    print 'Job', arg[0], arg[1], 'finished!'
    return 1


if __name__ == "__main__":
    main()
