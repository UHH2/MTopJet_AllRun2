#!/usr/bin/python
from datetime import datetime
import sys, multiprocessing, time, subprocess
from ROOT import *

n_jobs = 6


def main():
    global n_jobs                                                       # Needed to modify global copy of n_jobs
    pool = multiprocessing.Pool(processes = n_jobs)
    job_list = ['data', 'pseudo1', 'pseudo1695', 'pseudo1715', 'pseudo1735', 'pseudo1755'] # Create a List with [1,2,...,n_jobs]
    result = pool.map_async(submit_job, job_list)
    pool.close()
    pool.join()

def submit_job(arg):
    print 'Job', arg, 'submitted'
    channel = 'muon'
    subprocess.call(['./do_unfolding',arg,channel], stdout=open("/dev/null","w"), stderr=subprocess.STDOUT)
    print 'Job', arg, 'finished!'
    return 1


if __name__ == "__main__":
    main()
