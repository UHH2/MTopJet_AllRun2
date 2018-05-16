#!/usr/bin/python
from datetime import datetime
import sys, multiprocessing, time, subprocess
from ROOT import *

n_jobs = 1
fine = False


def main():
    #try:
    print "number of jobs: ", sys.argv[1]
    global n_jobs                                              # Needed to modify global copy of n_jobs
    n_jobs = int(sys.argv[1])
    global fine
    if sys.argv[2] == "fine":
           fine = True
    else:
           fine = False
    pool = multiprocessing.Pool(processes = n_jobs)
    job_list = [x+1 for x in range(n_jobs)]                    # Create a List with [1,2,...,n_jobs]

    result = pool.map_async(submit_job,job_list)
    pool.close()
    pool.join()

    print 'now sum up root files'    
    merge_files(n_jobs)

def submit_job(i):
    job = str(i)
    njob = str(n_jobs)
    print 'Job ',job, 'submitted'
    if fine:
        subprocess.call(['./hist_filler',njob,job,'fine','ttw'])
    else:
        subprocess.call(['./hist_filler',njob,job,'large','ttw'])
    print 'Job ', i, ' finished!'
    return 1

def merge_files(n_jobs):
    all_hists = range(n_jobs)
    for i in range(0, n_jobs):
        hist_string = "Histograms_"
        hist_string += str(i+1)
        hist_string += ".root "
        all_hists[i]= hist_string
    print all_hists
    if fine:
        subprocess.call(['rm','Histograms_fine.root'])
        subprocess.call(['hadd','Histograms_fine.root']+ all_hists)
    else:
        subprocess.call(['rm','Histograms.root'])
        subprocess.call(['hadd','Histograms.root']+ all_hists)

if __name__ == "__main__":
    main()
