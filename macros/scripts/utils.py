import sys
import os
import subprocess
import time
import glob

# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts))
        else:
            print '%r  %2.2f s' % \
                  (method.__name__, (te - ts))
        return result
    return timed

# ---------------------------------------------------------------------------
#
# ---------------------------------------------------------------------------

def parallelise(list_processes, MaxProcess=10, list_logfiles=[], cwd=None, time_=2):
  ntotal = len(list_processes)
  processes = []
  logfiles = []
  condition = len(list_logfiles)>0 and len(list_logfiles)==len(list_processes)
  for index, process in enumerate(list_processes):
    wait = True
    while wait:
      nrunning = 0
      ncompleted = 0
      idx=0
      for proc in processes:
        if proc.poll() == None :
          nrunning += 1
        else:
          ncompleted += 1
          if not logfiles[idx].closed:
            logfiles[idx].close()
            print 'Job "%s" has finished.' % logfiles[idx].name
        idx += 1
      if nrunning >= MaxProcess:
        percentage = float(ncompleted)/float(ntotal)*100
        sys.stdout.write( 'Already completed '+str(ncompleted)+' out of '+str(ntotal)+' jobs --> '+str(percentage)+'%. Currently running: '+str(nrunning)+' \r')
        sys.stdout.flush()
        time.sleep(5)
      else:
        print 'only %i jobs are running, going to spawn new ones.' % nrunning
        wait = False
    if condition:
      f = open(list_logfiles[index],'w')
    else:
      f = open("log_"+str(index)+".txt",'w')
    logfiles.append(f)
    if cwd:
      processes.append(subprocess.Popen(process[1:], stdout=f, cwd=process[0]))
      time.sleep(time_)
    else:
      processes.append(subprocess.Popen(process, stdout=f))

  for proc in processes:
    proc.wait()
  for file in logfiles:
    file.close()
  if not condition:
    # os.remove("log.txt")
    a = map(os.remove, glob.glob("log_*.txt"))
