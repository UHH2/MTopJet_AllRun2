#!/usr/bin/env python

import argparse, sys, os, json
from utils import parallelise

##################################################
#                                                #
#                   MAIN Program                 #
#                                                #
##################################################

masswindows = ['140toInf', '140to200', '140to160', '160to170', '170to180', '180to190', '190to200', '200toInf']
# ptwindows   = [  '0toInf', '400to450', '450to500', '500to600', '600toInf'] # default
# ptwindows   = [  '0toInf', '420to450', '450to500', '500to600', '600toInf'] # default shift
# ptwindows   = [  '0toInf', '400to480', '480to580', '580toInf'] # wide
# ptwindows   = [  '0toInf', '400to430', '430to480', '480to550', '550to650', '650toInf'] # fine
# ptwindows   = [  '0toInf', '400to440', '440to480', '480to550', '550to610', '610to700', '700toInf'] # fine
ptwindows   = [
    '0toInf', '400to450', '450to500', '500to600', '600toInf',
    '420to450',
    '400to480', '480to580', '580toInf',
    '400to440', '440to480', '480to550', '550to610', '610to700', '700toInf'
] # full

study = sys.argv[1]

windows = ptwindows if study=='pt' else masswindows

commands = [['./CreateFSRPlots', study, win] for win in windows]
logs = ['log/'+study+'_'+win+'.txt' for win in windows]
nProcess = len(commands)

for command in commands:
    print command
print logs

os.system('mkdir -p log')
parallelise(commands, nProcess, logs)
