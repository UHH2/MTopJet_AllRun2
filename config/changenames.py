#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, time, subprocess
import os

def main():
    path = os.getcwd()
    folder = sys.argv[1]
    os.chdir('./'+folder)
    print folder
    print os.getcwd()
    for file in sorted(os.listdir(os.getcwd())): #lists all files in dir
        if not "Alex" in file:
            print 'mv '+file+' /Dennis'
            os.system('mv '+file+' ./Dennis')
        if "Alex" in file:
            newname = file.replace("Alex_","")
            os.rename(file,newname)
    os.chdir(path)

main()
