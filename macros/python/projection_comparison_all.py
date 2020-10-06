#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess, os, sys, math

################################################################################
################################################################################
def write_latex(list_JEC, list_XC, number_files, ptbin):
    path      = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/uncert/pt_bins/combined/"

    outfile = open(path+"compare_all_projections"+ptbin+".tex","w")

    print ""
    print "Write Latex"

    # begin document
    outfile.write("\\documentclass{article}\n")
    outfile.write("\\usepackage[margin=0pt]{geometry}\n")
    outfile.write("\\usepackage{graphicx}\n")
    outfile.write("\\usepackage[english,german]{babel}\n")
    outfile.write("\\usepackage{multicol}\n")
    outfile.write("\\usepackage{graphicx}\n")
    outfile.write("\\pagestyle{empty}\n")
    outfile.write("\\begin{document}\n")

    line_count  = 0
    files_count = 0
    print number_files
    while files_count<number_files:
        print files_count
        if line_count == 0:
            outfile.write("\\begin{center}\n")
            outfile.write("\\begin{tabular}{c||c}\n")
            outfile.write("\\begin{large}\n")
            outfile.write("JEC\n")
            outfile.write("\\end{large}&\n")
            outfile.write("\\begin{large}\n")
            outfile.write("XCone\n")
            outfile.write("\\end{large}\\\\\n")
            outfile.write("\\hline\n")
        outfile.write("\\includegraphics[width=0.3\linewidth]{"+ptbin+"/"+list_JEC[files_count]+"}&\n")
        outfile.write("\\includegraphics[width=0.3\linewidth]{"+ptbin+"/"+list_XC[files_count]+"}\\\\\n")
        line_count+=1
        if line_count==4:
           line_count=0
           outfile.write("\\end{tabular}\n")
           outfile.write("\\end{center}\n")
           outfile.write("\\newpage\n")
        files_count+=1

    if line_count < 4:
       outfile.write("\\end{tabular}\n")
       outfile.write("\\end{center}\n")

    outfile.write("\\end{document}\n")
    outfile.close()

################################################################################
################################################################################
def insert_dash(string, insert, index):
    return string[:index] + insert + string[index:]

################################################################################
################################################################################
def insert_zero(list):
    for file in list:
        if len(file) == 11:
            new_file = insert_dash(file, "0", 3)
            index = list.index(file)
            list[index]=new_file
    return list

################################################################################
################################################################################
def switch_element(list):
    i=0
    len_switch = 12
    if isXC:
        len_switch = 11
    for file in sorted(list):
        if len(file)==len_switch:
            new_file = file
            index = list.index(file)
            list.remove(file)
            list.insert(i, new_file)
            i = i+1
    return list

##################################################################################################################### MAIN
path      = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/uncert/pt_bins/combined/"
ptbins    = ["_hh", "_hl", "_lh", "_ll"]
paths     = []
for ptbin in ptbins:
    paths.append(path+ptbin)

variation = ["JEC", "XC"]

###################################################################################################
## Getting List ###################################################################################
bins = []
bin=60
while bin<110:
    bins.append(str(bin))
    bin+=1
print len(bins)

os.chdir(path)

list_JEC = []
list_XC  = []

for ptbin in ptbins:
    path_full = path+ptbin
    print path_full
    list_files = sorted(os.listdir(path_full))
    for file in list_files:
       if "JEC" in file:
          list_JEC.append(file)
       if "XC" in file:
          list_XC.append(file)
    print len(list_JEC)
    print len(list_XC)

    write_latex(list_JEC, list_XC, len(list_JEC), ptbin)
    os.system("pdflatex "+path+"compare_all_projections"+ptbin+".tex")
    list_JEC = []
    list_XC = []


os.system("rm *.gz *.rar *.aux *.tex *.log")

# os.system("pdflatex .tex")
