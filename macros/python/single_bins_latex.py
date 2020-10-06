#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess, os, sys, math

################################################################################
################################################################################
def write_latex(PathPlots,list_files,variation,number_pages,addition):
    print ""
    print "Write Latex"
    print PathPlots+"\n"
    outfile = open(PathPlots+"/projections_"+variation+addition+".tex","w")
    outfile.write("\\documentclass{article}\n")
    outfile.write("\usepackage[margin=0pt]{geometry}")
    outfile.write("\\usepackage{graphicx} \n")
    outfile.write("\pagestyle{empty}\n")
    outfile.write("\n\\begin{document}\n")

    file=0
    page=0
    vertical_count = 0
    horizontal_count = 0
    plotwidth = 0.19
    number_files=len(list_files)
    while page<number_pages:
        while vertical_count<7:
            while horizontal_count<5:
                outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
                outfile.write("\\includegraphics[width=\\textwidth]{projection/")
                outfile.write(list_files[file])
                outfile.write("}\n\\end{minipage}\n")
                file = file+1
                if file == number_files:
                    break
                horizontal_count = horizontal_count+1

            if file == number_files:
                break
            outfile.write("\\vfill\n")
            horizontal_count = 0
            vertical_count = vertical_count+1

        if file == number_files:
            break
        vertical_count = 0
        page = page+1
    #end document
    outfile.write("\\end{document}")
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
def switch_element(list, isXC):
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
path_full = sys.argv[1]
addition  = sys.argv[2]
folder    = sys.argv[2]
print path_full
os.chdir(path_full)

###################################################################################################
## Getting List ###################################################################################
bin=1
bin_exists = False
list_files = sorted(os.listdir("./projection"))
list_JEC = []
list_XC = []
print len(list_files)

list_corrections_str = ["JEC", "XC"]
list_fits_str = ["xy", "xy2", "xyy2", "x2y2", "xx2yy2"]

for file in list_files:
   if "_XC"+addition+".pdf" in file:
       list_XC.append(file)
   if "_JEC"+addition+".pdf" in file:
       list_JEC.append(file)

print len(list_XC)
print len(list_JEC)

print("switch element")
list_XC  = switch_element(list_XC, True)
list_JEC = switch_element(list_JEC, False)

number_pages = math.ceil(len(list_XC)/35)+1 #JEC same size

write_latex(path_full, list_XC, "XC", number_pages, addition)
write_latex(path_full, list_JEC, "JEC", number_pages, addition)

os.system("pdflatex projections_XC"+addition+".tex")
os.system("pdflatex projections_JEC"+addition+".tex")
os.system("rm *.aux *.log *.nav *.gz *.out *.snm *.tex *.toc")
