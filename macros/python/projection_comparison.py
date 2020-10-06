import subprocess, os, sys, math

################################################################################
################################################################################
def write_latex(bins,number_pages,ptbin):
    path      = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/mc_data_uncert/pt_bins/"
    path_add  = "/btag/rebin180/masspeak/projection"
    path_nom  = "combined"+path_add
    path_1695 = "1695/combined"+path_add
    path_1755 = "1755/combined"+path_add

    outfile = open(path+"compare_all_projections_"+ptbin+".tex","w")

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

    line_count=0
    for bin in bins:
        if line_count == 0:
            outfile.write("\\begin{center}\n")
            outfile.write("\\begin{tabular}{c|c|c|||c|c|c}\n")
            outfile.write("&\n")
            outfile.write("\\begin{large}\n")
            outfile.write("JEC\n")
            outfile.write("\\end{large}&\n")
            outfile.write("&\n")
            outfile.write("&\n")
            outfile.write("\\begin{large}\n")
            outfile.write("XCone\n")
            outfile.write("\\end{large}&\n")
            outfile.write("\\\\\n")
            outfile.write("\\hline\n")
            outfile.write("1695&\n")
            outfile.write("nominal&\n")
            outfile.write("1755&\n")
            outfile.write("1695&\n")
            outfile.write("nominal&\n")
            outfile.write("1755\\\\\n")
            outfile.write("\\hline\n")
            outfile.write("\\hline\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_1695+"/Bin"+bin+"_JEC_"+ptbin+".pdf}&\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_nom+"/Bin"+bin+"_JEC_"+ptbin+".pdf}&\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_1755+"/Bin"+bin+"_JEC_"+ptbin+".pdf}&\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_1695+"/Bin"+bin+"_XC_"+ptbin+".pdf}&\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_nom+"/Bin"+bin+"_XC_"+ptbin+".pdf}&\n")
        outfile.write("\\includegraphics[width=0.15\linewidth]{"+path_1755+"/Bin"+bin+"_XC_"+ptbin+".pdf}\\\\\n")
        line_count+=1
        if line_count==8:
            line_count=0
            outfile.write("\\end{tabular}\n")
            outfile.write("\\end{center}\n")
            outfile.write("\\newpage\n")

    if line_count < 8:
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
path      = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/mc_data_uncert/pt_bins/"
path_add  = "/btag/rebin180/masspeak/projection"
variation = ["1695/combined", "combined", "1755/combined"]
paths     = []
for mass in variation:
    paths.append(path+mass+path_add)

variation = ["JEC", "XC"]
ptbins    = ["hh", "hl", "lh", "ll"]

###################################################################################################
## Getting List ###################################################################################
bins = []
bin=60
while bin<110:
    bins.append(str(bin))
    bin+=1
print len(bins)

os.chdir(path)

list_bin = []
for ptbin in ptbins:
    list_files = sorted(os.listdir(paths[1]))
    for file in list_files:
        for bin in bins:
            if bin in file:
                if "JEC" in file and ptbin in file: # Only once, XC is equal
                    list_bin.append(str(bin))
    write_latex(list_bin, 4, ptbin)
    os.system("pdflatex "+path+"compare_all_projections_"+ptbin+".tex")
    list_bin = []

os.system("rm *.gz *.rar *.aux *.tex *.log")

# os.system("pdflatex .tex")
