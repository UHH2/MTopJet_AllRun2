import subprocess, os, sys



################################################################################
################################################################################
def write_latex(PathPlots,list_files,variation,number_frames):
    print PathPlots
    outfile = open(PathPlots+"projections_"+variation+".tex","w")
    outfile.write("\\documentclass[aspectratio=169]{beamer}\n")
    outfile.write("\\usepackage[english]{babel}\n")
    outfile.write("\\usepackage{graphicx} \n")
    outfile.write("\\usecolortheme{beaver}\n")
    outfile.write("\n\\begin{document}\n")

    frame=0
    file=0
    number_files=len(list_files)
    while frame<number_frames:
        # begin new frame
        outfile.write("%--------------------------------------------\n")
        outfile.write("\\begin{frame}\n")
        plotwidth = 0.325
        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")

        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")

        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")

        outfile.write("\\hfill\n") #############################################
        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")

        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")

        outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
        outfile.write("\\begin{figure}\n")
        outfile.write("\\includegraphics[width=\\textwidth]{")
        outfile.write(list_files[file])
        file = file+1
        if file == number_files:
            outfile.write("} \n")
            outfile.write("\\end{figure}\n")
            outfile.write("\\end{minipage}\n")
            # end frame
            outfile.write("\\end{frame}\n")
            break

        outfile.write("} \n")
        outfile.write("\\end{figure}\n")
        outfile.write("\\end{minipage}\n")
        # end frame
        outfile.write("\\end{frame}\n")

        frame = frame +1


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
def switch_element(list):
    i=0
    for file in sorted(list):
        if len(file) == 11:
            new_file = file
            index = list.index(file)
            list.remove(file)
            list.insert(i, new_file)
            i = i+1
    return list

##################################################################################################################### MAIN
path = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"
# year = input("Which year? ")
year=sys.argv[1]
binning=sys.argv[2]
reconst = "/btag/"
# binning = input("Which binning (180, 90, 60, 45, 36,...)? ")
path_full = path+year+reconst+"rebin"+binning+"/projection/"

os.chdir(path_full)

###################################################################################################
## Getting List ###################################################################################
bin=1
bin_exists = False
list_files = sorted(os.listdir("."))
list_JEC = []
list_XC = []
print len(list_files)

list_corrections_str = ["JEC", "XC"]
list_fits_str = ["xy", "xy2", "xyy2", "x2y2", "xx2yy2"]

for file in list_files:
   if "_XC.pdf" in file and not "projections_" in file:
       list_XC.append(file)
   if "_JEC.pdf" in file and not "projections_" in file:
       list_JEC.append(file)

print len(list_XC)
print len(list_JEC)

list_XC = switch_element(list_XC)
list_JEC = switch_element(list_JEC)


number_frames = len(list_XC)/6+1 #JEC same size
print number_frames

write_latex(path_full, list_XC, "XC", number_frames)
write_latex(path_full, list_JEC, "JEC", number_frames)

os.system("pdflatex projections_XC.tex")
os.system("pdflatex projections_JEC.tex")
os.system("rm *.aux *.log *.nav *.out *.snm *.gz *.tex *.toc")
