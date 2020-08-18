import subprocess, os, sys, math

def main():
    path = "/afs/desy.de/user/p/paaschal/Plots/JEC_SYS/"
    path_latex = path+"Summary/"

    os.chdir(path_latex)

    list_years    = ["2016", "2017", "2018", "all"]
    list_binning  = ["180", "90", "60", "45", "36"]
    list_peak     = ["0", "20", "40"]
    # list_binning  = ["180", "90", "60", "45", "36"]
    # list_peak     = ["0", "20", "30", "40", "50"]
    list_fits     = ["", "_lin"]
    list_fits_str = ["", " only lin"]

    for year in list_years:
        write_latex_frame(path, path_latex, year, list_binning, list_fits, list_fits_str, list_peak)
        os.system("pdflatex Summary_"+year+".tex")
############################################################################################################################
############################################################################################################################
def write_latex_frame(PathPlots, PathLatex, year, list_binning, list_fits, list_fits_str, list_peak):
    print ""
    print "Write Latex"
    print PathPlots+"\n"

    plotwidth = 0.49
    outfile = open(PathLatex+"Summary_"+year+".tex","w")
    outfile.write("\\documentclass[aspectratio=169]{beamer}\n")
    outfile.write("\\usepackage[english]{babel}\n")
    outfile.write("\\usepackage{graphicx} \n")
    outfile.write("\\usepackage{verbatim}\n")
    outfile.write("\\usetheme{Copenhagen}\n")
    outfile.write("\\usecolortheme{beaver}\n")
    outfile.write("\\begin{document}\n")

    ################################################################################# FRAME
    for binning in list_binning:
        for fit in list_fits:
            for peak in list_peak:
                path_peak = "/masspeak/width"+peak
                index = list_fits.index(fit)
                outfile.write("%------------------------------------------------------------\n")
                outfile.write("\\begin{frame}\n")

                ####################################################### TITLE
                outfile.write("\\frametitle{")
                outfile.write(year)
                outfile.write(" - ")
                outfile.write(binning)
                outfile.write(" Bins")
                outfile.write(list_fits_str[index])
                if not peak=="0":
                    outfile.write(" peak "+peak)
                outfile.write("}\n")

                ####################################################### Plots
                outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
                # outfile.write("\\centering\n")
                with open(PathPlots+year+"/btag/rebin"+binning+"/jec_factor"+fit+".txt") as text:
                    for line in text:
                        if "GeV" in line:
                            continue
                        outfile.write(line)
                        outfile.write("\n")
                outfile.write("\n\\end{minipage}\n")
                outfile.write("\hfill\n")
                outfile.write("\\begin{minipage}{"+str(plotwidth)+"\\textwidth}\n")
                outfile.write("\\includegraphics[width=\\textwidth]{")
                outfile.write("../"+year+"/btag/rebin"+binning)
                if peak=="0":
                    outfile.write("/chi2_cont4z_points"+fit+".pdf")
                else:
                    outfile.write(path_peak+"/chi2_cont4z_points"+fit+".pdf")
                outfile.write("}\n\\end{minipage}\n")

                ####################################################### End Frame
                outfile.write("\\end{frame}\n")

    #end document
    outfile.write("\\end{document}")
    outfile.close()

############################################################################################################################
############################################################################################################################
main()
