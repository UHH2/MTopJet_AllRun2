import subprocess, os, sys, math, contextlib

path = os.path.abspath(os.getcwd())
print path


list_mass    = ["1695 ", "1755 "]
list_uncert  = ["mc", "factor"]
list_tf      = ["0 ", "1 "]

for tf in list_tf:
    for tf in list_tf:
        print "\n ---------------------------------------------------------- "+tf+" "+tf
        os.system("./JEC_SYS combined "+tf+tf)

for mass in list_mass:
    for uncert in list_uncert:
        for tf in list_tf:
            print "\n ---------------------------------------------------------- "+mass+" "+tf+" "+uncert
            # os.system("./JEC_SYS combined "+mass+tf+uncert+" &> /dev/null")
            os.system("./JEC_SYS_mTop combined "+mass+tf+uncert)
