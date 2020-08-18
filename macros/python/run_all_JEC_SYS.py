import subprocess, os, sys, math, contextlib

path = os.path.abspath(os.getcwd())
print path

# list_only    = ["0", "1"]
# list_binning = ["1", "2", "3", "4", "5"]
list_binning = ["1", "2", "5"]
list_years   = ["2016", "2017", "2018", "all"]
# list_peak    = ["0", "20", "40"]


# for year in list_years:
#     for binning in list_binning:
#         for only in list_only:
#             for peak in list_peak:
#                 print "\n ---------------------------------------------------------- "+year+" "+binning+" "+only+" "+peak
#                 os.system("./JEC_SYS "+year+" "+binning+" "+peak+" "+only+" &> /dev/null")

# for year in list_years:
#     for binning in list_binning:
#         for only in list_only:
#                 print "\n ---------------------------------------------------------- "+year+" "+binning+" "+only
#                 os.system("./JEC_SYS "+year+" "+binning+" 0 "+only+" &> /dev/null")

for year in list_years:
    for binning in list_binning:
        print "\n ---------------------------------------------------------- "+year+" "+binning
        os.system("./JEC_SYS "+year+" "+binning+" 0 0 &> /dev/null")
