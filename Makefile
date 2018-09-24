LIBRARY := SUHH2MTopJet
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -lSUHH2JetMETObjects

# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
# include ../Makefile.local
include ../Makefile.common
