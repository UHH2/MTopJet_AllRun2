LIBRARY := SUHH2MTopJet
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -lSUHH2JetMETObjects

FJINC=$(shell scram tool tag FASTJET INCLUDE)
FJLIB=$(shell scram tool tag FASTJET LIBDIR)

USERLDFLAGS += -L${FJLIB} -lfastjettools -lfastjet -lHOTVR -lNsubjettiness
USERCXXFLAGS := -I${FJINC}

# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
# include ../Makefile.local
include ../Makefile.common
