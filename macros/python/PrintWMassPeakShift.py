import glob, os, ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
from array import array
import numpy as np

dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/"

muon = "muon"
elec = "elec"

data_f   = "uhh2.AnalysisModuleRunner.DATA.DATA";
ttbar_f  = "uhh2.AnalysisModuleRunner.MC.TTbar";
wjets_f  = "uhh2.AnalysisModuleRunner.MC.WJets";
st_f     = "uhh2.AnalysisModuleRunner.MC.SingleTop";
other_f  = "uhh2.AnalysisModuleRunner.MC.other";
uhh_mc   = "uhh2.AnalysisModuleRunner.MC.";

jec_up   = "JEC_up";
jec_down = "JEC_down";
cor_up   = "COR_up";
cor_down = "COR_down";

y16 = "_2016v3";
y17 = "_2017v2";
y18 = "_2018";

f_data = {
    "muon": {
        "2016": dir+muon+"/"+data_f+y16+".root",
        "2017": dir+muon+"/"+data_f+y17+".root",
        "2018": dir+muon+"/"+data_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+data_f+y16+".root",
        "2017": dir+elec+"/"+data_f+y17+".root",
        "2018": dir+elec+"/"+data_f+y18+".root",
    },
}

f_ttbar = {
    "muon": {
        "2016": dir+muon+"/"+ttbar_f+y16+".root",
        "2017": dir+muon+"/"+ttbar_f+y17+".root",
        "2018": dir+muon+"/"+ttbar_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+ttbar_f+y16+".root",
        "2017": dir+elec+"/"+ttbar_f+y17+".root",
        "2018": dir+elec+"/"+ttbar_f+y18+".root",
    },
}

f_wjets = {
    "muon": {
        "2016": dir+muon+"/"+wjets_f+y16+".root",
        "2017": dir+muon+"/"+wjets_f+y17+".root",
        "2018": dir+muon+"/"+wjets_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+wjets_f+y16+".root",
        "2017": dir+elec+"/"+wjets_f+y17+".root",
        "2018": dir+elec+"/"+wjets_f+y18+".root",
    },
}

f_st = {
    "muon": {
        "2016": dir+muon+"/"+st_f+y16+".root",
        "2017": dir+muon+"/"+st_f+y17+".root",
        "2018": dir+muon+"/"+st_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+st_f+y16+".root",
        "2017": dir+elec+"/"+st_f+y17+".root",
        "2018": dir+elec+"/"+st_f+y18+".root",
    },
}

f_other = {
    "muon": {
        "2016": dir+muon+"/"+other_f+y16+".root",
        "2017": dir+muon+"/"+other_f+y17+".root",
        "2018": dir+muon+"/"+other_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+other_f+y16+".root",
        "2017": dir+elec+"/"+other_f+y17+".root",
        "2018": dir+elec+"/"+other_f+y18+".root",
    },
}

f_jec_up = {
    "muon": {
        "2016": dir+muon+"/"+jec_up+"/"+ttbar_f+y16+".root",
        "2017": dir+muon+"/"+jec_up+"/"+ttbar_f+y17+".root",
        "2018": dir+muon+"/"+jec_up+"/"+ttbar_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+jec_up+"/"+ttbar_f+y16+".root",
        "2017": dir+elec+"/"+jec_up+"/"+ttbar_f+y17+".root",
        "2018": dir+elec+"/"+jec_up+"/"+ttbar_f+y18+".root",
    },
}

f_jec_down = {
    "muon": {
        "2016": dir+muon+"/"+jec_down+"/"+ttbar_f+y16+".root",
        "2017": dir+muon+"/"+jec_down+"/"+ttbar_f+y17+".root",
        "2018": dir+muon+"/"+jec_down+"/"+ttbar_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+jec_down+"/"+ttbar_f+y16+".root",
        "2017": dir+elec+"/"+jec_down+"/"+ttbar_f+y17+".root",
        "2018": dir+elec+"/"+jec_down+"/"+ttbar_f+y18+".root",
    },
}

f_cor_up = {
    "muon": {
        "2016": dir+muon+"/"+cor_up+"/"+ttbar_f+y16+".root",
        "2017": dir+muon+"/"+cor_up+"/"+ttbar_f+y17+".root",
        "2018": dir+muon+"/"+cor_up+"/"+ttbar_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+cor_up+"/"+ttbar_f+y16+".root",
        "2017": dir+elec+"/"+cor_up+"/"+ttbar_f+y17+".root",
        "2018": dir+elec+"/"+cor_up+"/"+ttbar_f+y18+".root",
    },
}

f_cor_down = {
    "muon": {
        "2016": dir+muon+"/"+cor_down+"/"+ttbar_f+y16+".root",
        "2017": dir+muon+"/"+cor_down+"/"+ttbar_f+y17+".root",
        "2018": dir+muon+"/"+cor_down+"/"+ttbar_f+y18+".root",
    },
    "elec": {
        "2016": dir+elec+"/"+cor_down+"/"+ttbar_f+y16+".root",
        "2017": dir+elec+"/"+cor_down+"/"+ttbar_f+y17+".root",
        "2018": dir+elec+"/"+cor_down+"/"+ttbar_f+y18+".root",
    },
}
