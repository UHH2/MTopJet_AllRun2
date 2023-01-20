import glob, os, json, sys, math
import ROOT as rt
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)
from array import array
import numpy as np
from tdrstyle_all import *
import tdrstyle_all as TDR

dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/"

muon = "muon"
elec = "elec"

channels = [muon, elec]

processes = {
    'data' : "uhh2.AnalysisModuleRunner.DATA.DATA",
    'ttbar' : "uhh2.AnalysisModuleRunner.MC.TTbar",
    'wjets' : "uhh2.AnalysisModuleRunner.MC.WJets",
    'st' : "uhh2.AnalysisModuleRunner.MC.SingleTop",
    'other' : "uhh2.AnalysisModuleRunner.MC.other",
}

years = {
    '2016': '_2016v3',
    '2017': '_2017v2',
    '2018': '_2018',
}


print '+'+'-'*50+'+'
print '| Usage: python EfficiencyEPJC.py <year> <channel> |'
if len(sys.argv)!=3:
    print '<E> I need a year an a channel!'
    print '+'+'-'*50+'+'
    sys.exit()
else:
    year = sys.argv[1]
    print '| --year:    {:<37} |'.format(sys.argv[1])
    print '| --channel: {:<37} |'.format(sys.argv[2])
    print '+'+'-'*50+'+'

year = sys.argv[1]
channel = sys.argv[2]

files = {
    p: {
        c: {
            y: rt.TFile(dir+c+"/"+processes[p]+years[y]+".root") for y in years
        } for c in channels
    } for p in processes
}

hists = {
    'all_gen': 'efficiency_all_gen/n_events',
    'masscut_gen': 'efficiency_masscut_gen/n_events',
    'lep_pt_gen': 'efficiency_lep_pt_gen/n_events',
    'all': 'efficiency_all/n_events',
    'masscut': 'efficiency_masscut/n_events',
    'lep_pt': 'efficiency_lep_pt/n_events',
}

h_efficiency = {p: {c: {y: {} for y in years} for c in channels} for p in files}
efficiency = {p: {c: {y: {} for y in years} for c in channels} for p in files}
for p in files:
    for c in channels:
        for y in years:
            for h in hists:
                if p!='ttbar' and 'gen' in h: continue
                h_efficiency[p][c][y][h] = rt.TH1F(files[p][c][y].Get(hists[h]))
                h_efficiency[p][c][y][h].SetDirectory(0)
                efficiency[p][c][y][h] = h_efficiency[p][c][y][h].GetBinContent(1)

for p in files:
    for c in channels:
        efficiency[p][c]['combine'] = {}
        for h in hists:
            if 'ttbar'!=p and 'gen' in h: continue
            eff = efficiency[p][c]['2016'][h]+efficiency[p][c]['2017'][h]+efficiency[p][c]['2018'][h]
            efficiency[p][c]['combine'][h] = eff

years['combine'] = 'combine'
channels.append('combine')
for p in files:
    efficiency[p]['combine'] = {}
    for y in years:
        efficiency[p]['combine'][y] = {}
        for h in hists:
            if 'ttbar'!=p and 'gen' in h: continue
            eff = efficiency[p]['muon'][y][h]+efficiency[p]['elec'][y][h]
            efficiency[p]['combine'][y][h] = eff


for y in years:
    if not 'combine' in y: continue
    print ''
    print 'channels combine - year',y,' - Reco Level'
    print '{:<15}|{:8}|{:>17}|{:>17}'.format('','All cuts', 'no mass cut', 'no lep pt cut')
    print '-'*61
    for p in files:
        print '{:<15}|{:8.0f}|{:8.0f} ({:3.2f}%)|{:8.0f} ({:3.2f}%)'.format(p, efficiency[p]['combine'][y]['all'], efficiency[p]['combine'][y]['masscut'], (efficiency[p]['combine'][y]['all']/efficiency[p]['combine'][y]['masscut']*100), efficiency[p]['combine'][y]['lep_pt'], (efficiency[p]['combine'][y]['all']/efficiency[p]['combine'][y]['lep_pt']*100))

print ''
print 'channels combine - only ttbar - Gen Level'
print '{:<15}|{:8}|{:>17}|{:>17}'.format('','All cuts', 'no mass cut', 'no lep pt cut')
print '-'*61
for y in years:
    print '{:<15}|{:8.0f}|{:8.0f} ({:3.2f}%)|{:8.0f} ({:3.2f}%)'.format(y, efficiency['ttbar']['combine'][y]['all_gen'], efficiency['ttbar']['combine'][y]['masscut_gen'], (efficiency['ttbar']['combine'][y]['all_gen']/efficiency['ttbar']['combine'][y]['masscut_gen']*100), efficiency['ttbar']['combine'][y]['lep_pt_gen'], (efficiency['ttbar']['combine'][y]['all_gen']/efficiency['ttbar']['combine'][y]['lep_pt_gen']*100))
