import glob, os, json, sys, math
import ROOT as rt
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)
from array import array
import numpy as np
from tdrstyle_all import *
import tdrstyle_all as TDR

class Ratios:
    def __init__(self, g1, g2, equal):
        self.g1 = g1
        self.g2 = g2
        self.equal = equal

    def GetRatioGraph(self):
        ratio = self.g1.Clone()
        ratio.SetDirectory(0)
        nbins = self.g1.GetNbinsX()
        nbins2 = self.g2.GetNbinsX()
        for i in range(nbins):
            pt = 0
            N1 = self.g1.GetBinContent(i+1)
            N2 = self.g2.GetBinContent(i+1) if i<nbins2 else 0 # TODO: more elegant
            E1 = self.g1.GetBinError(i+1)
            E2 = self.g2.GetBinError(i+1)
            if N1==0 or N2==0:
                ratio.SetBinContent(i+1, 1 if self.equal else 0)
                ratio.SetBinError(i+1,0)
            else:
                r = N1/N2
                error = math.sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2))
                ratio.SetBinContent(i+1, r)
                ratio.SetBinError(i+1, error)
        return ratio

dir = "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/"

muon = "muon"
elec = "elec"

channels = [muon, elec]
# mtop = ['1665','1695','1715','1735','1755','1785']
mtop = ['1715','1735']

processes = {
    # 'data' : "uhh2.AnalysisModuleRunner.DATA.DATA",
    'ttbar' : "uhh2.AnalysisModuleRunner.MC.TTbar",
    # 'wjets' : "uhh2.AnalysisModuleRunner.MC.WJets",
    # 'st' : "uhh2.AnalysisModuleRunner.MC.SingleTop",
    # 'other' : "uhh2.AnalysisModuleRunner.MC.other",
}
uhh_mtop = "uhh2.AnalysisModuleRunner.MC.TTbar_mtop"

years = {
    '2016': '_2016v3',
    '2017': '_2017v2',
    '2018': '_2018',
}


print '+---------------------------------------------+'
print '| Usage: python Tau32_top.py <year> <channel> |'
if len(sys.argv)!=3:
    print '<E> I need a year an a channel!'
    print '+'+'-'*45+'+'
    sys.exit()
else:
    year = sys.argv[1]
    print '| --year:    {:<32} |'.format(sys.argv[1])
    print '| --channel: {:<32} |'.format(sys.argv[2])
    print '+'+'-'*45+'+'

year = sys.argv[1]
channel = sys.argv[2]

files = {
    p: {
        c: {
            y: rt.TFile(dir+c+"/"+processes[p]+years[y]+".root") for y in years
        } for c in channels
    } for p in processes
}

for m in mtop:
    files[m] = {
        c:{
            y: rt.TFile(dir+c+"/"+uhh_mtop+m+years[y]+".root") for y in years
        } for c in channels
    }

hist = 'comparison_topjet_xcone_pass_rec_masscut_140/ak8_hadjet_tau32'
hists = {}
for p in files:
    hists[p] = {}
    for c in channels:
        hists[p][c] = {}
        for y in years:
            hists[p][c][y] = rt.TH1F(files[p][c][y].Get(hist).Rebin(5))
            hists[p][c][y].SetDirectory(0)

for p in files:
    for c in channels:
        hist = hists[p][c]['2017'].Clone()
        hist.SetDirectory(0)
        hist.Add(hists[p][c]['2018'], 1)
        hists[p][c]['combine'] = hist

years['combine'] = 'combine'
channels.append('combine')
for p in files:
    hists[p]['combine'] = {}
    for y in years:
        hist = hists[p]['muon'][y].Clone()
        hist.SetDirectory(0)
        hist.Add(hists[p]['elec'][y], 1)
        hists[p]['combine'][y] = hist

# hists['bkg'] = {}
# hists['data-bkg'] = {}
# for c in channels:
#     hists['bkg'][c] = {}
#     hists['data-bkg'][c] = {}
#     for y in years:
#         hb = hists['other'][c][y]
#         hw = hists['wjets'][c][y]
#         hs = hists['st'][c][y]
#         hw.Add(hs, 1)
#         hb.Add(hw, 1)
#         hists['bkg'][c][y] = hb
#
#         hd = hists['data'][c][y]
#         hd.Add(hists['bkg'][c][y], -1)
#         hists['data-bkg'][c][y] = hd

for p in files:
    for c in channels:
        for y in years:
            hn = hists[p][c][y].Clone()
            hn.SetDirectory(0)
            hn.Scale(1/hn.Integral())
            hists[p][c][y] = hn

# for p in files:
#     for c in channels:
#         for y in years:
#             print p,c,y,hists[p][c][y].GetMaximum()

ratios = {}
for m in mtop:
    ratios[m] = Ratios(hists[m][channel][year], hists['ttbar'][channel][year], False).GetRatioGraph()

if "18" in year:
    TDR.cms_lumi = "59.7 fb^{-1}"
elif "17" in year:
    TDR.cms_lumi = "41.5 fb^{-1}"
elif "combine" in year:
    TDR.cms_lumi = "101 fb^{-1}"
elif "16" in year:
    TDR.cms_lumi = "35.9 fb^{-1}"
else:
    TDR.cms_lumi = "137.2 fb^{-1}"
TDR.extraText  = "Work in progress"

canv = tdrDiCanvas('tau32', 0, 1, 0.001, hists['ttbar'][channel][year].GetMaximum()*1.2, 0.81, 1.21, "#tau_{32}", "a.u.", 'Ratio', square=kRectangular, iPeriod=4, iPos=11)

canv.cd(1);
leg = tdrLeg(0.60, 0.10, 0.80, 0.34, 0.035, 42, rt.kBlack)

hists['ttbar'][channel][year].SetLineWidth(2)
tdrDraw(hists['ttbar'][channel][year], 'hist', rt.kFullCircle, rt.kGray+2, rt.kSolid, rt.kGray+2, 0, rt.kGray+2, 1)
leg.AddEntry(hists['ttbar'][channel][year], '1725', 'l')

colors = {'1665':rt.kAzure-7, '1695':rt.kGreen+2, '1715':rt.kYellow+2, '1735':rt.kRed+2, '1755':rt.kViolet-2, '1785':rt.kOrange+2}
for m in mtop:
    hists[m][channel][year].SetLineWidth(2)
    tdrDraw(hists[m][channel][year], 'hist', rt.kFullCircle, colors[m], rt.kSolid, colors[m], 0, colors[m], 0)
    leg.AddEntry(hists[m][channel][year], m, 'l')

canv.cd(2)
line = rt.TLine(0, 1, 1, 1)
tdrDrawLine(line, rt.kGray+2, lstyle=1, lwidth=2)
for m in mtop:
    ratios[m].SetLineWidth(2)
    tdrDraw(ratios[m], 'hist', rt.kFullCircle, colors[m], rt.kSolid, colors[m], 0, colors[m], 0)

canv.SaveAs("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/FSR/mtop/Tau32_mtop_"+year+'_'+channel+".pdf")
canv.SaveAs("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/FSR/mtop/Tau32_mtop_"+year+'_'+channel+".png")
