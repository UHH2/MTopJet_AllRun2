import glob, os, json, sys, math
import ROOT as rt
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)
rt.gErrorIgnoreLevel = rt.kWarning
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

print '+---------------------------------------------+'
print '| Usage: python Tau32_top.py <year> <channel> |'
if len(sys.argv)!=2:
    print '| <E> I need a year!                          |'
    print '+'+'-'*45+'+'
    sys.exit()
else:
    print '| --year:    {:<32} |'.format(sys.argv[1])
    # print '| --channel: {:<32} |'.format(sys.argv[2])
    print '+'+'-'*45+'+'

year = sys.argv[1]
channel = 'combine'

plots = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/WMass/SYS/"+year+"/"
os.system('mkdir -p '+plots)

file = rt.TFile("/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/files/WMassPlots.root")

bins = ['hh', 'hl', 'lh', 'll']
process = ['DATA','TTbar','TTbarJECup','TTbarJECdown','TTbarCORup','TTbarCORdown','TTbarJERup','TTbarJERdown','TTbarhdampup','TTbarhdampdown','TTbartuneup','TTbartunedown','TTbargluonmove','TTbarQCDbased','TTbarFSRup_sqrt2','TTbarFSRdown_sqrt2','TTbarFSRup_2','TTbarFSRdown_2','TTbarFSRup_4','TTbarFSRdown_4']
# "wmass"+addition+"__DATA__"+sYear
sys = ['JEC','COR','JER','hdamp','tune', 'CR','FSR_sqrt2','FSR_2','FSR_4']

years = ['2016','2017','2018','combine']
hists = {b: {p: {} for p in process} for b in bins}
for b in bins:
    for p in process:
        for y in years:
            hists[b][p][y] = file.Get('wmass_'+b+'__'+p+'__'+y)
            hists[b][p][y].SetDirectory(0)

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

for b in bins:
    for s in sys:
        draws = []
        if 'JEC'==s: draws=['TTbarJECup', 'TTbarJECdown']
        elif 'tune'==s: draws=['TTbartuneup', 'TTbartunedown']
        elif 'hdamp'==s: draws=['TTbarhdampup', 'TTbarhdampdown']
        elif 'COR'==s: draws=['TTbarCORup', 'TTbarCORdown']
        elif 'JER'==s: draws=['TTbarJERup', 'TTbarJERdown']
        elif 'CR'==s: draws=['TTbargluonmove', 'TTbarQCDbased']
        elif 'FSR_sqrt2'==s: draws=['TTbarFSRup_sqrt2', 'TTbarFSRdown_sqrt2']
        elif 'FSR_2'==s: draws=['TTbarFSRup_2', 'TTbarFSRdown_2']
        elif 'FSR_4'==s: draws=['TTbarFSRup_4', 'TTbarFSRdown_4']

        canv = tdrDiCanvas('wmass'+s+b, 60, 120, 0.001, hists[b]['TTbar'][year].GetMaximum()*1.2, 0.81, 1.21, "m_\mathrm(W)", "a.u.", 'Ratio', square=kRectangular, iPeriod=4, iPos=11)

        ratios = [Ratios(hists[b][d][year], hists[b]['TTbar'][year], False).GetRatioGraph() for d in draws]

        canv.cd(1);
        leg = tdrLeg(0.60, 0.70, 0.80, 0.89, 0.035, 42, rt.kBlack)

        hists[b]['TTbar'][year].SetLineWidth(2)
        tdrDraw(hists[b]['TTbar'][year], 'hist', rt.kFullCircle, rt.kGray+2, rt.kSolid, rt.kGray+2, 0, rt.kGray+2, 1)
        leg.AddEntry(hists[b]['TTbar'][year], '1725', 'l')

        styles = [rt.kSolid, rt.kDashed]
        for i in range(2):
            d = draws[i]
            hists[b][d][year].SetLineWidth(2)
            tdrDraw(hists[b][d][year], 'hist', rt.kFullCircle, rt.kRed+1, styles[i], rt.kRed+1, 0, rt.kRed+1, 0)
            leg.AddEntry(hists[b][d][year], d.replace('TTbar',''), 'l')

        canv.cd(2)
        line = rt.TLine(60, 1, 120, 1)
        tdrDrawLine(line, rt.kGray+2, lstyle=1, lwidth=2)
        for i in range(2):
            r = ratios[i]
            r.SetLineWidth(2)
            tdrDraw(r, 'hist', rt.kFullCircle, rt.kRed+1, styles[i], rt.kRed+1, 0, rt.kRed+1, 0)

        save = plots+'WMass_'+b+'_'+s+'_'+year+'_'+channel
        canv.SaveAs(save+'.pdf')
        canv.SaveAs(save+'.png')
