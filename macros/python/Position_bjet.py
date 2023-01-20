import sys, os, ROOT
ROOT.gROOT.SetBatch(1) # Batch mode, e.g. suppress TCanvas on screen
ROOT.gStyle.SetOptTitle(0);
ROOT.gROOT.ForceStyle()

file = '/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/muon/uhh2.AnalysisModuleRunner.MC.TTbar_2018_btag.root'
# save_afs = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/JMS/Projections/";
save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/BTag/";

entries = {
    '': [],
    'single': ['0: in had', '1: in lep', '2: in non'],
    'both':   ['-1: no option', ' 0: in non', ' 1: both had', ' 2: both lep', ' 3: in both', ' 4: 1 in had', ' 5: 1 in lep', ' 6: only 1 in had', ' 7: only 1 in lep'],
}

def Cosmetics(hist, color=ROOT.kBlack, style=1, width=1, xmin=0, xmax=3, xname='Location'):
    hist.SetStats(0)
    hist.SetLineColor(color)
    hist.SetLineStyle(style)
    hist.SetLineWidth(width)
    hist.SetFillColorAlpha(color, 0.4)
    hist.GetXaxis().SetRangeUser(xmin, xmax)
    hist.GetXaxis().SetTitle(xname)
    hist.GetYaxis().SetTitle('Events')

def Draw(hists, name='', option='hist', nleg=''):
    canvas = ROOT.TCanvas('canvas'+name)
    i=0
    while i < len(hists):
        add = ' same' if i>0 else ''
        hists[i].Draw(option+add)
        i += 1

    leg = ROOT.TLegend()
    leg = ROOT.TLegend(0.7, 0.5 if 'both' in name else 0.7, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    i=0
    while i < len(entries[nleg]):
        leg.AddEntry(0, entries[nleg][i], '')
        i += 1
    leg.Draw()
    canvas.SaveAs(save_nfs+name+'.pdf')
    leg.Clear()
    canvas.Close()

def Check(object, name):
    if not object:
        print " Failed to get "+name
        sys.exit (1)

def RelativeShare(hist):
    bin0 = hist.GetBinContent(1)
    bin1 = hist.GetBinContent(2)
    bin2 = hist.GetBinContent(3)
    combine = bin0+bin1+bin2
    integral = hist.GetIntegral()
    rbin0 = bin0/combine
    rbin1 = bin1/combine
    rbin2 = bin2/combine
    print("%.2f" % (round(rbin0, 4)*100))
    print("%.2f" % (round(rbin1, 4)*100))
    print("%.2f" % (round(rbin2, 4)*100))
    print bin0, bin1, bin2, combine, integral
    print bin0/combine, bin1/combine, bin2/combine

def Plot():
    inFile = ROOT.TFile.Open(file, 'READ')
    Check(inFile, 'file')

    btag1 = inFile.Get('PositionBTagHists/Position1stBtag')
    btag2 = inFile.Get('PositionBTagHists/Position2ndBtag')
    btagb = inFile.Get('PositionBTagHists/PositionBoth')
    Check(btag1, '1st bTag')
    Check(btag2, '2nd bTag')
    Check(btagb, 'both bTags')

    RelativeShare(btag1)

    Cosmetics(btag1, ROOT.kRed, 1, 1)
    Cosmetics(btag2, color=ROOT.kBlue, style=1, width=1)
    Cosmetics(btagb, color=ROOT.kBlue, style=1, width=1, xmin=-1, xmax=8)

    Draw([btag1], name='btag_1st', option='hist', nleg='single')
    Draw([btag2], name='btag_2nd', option='hist', nleg='single')
    Draw([btagb], name='btag_both', option='hist', nleg='both')

if __name__=="__main__":
    os.system('mkdir -p '+save_nfs)
    Plot()
