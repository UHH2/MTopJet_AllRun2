import sys, os, ROOT, math
from prettytable import PrettyTable
ROOT.gROOT.SetBatch(1) # Batch mode, e.g. suppress TCanvas on screen
ROOT.gStyle.SetOptTitle(0);
ROOT.gROOT.ForceStyle()

file = 'files/PaperPlots_Peak.root'
save_afs = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/JMS/Projections/";
save_nfs = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/JMS/Projections/";

def Cosmetics(hist, color=ROOT.kBlack, style=1, width=1, xmin=-2, xmax=2):
    hist.SetStats(0)
    hist.SetLineColor(color)
    hist.SetLineStyle(style)
    hist.SetLineWidth(width)
    hist.SetMarkerStyle(1)
    hist.SetFillColorAlpha(color, 0.4)
    hist.GetXaxis().SetRangeUser(xmin, xmax)
    hist.GetYaxis().SetRangeUser(128, 180)

def CosmeticsFit(fit, name='JEC', color=ROOT.kBlack, style=1, width=1, xmin=-2, xmax=2):
    fit.SetLineColor(color)
    fit.SetLineStyle(style)
    fit.SetLineWidth(width)
    fit.SetMarkerStyle(1)
    fit.SetFillColorAlpha(color, 0.4)
    fit.GetXaxis().SetRangeUser(xmin, xmax)
    ymin = 120
    ymax = 180 if 'JEC' in name else 140
    fit.GetYaxis().SetRangeUser(ymin, ymax)

def Draw(hists, name='', option=['hist']):
    canvas = ROOT.TCanvas('canvas'+name)
    i=0
    while i < len(hists):
        add = ' same' if i>0 else ''
        hists[i].Draw(option+add)
        i += 1
    canvas.Print(save_nfs+name+'.pdf')

def Draw1DChi2(hist, name='', option='hist', xmin=0, xmax=0, min_x=0, min_y=0, eu=0, ed=0):
    canvas = ROOT.TCanvas('canvas'+name)
    hist.Draw(option)
    lmin=ROOT.TLine(min_x, 128, min_x, min_y)
    lerr=ROOT.TLine(xmin, min_y+1, xmax, min_y+1)
    lerrU=ROOT.TLine(min_x+eu, 128, min_x+eu, min_y+1)
    lerrD=ROOT.TLine(min_x-ed, 128, min_x-ed, min_y+1)
    lmin.SetLineColor(ROOT.kGray+2)
    lerr.SetLineColor(ROOT.kGray+2)
    lerrU.SetLineColor(ROOT.kGray+2)
    lerrD.SetLineColor(ROOT.kGray+2)
    lmin.SetLineStyle(2)
    lerr.SetLineStyle(1)
    lerrU.SetLineStyle(1)
    lerrD.SetLineStyle(1)
    lmin.Draw('same')
    lerr.Draw('same')
    lerrU.Draw('same')
    lerrD.Draw('same')
    leg = ROOT.TLegend(0.3, 0.7, 0.7, 0.8)
    leg.AddEntry(lmin, 'Minimum ('+str("%.2f"%min_x)+','+str("%.2f"%min_y)+')', 'l')
    leg.AddEntry(lerr, '1#sigma ('+str("%.2f"%(min_x-ed))+','+str("%.2f"%(min_x+eu))+')', 'l')
    leg.Draw()
    canvas.Print(save_nfs+name+'.pdf')

def Get1DJMS(func, xmin=-3, xmax=3):
    # func = ROOT.TF1("func", f, -1.5, 1.5)
    minimum = func.GetMinimum()
    ymin = func.GetX(minimum, xmin, xmax)
    eu = func.GetX(minimum+1, ymin, xmax)-ymin
    ed = ymin-func.GetX(minimum+1, xmin, ymin)
    return [ymin, eu, ed, minimum]

def CreateHist(formula, axis='JEC'):
    isXCone = True if axis=='XCone' else False
    hold = 'y' if isXCone else 'x'
    YtoX = True if not isXCone else False
    range = [-1.0, 1.0] if isXCone else [-0.5, 1.7]
    step = 0.1
    nbins = int(abs(range[0]-range[1])/step)
    print nbins, int(nbins)
    hist = ROOT.TH1F('hist', axis, nbins, range[0], range[1])
    start = range[0]
    print hist.GetNbinsX()
    c = 1
    while c<=hist.GetNbinsX():
        print c,hist.GetXaxis().GetBinCenter(c),hist.GetXaxis().GetBinLowEdge(c),hist.GetXaxis().GetBinUpEdge(c)
        c+=1
    print "-----------------------"
    table = PrettyTable()
    table.field_names = ["Value", "bin", "fmin", "fup", "fdown", "chi2 min"]
    while start <= range[1]:
        value = str(start+step/2)
        f_tmp = formula.replace(hold, value)
        f = f_tmp
        if YtoX:
            f = f_tmp.replace('y', 'x')
        # func = ROOT.TF1("func", f, range[0], range[1])
        func = ROOT.TF1("func", f, range[0], range[1])
        CosmeticsFit(func, name=axis,xmin=range[0], xmax=range[1])
        Draw([func], 'fits/'+axis+'_'+str(start),'')
        min = Get1DJMS(func)
        bin = hist.FindBin(start+step/2)
        print value,bin,hist.GetXaxis().GetBinLowEdge(bin),hist.GetXaxis().GetBinUpEdge(bin)
        hist.SetBinContent(bin, min[3])
        hist.SetBinError(bin, (min[1]+min[2])/2)
        table.add_row([value,bin,"%.2f"%min[0],"%.2f"%min[1],"%.2f"%min[2],"%.2f"%min[3]])
        start += step

    print table
    fit = ROOT.TF1('parabel', '[0]*x*x+[1]*x+[2]', range[0], range[1])
    hist.Fit(fit, 'MQ', '', range[0], range[1])
    Cosmetics(hist, xmin=range[0], xmax=range[1])
    # Draw([hist], '1D_'+axis, 'PE')
    # def Draw1DChi2(hist, name='', option='hist', xmin=0, xmax=0, min_x=0, min_y=0, eu=0, ed=0):

    print fit.Eval(0.6), fit.GetParameter(0)
    fit2 = hist.GetFunction('parabel')
    factor = Get1DJMS(fit, range[0], range[1])
    Draw1DChi2(hist, name='1D_'+axis, option='PE', xmin=range[0], xmax=range[1], min_x=factor[0], min_y=factor[3], eu=factor[1], ed=factor[2])
    # Draw([fit], '1D_fit_'+axis, '')
    print axis,':\t', "%.2f"%factor[0],'+',"%.2f"%factor[1],'-',"%.2f"%factor[2],'('+str("%.2f"%factor[3])+')'

    return hist

def GetProjections():
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    inFile = ROOT.TFile.Open(file, 'READ')
    if not inFile :
        print "<E> Failed to get file "
        sys.exit (1)
    else:
        print 'Open '+file
    chi2 = inFile.Get('Functions/JMS_Chi2')
    if not chi2 :
        print "<E> Failed to get function "
        sys.exit (1)
    else:
        print 'Get function ... '
    print chi2.Eval(0.603918,  -0.060578)
    hist = chi2.GetHistogram()
    formula = chi2.GetExpFormula()
    print 'Start to create hists'
    print '\t ... JEC'
    CreateHist(str(formula), 'JEC')
    print '\t ... XCone'
    CreateHist(str(formula), 'XCone')

if __name__=="__main__":
    print '[Usage] Run from macro directory, not macro/python'
    os.system('mkdir -p '+save_nfs+'fits')
    GetProjections()
