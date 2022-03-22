import sys, os, ROOT
ROOT.gROOT.SetBatch(1) # Batch mode, e.g. suppress TCanvas on screen
ROOT.gStyle.SetOptTitle(0);
ROOT.gROOT.ForceStyle()


bins = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]

def GetValues(year):
    path = "/afs/desy.de/user/p/paaschal/WorkingArea/Plots/MTopJet/FSRuncertainty/"+year+"/combine/"
    values = []
    for bin in bins:
        dir = "Remove_"+bin
        file = open(path+dir+"/pt_0toInf/fsr_factor_cor.txt", 'r')
        lines = file.readlines()
        for l in lines:
            if "FSR" in l:
                string1 = l.replace("f_FSR  = ", "")
                string2 = string1.replace(" + ", " ")
                string3 = string2.replace(" - ", " ")
                values.append(string3.split(' '))
                # print(values)
    return values

def PlotValues(values, year):
    bins = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
    xvalues = []
    xup = []
    xdown = []
    central = []
    up = []
    down = []
    for v in values:
        central.append(float(v[0]));
        up.append(float(v[1]));
        down.append(float(v[2]));
        xup.append(0.);
        xdown.append(0.);
    dummy = ROOT.TH1F("dummy", "dummy", len(values), 0.5, 0.5+len(values))
    canvas = ROOT.TCanvas("canvas")
    i = 0
    while i<len(bins):
        dummy.GetXaxis().SetBinLabel(i+1, bins[i])
        xvalues.append(i+1);
        i += 1
    dummy.SetTitle("");
    dummy.SetStats(False);
    dummy.GetYaxis().SetTitle("1/#it{f}^{FSR}");
    dummy.GetXaxis().SetTitle("Removed Bin");
    ymin = 2 if "combine" in year else 0
    ymax = 4 if "combine" in year else 2
    dummy.GetYaxis().SetRangeUser(ymin, ymax);
    dummy.GetYaxis().SetTitleOffset(1.3);
    dummy.Draw();

    results = ROOT.TGraphAsymmErrors(len(values));
    i = 0
    while i<len(bins):
        results.SetPoint(i, xvalues[i], central[i])
        results.SetPointError(i, xdown[i], xup[i], down[i], up[i])
        i += 1
    results.SetMarkerStyle(8)
    results.Draw("P")
    errbandup = ROOT.TLine(0.5, central[0]+up[0], 0.5+len(central), central[0]+up[0])
    errbanddown = ROOT.TLine(0.5, central[0]-down[0], 0.5+len(central), central[0]-down[0])
    errbandup.SetLineColor(13);
    errbanddown.SetLineColor(13);
    errbandup.SetLineStyle(2);
    errbanddown.SetLineStyle(2);
    errbandup.Draw("SAME");
    errbanddown.Draw("SAME");
    results.Draw("P SAME");
    canvas.SaveAs("plots/RemoveBins_"+year+".pdf")

if __name__=="__main__":
    PlotValues(GetValues("combine"), "combine")
    PlotValues(GetValues("2016"), "2016")
