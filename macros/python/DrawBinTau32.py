import sys, ROOT, os
ROOT.gROOT.SetBatch(1) # Batch mode, e.g. suppress TCanvas on screen
# ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.ForceStyle()
ROOT.gErrorIgnoreLevel = ROOT.kWarning;

file = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/files/Tau32.root'
save = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/FSR/'

ptwindows   = [
    ['0toInf', '420to450', '400to450', '450to500', '500to600', '600toInf'],
    ['400to480', '480to580', '580toInf'],
    ['400to440', '440to480', '480to550', '550to610', '610to700', '700toInf']
]

names = ['DATA', 'TTbar', 'TTbarFSRupsqrt2', 'TTbarFSRdownsqrt2', 'TTbarFSRup2', 'TTbarFSRdown2', 'TTbarFSRup4', 'TTbarFSRdown4']
settings = { # [color, line style, marker, marker size]
    'DATA':              [ROOT.kBlack,    1, 8, 1, 'PE'],
    'TTbar':             [ROOT.kGray+2,   1, 0, 0, 'SAME HIST'],
    'TTbarFSRupsqrt2':   [ROOT.kAzure+7,  1, 0, 0, 'SAME HIST'],
    'TTbarFSRdownsqrt2': [ROOT.kAzure+7,  2, 0, 0, 'SAME HIST'],
    'TTbarFSRup2':       [ROOT.kGreen-4,  1, 0, 0, 'SAME HIST'],
    'TTbarFSRdown2':     [ROOT.kGreen-4,  2, 0, 0, 'SAME HIST'],
    'TTbarFSRup4':       [ROOT.kOrange+7, 1, 0, 0, 'SAME HIST'],
    'TTbarFSRdown4':     [ROOT.kOrange+7, 2, 0, 0, 'SAME HIST'],
}

# var+"_"+win+"__DATA__norm__"+year
def DrawPlots(study,year):
    inFile = ROOT.TFile.Open(file, 'READ')
    if not inFile :
        print " Failed to get file "
        sys.exit (1)
    for scale in ['', '__norm']:
        loop_list = 0
        outFile = 'Tau32_'+study+scale.replace('__', '_')+'_'+year+'.pdf'
        for windows in ptwindows:
            canvas = ROOT.TCanvas(windows[0]+scale+year, "", 800, 800)
            canvas.cd()
            canvas.Divide(3,2)

            loop_win = 1

            for win in windows:
                canvas.cd(loop_win)

                for name in names:
                    if year == 'combine':
                        if 'TTbarFSRupsqrt2' == name:   continue
                        if 'TTbarFSRdownsqrt2' == name: continue
                        if 'TTbarFSRup2' == name:       continue
                        if 'TTbarFSRdown2' == name:     continue
                    if year == '2016':
                        if 'TTbarFSRupsqrt2' == name:   continue
                        if 'TTbarFSRdownsqrt2' == name: continue
                        if 'TTbarFSRup4' == name:       continue
                        if 'TTbarFSRdown4' == name:     continue
                    hist = inFile.Get(study+'_'+win+'__'+name+scale+'__'+year)
                    if not hist :
                        print " Failed to get data histogram "
                        sys.exit (1)
                    hist.SetLineColor(settings[name][0])
                    hist.SetTitle(study+' '+win+' GeV')
                    hist.SetMarkerColor(settings[name][0])
                    hist.SetLineStyle(settings[name][1])
                    hist.SetMarkerStyle(settings[name][2])
                    hist.SetMarkerSize(settings[name][3])
                    hist.GetYaxis().SetRangeUser(0, hist.GetMaximum()*1.5)
                    hist.GetYaxis().SetTitle('Events')
                    hist.GetXaxis().SetTitle('#tau_{32}')
                    # print name, hist.GetMaximum(), settings[name][4]
                    hist.Draw(settings[name][4])
                loop_win += 1

            if   loop_list == 0:
                canvas.Print(save+outFile+'(', 'pdf')
            elif loop_list<len(ptwindows)-1:
                canvas.Print(save+outFile)
            elif loop_list==len(ptwindows)-1:
                canvas.Print(save+outFile+')', 'pdf')
            else:
                print 'You should not be here o.O'
            loop_list += 1

if __name__ == '__main__':
    DrawPlots('pt','combine')
    DrawPlots('pt','2016')
