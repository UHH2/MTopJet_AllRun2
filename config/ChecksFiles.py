import glob, os, sys, ROOT

# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT. kFatal) + ";")

class CheckFileNumbers():
    def __init__(self, year="", study="", jecversion = "", jetLabel="", systematic=""):
        USER = os.environ["USER"]
        self.sframeDir = "/nfs/dust/cms/user/"+USER+"/MTopJet_Run2/PostSel/"
        self.xmlDir = os.environ["CMSSW_BASE"]+"/src/UHH2/MTopJet/config/"
        self.systematic = systematic.replace("_","")
        if systematic == "":
            self.folder = year+"/"+study+"/"+jecversion+"/"+jetLabel+"/"
        else:
            self.folder = year+"/"+study+"/"+jecversion+"/"+jetLabel+"/"+systematic.split("_")[0]+"/"+systematic.split("_")[1]+"/"
        self.module = "DiJetJERC_DiJetHLT"
        self.PrefixrootFile = "uhh2.AnalysisModuleRunner."

    def Count(self):
        config_path  = self.xmlDir+self.folder
        storage_path = self.sframeDir+self.folder
        # print config_path
        tot_xml = 0
        tot_root = 0
        for workdir in glob.glob(config_path+"workdir_"+self.module+"_*"):
            file_dir = workdir.replace(config_path,storage_path)
            workdir_name = workdir.replace(config_path,"").replace("workdir_"+self.module+"_","").replace(self.systematic,"")
            # print "\t", file_dir
            # print "\t", workdir_name
            for xml in glob.glob(workdir+"/*xml"):
                xml_name = xml.replace(workdir,"")
                if "Result" in xml_name: continue
                if self.module in xml_name: continue
                tot_xml += 1
                number = xml_name.replace(workdir_name,"").replace(".xml","").replace("/","").replace("_","")
                # print xml_name, number
                number = int(number)
                root_file = file_dir+"/"+self.PrefixrootFile+("MC." if "QCD" in xml_name else "DATA.")+workdir_name+"_"+str(number-1)+".root"
                root_file = root_file.replace("__","_")
                # print number, xml, root_file
                if not os.path.isfile(root_file):
                    print "FILE doesn't exist", root_file
                    print "mkdir -p", file_dir.replace(self.xmlDir,self.sframeDir), "; sframe_main",xml
                else:
                    ntuple = ROOT.TFile(str(root_file))
                    if ntuple.IsZombie() or ntuple.ReadKeys()==0:
                        # print root_file
                        print "sframe_main",xml
                    else: tot_root += 1
        print self.folder, " "*(70-len(self.folder)) , "DONE. Counted: XML=", tot_xml, " "*(5-len(str(tot_xml))), "ROOT:", tot_root


if __name__ == '__main__':

    # years       = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]
    years    = ["2016", "2017", "2018"]
    Labels   = ["", "SYS"]
    channels = ["muon", "elec", "muon_elecJMS", "elec_muonJMS"]
    # Systematics = [""]

    Workdirs = {}
    for y in years:
        for c in channels:
            ch = "Mu" if c == 'muon' or 'muon_' in c else "El"
            isJMS = True if 'JMS' in c else False
            jms = "_JMS" if isJMS else ""
            chjms = ""
            if isJMS: chjms = "muon" if "_muon" in c else "elec"
            print chjms
            Workdirs[y] = [{c:{"":"Workdir_PostSel"+ch+"_"+y+jms}}]

   # SFrameUncerts_2017MTopJetelec         Workdir_PostSelEl_2016_JMSmuon        Workdir_PostSelMu_2017_JMSelec
   # SFrameUncerts_2017MTopJetmuon         Workdir_PostSelEl_2017                Workdir_PostSelMu_2018
   # SFrameUncerts_2016MTopJetJMSelecmuon  SFrameUncerts_2018MTopJetJMSelecmuon  Workdir_PostSelEl_2017_JMSmuon
   # SFrameUncerts_2016MTopJetJMSmuonelec  SFrameUncerts_2018MTopJetJMSmuonelec  Workdir_PostSelEl_2018
   # SFrameUncerts_2016MTopJetelec         SFrameUncerts_2018MTopJetelec         Workdir_PostSelEl_2018_JMSmuon
   # SFrameUncerts_2016MTopJetmuon         SFrameUncerts_2018MTopJetmuon         Workdir_PostSelMu_2016
   # SFrameUncerts_2017MTopJetJMSelecmuon  Workdir_PostSelMu_2016_JMSelec        Workdir_PostSelMu_2018_JMSelec
   # SFrameUncerts_2017MTopJetJMSmuonelec  Workdir_PostSelMu_2017

    # for year in years:
    #     for study in studies:
    #         for jecversion in JECVersions[year]:
    #             for jetLabel in JetLabels:
    #                 for systematic in Systematics:
    #                     CFN = CheckFileNumbers(year=year, study=study, jecversion=jecversion, jetLabel=jetLabel, systematic=systematic)
    #                     CFN.Count()
