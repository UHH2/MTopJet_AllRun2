#!/usr/bin/env python

import sys

# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py

#simple script that runs several sframe batch jobs and creates everything you might need
if __name__ == "__main__":
    debug = False
    remove = True #remove directories with old results

    #put your local sfram_batch dir in search path
    sys.path.append('/nfs/dust/cms/user/schwarzd/SFrameBatch/')
    #import the main function
    from sframe_batch import SFrameBatchMain

    #=========================================
    #important variables!
    #=========================================
    xmlfile = "MTopJetSelection_muon.xml"
    outputDir = "/nfs/dust/cms/user/schwarzd/MTopJet/Selection/muon/"
    option = "-slac "
    #option = " "
    #=========================================

    # the first value in systematics is the name of the variable in the xml file
    # the second value is the name of the directory
    # systematics = [ ['jecsmear_direction', 'JEC'], ['jersmear_direction', 'JER'], ['JetCorrection_direction', 'COR'] ]
    # variations = ['up','down']
    systematics = [ ['jersmear_direction', 'JER'], ['JetCorrection_direction', 'COR'] ]
    variations = ['up','down']
    #systematics = [ ['NonClosureUncertainty', 'NonClosure'] ]
    #variations = ['true']
    for sys in systematics:
        for value in variations:

            command_string = option+xmlfile+" -w workdirSelectionSYS."+sys[1]+"_"+value+" -o "+outputDir+sys[1]+"_"+value+" --ReplaceUserItem "+sys[0]+","+value

            if debug :print command_string.split(" ")
            else:
                #try:
                SFrameBatchMain(command_string.split(" "))
                """
                except:
                    print "SFrameBatch did crash during running:"
                    print command_string
                    sys.exit(1)
                """
