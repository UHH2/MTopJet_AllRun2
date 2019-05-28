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
    xmlfile = "MTopJetPostSelection_SYS_elec.xml"
    outputDir = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/elec/"
    option = "-slac "
    #option = " "
    #=========================================

    # the first value in systematics is the name of the variable in the xml file
    # the second value is the name of the directory
    systematics = [ ['ElID_variation', 'ELID'], ['ElReco_variation', 'ELRECO'], ['ElTrigger_variation', 'ELTR'], ['PU_variation', 'PU'] ]
    variations = ['up','down']

    systematics_btag = [ ['BTag_variation','BTAG'] ]
    variations_btag = ['up','down']

    # things are different for scale since the folder should always be named "SCALE_x"
    # and two values muR and muF have to be varied simultaneously
    systematics_scale = ['ScaleVariationMuR', 'ScaleVariationMuF']
    variations_scale = [ ['up','up'], ['up','none'], ['down','down'], ['down','none'], ['none','up'], ['none','down'] ]


    for sys in systematics:
        for value in variations:

            command_string = option+xmlfile+" -w workdirPostSelectionSYS."+sys[1]+"_"+value+" -o "+outputDir+sys[1]+"_"+value+" --ReplaceUserItem "+sys[0]+","+value

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


    for sys in systematics_btag:
        for value in variations_btag:

            command_string = option+xmlfile+" -w workdirPostSelectionSYS."+sys[1]+"_"+value+" -o "+outputDir+sys[1]+"_"+value+" --ReplaceUserItem "+sys[0]+","+value

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

    for value in variations_scale:

        command_string = option+xmlfile+" -w workdirPostSelectionSYS."+"SCALE"+"_"+value[0]+value[1]+" -o "+outputDir+"SCALE"+"_"+value[0]+value[1]+" --ReplaceUserItem "+variation_variables_scale[0]+","+value[0]+" --ReplaceUserItem "+variation_variables_scale[1]+","+value[1]

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
