#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
     #'name' : {'item name': 'item value', ...},
     # 'JEC_up'                : {'jecsmear_direction':'up'},
     # 'JEC_down'              : {'jecsmear_direction':'down'},
     # 'JER_up'                : {'jersmear_direction':'up'},
     # 'JER_down'              : {'jersmear_direction':'down'},
     # 'COR_up'                : {'JetCorrection_direction':'up'},
     # 'COR_down'              : {'JetCorrection_direction':'down'},
     'JMS_upup'              : {'JetMassScale_direction':'upup'},
     'JMS_updown'            : {'JetMassScale_direction':'updown'},
     'JMS_downup'            : {'JetMassScale_direction':'downup'},
     'JMS_downdown'          : {'JetMassScale_direction':'downdown'},
     'JMS_up'                : {'JetMassScale_direction':'up'},
     'JMS_down'              : {'JetMassScale_direction':'down'},
}
start_all_parallel = True

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude! ... do not forget the name!'
    exit(-1)

def insert_str(string, str_to_insert, index):
    return string[:index] + str_to_insert + string[index:]

def set_uncert_func(uncert_name):
    uncert = sys_uncerts[uncert_name]
    def do_set_uncert(element_tree):
        cycle = element_tree.getroot().find('Cycle')
        user_config = cycle.find('UserConfig')
        output_dir = cycle.get('OutputDirectory')

        cycle.set('OutputDirectory', os.path.join(output_dir, uncert_name+'/'))

        for name, value in uncert.iteritems():
            uc_item = list(i for i in user_config if i.get('Name') == name)
            assert uc_item, 'could not find item with name: %s' % name
            uc_item[0].set('Value', value)

    return do_set_uncert

def DirName(prefix, xmlname):
    newname = xmlname.replace('MTopJet', '')
    newname = xmlname.replace('Alex', '')
    newname = newname.replace('PostSelection', '')
    newname = newname.replace('SYS', '')
    newname = newname.replace('.xml', '')
    newname = newname.replace('_', '')
    newname = newname.replace('/', '')
    newname = insert_str(newname, "_" , 4)
    return prefix+'_'+newname
from varial.extensions.sframe import SFrame
from varial import tools
if start_all_parallel:
    ToolChain = tools.ToolChainParallel
else:
    ToolChain = tools.ToolChain

class MySFrameBatch(SFrame):
    def configure(self):
        self.xml_doctype = self.xml_doctype + """
<!--
   <ConfigParse NEventsBreak="100000" LastBreak="0" FileSplit="0"/>
   <ConfigSGE RAM="2" DISK="2" Mail="alexander.paasch@desy.de" Notification="as" Workdir="PostSelSYSEl_workdir"/>
-->
"""
        if os.path.exists(self.cwd + 'workdir'):
            opt = ' -rl --exitOnQuestion'
        else:
            opt = ' -sl --exitOnQuestion'

        self.exe = 'sframe_batch.py' + opt

sframe_tools = ToolChain(
    #'SFrameUncertsSR',
    DirName('SFrameUncerts',sys.argv[1]),
    list(
        SFrame(
            cfg_filename=sys.argv[1],
            xml_tree_callback=set_uncert_func(uncert),
            name='SFrame_' + uncert,
            halt_on_exception=False,
        )
        for uncert in sys_uncerts
    )
)

if __name__ == '__main__':
    varial.tools.Runner(sframe_tools)
