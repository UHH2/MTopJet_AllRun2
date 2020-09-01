#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')
# sys.path.append('/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/Varial')

##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
    # 'FSRup_sqrt2'              : {'PS_variation':'FSRup_sqrt2'},
    # 'FSRdown_sqrt2'            : {'PS_variation':'FSRdown_sqrt2'},
    # 'FSRup_2'                  : {'PS_variation':'FSRup_2'},
    # 'FSRdown_2'                : {'PS_variation':'FSRdown_2'},
    # 'FSRup_4'                  : {'PS_variation':'FSRup_4'},
    # 'FSRdown_4'                : {'PS_variation':'FSRdown_4'},
    # 'ISRup_sqrt2'              : {'PS_variation':'ISRup_sqrt2'},
    # 'ISRdown_sqrt2'            : {'PS_variation':'ISRdown_sqrt2'},
    # 'ISRup_2'                  : {'PS_variation':'ISRup_2'},
    # 'ISRdown_2'                : {'PS_variation':'ISRdown_2'},
    # 'ISRup_4'                  : {'PS_variation':'ISRup_4'},
    # 'ISRdown_4'                : {'PS_variation':'ISRdown_4'},
    'JMS_upup'                 : {'JetMassScale_direction':'upup'},
    'JMS_updown'               : {'JetMassScale_direction':'updown'},
    'JMS_downup'               : {'JetMassScale_direction':'downup'},
    'JMS_downdown'             : {'JetMassScale_direction':'downdown'},
}
start_all_parallel = False

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude!'
    exit(-1)

################################################################################
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

################################################################################
def insert_str(string, str_to_insert, index):
    return string[:index] + str_to_insert + string[index:]
def DirName(prefix, xmlname):
    newname = xmlname.replace('MTopJet', '')
    newname = newname.replace('PostSelection', '')
    newname = newname.replace('SYS', '')
    newname = newname.replace('muon', '')
    newname = newname.replace('.xml', '')
    newname = newname.replace('_', '')
    newname = newname.replace('/', '')
    # newname = insert_str(newname, "_" , 4)
    return prefix+'_'+newname

################################################################################
from varial.extensions.sframe import SFrame
from varial import tools
if start_all_parallel:
    ToolChain = tools.ToolChainParallel
else:
    ToolChain = tools.ToolChain

################################################################################
class MySFrameBatch(SFrame):

    def configure(self):
        self.xml_doctype = self.xml_doctype + """
<!--
   <ConfigParse NEventsBreak="100000" FileSplit="0" AutoResubmit="0" />
   <ConfigSGE RAM ="2" DISK ="2" Mail="alexander.paasch@desy.de" Notification="as" Workdir="workdir"/>
-->
"""
   	if os.path.exists(self.cwd + 'workdir'):
            opt = ' -f'
        else:
            opt = ' -s'

        self.exe = 'sframe_batch.py' + opt

sframe_tools = ToolChain(
    DirName('SFrameUncerts',sys.argv[1]),
    list(
	MySFrameBatch(
            cfg_filename=sys.argv[1],
            xml_tree_callback=set_uncert_func(uncert),
            name='SFrame_' + uncert,
            halt_on_exception=False,
        )
	for uncert in sys_uncerts
    )
)

################################################################################
if __name__ == '__main__':
    varial.tools.Runner(sframe_tools)
