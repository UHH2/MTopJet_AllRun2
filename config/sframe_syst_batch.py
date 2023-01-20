#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')
# sys.path.append('/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/Varial')

##################################### definition of UserConfig item changes ###

#'name' : {'item name': 'item value', ...},
sys_uncerts = {}
if "2017" in sys.argv[1] or "2018" in sys.argv[1]:
    sys_uncerts['FSRup_4'] = {'PS_variation':'FSRup_4'}
    sys_uncerts['FSRup_sqrt2'] = {'PS_variation':'FSRup_sqrt2'}
    sys_uncerts['FSRdown_sqrt2'] = {'PS_variation':'FSRdown_sqrt2'}
    sys_uncerts['FSRup_2'] = {'PS_variation':'FSRup_2'}
    sys_uncerts['FSRdown_2'] = {'PS_variation':'FSRdown_2'}
    sys_uncerts['FSRdown_4'] = {'PS_variation':'FSRdown_4'}
    sys_uncerts['ISRup_sqrt2'] = {'PS_variation':'ISRup_sqrt2'}
    sys_uncerts['ISRdown_sqrt2'] = {'PS_variation':'ISRdown_sqrt2'}
    sys_uncerts['ISRup_2'] = {'PS_variation':'ISRup_2'}
    sys_uncerts['ISRdown_2'] = {'PS_variation':'ISRdown_2'}
    sys_uncerts['ISRup_4'] = {'PS_variation':'ISRup_4'}
    sys_uncerts['ISRdown_4'] = {'PS_variation':'ISRdown_4'}
sys_uncerts['JMS_downdown'] = {'JetMassScale_direction':'downdown'}
sys_uncerts['JMS_upup'] = {'JetMassScale_direction':'upup'}
sys_uncerts['JEC_up'] = {'jecsmear_direction':'up'}
sys_uncerts['JEC_down'] = {'jecsmear_direction':'down'}
sys_uncerts['JER_up'] = {'jersmear_direction':'up'}
sys_uncerts['JER_down'] = {'jersmear_direction':'down'}
sys_uncerts['COR_up'] = {'JetCorrection_direction':'up'}
sys_uncerts['COR_down'] = {'JetCorrection_direction':'down'}
sys_uncerts['JMS_flavor_up'] = {'JetMassScale_Flavor':'up'}
sys_uncerts['JMS_flavor_down'] = {'JetMassScale_Flavor':'down'}

# sys_uncerts['JMS_updown'] = {'JetMassScale_direction':'updown'}
# sys_uncerts['JMS_downup'] = {'JetMassScale_direction':'downup'}
# sys_uncerts['JMS_up'] = {'JetMassScale_direction':'up'}
# sys_uncerts['JMS_down'] = {'JetMassScale_direction':'down'}

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
    newname = xmlname.replace('Alex', '')
    newname = newname.replace('PostSelection', '')
    newname = newname.replace('SYS', '')
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
   <ConfigParse NEventsBreak="50000" FileSplit="0" AutoResubmit="0" />
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
