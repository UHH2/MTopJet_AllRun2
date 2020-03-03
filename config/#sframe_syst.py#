#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     'TagEffi_3prong_up'     : {'TagEffi_variation':'up_3prong'},
     'TagEffi_3prong_down'   : {'TagEffi_variation':'down_3prong'},
     'TagEffi_2prong_up'     : {'TagEffi_variation':'up_2prong'},
     'TagEffi_2prong_down'   : {'TagEffi_variation':'down_2prong'},
     'TagEffi_1prong_up'     : {'TagEffi_variation':'up_1prong'},
     'TagEffi_1prong_down'   : {'TagEffi_variation':'down_1prong'},
     'ScaleVariationMuF_up'  : {'ScaleVariationMuF':'up'},
     'ScaleVariationMuF_down': {'ScaleVariationMuF':'down'},
     'ScaleVariationMuR_up'  : {'ScaleVariationMuR':'up'},
     'ScaleVariationMuR_down': {'ScaleVariationMuR':'down'},
     'PU_up'                 : {'PU_variation':'up'},
     'PU_down'               : {'PU_variation':'down'},
     'MUID_up'               : {'MuonID_variation':'up'},
     'MUID_down'             : {'MuonID_variation':'down'},
     'MUTR_up'               : {'MuonTrigger_variation':'up'},
     'MUTR_down'             : {'MuonTrigger_variation':'down'},
     'FSR_up_2'              : {'PS_variation':'FSRup_2'},
     'FSR_down_2'            : {'PS_variation':'FSRdown_2'},
     'BTAG_bc_up'            : {'BTag_variation':'up_bc'},
     'BTAG_bc_down'          : {'BTag_variation':'down_bc'},
     'BTAG_udsg_up'          : {'BTag_variation':'up_udsg'},
     'BTAG_udsg_down'        : {'BTag_variation':'down_udsg'},
}
start_all_parallel = True

############################################################### script code ###
import varial
import sys
import os

if len(sys.argv) != 2:
    print 'Plz. give me da name of da sframe-config! ... dude!'
    exit(-1)


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
    newname = xmlname.replace('TopTaggingSF', '')
    newname = newname.replace('PostSelection', '')
    newname = newname.replace('SYS', '')
    newname = newname.replace('.xml', '')
    newname = newname.replace('_', '')
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
   <ConfigSGE RAM="2" DISK="2" Mail="dennis.schwarz@desy.de" Notification="as" Workdir="/nfs/dust/cms/user/schwarzd/PostSelMu_workdir"/>
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
