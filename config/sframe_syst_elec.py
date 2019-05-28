#!/usr/bin/env python

import os, sys, itertools

sys.path.append('/nfs/dust/cms/user/tholenhe/installs/varial-stable/Varial')


##################################### definition of UserConfig item changes ###

sys_uncerts = {
    # 'name' : {'item name': 'item value', ...},
     'SCALE_upup'            : {'ScaleVariationMuR':'up','ScaleVariationMuF':'up'},
     'SCALE_upnone'          : {'ScaleVariationMuR':'up','ScaleVariationMuF':'none'},
     'SCALE_noneup'          : {'ScaleVariationMuR':'none','ScaleVariationMuF':'up'},
     'SCALE_nonedown'        : {'ScaleVariationMuR':'none','ScaleVariationMuF':'down'},
     'SCALE_downnone'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'none'},
     'SCALE_downdown'        : {'ScaleVariationMuR':'down','ScaleVariationMuF':'down'},
     'PU_up'                 : {'PU_variation':'up'},
     'PU_down'               : {'PU_variation':'down'},
     'ELID_up'               : {'ElID_variation':'up'},
     'ELID_down'             : {'ElID_variation':'down'},
     'ELTR_up'               : {'ElTrigger_variation':'up'},
     'ELTR_down'             : {'ElTrigger_variation':'down'},
     'ELRECO_up'             : {'ElReco_variation':'up'},
     'ELRECO_down'           : {'ElReco_variation':'down'},
     'BTAG_up'               : {'BTag_variation':'up'},
     'BTAG_down'             : {'BTag_variation':'down'},
     # 'BTAG_bc_up'            : {'BTag_variation':'up_bc'},
     # 'BTAG_bc_down'          : {'BTag_variation':'down_bc'},
     # 'BTAG_udsg_up'          : {'BTag_variation':'up_udsg'},
     # 'BTAG_udsg_down'        : {'BTag_variation':'down_udsg'},
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
   <ConfigSGE RAM="2" DISK="2" Mail="dennis.schwarz@desy.de" Notification="as" Workdir="/nfs/dust/cms/user/schwarzd/PostSelEl_workdir"/>
-->
"""
        if os.path.exists(self.cwd + 'workdir'):
            opt = ' -rl --exitOnQuestion'
        else:
            opt = ' -sl --exitOnQuestion'

        self.exe = 'sframe_batch.py' + opt


sframe_tools = ToolChain(
    'SFrameUncertsSR',
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
