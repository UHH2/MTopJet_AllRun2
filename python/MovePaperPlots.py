import sys, os, ROOT

def MovePaperPlots():
    path_in  = '/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/macros/plots/PaperPlots/'
    path_out = '/afs/desy.de/user/p/paaschal/WorkingArea/CMSNotes/Papers/TOP-21-012/figs/'
    files = {
        'chi2_JMS.pdf':      'chi2_JMS.pdf',
        'tau32_2016.pdf':    'fsr_tau32_2016.pdf',
        'tau32_combine.pdf': 'fsr_tau32_combine.pdf',
        'wmass_hl.pdf':      'Wjet_mass_all_hl_combine_norm.pdf',
        'wmass_ll.pdf':      'Wjet_mass_all_ll_combine_norm.pdf',
        'wmass_hh.pdf':      'Wjet_mass_all_hh_combine_norm.pdf',
        'wmass_lh.pdf':      'Wjet_mass_all_lh_combine_norm.pdf'
        }

    for key in files:
        file_in  = path_in+key
        file_out = path_out+files[key]
        command = 'cp '+file_in+' '+file_out
        print command
        os.system(command)


if __name__=='__main__':
    MovePaperPlots()
