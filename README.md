# MTopJet
Measurement of the jet mass distribution in boosted top quark decays at 13 TeV


1) MTopJetPreSelectionModule

   just require 2 jets

   verly loose cuts on MET, lepton pT


2) MTopJetSelectionModule

   basic lepton+jets ttbar selection with MET, Lepton pT, lepton isolation ...


3) MTopJetPostSelectionModule

   select boosted ttbar

   define measurement phase-space

   define sideband regions

# Unfolding
to compile scripts do:

'make lib'

'make mtopjet'

To Run Unfolding:

1) ./binning fine

2) ./hist_filler fine true/false

3) ./SmearInput

4) ./do_unfolding data/pseudo1/pseudo2 fine
