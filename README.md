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