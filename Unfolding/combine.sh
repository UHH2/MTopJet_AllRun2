#!/bin/bash

rm Histograms_combine.root

hadd Histograms_combine.root Histograms_elec.root Histograms_muon.root
