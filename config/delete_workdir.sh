#!/bin/bash

if [ -z "$1" ]; then # Check if first argument is not empty !!!!
    echo "No argument supplied"
    return # exit closes zsh shell
fi

if [[ $1 = "workdir" ]]; then # Check if first argument is equal to

  if [ \! -z "$2" ]; then # negation
    echo "Deleting all SYS workdirs for $2 channel"
    rm -r /nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/$2/*/workdir
  else
    echo "I am supposed to delete all SYS workdirs, but no channel is given"
    return
  fi

else # Workdirs and SFrameUncerts

  if [[ $1 = *"ElecSF"* ]]; then
    echo "Deleting /nfs/dust/cms/user/paaschal/MTopJet_Run2/ElecSF/$1"
    rm -r /nfs/dust/cms/user/paaschal/MTopJet_Run2/ElecSF/$1
    echo "Deleting ./$1"
    rm -r /nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/config/$1
    return
  fi

  channel="WRONG CHANNEL"
  if [[ $1 = *"Mu"* ]] || [[ $2 = "muon" ]]; then
    channel="muon"
    if [[ $1 = *"JMSelec"* ]]; then
      channel="muon_elecJMS"
    fi
  elif [[ $1 = *"El"* ]] || [[ $2 = "elec" ]]; then
    channel="elec"
    if [[ $1 = *"JMSmuon"* ]]; then
      channel="elec_muonJMS"
    fi
  fi

  year="WRONG YEAR"
  if [[ $1 = *"2016"* ]]; then
    year="2016"
  elif [[ $1 = *"2017"* ]]; then
    year="2017"
  elif [[ $1 = *"2018"* ]]; then
    year="2018"
  fi

  echo "$year/$channel"

  if [[ ! $channel = "WRONG CHANNEL" ]] || [[ ! $year = "WRONG YEAR" ]]; then # negation; different structur of [<...>] and [[<...>]]

    if [ -d "./$1" ]; then # Check if folder exists
      echo "Deleting ./$1"
      rm -r /nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/config/$1
    fi

    if [ -d "/nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/$channel/$1" ]; then # Check if folder exists
      echo "Deleting /nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/$channel/$1"
      rm -r /nfs/dust/cms/user/paaschal/MTopJet_Run2/PostSel/$channel/$1
    fi

  fi
fi

return
