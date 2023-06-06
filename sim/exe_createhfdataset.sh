#!/bin/bash

nEvents=$1
maxRunTime=$2
energy=$3
trigger=$4
tune=$5
process=$6
hfPdg=$7
kinem=$8

oDir=$PWD
seed=${PWD##*/}

echo "../CreateHFDataset.cc($nEvents, $maxRunTime, $energy, $trigger, $tune, $process, $hfPdg, $kinem, $seed, \"$oDir\")"
time root -l -q -b "../CreateHFDataset.cc($nEvents, $maxRunTime, $energy, $trigger, $tune, $process, $hfPdg, $kinem, $seed, \"$oDir\")"
