#!/bin/bash

nEvents=$1
maxRunTime=$2
energy=$3
tune=$4
process=$5

oDir=$PWD
seed=${PWD##*/}

echo "../StudyMultiplicity.cc($nEvents, $maxRunTime, $energy, $tune, $process, $seed, \"$oDir\")"
time root -l -q -b "../StudyMultiplicity.cc($nEvents, $maxRunTime, $energy, $tune, $process, $seed, \"$oDir\")"
