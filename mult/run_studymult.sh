#!/bin/bash

currentDir=$PWD

append='false'

lastJob=0
nJobs=0
oDir=''

while getopts 'a' flag; do
    case "${flag}" in
        a)
            shift
            append='true';;
        *)
            nJobs=$1
            exit 1 ;;
    esac
done



nJobs=$1
shift

oDir=$1
shift

cd $oDir
if [[ "$append" == 'true' ]]; then 
    lastJob=$(ls -d $oDir/*/ | wc -l) ;
fi
echo "Writing output to " $oDir

mkdir -p $oDir || exit 1
cp /home/ge86rim/phsw/fempy/mult/StudyMultiplicity.cc $oDir || exit 1
cp /home/ge86rim/phsw/fempy/mult/exe_studymult.sh $oDir/exe.sh || exit 1

for ((iJob = $((1 + lastJob)); iJob <= $((nJobs + lastJob)); iJob++)); do
    mkdir -p $((iJob)) || continue
    
    cd $iJob || continue
    echo "job: $iJob sbatch --time=3:00:00 --exclude=kane,lambert --mem 2000 ../exe.sh $@"
    sbatch --time=3:00:00 --exclude=kane,lambert --mem 2000 ../exe.sh $@
    cd ..
done

cd $currentDir