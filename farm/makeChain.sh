#!/bin/bash

datasets=($(ls | grep 'seed-1$'))

for dataset in "${datasets[@]}"; do
    pattern=${dataset[@]::-7}
    echo $pattern

    find . -wholename "./${pattern}*/tracksummary_ambi.root" -exec realpath {} \; > ${pattern}.chain
done
