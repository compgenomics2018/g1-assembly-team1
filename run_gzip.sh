#!/bin/bash

echo "Run FasQC"

for f in $(ls /projects/data/team1_genomeAssembly/trimming)
do
    echo $f
    gzip /projects/data/team1_genomeAssembly/trimming/${f}
done
