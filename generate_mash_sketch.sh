#!/bin/bash

fold=/projects/data/team1_genomeAssembly/trimming
mkdir ${fold}/mash_sketches
cd $fold/mash_sketches


for f in $(cat /projects/home/shladyshau3/g1-assembly-team1/IDs.txt)
do
    echo $f
    cat ${fold}/fastq/${f}* > ${fold}/mash_sketches/${f}.fastq
    mash sketch -m 2 ${fold}/mash_sketches/${f}.fastq
    rm ${fold}/mash_sketches/${f}.fastq
done

