#!/bin/bash

echo "Run FasQC"

for f in $(ls /projects/data/team1_genomeAssembly/trimming/fastq)
do
    echo $f
    /projects/home/shladyshau3/FastQC/fastqc /projects/data/team1_genomeAssembly/trimming/fastq/${f} -o /projects/data/team1_genomeAssembly/trimming/fastqc/
done
