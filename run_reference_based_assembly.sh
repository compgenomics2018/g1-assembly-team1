#!/bin/bash

out_fold=/projects/data/team1_genomeAssembly/ref_assembly
mkdir $out_fold


echo "Run reference based assembly"


for f in $(cat /projects/home/shladyshau3/g1-assembly-team1/IDs_short.txt)
do
    echo $f
    perl reference_based_assembly.pl -1 /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_forward_paired.fq -2 /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_reverse_paired.fq -r /projects/home/shladyshau3/reference_genomes/GCA_002893825.1_ASM289382v1_genomic.fna -d /projects/data/team1_genomeAssembly/ref_assembly/${f} -o ${f}_ref_assembly.fastq
done

