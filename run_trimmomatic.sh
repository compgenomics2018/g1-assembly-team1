#!/bin/bash

echo "Run trimmomatic for all files"

for f in $(ls /projects/data/team1_genomeAssembly/downloads/fastq)
do
    bn=$(basename ${f} | cut -d . -f 1)
    echo $bn
    java -jar /projects/home/shladyshau3/Trimmomatic-0.36/trimmomatic-0.36.jar PE /projects/data/team1_genomeAssembly/downloads/split_fastq/${bn}_1.fastq /projects/data/team1_genomeAssembly/downloads/split_fastq/${bn}_2.fastq /projects/data/team1_genomeAssembly/trimming/${bn}_forward_paired.fq /projects/data/team1_genomeAssembly/trimming/${bn}_forward_unpaired.fq /projects/data/team1_genomeAssembly/trimming/${bn}_reverse_paired.fq /projects/data/team1_genomeAssembly/trimming/${bn}_reverse_unpaired.fq ILLUMINACLIP:/projects/home/shladyshau3/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:20 CROP:245 HEADCROP:20 
done
