#!/bin/bash


echo "Run trimmomatic for all files"

for f in $(cat /projects/home/shladyshau3/g1-assembly-team1/IDs.txt)
do
    echo $f
    java -jar /projects/home/shladyshau3/Trimmomatic-0.36/trimmomatic-0.36.jar PE -trimlog /projects/data/team1_genomeAssembly/trimming2/log/${f}.log /projects/data/team1_genomeAssembly/downloads/split_fastq/${f}_1.fastq /projects/data/team1_genomeAssembly/downloads/split_fastq/${f}_2.fastq /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_forward_paired.fq /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_forward_unpaired.fq /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_reverse_paired.fq /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_reverse_unpaired.fq ILLUMINACLIP:/projects/home/shladyshau3/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 CROP:245 HEADCROP:15 SLIDINGWINDOW:4:20 MINLEN:20 
done



echo "Run FastQC"

for f in $(ls /projects/data/team1_genomeAssembly/trimming)
do
    echo $f
    /projects/home/shladyshau3/FastQC/fastqc /projects/data/team1_genomeAssembly/trimming2/fastq/${f} -o /projects/data/team1_genomeAssembly/trimming2/fastqc
done
