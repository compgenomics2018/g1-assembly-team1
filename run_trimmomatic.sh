#!/bin/bash

out_fold=/projects/data/team1_genomeAssembly/trimming2
mkdir $out_fold


echo "Run trimmomatic for all files"

counter = 0

mkdir $out_fold/fastq

for f in $(cat /projects/home/shladyshau3/g1-assembly-team1/IDs.txt)
do
    if [$counter < 101]
    then
        echo $f
        java -jar /projects/home/shladyshau3/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
        /projects/data/team1_genomeAssembly/downloads/split_fastq/${f}_1.fastq \
        /projects/data/team1_genomeAssembly/downloads/split_fastq/${f}_2.fastq \
        /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_forward_paired.fq \
        /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_forward_unpaired.fq \
        /projects/data/team1_genomeAssembly/trimming2/fastq/${f}_reverse_paired.fq \
        /projects/data/team1_genomeAssembly/trimming/fastq/${f}_reverse_unpaired.fq \
        ILLUMINACLIP:/projects/home/shladyshau3/map_reads_to_adapters/adapters.fa:5:30:10 \
        SLIDINGWINDOW:4:20 MINLEN:20 >> $out_fold/trimmomatic.log
        count=`expr $count + 1`
    fi
done

echo "Running FastQC"

mkdir $out_fold/fastqc

for f in $(ls /projects/data/team1_genomeAssembly/trimming2/fastq)
do
    echo $f
    /projects/home/shladyshau3/FastQC/fastqc /projects/data/team1_genomeAssembly/trimming2/fastq/${f} \
    -o /projects/data/team1_genomeAssembly/trimming2/fastqc
done


echo "Run MultiQC"

mkdir $out_fold/multiqc

multiqc $out_fold/fastqc/*forward_paired* -O $out_fold/multiqc/forward_paired
multiqc $out_fold/fastqc/*reverse_paired* -O $out_fold/multiqc/reverse_paired
multiqc $out_fold/trimmomatic.log -O $out_fold/multiqc/trimming2


