#!/bin/bash

wget ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

awk -F "\t" '/Klebsiella/ && /Complete/ && !/phage/ { split($20,a,"/"); print $20 "/" a[10] "_genomic.fna.gz" }' assembly_summary_refseq.txt > download_list.txt

mkdir genomes_download
wget -i ./download_list.txt -P ./genomes_download


for f in `ls genomes_download`;
do
    nf=`echo $f $f | cut -f 1,2,3 -d '.'`
    echo $nf
    gunzip genomes_download/$f
done 
