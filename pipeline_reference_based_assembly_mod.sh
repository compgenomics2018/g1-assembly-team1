#!/bin/bash

outp_fold=/projects/home/shladyshau3/pipeline_ref_based_assembly/output_mod2
ref_file=/projects/home/shladyshau3/reference_genomes/GCA_002893825.1_ASM289382v1_genomic.fna
inp_fold=/projects/data/team1_genomeAssembly/trimming/fastq
id=SRR3467249

mkdir $outp_fold
	

bwa index -a is $ref_file

bwa mem $ref_file ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/aln.sam

samtools sort ${outp_fold}/aln.sam > ${outp_fold}/aln_sorted.bam

samtools index ${outp_fold}/aln_sorted.bam

samtools view -b -f 4 ${outp_fold}/aln_sorted.bam > ${outp_fold}/aln_unmapped.bam

samtools bam2fq ${outp_fold}/aln_unmapped.bam > ${outp_fold}/unmapped_reads.fastq

samtools mpileup -v --no-BAQ -f $ref_file ${outp_fold}/aln_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A > ${outp_fold}/assembly.fasta 

quast.py -t 6 -m 1000 ${outp_fold}/assembly.fasta -o ${outp_fold}/quast_report

rm ${outp_fold}/aln.sam

rm ${outp_fold}/aln_sorted.bam

rm ${outp_fold}/aln_sorted.bam.bai

rm ${outp_fold}/aln_unmapped.bam


