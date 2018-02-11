#!/bin/bash

outp_fold=/projects/home/shladyshau3/pipeline_ref_based_assembly/output
ref_file=/projects/home/shladyshau3/reference_genomes/GCA_002893825.1_ASM289382v1_genomic.fna
inp_fold=/projects/home/shladyshau3/pipeline_ref_based_assembly/fastq
id=SRR3467249

mkdir $outp_fold
mkdir ${outp_fold}/alignment
mkdir ${outp_fold}/assembly
	

bwa index -a is $ref_file

bwa aln $ref_file ${inp_fold}/${id}_forward_paired.fq > ${outp_fold}/alignment/aln_1.sai

bwa aln $ref_file ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/alignment/aln_2.sai

bwa sampe $ref_file ${outp_fold}/alignment/aln_1.sai ${outp_fold}/alignment/aln_2.sai ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/alignment/aln.sam

samtools sort ${outp_fold}/alignment/aln.sam > ${outp_fold}/alignment/aln_sorted.bam

samtools index ${outp_fold}/alignment/aln_sorted.bam

samtools view -b -f 4 ${outp_fold}/alignment/aln_sorted.bam > ${outp_fold}/alignment/aln_unmapped.bam

samtools bam2fq ${outp_fold}/alignment/aln_unmapped.bam > ${outp_fold}/alignment/unmapped_reads.fastq

samtools mpileup -v --no-BAQ -f $ref_file ${outp_fold}/alignment/aln_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A > ${outp_fold}/assembly/assembly.fasta 

quast.py -t 6 -m 1000 ${outp_fold}/assembly/assembly.fasta -o ${outp_fold}/assembly/quast_report

rm ${outp_fold}/alignment/aln_1.sai

rm ${outp_fold}/alignment/aln_2.sai

rm ${outp_fold}/alignment/aln.sam

rm ${outp_fold}/alignment/aln_unmapped.bam


