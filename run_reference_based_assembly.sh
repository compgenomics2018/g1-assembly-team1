#!/bin/bash

inp_fold=/projects/data/team1_genomeAssembly/trimming2/fastq/trimmed
outp_fold=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all
mkdir $outp_fold
best_gen_file=/projects/data/team1_genomeAssembly/reference_based_assembly/best_genomes_all.txt
	
mkdir $outp_fold
	
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo $id
    #mkdir ${outp_fold}/${id}
    #ref_file=`grep $id ${best_gen_file} | awk '{print $2}'`
    #echo $ref_file
    #bwa index -a is $ref_file
    #bwa mem $ref_file ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/${id}/aln.sam
    #samtools sort ${outp_fold}/${id}/aln.sam > ${outp_fold}/${id}/aln_sorted.bam
    #samtools index ${outp_fold}/${id}/aln_sorted.bam
    #samtools view -b -f 4 ${outp_fold}/${id}/aln_sorted.bam > ${outp_fold}/${id}/aln_unmapped.bam
    #samtools bam2fq ${outp_fold}/${id}/aln_unmapped.bam > ${outp_fold}/${id}/unmapped_reads.fastq
    #samtools mpileup -v --no-BAQ -f $ref_file ${outp_fold}/${id}/aln_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A > ${outp_fold}/${id}/assembly.fasta 
    #quast.py -t 6 -m 1000 ${outp_fold}/${id}/assembly.fasta -o ${outp_fold}/${id}/quast_report

    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${outp_fold}/${id}/assembly.fasta  > ${outp_fold}/${id}/assembly_upper.fasta
    fasta_split_scaffolds.py ${outp_fold}/${id}/assembly_upper.fasta 20 ${outp_fold}/${id}/assembly_upper_splitted.fasta
    echo "lib ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 250 0.75 FR" > ${outp_fold}/${id}/lib.txt
    cd ${outp_fold}/${id}
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib.txt -s ${outp_fold}/${id}/assembly_upper_splitted.fasta -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_extend

    #bwa index -a is ${outp_fold}/${id}/sspace_extend.final.scaffolds.fasta
    #bwa mem ${outp_fold}/${id}/sspace_extend.final.scaffolds.fasta ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/${id}/aln_scaffold.sam
    #samtools sort ${outp_fold}/${id}/aln_scaffold.sam >${outp_fold}/${id}/aln_scaffold_sorted.bam
    #samtools index ${outp_fold}/${id}/aln_scaffold_sorted.bam
    #samtools view -b -f 4 ${outp_fold}/${id}/aln_scaffold_sorted.bam > ${outp_fold}/${id}/aln_unmapped_scaffold.bam
    #samtools bam2fq ${outp_fold}/${id}/aln_unmapped_scaffold.bam > ${outp_fold}/${id}/unmapped_scaffold_reads.fastq

    #rm ${outp_fold}/${id}/aln.sam
    #rm ${outp_fold}/${id}/aln_sorted.bam
    #rm ${outp_fold}/${id}/aln_sorted.bam.bai
    #rm ${outp_fold}/${id}/aln_unmapped.bam
    #rm ${ref_file}.amb
    #rm ${ref_file}.ann
    #rm ${ref_file}.bwt
    #rm ${ref_file}.fai
    #rm ${ref_file}.pac
    #rm ${ref_file}.sa

done
