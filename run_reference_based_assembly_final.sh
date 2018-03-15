#!/bin/bash

inp_fold=/projects/data/team1_genomeAssembly/trimming2/fastq/trimmed
outp_fold=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all
mkdir $outp_fold
best_gen_file=/projects/data/team1_genomeAssembly/reference_based_assembly/best_genomes_all.txt
	
mkdir $outp_fold

mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/assembly_selected
mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final
mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final_contigs

for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo $id
    mkdir ${outp_fold}/${id}
    ref_file=`grep $id ${best_gen_file} | awk '{print $2}'`
    echo $ref_file
    bwa index -a is $ref_file
    bwa mem $ref_file ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq > ${outp_fold}/${id}/aln.sam
    samtools sort ${outp_fold}/${id}/aln.sam > ${outp_fold}/${id}/aln_sorted.bam
    samtools index ${outp_fold}/${id}/aln_sorted.bam
    samtools view -b -f 4 ${outp_fold}/${id}/aln_sorted.bam > ${outp_fold}/${id}/aln_unmapped.bam
    samtools bam2fq ${outp_fold}/${id}/aln_unmapped.bam > ${outp_fold}/${id}/unmapped_reads.fastq
    samtools mpileup -v --no-BAQ -f $ref_file ${outp_fold}/${id}/aln_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A > ${outp_fold}/${id}/assembly.fasta 
    quast.py -t 6 -m 1000 ${outp_fold}/${id}/assembly.fasta -o ${outp_fold}/${id}/quast_report

    awk '{if(!/>/){gsub("n","N"); print $0}else{print $0}}' ${outp_fold}/${id}/assembly.fasta > ${outp_fold}/${id}/assembly_n2N.fasta
    cat ${outp_fold}/${id}/assembly.fasta | awk '!/>/ {print}' |  grep -o 'N*' | awk '{print length($0)}'> ${outp_fold}/${id}/assembly_n2N_length_of_Ns.txt
    fasta_split_scaffolds.py ${outp_fold}/${id}/assembly_n2N.fasta 10 ${outp_fold}/${id}/assembly_n2N_splitted.fasta
    cat ${inp_fold}/${id}_forward_unpaired.fq ${inp_fold}/${id}_reverse_unpaired.fq > ${outp_fold}/${id}/unpaired_reads.fq
    cd ${outp_fold}/${id}
    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 500 0.75 FR" > ${outp_fold}/${id}/lib_500.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_500.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_500

    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 100 0.75 FR" > ${outp_fold}/${id}/lib_100.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_100.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_100

    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 300 0.75 FR" > ${outp_fold}/${id}/lib_300.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_300.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_300

    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 700 0.75 FR" > ${outp_fold}/${id}/lib_700.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_700.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_700

    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 900 0.75 FR" > ${outp_fold}/${id}/lib_900.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_900.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_900

    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 1100 0.75 FR" > ${outp_fold}/${id}/lib_1100.txt
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_1100.txt -s ${outp_fold}/${id}/assembly_n2N_splitted.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b sspace_insert_1100

    
    fasta_split_scaffolds.py ${outp_fold}/${id}/assembly_n2N.fasta 4 ${outp_fold}/${id}/assembly_n2N_splitted_by_4.fasta
    skesa --fastq ${outp_fold}/${id}/unmapped_reads.fastq --contigs_out ${outp_fold}/${id}/unmapped_skesa_assembly.fasta
    
    cat ${outp_fold}/${id}/assembly_n2N_splitted_by_4.fasta ${outp_fold}/${id}/unmapped_skesa_assembly.fasta > ${outp_fold}/${id}/ref_and_denovo_contigs.fasta
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_500.txt -s ${outp_fold}/${id}/ref_and_denovo_contigs.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b ref_and_denovo

    cat ${outp_fold}/${id}/assembly_n2N_splitted.fasta ${outp_fold}/${id}/unmapped_skesa_assembly.fasta > ${outp_fold}/${id}/ref_10_and_denovo_contigs.fasta
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_500.txt -s ${outp_fold}/${id}/ref_10_and_denovo_contigs.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b ref_10_and_denovo

    quast.py -t 6 -s -m 500 -s ${outp_fold}/${id}/ref_10_and_denovo.final.scaffolds.fasta -o ${outp_fold}/${id}/quast_report_ref_10_and_denovo
    mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final/${id}
    cp -r ${outp_fold}/${id}/quast_report_ref_10_and_denovo/* /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final/${id}/
    
    quast.py -t 6 -m 500 -s ${outp_fold}/${id}/ref_10_and_denovo.final.scaffolds.fasta -o ${outp_fold}/${id}/quast_report_ref_10_and_denovo
    mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final/${id}
    cp -r ${outp_fold}/${id}/quast_report_ref_10_and_denovo/* /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final/${id}/

    quast.py -t 6 -m 500 ${outp_fold}/${id}/ref_and_denovo_contigs.fasta -o ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs
    mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final_contigs/${id}
    cp -r ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/* /projects/data/team1_genomeAssembly/reference_based_assembly/quast_final_contigs/${id}/

    rm -r ${outp_fold}/${id}/reads/*
    rm -r ${outp_fold}/${id}/bowtieoutput/*
    rm -r ${outp_fold}/${id}/intermediate_results/*

    rm ${outp_fold}/${id}/aln.sam
    rm ${outp_fold}/${id}/aln_sorted.bam
    rm ${outp_fold}/${id}/aln_sorted.bam.bai
    rm ${outp_fold}/${id}/aln_unmapped.bam
    rm ${ref_file}.amb
    rm ${ref_file}.ann
    rm ${ref_file}.bwt
    rm ${ref_file}.fai
    rm ${ref_file}.pac
    rm ${ref_file}.sa

    mkdir /projects/data/team1_genomeAssembly/reference_based_assembly/assembly_selected/${id}
    cp ${outp_fold}/${id}/ref_10_and_denovo* /projects/data/team1_genomeAssembly/reference_based_assembly/assembly_selected/${id}/

done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/N50_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'N50' | awk '{print $3}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/N50_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_scaffolds_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of scaffolds' | awk '{print $6}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_scaffolds_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_Ns_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of N' | awk '{print $6}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_Ns_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/max_scaffold_length_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Max scaffold size' | awk '{print $5}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/max_scaffold_length_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/GC_content_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/report.txt | grep 'GC (%)' | awk '{print $3}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/GC_content_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/total_length_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/report.txt | grep 'Total length (>= 0 bp)'| awk '!/All statistics/'  | awk '{print $6}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/total_length_summary.txt
done


