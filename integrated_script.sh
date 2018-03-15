#!/bin/bash

##Script for Trimming
mkdir temp_trimming2
mypwd=`pwd`
input_fold=$mypwd/fastqRaw # fastqRaw is a folder created before storing raw fastq files
out_fold=$mypwd/temp_trimming2
trimmomatic_path=/projects/data/team1_genomeAssembly/local/bin
adapters_path=/projects/data/team1_genomeAssembly/trimming2
echo "Run trimmomatic for all files"
mkdir $out_fold/fastq


for f in $(cat $mypwd/ID.txt)
do
    java -jar ${trimmomatic_path}/trimmomatic-0.36.jar PE \
    ${input_fold}/${f}_1.fastq \
    ${input_fold}/${f}_2.fastq \
    ${out_fold}/fastq/${f}_forward_paired.fq \
    ${out_fold}/fastq/${f}_forward_unpaired.fq \
    ${out_fold}/fastq/${f}_reverse_paired.fq \
    ${out_fold}/fastq/${f}_reverse_unpaired.fq \
    ILLUMINACLIP:${adapters_path}/adapters.fa:5:30:10 \
    SLIDINGWINDOW:4:20 MINLEN:20
done

echo "Running FastQC"

mkdir $out_fold/fastq/fastqc

for f in $(ls ${out_fold}/fastq)
do
    echo $f
    fastqc ${out_fold}/fastq/${f} \
    -o ${out_fold}/fastq/fastqc
done

##End of script for trimming

###Script for Generating MASH tree
mkdir mash_sketches_samples
inp_fold1=${out_fold}/fastq
outp_fold1=$mypwd/mash_sketches_samples
ref_gen_path=$mypwd/reference_genomes # may need attention for test

for f in `ls $inp_fold1 | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo $f
    cat ${inp_fold1}/${f}* > ${outp_fold1}/${f}.fastq
    mash sketch -m 2 ${outp_fold1}/${f}.fastq
    rm ${outp_fold1}/${f}.fastq
done

mash paste ${outp_fold1}/all ${outp_fold1}/*

mash sketch ${ref_gen_path}/*.fna -o $mypwd/all_genomes

mash dist ${outp_fold1}/all.msh $mypwd/all_genomes.msh > $mypwd/ref_gen_vs_sample_all.dist


cat "" > $mypwd/best_genomes_all.txt
for f in `ls $inp_fold1 | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    g=`cat $mypwd/ref_gen_vs_sample_all.dist | awk '/'${f}'/ {print}' | sort -k3 -n | head -n1 | awk '{print $2}'`

    echo ${f} ${g}
    echo ${f} ${g} >> $mypwd/best_genomes_all.txt
   # echo $f "\t" cat ref_gen_vs_sample.dist | awk '/${f}/ {print}' | awk 'BEGIN{a=1000}{if ($3<0+a) a=$3} END{print $2}' | awk '{split($1,a,"/"); print a[3]}'
done

##End of script for mash tree

##Script for Running Reference Assembly
inp_fold=${out_fold}/fastq
outp_fold=$mypwd/assembly_all
mkdir $outp_fold
best_gen_file=$mypwd/best_genomes_all.txt


mkdir $mypwd/quast_final
mkdir $mypwd/quast_final_contigs

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
    fasta_split_scaffolds.py ${outp_fold}/${id}/assembly_n2N.fasta 10 ${outp_fold}/${id}/assembly_n2N_splitted.fasta
    cat ${inp_fold}/${id}_forward_unpaired.fq ${inp_fold}/${id}_reverse_unpaired.fq > ${outp_fold}/${id}/unpaired_reads.fq

 cd ${outp_fold}/${id}
    echo "test path: " `pwd`
    echo "lib1 ${inp_fold}/${id}_forward_paired.fq ${inp_fold}/${id}_reverse_paired.fq 500 0.75 FR" > ${outp_fold}/${id}/lib_500.txt
    skesa --fastq ${outp_fold}/${id}/unmapped_reads.fastq --contigs_out ${outp_fold}/${id}/unmapped_skesa_assembly.fasta
    cat ${outp_fold}/${id}/assembly_n2N_splitted.fasta ${outp_fold}/${id}/unmapped_skesa_assembly.fasta > ${outp_fold}/${id}/ref_10_and_denovo_contigs.fasta
    SSPACE_Basic.pl -l ${outp_fold}/${id}/lib_500.txt -s ${outp_fold}/${id}/ref_10_and_denovo_contigs.fasta -u ${outp_fold}/${id}/unpaired_reads.fq -x 1 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 4 -b ref_10_and_denovo

    quast.py -t 6 -m 500 -s ${outp_fold}/${id}/ref_10_and_denovo.final.scaffolds.fasta -o ${outp_fold}/${id}/quast_report_ref_10_and_denovo
    mkdir $mypwd/quast_final/${id}
    cp -r ${outp_fold}/${id}/quast_report_ref_10_and_denovo/* $mypwd/quast_final/${id}/

    quast.py -t 6 -m 500 ${outp_fold}/${id}/ref_10_and_denovo_contigs.fasta -o ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs
    mkdir $mypwd/quast_final_contigs/${id}
    cp -r ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/* $mypwd/quast_final_contigs/${id}/

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

done

cat '' > $mypwd/N50_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'N50' | awk '{print $3}'` >> $mypwd/N50_summary.txt
done

cat '' > $mypwd/number_of_scaffolds_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of scaffolds' | awk '{print $6}'` >> $mypwd/number_of_scaffolds_summary.txt
done

cat '' > $mypwd/number_of_Ns_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of N' | awk '{print $6}'` >> $mypwd/number_of_Ns_summary.txt
done

cat '' > $mypwd/max_scaffold_length_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Max scaffold size' | awk '{print $5}'` >> $mypwd/max_scaffold_length_summary.txt
done

cat '' > $mypwd/GC_content_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/report.txt | grep 'GC (%)' | awk '{print $3}'` >> $mypwd/GC_content_summary.txt
done

cat '' > $mypwd/total_length_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/quast_report_ref_10_and_denovo_contigs/report.txt | grep 'Total length (>= 0 bp)'| awk '!/All statistics/'  | awk '{print $6}'` >> $mypwd/total_length_summary.txt
done
